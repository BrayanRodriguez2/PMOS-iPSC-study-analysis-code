#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/03_sample_qc.py   (original step script: 09_sample_qc_decide.py)
# Produces: Sample QC (call rate, heterozygosity, sex, kinship) and the list of excluded samples.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 09 (decision) - turn the per-sample plink2 metrics into keep / exclude decisions.

Thresholds (standard for array genotype QC, Anderson et al. 2010 Nat Protoc 5:1564;
KING: Manichaikul et al. 2010 Bioinformatics 26:2867):
  call rate          >= 0.98        sample call rate on autosomal variants
  heterozygosity F   within mean +/- 3 SD of the cohort (autosomal, LD-pruned, MAF >= 5 %).
                     Strongly positive F = excess homozygosity / inbreeding / low DNA quality;
                     strongly negative F = excess heterozygosity / sample contamination.
  sex                the cohort is female (endometrial stromal fibroblast donors, .fam sex=2).
                     A sample is flagged if its chrX F statistic does NOT call it female
                     (SNPSEX != 2), i.e. male or ambiguous (F between 0.2 and 0.8).
                     A sample with sex=0 (unknown) in the .fam that is genetically female is
                     NOT excluded: only the pedigree field was empty.
  relatedness        KING kinship >= 0.177 (first degree; >= 0.354 = duplicate / MZ twin):
                     one of each pair is removed (the one with the LOWER call rate; ties ->
                     larger IID). Pairs with 0.0884 <= kinship < 0.177 (second degree) are
                     only REPORTED: for a descriptive PRS in 138 donors they do not bias anything.

Inputs : .fam, plink2 .smiss (autosomal call rate), .smiss for chrY, .het, .sexcheck, .kin0,
         donor_key.tsv (the 14 donors with iPSC lines; used only to report overlaps).
Outputs: sample_qc_table.tsv (every sample, all metrics and flags), samples_excluded.tsv
         (excluded samples and reason), related_pairs.tsv, keep_samples.txt (#FID IID, for
         plink2 --keep) and qc_samples_report.txt.
Standard library only.
"""
import argparse
import os
import statistics
import sys

MIN_CALL_RATE = 0.98
HET_SD = 3.0
KIN_FIRST_DEGREE = 0.177
KIN_REPORT = 0.0884


def read_plink_table(path):
    """Read a plink2 tab-delimited table whose header starts with '#'. Returns list of dicts."""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, f)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("fam", "smiss", "ysmiss", "het", "sexcheck", "kin0", "donor_key", "outdir"):
        ap.add_argument("--" + a.replace("_", "-"), required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rep = []

    def log(msg=""):
        print(msg, flush=True)
        rep.append(msg)

    # ---- load -------------------------------------------------------------
    samples = {}  # IID -> dict
    with open(args.fam) as fh:
        for line in fh:
            f = line.split()
            samples[f[1]] = {"FID": f[0], "IID": f[1], "pedsex": f[4], "flags": []}
    for r in read_plink_table(args.smiss):
        s = samples[r["IID"]]
        s["call_rate"] = 1.0 - float(r["F_MISS"])
    for r in read_plink_table(args.ysmiss):
        samples[r["IID"]]["chrY_call_rate"] = 1.0 - float(r["F_MISS"])
    het = {r["IID"]: float(r["F"]) for r in read_plink_table(args.het)}
    mu, sd = statistics.mean(het.values()), statistics.stdev(het.values())
    for iid, s in samples.items():
        s["het_F"] = het[iid]
        s["het_z"] = (het[iid] - mu) / sd
    for r in read_plink_table(args.sexcheck):
        s = samples[r["IID"]]
        s["snpsex"] = r["SNPSEX"]
        s["xF"] = float(r["F"]) if r["F"] not in ("nan", "NA") else float("nan")
    donors = {}
    with open(args.donor_key) as fh:
        next(fh)
        for line in fh:
            line_f = line.rstrip("\n").split("\t")
            donors[line_f[1]] = (line_f[0], line_f[2])

    log(f"[samples] {len(samples)} in .fam | pedsex codes: "
        f"{ {k: sum(1 for s in samples.values() if s['pedsex'] == k) for k in sorted({s['pedsex'] for s in samples.values()}) } }")
    log(f"[thresholds] call rate >= {MIN_CALL_RATE} | |het F - mean| <= {HET_SD} SD (mean {mu:.4f}, SD {sd:.4f}) | "
        f"female by chrX F | KING first degree >= {KIN_FIRST_DEGREE}")

    # ---- per-sample filters ------------------------------------------------
    for s in samples.values():
        if s["call_rate"] < MIN_CALL_RATE:
            s["flags"].append("low_call_rate")
        if abs(s["het_z"]) > HET_SD:
            s["flags"].append("heterozygosity_outlier_" + ("high" if s["het_z"] > 0 else "low"))
        if s["snpsex"] != "2":
            s["flags"].append("not_female_by_chrX_F")
        elif s["pedsex"] != "2":
            s.setdefault("notes", []).append("pedsex_missing_but_genetically_female")
        if s["pedsex"] == "2" and s["snpsex"] == "1":
            s["flags"].append("sex_mismatch_ped_female_snp_male")

    # ---- relatedness ---------------------------------------------------------
    pairs = []
    for r in read_plink_table(args.kin0):
        pairs.append((float(r["KINSHIP"]), r["IID1"], r["IID2"]))
    pairs.sort(reverse=True)
    with open(os.path.join(args.outdir, "related_pairs.tsv"), "w") as out:
        out.write("kinship\tIID1\tIID2\tcategory\n")
        for k, a, b in pairs:
            cat = "duplicate_or_MZ" if k >= 0.354 else "first_degree" if k >= KIN_FIRST_DEGREE else "second_degree_reported_only"
            out.write(f"{k:.4f}\t{a}\t{b}\t{cat}\n")
    dropped_by_kin = set()
    for k, a, b in pairs:
        if k < KIN_FIRST_DEGREE:
            break
        if a in dropped_by_kin or b in dropped_by_kin:
            continue
        # drop the sample with the lower call rate; ties -> larger IID (deterministic)
        ca, cb = samples[a]["call_rate"], samples[b]["call_rate"]
        drop = a if (ca, -int(a)) < (cb, -int(b)) else b
        samples[drop]["flags"].append(f"related_to_{b if drop == a else a}_kinship_{k:.3f}")
        dropped_by_kin.add(drop)
    log(f"[kinship] pairs >= {KIN_REPORT}: {len(pairs)} | >= {KIN_FIRST_DEGREE}: {sum(1 for p in pairs if p[0] >= KIN_FIRST_DEGREE)} "
        f"| samples dropped for relatedness: {len(dropped_by_kin)}")

    # ---- outputs -------------------------------------------------------------
    cols = ["FID", "IID", "pedsex", "snpsex", "xF", "call_rate", "chrY_call_rate", "het_F", "het_z"]
    with open(os.path.join(args.outdir, "sample_qc_table.tsv"), "w") as out:
        out.write("\t".join(cols + ["is_donor_of_ipsc_line", "excluded", "reasons", "notes"]) + "\n")
        for iid in sorted(samples, key=int):
            s = samples[iid]
            d = f"{donors[iid][0]}({donors[iid][1]})" if iid in donors else ""
            vals = [f"{s.get(c, float('nan')):.4f}" if isinstance(s.get(c), float) else str(s.get(c, "NA")) for c in cols]
            out.write("\t".join(vals + [d, "yes" if s["flags"] else "no", ";".join(s["flags"]), ";".join(s.get("notes", []))]) + "\n")
    excluded = [s for s in samples.values() if s["flags"]]
    kept = [s for s in samples.values() if not s["flags"]]
    with open(os.path.join(args.outdir, "samples_excluded.tsv"), "w") as out:
        out.write("FID\tIID\tpedsex\tsnpsex\txF\tcall_rate\thet_F\thet_z\tdonor_of_ipsc_line\treasons\n")
        for s in sorted(excluded, key=lambda x: int(x["IID"])):
            d = f"{donors[s['IID']][0]}({donors[s['IID']][1]})" if s["IID"] in donors else ""
            out.write(f"{s['FID']}\t{s['IID']}\t{s['pedsex']}\t{s['snpsex']}\t{s['xF']:.3f}\t{s['call_rate']:.4f}\t"
                      f"{s['het_F']:.4f}\t{s['het_z']:.2f}\t{d}\t{';'.join(s['flags'])}\n")
    with open(os.path.join(args.outdir, "keep_samples.txt"), "w") as out:
        out.write("#FID\tIID\n")
        for s in sorted(kept, key=lambda x: int(x["IID"])):
            out.write(f"{s['FID']}\t{s['IID']}\n")

    log("")
    log(f"[result] samples kept: {len(kept)} | excluded: {len(excluded)}")
    for s in sorted(excluded, key=lambda x: int(x["IID"])):
        log(f"  excluded IID {s['IID']} (FID {s['FID']}): {';'.join(s['flags'])} | call rate {s['call_rate']:.4f} | "
            f"het F {s['het_F']:.3f} (z {s['het_z']:.2f}) | chrX F {s['xF']:.3f} | snpsex {s['snpsex']} | pedsex {s['pedsex']}"
            + (f" | donor of {donors[s['IID']][0]}" if s["IID"] in donors else ""))
    amb = sorted(i for i, s in samples.items() if s["pedsex"] == "0")
    exc_ids = sorted(s["IID"] for s in excluded)
    log("")
    log(f"[check vs the manuscript authors (138 valid of 140)] samples with sex=0 in .fam: {amb} | excluded by this QC: {exc_ids} | "
        f"same set: {'YES' if set(amb) == set(exc_ids) else 'NO'}")
    for i in amb:
        s = samples[i]
        log(f"  sex=0 sample IID {i}: chrX F {s['xF']:.3f} -> SNPSEX {s['snpsex']} | call rate {s['call_rate']:.4f} | "
            f"het F {s['het_F']:.3f} | chrY call rate {s.get('chrY_call_rate', float('nan')):.3f}")
    log(f"[donors of the 14 iPSC lines] excluded: {[d for d in exc_ids if d in donors]} | present in the keep list: "
        f"{sum(1 for s in kept if s['IID'] in donors)} of {len(donors)}")
    ys = sorted(s.get("chrY_call_rate", float("nan")) for s in samples.values())
    log(f"[chrY call rate across samples] min {ys[0]:.3f} | median {statistics.median(ys):.3f} | max {ys[-1]:.3f}")
    with open(os.path.join(args.outdir, "qc_samples_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
