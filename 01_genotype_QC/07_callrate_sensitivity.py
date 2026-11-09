#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/07_callrate_sensitivity.py   (original step script: 09c_callrate_sensitivity.py)
# Produces: Sensitivity of the sample set to the call-rate threshold.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 09c - How many samples would step 09 exclude with a different sample call-rate threshold?

Nothing is recomputed with plink2. The per-sample metrics of step 09 (sample_qc_table.tsv) do not depend
on the call-rate threshold: heterozygosity F uses mean and SD of ALL 140 samples (recovered from het_z) and
the relatedness filter only drops the lower-call-rate member of the pairs with kinship >= 0.177 (the one pair
of this kind, whose members are far from any call-rate threshold considered here). Therefore changing the call-rate
threshold while KEEPING the heterozygosity and relatedness filters is exactly:
    excluded = (call rate < threshold) OR (|het z| > SD limit) OR (dropped for relatedness)
The "related" flag is read from the reasons column of step 09.
Grid: call-rate thresholds 0.95, 0.96, 0.97, 0.975, 0.98 (current) x heterozygosity limit 3 SD (current) and 4 SD
(the second only to show how much the heterozygosity limit matters: one sample has z < -3).
Outputs (--outdir): callrate_sensitivity.tsv (one row per scenario) and callrate_sensitivity_report.txt.
Standard library only.
"""
import argparse
import os

THRESHOLDS = [0.95, 0.96, 0.97, 0.975, 0.98]
HET_LIMITS = [3.0, 4.0]


def read_tsv(path):
    with open(path) as fh:
        h = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in fh]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-table", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    S = []
    for r in read_tsv(args.qc_table):
        S.append(dict(iid=r["IID"], call=float(r["call_rate"]), z=float(r["het_z"]), pedsex=r["pedsex"],
                      donor=r["is_donor_of_ipsc_line"], related="related_to" in r["reasons"], sexflag=("sex" in r["reasons"] or "female" in r["reasons"])))
    rep, rows = [], []
    for het in HET_LIMITS:
        for thr in THRESHOLDS:
            ex = []
            for s in S:
                why = []
                if s["call"] < thr:
                    why.append("call rate")
                if abs(s["z"]) > het:
                    why.append("het")
                if s["related"]:
                    why.append("kinship")
                if s["sexflag"]:
                    why.append("sex")
                if why:
                    ex.append((s["iid"], "+".join(why)))
            ex.sort(key=lambda x: int(x[0]))
            donors_ex = [f"{s['donor']}" for s in S if s["iid"] in {i for i, _ in ex} and s["donor"]]
            amb_ex = [s["iid"] for s in S if s["pedsex"] == "0" and s["iid"] in {i for i, _ in ex}]
            rows.append(dict(het_sd=het, thr=thr, n_ex=len(ex), n_keep=len(S) - len(ex),
                             ex="; ".join(f"{i}({w})" for i, w in ex), donors=", ".join(donors_ex) or "none",
                             amb=", ".join(amb_ex) or "none"))
    with open(os.path.join(args.outdir, "callrate_sensitivity.tsv"), "w") as out:
        out.write("het_limit_SD\tcall_rate_threshold\tn_excluded\tn_kept\texcluded_IIDs(reason)\tiPSC_donors_excluded\tsex0_samples_excluded\n")
        for r in rows:
            out.write(f"{r['het_sd']:g}\t{r['thr']}\t{r['n_ex']}\t{r['n_keep']}\t{r['ex']}\t{r['donors']}\t{r['amb']}\n")
    for r in rows:
        mark = "  <-- current (step 09)" if (r["het_sd"], r["thr"]) == (3.0, 0.98) else ("  <-- requested scenario" if (r["het_sd"], r["thr"]) == (3.0, 0.97) else "")
        rep.append(f"het limit {r['het_sd']:g} SD | call rate >= {r['thr']}: excluded {r['n_ex']}, kept {r['n_keep']} | "
                   f"excluded: {r['ex']} | iPSC donors excluded: {r['donors']} | sex=0 samples excluded: {r['amb']}{mark}")
    # call-rate-only counts, for reference
    rep.append("call rate alone (no other filter): " + ", ".join(f"< {t}: {sum(1 for s in S if s['call'] < t)}" for t in THRESHOLDS))
    with open(os.path.join(args.outdir, "callrate_sensitivity_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
