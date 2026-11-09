#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/04_lead_snps_lookup.py   (original step script: 11_lead_snps_lookup.py)
# Produces: The 29 lead SNPs in the summary statistics and genome-wide significant loci.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 11 - look up the 29 lead SNPs of the PCOS GWAS in the deposited summary statistics
(GCST90570579, harmonised, GRCh38) and cross-check them against the loci that are
genome-wide significant in THAT file.

Source of the 29 lead SNPs
  Moolhuijsen LME ... Karaderi T, Nat Genet 2026;58:1040-1050 (doi 10.1038/s41588-026-02543-9).
  The 29 rsIDs below are the column labels of the heat-map of Fig. 2 ("29 PCOS loci"), extracted
  from the PDF text with pdftotext (step 08). The Supplementary Tables (2 and 4), which would give
  the official lead-SNP list, are NOT in the PDF, so the list was NOT checked against them.
  Sanity checks done here: exactly 29 rsIDs, every one present in the sumstats, and the text of the
  paper itself assigns rs732310 to AMH (consistent with it being a lead SNP).

IMPORTANT caveat (from the paper's Methods): the 29 loci come from the "three-strata" meta-analysis
(544,513 women, WITH 23andMe) whereas the deposited GCST90570579/80 are the "two-strata" model
WITHOUT 23andMe (456,570; the file used here is the European, age-adjusted version, n=420,703).
Some leads may therefore be weaker (P > 5e-8) in the deposited data.

Cross-check by distance: variants with P < 5e-8 in the deposited file are grouped into loci by
greedy distance pruning (best P first, +/- 500 kb) -- NO LD is used, so close independent signals
merge into one locus. This is only meant to show whether the deposited data recover ~29 loci.

Outputs (--outdir): lead_snps_lookup.tsv, sumstats_significant_loci.tsv, lead_snps_report.txt
Standard library only.
"""
import argparse
import gzip
import os

LEAD_RSIDS = [
    "rs1641518", "rs8047587", "rs76658602", "rs3945628", "rs13389190", "rs11075468", "rs6712151",
    "rs17880096", "rs2070729", "rs853854", "rs10876864", "rs12488257", "rs12149203", "rs804263",
    "rs74690816", "rs12526762", "rs2183958", "rs13004711", "rs9312937", "rs670441", "rs1784692",
    "rs117568227", "rs1824982", "rs4944653", "rs7925543", "rs148982377", "rs732310", "rs11031005",
    "rs7589428",
]
P_GWS = 5e-8
WINDOW = 500_000


def chr_key(c):
    return (0, int(c)) if c.isdigit() else (1, {"X": 23, "Y": 24, "MT": 25}.get(c, 99))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sumstats", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rep = []

    def log(m=""):
        print(m, flush=True)
        rep.append(m)

    assert len(set(LEAD_RSIDS)) == 29, "expected exactly 29 distinct lead rsIDs"
    lead = set(LEAD_RSIDS)
    found, sig = {}, []
    with gzip.open(args.sumstats, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                p = float(f[ix["p_value"]])
            except ValueError:
                continue
            rs = f[ix["rsid"]]
            rec = (f[ix["chromosome"]], int(f[ix["base_pair_location"]]), f[ix["effect_allele"]], f[ix["other_allele"]],
                   f[ix["beta"]], f[ix["standard_error"]], f[ix["effect_allele_frequency"]], p, f[ix["n"]], rs)
            if rs in lead:
                found[rs] = rec
            if p < P_GWS:
                sig.append(rec)
    log(f"[lead SNPs] found in sumstats: {len(found)} of {len(lead)} | missing: {sorted(lead - set(found))}")
    log(f"[sumstats] variants with P < {P_GWS:g}: {len(sig):,}")

    # ---- distance-pruned loci from the deposited data ----------------------------------
    sig.sort(key=lambda r: r[7])
    loci = []  # (chr, pos, p, rsid, n_sig_variants)
    assigned = set()
    for i, r in enumerate(sig):
        if i in assigned:
            continue
        members = [j for j, q in enumerate(sig) if j not in assigned and q[0] == r[0] and abs(q[1] - r[1]) <= WINDOW]
        assigned.update(members)
        loci.append((r[0], r[1], r[7], r[9], len(members)))
    loci.sort(key=lambda x: (chr_key(x[0]), x[1]))
    log(f"[loci] distance-pruned (+/-{WINDOW // 1000} kb) loci with P < {P_GWS:g} in the deposited file: {len(loci)}")

    def nearest_locus(chrom, pos):
        c = [(abs(l[1] - pos), l) for l in loci if l[0] == chrom]
        return min(c, key=lambda x: x[0]) if c else (None, None)

    with open(os.path.join(args.outdir, "lead_snps_lookup.tsv"), "w") as out:
        out.write("rsid\tchr\tpos_grch38\teffect_allele\tother_allele\tbeta\tse\teaf\tp_value\tn\tgws_in_deposited\t"
                  "nearest_deposited_locus_lead\tdistance_bp\n")
        n_gws = n_near = 0
        for rs in LEAD_RSIDS:
            r = found.get(rs)
            if r is None:
                out.write(f"{rs}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tno\tNA\tNA\n")
                continue
            d, loc = nearest_locus(r[0], r[1])
            gws = r[7] < P_GWS
            n_gws += gws
            n_near += (d is not None and d <= WINDOW)
            out.write(f"{rs}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]:.3g}\t{r[8]}\t"
                      f"{'yes' if gws else 'no'}\t{loc[3] if loc else 'NA'}\t{d if d is not None else 'NA'}\n")
    log(f"[cross-check] paper leads with P < {P_GWS:g} in the deposited file: {n_gws} of {len(found)} | "
        f"leads within {WINDOW // 1000} kb of a deposited-data locus: {n_near} of {len(found)}")
    lead_pos = [(r[0], r[1]) for r in found.values()]
    extra = [l for l in loci if not any(c == l[0] and abs(p - l[1]) <= WINDOW for c, p in lead_pos)]
    log(f"[cross-check] deposited-data loci with no paper lead within {WINDOW // 1000} kb: {len(extra)} "
        f"{[(l[0], l[1], l[3]) for l in extra]}")
    with open(os.path.join(args.outdir, "sumstats_significant_loci.tsv"), "w") as out:
        out.write("chr\tpos_grch38\tp_value\trsid\tn_significant_variants_in_window\n")
        for l in loci:
            out.write(f"{l[0]}\t{l[1]}\t{l[2]:.3g}\t{l[3]}\t{l[4]}\n")
    with open(os.path.join(args.outdir, "lead_snps_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
