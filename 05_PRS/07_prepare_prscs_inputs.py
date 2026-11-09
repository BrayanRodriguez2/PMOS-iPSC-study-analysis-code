#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/07_prepare_prscs_inputs.py   (original step script: 14_prepare_sumstats.py)
# Produces: Summary statistics in PRS-CS format, n_gwas and the whole-panel .bim.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 14 - convert the harmonised GWAS Catalog file (GCST90570579, GWAS-SSF, GRCh38) into the
input format of PRS-CS:  header "SNP A1 A2 BETA SE" (whitespace separated).

Mapping:   SNP = rsid | A1 = effect_allele | A2 = other_allele | BETA = beta | SE = standard_error
           (in GWAS-SSF "beta" is the effect of the effect_allele; PRS-CS wants BETA for A1.)
Filters (each counted in the report; nothing is dropped silently):
  * rsID must look like rs<digits>                 (PRS-CS matches SNPs by rsID)
  * both alleles must be a single A/C/G/T          (indels cannot be handled by PRS-CS)
  * BETA and SE must be numeric and SE > 0
  * duplicated rsIDs: the FIRST occurrence is kept (file is sorted by position)
--n_gwas for PRS-CS = MEDIAN of the column `n` over ALL rows of the file (as requested), written to
n_gwas.txt. The median over the rows kept is also reported for comparison. Note that `n` varies per
variant (variants present in fewer studies have smaller n); the median is a summary, not the total
sample size of the study (n = 420,703 in the file metadata).
Outputs (--outdir): sumstats_prscs.txt, n_gwas.txt, sumstats_prep_report.txt.
Standard library only.
"""
import argparse
import gzip
import math
import os
import re
import statistics

RS_RE = re.compile(r"^rs\d+$")


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

    seen = set()
    n_all, n_kept_vals = [], []
    c = dict(rows=0, no_rsid=0, not_snp=0, bad_beta_se=0, duplicate=0, kept=0)
    with gzip.open(args.sumstats, "rt") as fh, open(os.path.join(args.outdir, "sumstats_prscs.txt"), "w") as out:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(h)}
        out.write("SNP\tA1\tA2\tBETA\tSE\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            c["rows"] += 1
            try:
                n_all.append(int(float(f[ix["n"]])))
                nval = n_all[-1]
            except ValueError:
                nval = None
            rs, a1, a2 = f[ix["rsid"]], f[ix["effect_allele"]].upper(), f[ix["other_allele"]].upper()
            if not RS_RE.match(rs):
                c["no_rsid"] += 1
                continue
            if len(a1) != 1 or len(a2) != 1 or a1 not in "ACGT" or a2 not in "ACGT":
                c["not_snp"] += 1
                continue
            try:
                beta, se = float(f[ix["beta"]]), float(f[ix["standard_error"]])
            except ValueError:
                c["bad_beta_se"] += 1
                continue
            if not (math.isfinite(beta) and math.isfinite(se) and se > 0):
                c["bad_beta_se"] += 1
                continue
            if rs in seen:
                c["duplicate"] += 1
                continue
            seen.add(rs)
            out.write(f"{rs}\t{a1}\t{a2}\t{f[ix['beta']]}\t{f[ix['standard_error']]}\n")
            c["kept"] += 1
            if nval is not None:
                n_kept_vals.append(nval)
    med_all = int(round(statistics.median(n_all)))
    med_kept = int(round(statistics.median(n_kept_vals)))
    with open(os.path.join(args.outdir, "n_gwas.txt"), "w") as out:
        out.write(f"{med_all}\n")
    log(f"[rows] {c['rows']:,} | dropped: no valid rsID {c['no_rsid']:,}, not a single-nucleotide SNP {c['not_snp']:,}, "
        f"bad BETA/SE {c['bad_beta_se']:,}, duplicated rsID {c['duplicate']:,} | KEPT {c['kept']:,}")
    log(f"[n] median of `n` over all rows: {med_all:,} (used as --n_gwas) | over kept rows: {med_kept:,} | "
        f"min {min(n_all):,} | max {max(n_all):,}")
    with open(os.path.join(args.outdir, "sumstats_prep_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
