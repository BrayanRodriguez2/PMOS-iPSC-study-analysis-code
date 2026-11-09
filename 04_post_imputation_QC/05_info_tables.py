#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/05_info_tables.py   (original step script: 40_info_tables.py)
# Produces: Per-chromosome tables from the server info files.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 40 - Per-chromosome table of the TOPMed imputation "info" file (post-imputation QC, step 1 of 6).

The info.gz of this server version (Minimac 4.1.6, TOPMed Imputation Server 2) is NOT a column table: it is a sites-only VCF whose INFO
field carries  AF (estimated alt allele frequency), MAF, AVG_CS (average call score), R2 (estimated imputation accuracy), ER2 (empirical
leave-one-out R2, only for genotyped variants) and the flags IMPUTED / TYPED. This script turns it into a plain table:
  CHROM POS ID REF ALT AF MAF R2 ER2 AVG_CS TYPE      with TYPE in {Imputed, Typed, Typed_only}
(Typed = genotyped on the array AND present in the panel; Typed_only is only reported if a TYPED_ONLY flag appears; the server's list of
"typed only" sites is in statistics/snps-typed-only.txt.) AF/MAF/R2 are computed by the server over the 140 uploaded samples.
Input : unzipped/chrN.info.gz (step 39)       Output: info_tables/chrN.tsv.gz  (+ one line of totals printed)
"""
import argparse, gzip, io, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--info", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
with gzip.open(a.info, "rt") as fh:
    lines = [l for l in fh if not l.startswith("##")]
df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
df.columns = [c.lstrip("#") for c in df.columns]
for key in ["AF", "MAF", "AVG_CS", "R2", "ER2"]:
    df[key] = pd.to_numeric(df["INFO"].str.extract(rf"(?:^|;){key}=([^;]+)")[0], errors="coerce")
typed = df["INFO"].str.contains(r"(?:^|;)TYPED(?:;|$)", regex=True)
typed_only = df["INFO"].str.contains("TYPED_ONLY", regex=False)
imputed = df["INFO"].str.contains(r"(?:^|;)IMPUTED(?:;|$)", regex=True)
df["TYPE"] = "Imputed"
df.loc[typed, "TYPE"] = "Typed"
df.loc[typed_only, "TYPE"] = "Typed_only"
combo = pd.Series(list(zip(typed, imputed, typed_only))).value_counts().to_dict()    # (TYPED, IMPUTED, TYPED_ONLY) -> count
# Minimac4 marks genotyped variants that are also in the panel with TYPED (and, in this server version, some also with IMPUTED:
# TYPED wins = the variant was genotyped on the array; "Imputed" = IMPUTED flag without TYPED)
assert (typed | typed_only | imputed).all(), "a variant has none of the flags"
df["POS"] = df["POS"].astype(int)
os.makedirs(os.path.dirname(a.out), exist_ok=True)
df[["CHROM", "POS", "ID", "REF", "ALT", "AF", "MAF", "R2", "ER2", "AVG_CS", "TYPE"]].to_csv(a.out, sep="\t", index=False, compression="gzip", float_format="%.6g")
print(f"flag combinations (TYPED, IMPUTED, TYPED_ONLY): {combo}")
print(f"{os.path.basename(a.info)}: {len(df):,} variants | {df.TYPE.value_counts().to_dict()} | R2 min {df.R2.min():.3f} (Imputed {df[df.TYPE=='Imputed'].R2.min():.3f}) | "
      f"ER2 available for {int(df.ER2.notna().sum()):,} | duplicated positions {int(df.duplicated(['POS']).sum()):,}")
