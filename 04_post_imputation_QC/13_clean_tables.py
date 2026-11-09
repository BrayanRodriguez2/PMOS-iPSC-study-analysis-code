#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/13_clean_tables.py   (original step script: 45_clean_tables.py)
# Produces: Variant table and counts of the clean dataset.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 45b (python part) - Variant table and filter-count table of the clean imputed dataset (v2_imputed_clean.pgen/pvar/psam).

variant_table.tsv.gz (one row per variant of the clean file, so that R2 >= 0.8 or ER2 filters can be applied later, e.g. for the PRS sensitivity):
  ID (chr:pos:ref:alt, as in the pvar), rsID (original ID of the TOPMed panel), CHROM, POS, REF, ALT, R2, ER2 (genotyped only), MAF_server (140 samples),
  TYPE (Typed/Imputed), AF129 and HWE_p129 (plink2 --freq/--hardy over the 129 QC samples, ALT allele = ALT of the pvar).
counts_clean.tsv: variants after each filter per chromosome (written by 45a) + totals, and the final composition.
"""
import argparse, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--clean-dir", required=True); ap.add_argument("--tables-dir", required=True); ap.add_argument("--unzipped-dir", required=True); ap.add_argument("--qc-only-dir", default=None)
a = ap.parse_args()
C = a.clean_dir
pv = pd.read_csv(os.path.join(C, "v2_imputed_clean.pvar"), sep="\t", comment="#", header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "INFO"], dtype={"CHROM": str, "ID": str})
print(f"final pvar: {len(pv):,} variants; unique IDs: {pv.ID.is_unique}", flush=True)
parts = []
for c in range(1, 23):
    d = pd.read_csv(os.path.join(a.tables_dir, f"chr{c}.tsv.gz"), sep="\t", usecols=["POS", "ID", "REF", "ALT", "MAF", "R2", "ER2", "TYPE"])
    d["CHROM"] = str(c); parts.append(d.rename(columns={"ID": "rsID", "MAF": "MAF_server"}))
info = pd.concat(parts, ignore_index=True)
info["ID"] = info.CHROM + ":" + info.POS.astype(str) + ":" + info.REF + ":" + info.ALT
t = pv[["CHROM", "POS", "ID", "REF", "ALT"]].merge(info[["ID", "rsID", "R2", "ER2", "MAF_server", "TYPE"]].drop_duplicates("ID"), on="ID", how="left")
af = pd.read_csv(os.path.join(C, "stats129.afreq"), sep="\t").rename(columns={"#CHROM": "CHROM"})[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AF129"})
hw = pd.read_csv(os.path.join(C, "stats129.hardy"), sep="\t").rename(columns={"#CHROM": "CHROM"})[["ID", "P"]].rename(columns={"P": "HWE_p129"})
t = t.merge(af, on="ID", how="left").merge(hw, on="ID", how="left")
t.to_csv(os.path.join(C, "variant_table.tsv.gz"), sep="\t", index=False, compression="gzip", float_format="%.6g")
miss = int(t.R2.isna().sum())
cnt = pd.concat([pd.read_csv(os.path.join(C, f"counts_chr{c}.tsv"), sep="\t", header=None, names=["chrom", "delivered_R2_ge_0.3", "snps_acgt_biallelic", "maf_ge_0.01_in_129", "hwe_p_ge_1e-6_in_129", "final_unique_ids"]) for c in range(1, 23)])
# "delivered" could not be counted inside 45a (no bcftools in that env): take it from the step-39 check tables (variants of the dose VCF)
deliv = {r[0]: int(r[3]) for r in (l.rstrip("\n").split("\t") for f in os.listdir(a.unzipped_dir) if f.startswith("check_chr") for l in open(os.path.join(a.unzipped_dir, f)))}
cnt["delivered_R2_ge_0.3"] = cnt.chrom.map(deliv)
cnt.loc[len(cnt)] = ["total"] + [int(cnt[c].sum()) for c in cnt.columns[1:]]
cnt.to_csv(os.path.join(C, "counts_clean.tsv"), sep="\t", index=False)
print(cnt.to_string(index=False))
L = [f"clean file: {len(t):,} variants x 135 samples (129 main / 135 sensitivity via --keep); variants without a match in the info tables: {miss}",
     f"composition: Typed {int((t.TYPE == 'Typed').sum()):,} | Imputed {int((t.TYPE == 'Imputed').sum()):,}",
     f"R2 >= 0.8: {int((t.R2 >= 0.8).sum()):,} ({100 * (t.R2 >= 0.8).mean():.1f} %) | R2 >= 0.9: {int((t.R2 >= 0.9).sum()):,} | median R2 {t.R2.median():.3f}",
     f"MAF (129) quantiles: {t.AF129.clip(upper=1).apply(lambda x: min(x, 1 - x)).quantile([0.01, 0.25, 0.5, 0.75]).round(4).to_dict()}",
     f"genotyped variants with ER2 < 0.5 still in the clean file: {int(((t.TYPE == 'Typed') & (t.ER2 < 0.5)).sum())} | ER2 < 0.3: {int(((t.TYPE == 'Typed') & (t.ER2 < 0.3)).sum())}",
     f"HWE p (129) min {t.HWE_p129.min():.2e}"]
if a.qc_only_dir:
    k75 = pd.read_csv(os.path.join(a.qc_only_dir, "keep_fin_variants_88.tsv"), sep="\t", dtype={"CHROM": str})
    c446 = pd.read_csv(os.path.join(a.qc_only_dir, "compare_1000g_446.tsv"), sep="\t", dtype={"CHROM": str})
    key = lambda d: d.CHROM.str.replace("chr", "") + ":" + d.POS.astype(str) + ":" + d.REF + ":" + d.ALT
    s75 = set(key(k75)); s358 = set(key(c446[~c446.abs_diff_topmed_gt_0_2.astype(bool)]))
    lo = t[(t.TYPE == "Typed") & (t.ER2 < 0.5)]
    L.append(f"in the clean file: of the 75 kept FIN-supported variants {int(t.ID.isin(s75).sum())}; of the 358 kept discordant {int(t.ID.isin(s358).sum())}")
    L.append(f"genotyped with ER2 < 0.5 in the clean file: {len(lo)} (in the 75: {int(lo.ID.isin(s75).sum())}; in the 358: {int(lo.ID.isin(s358).sum())}; others: {int((~lo.ID.isin(s75 | s358)).sum())})")
L.append(f"variants with an rsID in the TOPMed panel ID: {int(t.rsID.astype(str).str.startswith('rs').sum()):,} of {len(t):,}")
open(os.path.join(C, "clean_report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
