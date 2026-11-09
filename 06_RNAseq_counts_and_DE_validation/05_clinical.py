#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/05_clinical.py   (original step script: 65a_clinical.py)
# Produces: Donor clinical table and group comparisons from the corrected workbook.
# Run: python 06_RNAseq_counts_and_DE_validation/05_clinical.py <args>   (or from a SLURM script; see the usage line below)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 65a - donor clinical table for the expression validation, from the CORRECTED clinical workbook (Clinical_Data_250526_corrected 061026, sheet
'Matching_Final sent from Oulu', one row per iPSC line; TSV made by 00_common/01_xlsx_to_tsv.py), compared with the original analysis' Table S1 (Additional_file_1, sheet S1_donors).
  * line names: the sheet gives 'HEL362' etc.; the full line name (HEL362.1) is taken from rna_samples.tsv (the FASTQ HEL362_3 is the clone HEL362.1).
  * FAI is recomputed as 100 x testosterone (nmol/L) / SHBG (nmol/L) (the sheet gives testosterone in pmol/L) and HOMA-IR as glucose x insulin / 22.5; both are compared
    with the sheet values (the corrected FAI of HEL359 = 6.38 is the one in the sheet).
  * every numeric column is compared with Table S1 (differences are listed; none is expected after the correction).
  * Mann-Whitney U (two-sided, scipy 'auto') PMOS vs control for AMH, testosterone, FAI, FG score, BMI, HOMA-IR (expected in Fig 2: AMH < 0.001, testosterone 0.004, FAI 0.002, FG 0.002, BMI 0.16, HOMA-IR 0.21).
Outputs (--out): donor_meta.tsv (analysis table), clinical_vs_S1.tsv, clinical_mannwhitney.tsv, clinical_report.txt
"""
import argparse, os
import numpy as np, pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--clinical", required=True); ap.add_argument("--s1", required=True); ap.add_argument("--rna-samples", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rep = []
log = lambda m="": (print(m), rep.append(str(m)))

raw = pd.read_csv(a.clinical, sep="\t", header=0)
raw = raw[raw["Sample_ID"].notna() & raw["iPSC line"].astype(str).str.startswith("HEL")].iloc[:14].copy()
rs = pd.read_csv(a.rna_samples, sep="\t")
full = {r.ipsc_line.split(".")[0]: r.ipsc_line for r in rs.itertuples()}
raw["line"] = raw["iPSC line"].map(full)
raw["esf_id"] = raw["Sample_ID"].astype(float).astype(int)
assert raw.line.notna().all() and len(raw) == 14, raw[["iPSC line", "line"]]
chk = rs.set_index("ipsc_line").esf_id.reindex(raw.line).values
assert (chk == raw.esf_id.values).all(), "donor ID of the clinical sheet does not match rna_samples.tsv"
d = pd.DataFrame({"line": raw.line, "esf_id": raw.esf_id, "group": np.where(raw.Group.str.upper().isin(["CTRL", "CONTROL"]), "CTRL", "PMOS"),
                  "fg_score": raw.FGscore.astype(float), "age": raw.AGE.astype(float), "bmi": raw.BMI.astype(float), "lh": raw.LH.astype(float), "shbg": raw.SHBG.astype(float),
                  "amh": raw.AMH.astype(float), "testosterone": raw.Testo.astype(float), "dht": raw.DHT.astype(float), "fai": raw.FAI.astype(float), "hba1c": raw.HbA1c.astype(float),
                  "glucose": raw.Gluk_mM.astype(float), "insulin": raw.Insuliini_mIU_L.astype(float), "homa_ir": raw.HOMA_IR.astype(float)})
d["fai_recalc"] = 100 * (d.testosterone / 1000.0) / d.shbg
d["homa_recalc"] = d.glucose * d.insulin / 22.5
d = d.set_index("line")
log(f"clinical sheet: {len(d)} lines ({(d.group == 'PMOS').sum()} PMOS, {(d.group == 'CTRL').sum()} CTRL); BMI categories: lean < 25 = {(d.bmi < 25).sum()} ({((d.bmi < 25) & (d.group == 'CTRL')).sum()} CTRL / {((d.bmi < 25) & (d.group == 'PMOS')).sum()} PMOS), obese > 30 = {(d.bmi > 30).sum()}, between = {((d.bmi >= 25) & (d.bmi <= 30)).sum()}")
log(f"FAI recalculated (100 x T[nmol/L] / SHBG) vs the FAI column: max |diff| = {(d.fai - d.fai_recalc).abs().max():.3f} (lines: {list(d.index[(d.fai - d.fai_recalc).abs() > 0.02])}); "
    f"HOMA-IR recalculated vs column: max |diff| = {(d.homa_ir - d.homa_recalc).abs().max():.3f} (lines: {list(d.index[(d.homa_ir - d.homa_recalc).abs() > 0.02])})")
d.to_csv(os.path.join(a.out, "donor_meta.tsv"), sep="\t", float_format="%.6g")

s1 = pd.read_csv(a.s1, sep="\t"); s1 = s1[s1.iPSC_line_RNAseq.astype(str).str.startswith("HEL")].set_index("iPSC_line_RNAseq")
mp = {"fg_score": "FGscore", "age": "AGE", "bmi": "BMI", "lh": "LH", "shbg": "SHBG", "amh": "AMH", "testosterone": "Testo", "dht": "DHT", "fai": "FAI", "hba1c": "HbA1c", "glucose": "Gluk_mM",
      "insulin": "Insuliini_mIU_L", "homa_ir": "HOMA_IR"}
rows = []
for k, v in mp.items():
    for l in d.index:
        o, t = d.loc[l, k], float(s1.loc[l, v])
        if abs(o - t) > max(0.006 * abs(t), 1e-6):
            rows.append(dict(line=l, variable=k, corrected_clinical=o, S1=t, diff=o - t))
grp_s1 = s1.Group.map({"Control": "CTRL", "PMOS": "PMOS"}); gdiff = [l for l in d.index if grp_s1[l] != d.loc[l, "group"]]
dv = pd.DataFrame(rows, columns=["line", "variable", "corrected_clinical", "S1", "diff"]); dv.to_csv(os.path.join(a.out, "clinical_vs_S1.tsv"), sep="\t", index=False, float_format="%.5g")
log(f"differences with the original Table S1 (> 0.6 % relative): {len(dv)}; group labels differing: {gdiff}")
if len(dv): log(dv.to_string(index=False))
mw = []
for k, exp in [("amh", "<0.001"), ("testosterone", "0.004"), ("fai", "0.002"), ("fg_score", "0.002"), ("bmi", "0.16"), ("homa_ir", "0.21")]:
    x, y = d[d.group == "PMOS"][k], d[d.group == "CTRL"][k]
    mw.append(dict(variable=k, median_PMOS=x.median(), median_CTRL=y.median(), p_mannwhitney=stats.mannwhitneyu(x, y, alternative="two-sided").pvalue, expected_p=exp))
mw = pd.DataFrame(mw); mw.to_csv(os.path.join(a.out, "clinical_mannwhitney.tsv"), sep="\t", index=False, float_format="%.4g")
log("Mann-Whitney PMOS vs control (Fig 2A-F):\n" + mw.round(4).to_string(index=False))
open(os.path.join(a.out, "clinical_report.txt"), "w").write("\n".join(rep) + "\n")
