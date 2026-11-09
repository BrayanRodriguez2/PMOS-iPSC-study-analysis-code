#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/07_update_S2_S3.py   (original step script: 83_update_S2_S3.py)
# Produces: Updated supplementary tables S2 (sample QC) and S3 (lead variants).
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 83 - updated versions of Tables S2 (genotype QC per sample, 140 women) and S3 (29 lead variants) of the original analysis' final package, filled with OUR results, in the original column format where possible
(the original columns first, same names; extra columns of ours after; the original original values for comparison in columns prefixed 'her_').
S2 sources: sample_qc_table.tsv (step 09: pedigree and genetic sex, call rate, heterozygosity F, X F, decisions), sample_intensity_metrics.tsv (step 24: contamination estimate alpha_hat),
  king.kin0 (step 09: kinship of all pairs), x_het_*.scount (step 83b: X heterozygosity), pca129.eigenvec (step 80: PC1-PC5 of the 129 QC samples; blank for the 11 excluded, which are not in that PCA),
  donor_key.tsv (iPSC line). Contamination status: alpha_hat >= 0.09 'contaminated' (excluded; the five contaminated samples: ~10-12 %), 0.035-0.09 'possible low-level' (kept; 11 samples: a set of samples, alpha 3.7-7.4 %; the alpha estimate has not been validated yet), else none.
S3 sources: grs29_orientation.tsv (step 51b: risk allele, risk-allele frequency in the 129 QC samples, TOPMed status, R2, ER2, allele orientation), leads_imputed_coverage.tsv (step 44: MAF of the server),
  raw_129.vmiss (step 10c: call rate of the on-array leads in the 129) and lead_qc_fate.txt (leads removed by the SNP QC); GRCh37 position, log OR, SE and P are COPIED from the original table S3 (ST2 of the paper; not recomputed).
Outputs (--out): S2_genotype_QC_updated.tsv, S3_lead_variants_updated.tsv, S2_S3_updated.xlsx (sheets S2_updated, S3_updated, README_columns), update_S2_S3_report.txt
Usage: python 08_preprint_tasks/07_update_S2_S3.py <PRS dir> <the original S2 tsv> <the original S3 tsv> <out_dir>
"""
import os, sys
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

P, s2f, s3f, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True); V = f"{P}/validation/results"
rep = []
log = lambda m="": (print(m), rep.append(str(m)))
her2 = pd.read_csv(s2f, sep="\t").set_index("Genotype_array_ID"); her2.index = her2.index.astype(int)
qc = pd.read_csv(f"{V}/qc/samples/sample_qc_table.tsv", sep="\t").set_index("IID")
cm = pd.read_csv(f"{V}/contamination/sample_intensity_metrics.tsv", sep="\t"); cm = cm[cm.IID.astype(str).str.fullmatch(r"\d+")].copy(); cm["IID"] = cm.IID.astype(int)   # positions without an assigned sample are 'NA_unassigned'
cm = cm.set_index("IID")   # index by sample ID (v1.0.1 fix: this statement was inside the comment above, so cm.loc[i] used the row position)
kin = pd.read_csv(f"{V}/qc/samples/king.kin0", sep="\t"); kin = kin.rename(columns={"#FID1": "FID1"})
xh = pd.read_csv(f"{V}/array_pca/x_het_nonPAR.scount", sep="\t").set_index("#IID")
pca = pd.read_csv(f"{V}/array_pca/pca129.eigenvec", sep="\t").set_index("IID")
key = pd.read_csv(f"{P}/metadata/donor_key.tsv", sep="\t").set_index("esf_id")
ids = sorted(qc.index)
assert len(ids) == 140 and set(ids) == set(her2.index), "sample sets differ"
# max kinship per sample over all pairs
mk = {}
for r in kin.itertuples():
    for a, b in ((r.IID1, r.IID2), (r.IID2, r.IID1)):
        if r.KINSHIP > mk.get(a, (-9, None))[0]: mk[a] = (r.KINSHIP, b)
rows = []
for i in ids:
    q = qc.loc[i]; a = cm.loc[i, "alpha_hat"] if i in cm.index else np.nan; az = cm.loc[i, "alpha_robust_z"] if i in cm.index else np.nan
    cstat = "contaminated (alpha ~%d %%)" % round(100 * a) if a >= 0.09 else ("possible low-level (alpha ~%.1f %%; kept, estimate not validated)" % (100 * a) if a >= 0.035 else "no evidence")
    reasons = [] if pd.isna(q.reasons) else str(q.reasons).split(";")
    txt = []
    for r_ in reasons:
        txt.append({"low_call_rate": "call rate < 0.98 (%.4f)" % q.call_rate, "heterozygosity_outlier_high": "excess homozygosity: heterozygosity F %.3f (z %.1f)" % (q.het_F, q.het_z),
                    "heterozygosity_outlier_low": "excess heterozygosity: F %.3f (z %.1f)" % (q.het_F, q.het_z)}.get(r_, "first-degree relative of another sample (see max_kinship_partner); the lower-quality sample removed" if r_.startswith("related") else r_))
    if a >= 0.09: txt.append("DNA contamination (alpha_hat %.3f), not a sample swap" % a)
    kv, kp = mk.get(i, (np.nan, None)); pcr = pca.loc[i, [f"PC{k}" for k in range(1, 6)]].values if i in pca.index else [np.nan] * 5
    xhet = xh.loc[i, "HET_CT"] / (xh.loc[i, "HET_CT"] + xh.loc[i, "HOM_CT"])
    rows.append(dict(Genotype_array_ID=i, iPSC_line=key.ipsc_line.get(i, ""), call_rate=q.call_rate, het_F=q.het_F, x_het=xhet, PC1=pcr[0], PC2=pcr[1], PC3=pcr[2], PC4=pcr[3], PC5=pcr[4],
                     flag="; ".join(txt + ([] if a < 0.035 or a >= 0.09 else ["contamination estimate alpha_hat %.3f: %s" % (a, cstat)])), standard_QC=("fail" if q.excluded == "yes" else "pass"), QC_reason="; ".join(txt),
                     sex_pedigree={2: "female", 1: "male", 0: "missing"}.get(int(q.pedsex), "?"), sex_genetic={2: "female", 1: "male", 0: "?"}.get(int(q.snpsex), "?"), het_F_z=q.het_z, chrX_F=q.xF,
                     max_kinship=kv, max_kinship_partner=kp, relatedness=("first-degree (%.3f) with %s" % (kv, kp) if kv >= 0.177 else ("second-degree-range (%.3f) with %s, reported only" % (kv, kp) if kv >= 0.0884 else "none >= 0.0884")),
                     contamination_alpha_hat=a, contamination_robust_z=az, contamination_status=cstat, final_decision=("excluded" if q.excluded == "yes" else "in the 129 QC set"),
                     in_135_sensitivity_set=("yes" if i not in {int(x) for x in __import__("os").environ.get("EXCLUDED_FROM_135_IDS", "").split(",") if x} else "no"),
                     her_call_rate=her2.loc[i, "call_rate"], her_het_F=her2.loc[i, "het_F"], her_x_het=her2.loc[i, "x_het"], her_flag=her2.loc[i, "flag"], her_standard_QC=her2.loc[i, "standard_QC_2Nov2026"], her_QC_reason=her2.loc[i, "QC_reason"]))
S2 = pd.DataFrame(rows)
S2.to_csv(f"{out}/S2_genotype_QC_updated.tsv", sep="\t", index=False, float_format="%.5g")
ex = S2[S2.standard_QC == "fail"]
log(f"S2: {len(S2)} samples; {int((S2.standard_QC == 'pass').sum())} pass, {len(ex)} excluded: " + ", ".join(f"{r.Genotype_array_ID}" for r in ex.itertuples()))
log(f"   contamination status: " + str(S2.contamination_status.str.split(" \\(").str[0].value_counts().to_dict()))
hp = her2.standard_QC_2Nov2026.astype(str).str.lower() if "standard_QC_2Nov2026" in her2 else None
her_fail = set(her2.index[her2.standard_QC_2Nov2026.astype(str).str.lower() != "pass"]); our_fail = set(ex.Genotype_array_ID)
log(f"   the original 'standard_QC_2Nov2026': {len(her_fail)} not 'pass' ({sorted(her_fail)}); ours {len(our_fail)} ({sorted(our_fail)}); only the original: {sorted(her_fail - our_fail)}; only ours: {sorted(our_fail - her_fail)}")
log(f"   call rate ours vs the original: median difference {(S2.call_rate - S2.her_call_rate).median():+.4f}, max |diff| {(S2.call_rate - S2.her_call_rate).abs().max():.4f}; het_F: max |diff| {(S2.het_F - S2.her_het_F).abs().max():.4f}; x_het ours {S2.x_het.median():.3f} (median) vs the original {S2.her_x_het.median():.3f}")
# ---------------------------------------------------------------- S3
her3 = pd.read_csv(s3f, sep="\t"); her3["rsid"] = her3.rsid.str.strip()
ori = pd.read_csv(f"{P}/prscs/weights_imputed/grs29_orientation.tsv", sep="\t").set_index("rsid")
cov = pd.read_csv(f"{V}/imputation_qc/leads_imputed_coverage.tsv", sep="\t").set_index("lead")
vm = pd.read_csv(f"{V}/qc/variants/raw_129.vmiss", sep="\t").set_index("ID")
rows = []
for r in her3.itertuples():
    o = ori.loc[r.rsid]; c = cov.loc[r.rsid]
    cr = 1 - vm.loc[r.rsid, "F_MISS"] if r.rsid in vm.index else np.nan
    on_raw = int(r.on_chip) == 1; typed = o.type == "Typed"
    note = []
    if on_raw and not typed: note.append("on the raw array but removed by the SNP QC (" + {"rs148982377": "monomorphic call: assay failure", "rs11031005": "call rate 75 % in the 129"}.get(r.rsid, "see lead_qc_fate") + "); recovered by imputation")
    if bool(o.palindromic): note.append("A/T or C/G: orientation checked by frequency")
    rows.append(dict(rsid=r.rsid, GRCh37_position=r.GRCh37_position, GRCh38_position=int(o.pos_grch38), nearest=r.nearest, consensus=r.consensus, risk_allele=o.risk_allele, other_allele=o.other_allele, RAF_GWAS=r.RAF_GWAS,
                     log_OR=r.log_OR, SE=r.SE, P_GWAS=r.P_GWAS, on_chip=int(on_raw), strand=r.strand if on_raw else "", RAF_Oulu_129_QC=o.cohort_risk_freq, MAF_Oulu_129_QC=min(o.cohort_risk_freq, 1 - o.cohort_risk_freq),
                     call_rate_129_QC=(cr if on_raw else np.nan), usable=1, note="; ".join(note), after_imputation=("genotyped (TOPMed Typed), R2 %.3f, ER2 %.3f" % (o.R2, o.ER2) if typed else "imputed (TOPMed r3), R2 %.3f" % o.R2),
                     status_now=o.type, TOPMed_variant_ID=o.clean_id, R2=o.R2, ER2=(o.ER2 if typed else np.nan), MAF_TOPMed_140=c.MAF, RAF_diff_vs_GWAS=o.freq_diff_cohort_vs_ST2, allele_match_to_clean_file=o.match,
                     sumstats_sign_agrees_with_risk_allele=(o.sumstats_sign_vs_ST2_risk == "+"), her_on_chip=r.on_chip, her_RAF_Oulu_cohort_n140=r.RAF_Oulu_cohort_n140, her_call_rate=r.call_rate, her_usable=r.usable, her_note=r.note))
S3 = pd.DataFrame(rows); S3.to_csv(f"{out}/S3_lead_variants_updated.tsv", sep="\t", index=False, float_format="%.5g")
log(f"S3: {len(S3)} leads; genotyped (Typed) {int((S3.status_now == 'Typed').sum())}, imputed {int((S3.status_now == 'Imputed').sum())}; R2 min {S3.R2.min():.3f}; risk-allele frequency difference vs GWAS: max |diff| {S3.RAF_diff_vs_GWAS.abs().max():.3f}, correlation {np.corrcoef(S3.RAF_Oulu_129_QC, S3.RAF_GWAS)[0, 1]:.3f}; sign of the deposited beta agrees in {int(S3.sumstats_sign_agrees_with_risk_allele.sum())}/29")
log("   on the raw array: %d (the original S3: %d); usable after our QC: typed on the array %d, the other %d recovered by imputation" % (int(S3.on_chip.sum()), int((her3.on_chip == 1).sum()), int((S3[S3.on_chip == 1].status_now == "Typed").sum()), 29 - int((S3[S3.on_chip == 1].status_now == "Typed").sum())))
# ---------------------------------------------------------------- xlsx
readme = [("sheet", "column", "content / source", "status"),
          ("S2_updated", "Genotype_array_ID ... QC_reason", "The original columns, same names and order. call_rate / het_F: ours (step 09, 129-sample QC; call rate over autosomal SNPs of the raw array). x_het: heterozygous fraction of X calls (non-PAR) of the raw array; the original definition is not documented and the original values are higher (~0.33-0.36 vs ~0.23), so the two are NOT comparable. PC1-PC5: PCA of the 129 QC samples (step 80); blank for the 11 excluded. flag, standard_QC, QC_reason: ours.", "filled (x_het with a different, documented definition)"),
          ("S2_updated", "sex_pedigree, sex_genetic, het_F_z, chrX_F, max_kinship, max_kinship_partner, relatedness", "Sex, heterozygosity z-score, X inbreeding F and KING kinship (max over all pairs; >= 0.177 first degree, 0.0884-0.177 second-degree range, reported only).", "new columns"),
          ("S2_updated", "contamination_alpha_hat, contamination_robust_z, contamination_status", "Contamination estimate from the BAF of the raw intensities (steps 23-24; alpha ~ fraction of foreign DNA): >= 0.09 contaminated (the five contaminated samples), 0.035-0.09 possible low-level (11 samples: a set of samples; alpha 3.7-7.4 %; all kept; the estimate is not validated yet, e.g. by in-silico mixtures).", "new columns"),
          ("S2_updated", "final_decision, in_135_sensitivity_set", "Final decision (129 QC set) and membership of the 135-sample sensitivity set (call rate >= 0.97; excludes only a few samples).", "new columns"),
          ("S2_updated", "her_* columns", "The original values of Table S2 for comparison.", "copied from the original table"),
          ("S3_updated", "rsid ... after_imputation", "The original columns in the original order. GRCh37_position, nearest, consensus, RAF_GWAS, log_OR, SE, P_GWAS: COPIED from the original table (ST2 of Moolhuijsen et al.), not recomputed. risk/other allele: ST2 (checked against the deposited summary statistics: sign agrees in 29/29). on_chip: raw array (7 leads). strand: copied from the original table (array alignment strand; not re-derived). usable: 1 for all 29 (typed or imputed).", "filled / copied"),
          ("S3_updated", "RAF_Oulu_129_QC, MAF_Oulu_129_QC, call_rate_129_QC", "Risk-allele frequency and MAF in the 129 QC samples (from the clean imputed file, dosage-based); call rate of the on-array leads in the 129 (raw array).", "replaces the original RAF_Oulu_cohort_n140 and call_rate (n = 140), which cannot be reproduced because the clean file has 135 samples (129 for frequencies)"),
          ("S3_updated", "status_now, TOPMed_variant_ID, R2, ER2, MAF_TOPMed_140, RAF_diff_vs_GWAS, allele_match_to_clean_file, sumstats_sign_agrees_with_risk_allele", "Genotyped/imputed status in the TOPMed r3 output, server R2 (ER2 = leave-one-out for genotyped variants), server MAF (140 samples), difference of the risk-allele frequency in the 129 vs the GWAS (ST2), orientation of the alleles in the clean file (direct = REF/ALT as in ST2).", "new columns"),
          ("S3_updated", "her_* columns", "The original on_chip, RAF_Oulu_cohort_n140, call_rate, usable, note.", "copied from the original table")]
wb = Workbook(); ws = wb.active; ws.title = "S2_updated"
for name, df in [("S2_updated", S2), ("S3_updated", S3)]:
    w = ws if name == "S2_updated" else wb.create_sheet(name)
    w.append(list(df.columns))
    for r in df.itertuples(index=False): w.append([None if (isinstance(x, float) and np.isnan(x)) else (x.item() if hasattr(x, "item") else x) for x in r])
    for c in w[1]: c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="DDEBF7"); c.alignment = Alignment(wrap_text=True, vertical="top")
    w.freeze_panes = "C2"
    for j, col in enumerate(df.columns, 1): w.column_dimensions[get_column_letter(j)].width = min(46, max(10, len(str(col)) + 2, int(df[col].astype(str).str.len().quantile(.9)) + 2))
w = wb.create_sheet("README_columns")
for r in readme: w.append(list(r))
for c in w[1]: c.font = Font(bold=True)
for col, wd in zip("ABCD", (14, 44, 120, 60)): w.column_dimensions[col].width = wd
wb.save(f"{out}/S2_S3_updated.xlsx")
open(f"{out}/update_S2_S3_report.txt", "w").write("\n".join(rep) + "\n")
