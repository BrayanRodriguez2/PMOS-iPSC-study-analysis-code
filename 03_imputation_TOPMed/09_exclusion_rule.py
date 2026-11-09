#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/09_exclusion_rule.py   (original step script: 35_exclusion_rule.py)
# Produces: Rule that removes problematic genotyped SNPs before the final upload.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 35 - Exclusion list for the v2 VCFs (decision taken by the analyst after the TOPMed QC-only review) and safety check against the
29 lead SNPs of the PCOS GWAS (ST2) and their LD proxies. Nothing is removed here: this script only WRITES the list; step 36 builds the VCFs.

Rule (all 140 samples are kept; the sample QC (129) is applied after imputation):
  R1  remove the 182 SNPs excluded by the server as "Allele mismatch" (alt allele of the array is not the panel's);
  R2  of the 88 SNPs whose |ALT freq here - ALT freq in TOPMed| > 0.2: REMOVE those that also differ by > 0.2 from 1000G FIN
      (99 Finns) or have no usable 1000G record (no record / different allele pair: cannot be checked, treated as suspect);
      KEEP those within 0.2 of 1000G FIN (plausible true Finnish variants) and list them separately;
  R3  keep the other 358 mismatched-frequency SNPs.
Lead check (before anything is removed): for every SNP to be removed -> (a) is it one of the 29 leads (position or rsID)?
(b) is it within +-500 kb of a lead with r2 >= 0.8 to that lead in 1000G EUR (503 samples; same definition as step 12)? (c) is it the
best proxy / best usable proxy recorded in locus_coverage.tsv (step 12)? Candidates with no 1000G record cannot be a step-12 proxy
(step 12 computed proxies from 1000G genotypes) and have r2 = NA; they are still listed when they lie within 500 kb of a lead.
Inputs : results/topmed_qc_only/{freq_mismatch_446_verified.tsv, allele_mismatch_182_verified.tsv, compare_1000g_446.tsv, g1000_gt.tsv,
         g1000_samples.txt}, results/loci/{lead_snps_lookup.tsv, locus_coverage.tsv}
Outputs: results/topmed_qc_only/{exclude_v2.tsv, keep_fin_variants_88.tsv, leads_check.tsv, exclusion_rule_report.txt}
"""
import argparse, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--outdir", required=True); ap.add_argument("--leads", required=True); ap.add_argument("--locus-cov", required=True)
a = ap.parse_args()
O = a.outdir
GT2D = {"0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2, "0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}
c = pd.read_csv(os.path.join(O, "compare_1000g_446.tsv"), sep="\t", dtype={"CHROM": str})
e = pd.read_csv(os.path.join(O, "allele_mismatch_182_verified.tsv"), sep="\t", dtype={"CHROM": str})
big = c.abs_diff_topmed_gt_0_2.astype(bool)
fin_close = (c.abs_diff_fin <= 0.2)                       # NaN -> False (no FIN frequency)
rm88 = c[big & ~fin_close].copy(); keep88 = c[big & fin_close].copy()
rm88["reason"] = np.where(rm88.fin_af.isna(), "freq_vs_TOPMed_gt_0.2_and_no_1000G_record_for_FIN", "freq_vs_TOPMed_gt_0.2_and_FIN_diff_gt_0.2")
rm182 = pd.DataFrame({"ID": e.ID, "CHROM": e.CHROM, "POS": e.POS, "REF": e.REF, "ALT": e.ALT, "reason": "TOPMed_allele_mismatch"})
ex = pd.concat([rm182, rm88[["ID", "CHROM", "POS", "REF", "ALT", "reason"]]], ignore_index=True)
assert ex.ID.is_unique, "an ID is in both lists"
ex.to_csv(os.path.join(O, "exclude_v2.tsv"), sep="\t", index=False)
keep88.to_csv(os.path.join(O, "keep_fin_variants_88.tsv"), sep="\t", index=False, float_format="%.5g")

# ---------------------------------------------------------------- lead / proxy check
leads = pd.read_csv(a.leads, sep="\t", dtype={"chr": str}).dropna(subset=["pos_grch38"])
leads["pos_grch38"] = leads.pos_grch38.astype(int)
cov = pd.read_csv(a.locus_cov, sep="\t", dtype={"chr": str})
samples = [l.strip() for l in open(os.path.join(O, "g1000_samples.txt")) if l.strip()]
g = pd.read_csv(os.path.join(O, "g1000_gt.tsv"), sep="\t", header=None, dtype=str); g.columns = ["CHROM", "POS", "REF", "ALT"] + samples
g["CHROM"] = g.CHROM.str.replace("chr", "", regex=False); g["POS"] = g.POS.astype(int)
G = np.vstack([g[s].map(GT2D).fillna(np.nan).to_numpy(dtype=float) for s in samples]).T
rec = {(r.CHROM, r.POS): i for i, r in enumerate(g.itertuples())}
cand = ex.assign(chrn=ex.CHROM.str.replace("chr", "", regex=False))
rows = []
for L in leads.itertuples():
    lead_i = rec.get((str(L.chr), L.pos_grch38))
    near = cand[(cand.chrn == str(L.chr)) & ((cand.POS - L.pos_grch38).abs() <= 500_000)]
    for r in near.itertuples():
        r2 = np.nan
        i = rec.get((r.chrn, r.POS))
        if lead_i is not None and i is not None:
            x, y = G[lead_i], G[i]; ok = ~np.isnan(x) & ~np.isnan(y)
            if ok.sum() > 10 and x[ok].std() > 0 and y[ok].std() > 0: r2 = float(np.corrcoef(x[ok], y[ok])[0, 1] ** 2)
        rows.append(dict(lead=L.rsid, lead_chr=L.chr, lead_pos=L.pos_grch38, candidate_id=r.ID, candidate_pos=r.POS, reason_removed=r.reason,
                         distance_bp=int(r.POS - L.pos_grch38), r2_with_lead_1000G_EUR=r2,
                         is_the_lead=(r.POS == L.pos_grch38 or r.ID == L.rsid), lead_in_1000G=lead_i is not None, candidate_in_1000G=i is not None))
lc = pd.DataFrame(rows)
cols = ["lead", "lead_chr", "lead_pos", "candidate_id", "candidate_pos", "reason_removed", "distance_bp", "r2_with_lead_1000G_EUR", "is_the_lead", "lead_in_1000G", "candidate_in_1000G"]
lc = lc[cols] if len(lc) else pd.DataFrame(columns=cols)
lc.to_csv(os.path.join(O, "leads_check.tsv"), sep="\t", index=False, float_format="%.4g")
direct_lead = ex[ex.ID.isin(leads.rsid) | ex.apply(lambda r: ((leads.chr == r.CHROM.replace("chr", "")) & (leads.pos_grch38 == r.POS)).any(), axis=1)]
named = set(cov.best_proxy_qc_array.dropna()) | set(cov.best_usable_proxy.dropna())
named_pos = set(zip(cov.chr.astype(str), cov.best_proxy_pos.dropna().astype(int)))
hit_named = ex[ex.ID.isin(named) | ex.apply(lambda r: (r.CHROM.replace("chr", ""), r.POS) in named_pos, axis=1)]
proxy_hits = lc[(lc.r2_with_lead_1000G_EUR >= 0.8)]
R = [f"TO REMOVE: {len(ex)} SNPs = {len(rm182)} allele mismatch + {len(rm88)} of the 88 (|diff vs TOPMed| > 0.2) that are not within 0.2 of 1000G FIN",
     f"  of the 88: removed {len(rm88)} ({rm88.reason.value_counts().to_dict()}); KEPT {len(keep88)} (|study - FIN| <= 0.2; listed in keep_fin_variants_88.tsv)",
     f"  the other 358 mismatched-frequency SNPs are kept: {int((~big).sum())}",
     f"  overlap between the allele-mismatch list and the 446: {len(set(e.ID) & set(c.ID))}",
     f"LEAD CHECK (29 leads; {len(leads)} with a GRCh38 position):",
     f"  candidates that ARE a lead (position or rsID): {len(direct_lead)} {direct_lead.ID.tolist()}",
     f"  candidates within +-500 kb of a lead: {len(lc)} (distinct SNPs {lc.candidate_id.nunique() if len(lc) else 0}); with r2 >= 0.8 to a lead in 1000G EUR: {len(proxy_hits)}",
     f"  candidates that are the best / best-usable proxy recorded in step 12: {len(hit_named)} {hit_named.ID.tolist()}"]
if len(lc): R.append("  near-lead candidates:\n" + lc.to_string(index=False))
open(os.path.join(O, "exclusion_rule_report.txt"), "w").write("\n".join(R) + "\n")
print("\n".join(R))
