#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/09_prscs_coverage.py   (original step script: 44_prscs_coverage.py)
# Produces: Coverage of the PRS-CS LD panel and of the 29 lead SNPs in the imputed data.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 44 - Post-imputation QC, part (f): coverage of the PRS-CS LD panel (UKBB EUR, HapMap3 SNPs) and of the 29 lead SNPs (ST2) in the imputed data,
compared with the array-only situation of step 17/18 (5 of 29 leads genotyped, 15/29 usable by PRS-CS, 169 k panel SNPs usable). Diagnostic only.

Method
  * PRS-CS panel = snpinfo_ukbb_hm3 (rsID, A1, A2; its BP is GRCh37, so positions are NOT taken from it).
  * GRCh38 chromosome/position/alleles of each panel SNP come from the harmonised GWAS summary statistics (GCST90570579, rsid column). A panel SNP
    is "potentially usable by PRS-CS" if it is in the summary statistics that were prepared for PRS-CS (prscs/inputs/sumstats_prscs.txt) with the same
    allele set (direct or strand-complement, as PRS-CS accepts both).
  * "Present in the imputed data" = a record at the same chromosome and GRCh38 position whose REF/ALT set equals the panel's allele pair (direct or
    complement; palindromic A/T-C/G SNPs only match directly, as the strand cannot be told apart). Presence is reported for R2 >= 0.3 (everything the server delivered
    passed it) and R2 >= 0.8, split into Typed (genotyped) / Imputed, and with the extra condition MAF >= 0.01 (server estimate over 140 samples).
  * Leads: position/alleles from lead_snps_lookup.tsv; status in the imputed data (Typed/Imputed/absent), R2, ER2, MAF; whether the lead is in the PRS-CS panel.
    A lead that is not in the HapMap3 panel cannot be used by PRS-CS itself (its proxies in the panel were not recomputed here).
Inputs : info_tables/chr*.tsv.gz (step 40), reference/ldblk_ukbb_eur/snpinfo_ukbb_hm3, gwas/GCST90570579_harmonised/*.h.tsv.gz,
         prscs/inputs/{sumstats_prscs.txt, V3_panel_ukbb.bim}, results/loci/{lead_snps_lookup.tsv, locus_coverage.tsv}
Outputs: results/imputation_qc/{prscs_panel_coverage.tsv, leads_imputed_coverage.tsv, coverage_report.txt}, figures/F7_coverage.png
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tables-dir", required=True); ap.add_argument("--snpinfo", required=True); ap.add_argument("--sumstats-gz", required=True)
ap.add_argument("--sumstats-prscs", required=True); ap.add_argument("--leads", required=True)
ap.add_argument("--locus-cov", required=True); ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(os.path.join(a.outdir, "figures"), exist_ok=True)
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
comp = lambda s: "".join(COMP.get(x, "N") for x in s)
pal = lambda s: len(s) == 2 and comp(s[0]) == s[1]

# imputed variants keyed by (chr, pos)
parts = []
for c in range(1, 23):
    d = pd.read_csv(os.path.join(a.tables_dir, f"chr{c}.tsv.gz"), sep="\t", usecols=["POS", "ID", "REF", "ALT", "MAF", "R2", "ER2", "TYPE"])
    d["chr"] = c; parts.append(d)
imp = pd.concat(parts, ignore_index=True)
imp["aset"] = [min(r, al) + max(r, al) if len(r) == 1 and len(al) == 1 else "NA" for r, al in zip(imp.REF, imp.ALT)]
imp = imp[imp.aset != "NA"]
print(f"imputed SNPs (single-base alleles): {len(imp):,}", flush=True)

# panel + sumstats
panel = pd.read_csv(a.snpinfo, sep="\t", dtype={"CHR": int})
ss = pd.read_csv(a.sumstats_gz, sep="\t", usecols=["chromosome", "base_pair_location", "effect_allele", "other_allele", "rsid"], dtype=str).dropna(subset=["rsid"])
ss = ss[ss.chromosome.isin([str(i) for i in range(1, 23)])].drop_duplicates("rsid")
ss = ss.rename(columns={"rsid": "SNP", "chromosome": "chr38", "base_pair_location": "pos38"}); ss["pos38"] = ss.pos38.astype(int); ss["chr38"] = ss.chr38.astype(int)
used = pd.read_csv(a.sumstats_prscs, sep=r"\s+", usecols=["SNP", "A1", "A2"])
pp = panel.merge(ss, on="SNP", how="inner").merge(used.rename(columns={"A1": "A1_ss", "A2": "A2_ss"}), on="SNP", how="inner")
pp["paset"] = [min(x, y) + max(x, y) for x, y in zip(pp.A1, pp.A2)]
pp["ssset"] = [min(x, y) + max(x, y) for x, y in zip(pp.A1_ss, pp.A2_ss)]
ok = (pp.paset == pp.ssset) | (pp.paset == pp.ssset.map(comp).map(lambda s: min(s) + max(s) if len(s) == 2 else s))
pp = pp[ok].copy()
L = [f"PRS-CS panel SNPs: {len(panel):,} | with GRCh38 position from the summary statistics and allele-compatible with the sumstats used for PRS-CS: {len(pp):,}"]


impr = imp[["chr", "POS", "aset", "ID", "MAF", "R2", "ER2", "TYPE"]].rename(columns={"MAF": "MAF_srv", "ID": "ID_srv", "chr": "chr_srv"})   # the panel also has a MAF column


def match(df):
    """presence of each panel SNP in the imputed data (direct allele set, or complement when not palindromic)"""
    direct = df.merge(impr, left_on=["chr38", "pos38", "paset"], right_on=["chr_srv", "POS", "aset"], how="left")
    comp_set = df.paset.map(lambda s: min(comp(s)) + max(comp(s)))
    nonpal = ~df.paset.map(lambda s: comp(s[0]) == s[1])
    dd = df[nonpal].assign(cset=comp_set[nonpal])
    cm = dd.merge(impr, left_on=["chr38", "pos38", "cset"], right_on=["chr_srv", "POS", "aset"], how="left")
    return direct, cm


direct, cm = match(pp)
hit = pd.concat([direct[direct.R2.notna()], cm[cm.R2.notna()]]).sort_values("R2", ascending=False).drop_duplicates("SNP")
L.append(f"panel SNPs present in the imputed data: {len(hit):,} ({100 * len(hit) / len(pp):.1f} % of the usable panel set)")
rows = []
for nm, m in [("present (R2 >= 0.3 = all delivered)", hit.R2 >= 0.3), ("R2 >= 0.8", hit.R2 >= 0.8), ("R2 >= 0.3 and MAF >= 0.01", (hit.R2 >= 0.3) & (hit.MAF_srv >= 0.01)),
              ("R2 >= 0.8 and MAF >= 0.01", (hit.R2 >= 0.8) & (hit.MAF_srv >= 0.01))]:
    s = hit[m]
    rows.append(dict(criterion=nm, n=len(s), pct_of_usable_panel=100 * len(s) / len(pp), pct_of_full_panel=100 * len(s) / len(panel),
                     n_typed=int((s.TYPE == "Typed").sum()), n_imputed=int((s.TYPE == "Imputed").sum())))
cov = pd.DataFrame(rows)
before = 169093      # V1 of step 15-17 (array QC SNPs in the UKBB panel and in the sumstats); V3_panel_ukbb.bim is the whole panel, not the array
cov.loc[len(cov)] = ["BEFORE imputation (V1): array QC SNPs usable by PRS-CS", before, 100 * before / len(pp), 100 * before / len(panel), before, 0]
cov.to_csv(os.path.join(a.outdir, "prscs_panel_coverage.tsv"), sep="\t", index=False, float_format="%.2f")
L.append(cov.round(2).to_string(index=False))
chr_cov = hit[hit.R2 >= 0.8].groupby("chr38").size().reindex(range(1, 23)).fillna(0).astype(int)

# leads
leads = pd.read_csv(a.leads, sep="\t", dtype={"chr": str}).dropna(subset=["pos_grch38"])
lc = pd.read_csv(a.locus_cov, sep="\t", dtype={"chr": str}).set_index("rsid")
lrows = []
for r in leads.itertuples():
    pa = min(r.effect_allele, r.other_allele) + max(r.effect_allele, r.other_allele)
    cand = imp[(imp.chr == int(r.chr)) & (imp.POS == int(r.pos_grch38)) & ((imp.aset == pa) | (imp.aset == min(comp(pa)) + max(comp(pa))) )]
    inpanel = r.rsid in set(panel.SNP)
    row = dict(lead=r.rsid, chr=r.chr, pos_grch38=int(r.pos_grch38), p_value=r.p_value, alleles=f"{r.effect_allele}/{r.other_allele}", in_prscs_panel=inpanel,
               before_on_qc_array=lc.loc[r.rsid, "lead_on_qc_array"], before_usable_by_prscs=lc.loc[r.rsid, "covered_usable_by_prscs"])
    if len(cand):
        b = cand.sort_values("R2", ascending=False).iloc[0]
        row.update(status=b.TYPE, imputed_id=b.ID, R2=b.R2, ER2=b.ER2, MAF=b.MAF, present_R2_ge_0_3=True, present_R2_ge_0_8=bool(b.R2 >= 0.8))
    else:
        row.update(status="absent (R2 < 0.3 or not in panel)", imputed_id="", R2=np.nan, ER2=np.nan, MAF=np.nan, present_R2_ge_0_3=False, present_R2_ge_0_8=False)
    row["usable_by_prscs_now_R2_ge_0_3"] = bool(inpanel and row["present_R2_ge_0_3"]); row["usable_by_prscs_now_R2_ge_0_8"] = bool(inpanel and row["present_R2_ge_0_8"])
    lrows.append(row)
lt = pd.DataFrame(lrows)
lt.to_csv(os.path.join(a.outdir, "leads_imputed_coverage.tsv"), sep="\t", index=False, float_format="%.4g")
L.append(f"\n29 leads: present in the imputed data (R2 >= 0.3): {int(lt.present_R2_ge_0_3.sum())} | R2 >= 0.8: {int(lt.present_R2_ge_0_8.sum())} | genotyped (Typed): {int((lt.status == 'Typed').sum())} | imputed: {int((lt.status == 'Imputed').sum())}")
L.append(f"  in the PRS-CS panel: {int(lt.in_prscs_panel.sum())} | usable by PRS-CS now (lead itself in panel and present): R2 >= 0.3 {int(lt.usable_by_prscs_now_R2_ge_0_3.sum())}, R2 >= 0.8 {int(lt.usable_by_prscs_now_R2_ge_0_8.sum())}"
         f" | BEFORE imputation (lead or proxy usable): {int((lt.before_usable_by_prscs.astype(str) == 'yes').sum())}/29; leads on the QC array: {int((lt.before_on_qc_array.astype(str) != 'no').sum())}/29")
L.append(lt[["lead", "chr", "p_value", "status", "R2", "ER2", "MAF", "in_prscs_panel", "before_on_qc_array", "before_usable_by_prscs"]].round(4).to_string(index=False))
open(os.path.join(a.outdir, "coverage_report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
lab = ["array QC\n(before)", "imputed\nR2 >= 0.3", "imputed\nR2 >= 0.8", "imputed\nR2 >= 0.8\nMAF >= 0.01"]
vals = [before, int(cov.n.iloc[0]), int(cov.n.iloc[1]), int(cov.n.iloc[3])]
ax[0].bar(lab, vals, color=["#7f8c8d", "#2c7fb8", "#2c7fb8", "#2c7fb8"]); ax[0].axhline(len(pp), color="k", ls=":", lw=.8); ax[0].set_ylabel("PRS-CS panel SNPs available")
for i, v in enumerate(vals): ax[0].text(i, v, f"{v:,}\n({100 * v / len(pp):.0f} %)", ha="center", va="bottom", fontsize=8)
ax[0].set_title(f"HapMap3 SNPs of the UKBB EUR panel usable with the sumstats (dotted = {len(pp):,})", fontsize=9)
ax[1].bar(range(1, 23), chr_cov.values, color="#2c7fb8"); ax[1].set_xlabel("chromosome"); ax[1].set_ylabel("panel SNPs with R2 >= 0.8"); ax[1].set_title("by chromosome", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "figures", "F7_coverage.png"), dpi=160)
