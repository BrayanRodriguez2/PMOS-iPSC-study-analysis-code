#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/14_er2_filter_list.py   (original step script: 47_er2_filter_list.py)
# Produces: List of genotyped variants with ER2 < 0.5 and their PRS-CS panel membership.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 47a - List of the genotyped variants with ER2 < 0.5 (leave-one-out R2 of the TOPMed imputation) that are removed from the clean imputed dataset
(decision of the analyst, 2026-11-04), and how many of them belong to the PRS-CS LD panel.

Input : imputation/clean/variant_table.tsv.gz (step 45), results/topmed_qc_only/{keep_fin_variants_88.tsv, compare_1000g_446.tsv} (subsets 75 / 358),
        PRS-CS panel snpinfo_ukbb_hm3, harmonised sumstats (GRCh38 positions), prscs/inputs/sumstats_prscs.txt
Output: imputation/clean/removed_er2_lt_0.5.tsv (the list: new ID, rsID, chr, pos, alleles, R2, ER2, MAF, subset, in_prscs_panel flags),
        removed_er2_lt_0.5_ids.txt (one ID per line, for plink2 --exclude), er2_filter_report.txt
PRS-CS panel membership, two ways (both reported): (1) by rsID in snpinfo_ukbb_hm3 (1,117,425 HapMap3 SNPs); (2) "usable by PRS-CS" = the panel SNP
is in the sumstats prepared for PRS-CS with compatible alleles and the same GRCh38 chromosome/position/allele pair (direct or strand-complement) as the
variant (same definition as step 44; 1,114,978 usable panel SNPs).
"""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--clean-dir", required=True); ap.add_argument("--qc-only-dir", required=True); ap.add_argument("--snpinfo", required=True)
ap.add_argument("--sumstats-gz", required=True); ap.add_argument("--sumstats-prscs", required=True)
a = ap.parse_args()
C = a.clean_dir
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
comp = lambda s: "".join(COMP.get(x, "N") for x in s)
srt = lambda s: min(s) + max(s)
t = pd.read_csv(os.path.join(C, "variant_table.tsv.gz"), sep="\t", dtype={"CHROM": str})
rm = t[(t.TYPE == "Typed") & (t.ER2 < 0.5)].copy()
k75 = pd.read_csv(os.path.join(a.qc_only_dir, "keep_fin_variants_88.tsv"), sep="\t", dtype={"CHROM": str})
c446 = pd.read_csv(os.path.join(a.qc_only_dir, "compare_1000g_446.tsv"), sep="\t", dtype={"CHROM": str})
key = lambda d: d.CHROM.str.replace("chr", "") + ":" + d.POS.astype(str) + ":" + d.REF + ":" + d.ALT
s75, s358 = set(key(k75)), set(key(c446[~c446.abs_diff_topmed_gt_0_2.astype(bool)]))
rm["subset"] = "other genotyped"
rm.loc[rm.ID.isin(s358), "subset"] = "358 discordant kept"
rm.loc[rm.ID.isin(s75), "subset"] = "75 kept (FIN-supported)"
# PRS-CS panel
panel = pd.read_csv(a.snpinfo, sep="\t")
ss = pd.read_csv(a.sumstats_gz, sep="\t", usecols=["chromosome", "base_pair_location", "rsid"], dtype=str).dropna(subset=["rsid"]).drop_duplicates("rsid")
ss = ss[ss.chromosome.isin([str(i) for i in range(1, 23)])].rename(columns={"rsid": "SNP", "chromosome": "chr38", "base_pair_location": "pos38"})
used = pd.read_csv(a.sumstats_prscs, sep=r"\s+", usecols=["SNP", "A1", "A2"]).rename(columns={"A1": "A1s", "A2": "A2s"})
pp = panel.merge(ss, on="SNP").merge(used, on="SNP")
pp["pa"] = [srt(x + y) for x, y in zip(pp.A1, pp.A2)]; pp["sa"] = [srt(x + y) for x, y in zip(pp.A1s, pp.A2s)]
pp = pp[(pp.pa == pp.sa) | (pp.pa == pp.sa.map(lambda s: srt(comp(s))))]
usable = set(zip(pp.chr38.astype(str), pp.pos38.astype(str), pp.pa))
rm["in_panel_by_rsID"] = rm.rsID.isin(set(panel.SNP))
rm["aset"] = [srt(r + al) for r, al in zip(rm.REF, rm.ALT)]
rm["in_prscs_usable_panel"] = [((str(c), str(p), s) in usable) or ((str(c), str(p), srt(comp(s))) in usable) for c, p, s in zip(rm.CHROM, rm.POS, rm.aset)]
rm = rm.sort_values(["subset", "ER2"])
cols = ["ID", "rsID", "CHROM", "POS", "REF", "ALT", "R2", "ER2", "MAF_server", "AF129", "subset", "in_panel_by_rsID", "in_prscs_usable_panel"]
rm[cols].to_csv(os.path.join(C, "removed_er2_lt_0.5.tsv"), sep="\t", index=False, float_format="%.5g")
rm[["ID"]].to_csv(os.path.join(C, "removed_er2_lt_0.5_ids.txt"), index=False, header=False)
L = [f"genotyped variants with ER2 < 0.5 in the clean file (to remove): {len(rm):,} | unique IDs {rm.ID.is_unique}",
     f"  ER2 < 0.3: {int((rm.ER2 < 0.3).sum()):,}; ER2 median {rm.ER2.median():.3f}; R2 (server) median {rm.R2.median():.3f}; MAF median {rm.MAF_server.median():.3f}",
     "  by subset: " + str(rm.subset.value_counts().to_dict()),
     f"IN THE PRS-CS PANEL: by rsID in snpinfo_ukbb_hm3: {int(rm.in_panel_by_rsID.sum())} of {len(rm)} ({100 * rm.in_panel_by_rsID.mean():.1f} %) | usable by PRS-CS (sumstats, alleles, GRCh38 position): {int(rm.in_prscs_usable_panel.sum())} ({100 * rm.in_prscs_usable_panel.mean():.1f} %)",
     "  in the usable panel, by subset: " + str(rm[rm.in_prscs_usable_panel].subset.value_counts().to_dict()),
     f"reference: the usable panel has {len(pp):,} SNPs; of the clean file's variants {int(t.rsID.isin(set(panel.SNP)).sum()):,} are in the panel by rsID",
     f"clean file after removal: {len(t) - len(rm):,} variants"]
open(os.path.join(C, "er2_filter_report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
