#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/15_prscs_bim_imputed.py   (original step script: 49_prscs_bim_imputed.py)
# Produces: PRS-CS .bim with panel rsIDs matched to the imputed data by position and alleles.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 49 - PRS-CS inputs for the IMPUTED clean dataset: a .bim carrying the RSID of the LD-panel SNP for every panel SNP that is present in the clean file
(PRS-CS matches by rsID; the clean PLINK2 file uses IDs chr:pos:ref:alt), and the map between the two identifiers, used later to return the weights to the
chr:pos:ref:alt IDs for scoring.

Method (same logic as step 44, but now with the CLEAN file = v2_imputed_clean after the MAF/HWE/ER2 filters, via clean/variant_table.tsv.gz)
  1. Panel = snpinfo of the LD panel (UKBB EUR for V1/V3, 1000G EUR for V2; HapMap3 SNPs; columns CHR SNP BP A1 A2 MAF; BP is GRCh37, NOT used).
  2. GRCh38 chromosome/position of each panel SNP = those of the same rsID in the harmonised summary statistics (GCST90570579, first occurrence).
  3. A panel SNP is kept only if the sumstats prepared for PRS-CS (prscs/inputs/sumstats_prscs.txt) have the same alleles (direct or strand complement).
  4. Match to the clean file by (chromosome, GRCh38 position, allele SET): direct (panel alleles = REF/ALT in either order) or complementary strand.
     Strand-ambiguous panel SNPs (A/T or C/G) are EXCLUDED (the strand cannot be told from the alleles); so are
       - panel SNPs with more than one clean variant at that position with the same allele set (duplicates/ambiguous),
       - clean variants claimed by more than one panel rsID (the clean variant would be scored twice).
  5. Output bim (PLINK .bim columns: CHR, rsID, 0, BP38, A1, A2) with the clean alleles (ALT as A1, REF as A2); PRS-CS accepts allele pairs on either strand.
Outputs (--outdir): V1_imputed.bim (UKBB), V2_imputed.bim (1000G), map_ukbb.tsv.gz, map_1kg.tsv.gz
  (rsID, CHR, BP38, panel_A1, panel_A2, ID (clean), REF, ALT, strand, R2, ER2, TYPE, MAF_server, AF129), bim_counts.tsv (SNPs after each step) and bim_report.txt.
The V3 bim (whole panel) already exists (prscs/inputs/V3_panel_ukbb.bim) and does not depend on the genotypes; its PRS-CS output is re-used as is.
"""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--variant-table", required=True); ap.add_argument("--snpinfo-ukbb", required=True); ap.add_argument("--snpinfo-1kg", required=True)
ap.add_argument("--sumstats-gz", required=True); ap.add_argument("--sumstats-prscs", required=True); ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
comp = lambda s: "".join(COMP.get(x, "N") for x in s)
aset = lambda x, y: min(x, y) + max(x, y)

vt = pd.read_csv(a.variant_table, sep="\t", dtype={"CHROM": int}, usecols=["CHROM", "POS", "ID", "REF", "ALT", "R2", "ER2", "MAF_server", "TYPE", "AF129"])
vt["aset"] = [aset(r, al) for r, al in zip(vt.REF, vt.ALT)]
print(f"clean variants: {len(vt):,}", flush=True)
grp = vt.groupby(["CHROM", "POS", "aset"]).size().rename("n_clean_same_site").reset_index()
vt = vt.merge(grp, on=["CHROM", "POS", "aset"], how="left")

ss = pd.read_csv(a.sumstats_gz, sep="\t", usecols=["chromosome", "base_pair_location", "rsid"], dtype=str).dropna(subset=["rsid"])
ss = ss[ss.chromosome.isin([str(i) for i in range(1, 23)])].drop_duplicates("rsid")
ss = ss.rename(columns={"rsid": "SNP", "chromosome": "CHR38", "base_pair_location": "BP38"}); ss["CHR38"] = ss.CHR38.astype(int); ss["BP38"] = ss.BP38.astype(int)
used = pd.read_csv(a.sumstats_prscs, sep=r"\s+", usecols=["SNP", "A1", "A2"]).rename(columns={"A1": "A1_ss", "A2": "A2_ss"})
used["ssset"] = [aset(x, y) for x, y in zip(used.A1_ss, used.A2_ss)]

L, counts = [], []


def one(name, snpinfo, bim_out, map_out):
    panel = pd.read_csv(snpinfo, sep="\t", dtype={"CHR": int}).rename(columns={"A1": "panel_A1", "A2": "panel_A2"})
    c = {"panel": name, "1_panel_snps": len(panel)}
    p = panel.merge(ss, on="SNP", how="inner"); c["2_with_grch38_position_from_sumstats"] = len(p)
    p = p.merge(used[["SNP", "ssset"]], on="SNP", how="inner")
    p["paset"] = [aset(x, y) for x, y in zip(p.panel_A1, p.panel_A2)]
    ok = (p.paset == p.ssset) | (p.paset == p.ssset.map(lambda s: aset(comp(s)[0], comp(s)[1]) if len(s) == 2 else s))
    p = p[ok].copy(); c["3_allele_compatible_with_prscs_sumstats"] = len(p)
    p["pal"] = [comp(x) == y for x, y in zip(p.panel_A1, p.panel_A2)]
    c["3b_of_which_strand_ambiguous_AT_CG"] = int(p.pal.sum())
    # ---- direct and complement matches against the clean file
    key = vt[["CHROM", "POS", "aset", "ID", "REF", "ALT", "R2", "ER2", "MAF_server", "TYPE", "AF129", "n_clean_same_site"]]
    d = p.merge(key, left_on=["CHR38", "BP38", "paset"], right_on=["CHROM", "POS", "aset"], how="left"); d["strand"] = "direct"
    pn = p[~p.pal].copy(); pn["cset"] = [aset(comp(x), comp(y)) for x, y in zip(pn.panel_A1, pn.panel_A2)]
    cm = pn.merge(key, left_on=["CHR38", "BP38", "cset"], right_on=["CHROM", "POS", "aset"], how="left"); cm["strand"] = "complement"
    cm = cm.drop(columns=["cset"])
    hit = pd.concat([d[d.ID.notna()], cm[cm.ID.notna()]], ignore_index=True)
    c["4_present_in_clean_file_position_and_alleles"] = int(hit.SNP.nunique())
    # a palindromic SNP can only have matched 'direct', but its strand is unknown -> excluded
    hit = hit[~hit.pal]; c["5_after_excluding_strand_ambiguous"] = int(hit.SNP.nunique())
    # panel SNP matching clean variants through BOTH strands or with several clean variants at the site -> ambiguous
    nm = hit.groupby("SNP").size(); amb = set(nm[nm > 1].index) | set(hit[hit.n_clean_same_site > 1].SNP)
    c["5b_dropped_ambiguous_site_or_two_matches"] = len(amb)
    hit = hit[~hit.SNP.isin(amb)]
    # a clean variant claimed by two different panel rsIDs -> both dropped
    nid = hit.groupby("ID").SNP.nunique(); dupid = set(nid[nid > 1].index)
    c["5c_dropped_clean_variant_claimed_by_two_rsids"] = int(hit.ID.isin(dupid).sum())
    hit = hit[~hit.ID.isin(dupid)]
    # duplicated panel rsID rows (should not exist)
    c["5d_duplicated_rsid_rows"] = int(hit.SNP.duplicated().sum()); hit = hit.drop_duplicates("SNP")
    hit = hit.sort_values(["CHR38", "BP38"])
    c["6_final_snps_in_prscs_bim_and_map"] = len(hit)
    c["6b_of_which_R2_ge_0.8"] = int((hit.R2 >= 0.8).sum()); c["6c_of_which_typed"] = int((hit.TYPE == "Typed").sum())
    c["6d_pct_of_usable_panel (step 3)"] = round(100 * len(hit) / c["3_allele_compatible_with_prscs_sumstats"], 2)
    c["6e_strand_complement"] = int((hit.strand == "complement").sum())
    with open(bim_out, "w") as out:
        for r in hit.itertuples():
            out.write(f"{r.CHR38}\t{r.SNP}\t0\t{r.BP38}\t{r.ALT}\t{r.REF}\n")
    m = hit[["SNP", "CHR38", "BP38", "panel_A1", "panel_A2", "ID", "REF", "ALT", "strand", "R2", "ER2", "TYPE", "MAF_server", "AF129"]].rename(columns={"SNP": "rsID", "CHR38": "CHR"})
    m.to_csv(map_out, sep="\t", index=False, compression="gzip", float_format="%.6g")
    counts.append(c)
    L.append(f"[{name}] " + " | ".join(f"{k}: {v:,}" if isinstance(v, int) else f"{k}: {v}" for k, v in c.items() if k != "panel"))


one("UKBB EUR (V1, V3)", a.snpinfo_ukbb, os.path.join(a.outdir, "V1_imputed.bim"), os.path.join(a.outdir, "map_ukbb.tsv.gz"))
one("1000G EUR (V2)", a.snpinfo_1kg, os.path.join(a.outdir, "V2_imputed.bim"), os.path.join(a.outdir, "map_1kg.tsv.gz"))
pd.DataFrame(counts).T.reset_index().to_csv(os.path.join(a.outdir, "bim_counts.tsv"), sep="\t", index=False, header=False)
open(os.path.join(a.outdir, "bim_report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
