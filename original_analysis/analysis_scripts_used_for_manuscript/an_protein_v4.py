"""v4 (4 Nov 2026): HipSci protein cross-reference rebuilt from the ORIGINAL HipSci data
(Mirauta et al. 2020, eLife 9:e57390), replacing the v3 version that used values copied
from the September supplement (about 5,500 genes; symbol matching only).
 
Inputs (all public, from the eLife supplementary zip elife-57390-supp-v1.zip):
  SD2  Figure 1 source data 2: HipSci RNA, Ensembl gene x 202 lines
  SD3  Figure 1 source data 3: HipSci protein (TMT-MS), UniProt x 202 lines + reference line
  QTLp Figure 2 source data 1: genes with a cis pQTL (FDR < 10%)
  QTLe Figure 2 source data 2: genes with a cis eQTL (FDR < 10%)
  HGNC protein-coding gene list (symbol, Ensembl ID, previous symbols; rest.genenames.org, 4 Nov 2026)
 
Gene matching: our symbol -> HGNC (current symbol, else unique previous symbol) -> Ensembl ID
-> HipSci. Where one gene maps to several HipSci protein groups, the group with most peptides
is used. Our original symbols are kept in all outputs.
 
Same expressed-gene set (16,489 genes: normalised count >= 10 in >= 7 of 14 lines), models and
deciles as v3 (taken from Table S13d v3 / S4 so nothing else changes).
 
Usage: python3 an_protein_v4.py <package_dir> <SD2> <SD3> <QTLp> <QTLe> <hgnc.tsv> <out_dir>
"""
import sys, os
import numpy as np, pandas as pd
from scipy import stats
 
PKG, SD2, SD3, QTLP, QTLE, HGNC, OUT = sys.argv[1:8]
os.makedirs(OUT, exist_ok=True)
R = PKG + "/4_Analysis_results"
log = open(f"{OUT}/log_protein_v4.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n")
 
# ---------- HGNC mapping ----------
h = pd.read_csv(HGNC, sep="\t", dtype=str).fillna("")
cur = dict(zip(h.symbol, h.ensembl))
prev = {}
for s, e, pv in zip(h.symbol, h.ensembl, h.prev):
    for ps in filter(None, pv.split("|")):
        prev.setdefault(ps, set()).add((s, e))
prev = {k: next(iter(v)) for k, v in prev.items() if len(v) == 1 and k not in cur}
def map_gene(sym):
    if sym in cur: return sym, cur[sym], True
    if sym in prev: return prev[sym][0], prev[sym][1], True
    return None, None, False
 
# ---------- v3 gene table (same genes, models, deciles) ----------
v3 = pd.read_csv(f"{R}/protein_crossref/S13d_protein_crossref_all_genes.tsv.gz", sep="\t")
base = v3[["gene", "mean_norm_count", "log2FC_PMOS_vs_control", "p", "q_BH", "p_adjXIST", "expression_decile"]].copy()
m = [map_gene(g) for g in base.gene]
base["hgnc_symbol"] = [x[0] for x in m]; base["ensembl_gene_id"] = [x[1] or None for x in m]
base["protein_coding"] = [x[2] for x in m]
say(f"expressed genes: {len(base)}; protein-coding (HGNC): {base.protein_coding.sum()}")
 
# ---------- HipSci protein ----------
p = pd.read_csv(SD3, sep="\t", index_col=0, low_memory=False)
pcols = [c for c in p.columns if "@" in c and not c.startswith("HPSI0314i-bubh_3")]  # drop reference-line replicates
p = p.sort_values("protein_npeptides", ascending=False).drop_duplicates("ensembl_gene_id")
pl = np.log10(p.set_index("ensembl_gene_id")[pcols].astype(float))
pl.columns = [c.split("@")[0] for c in pl.columns]
pl = pl.T.groupby(level=0).mean().T                     # one value per line
prot = pd.DataFrame({"uniprot_id": p.set_index("ensembl_gene_id").feature_id,
                     "hipsci_mean_protein": p.set_index("ensembl_gene_id")[pcols].astype(float).mean(1),
                     "hipsci_protein_nlines": p.set_index("ensembl_gene_id").protein_nlines,
                     "hipsci_protein_npeptides": p.set_index("ensembl_gene_id").protein_npeptides,
                     "cv_protein": p.set_index("ensembl_gene_id").protein_cv})
say(f"HipSci protein table: {len(prot)} genes (Ensembl) from {p.shape[0]} best protein groups; {len(pcols)} line columns")
 
# ---------- HipSci RNA and per-gene RNA-protein correlation across lines ----------
r = pd.read_csv(SD2, sep="\t", index_col=0)
common = sorted(set(r.columns) & set(pl.columns))
genes = sorted(set(r.index) & set(pl.index))
rr, pp = r.loc[genes, common].values, pl.loc[genes, common].values
corr = []
for a, b in zip(rr, pp):
    ok = ~np.isnan(a) & ~np.isnan(b)
    corr.append(stats.pearsonr(a[ok], b[ok]).statistic if ok.sum() >= 30 else np.nan)
prot["hipsci_mean_rna_SD2"] = r.mean(1)   # SD2 units (log scale), not comparable to v3 column
prot["hipsci_rna_protein_corr"] = pd.Series(corr, index=genes)
say(f"RNA-protein correlation: {prot.hipsci_rna_protein_corr.notna().sum()} genes, {len(common)} matched lines, "
    f"median r = {prot.hipsci_rna_protein_corr.median():.3f}")
 
# ---------- QTL flags ----------
pq = pd.read_csv(QTLP, sep="\t"); eq = pd.read_csv(QTLE, sep="\t")
pq_ids = set(pq.ensembl_gene_id); eq_ids = set(eq.iloc[:, 0])
prot["pqtl"] = prot.index.isin(pq_ids); prot["eqtl"] = prot.index.isin(eq_ids)
say(f"pQTL genes: {len(pq_ids)}, eQTL genes: {len(eq_ids)}")
 
# ---------- join ----------
d = base.merge(prot, left_on="ensembl_gene_id", right_index=True, how="left")
d["protein_detected_HipSci"] = d.hipsci_mean_protein.notna()
d.loc[~d.protein_detected_HipSci, ["pqtl", "eqtl"]] = np.nan
cols = ["gene", "hgnc_symbol", "ensembl_gene_id", "protein_coding", "mean_norm_count", "expression_decile",
        "log2FC_PMOS_vs_control", "p", "q_BH", "p_adjXIST", "protein_detected_HipSci", "uniprot_id",
        "hipsci_mean_rna_SD2", "hipsci_mean_protein", "cv_protein", "hipsci_protein_nlines",
        "hipsci_protein_npeptides", "hipsci_rna_protein_corr", "pqtl", "eqtl"]
d = d[cols].sort_values("p")
n_all, n_pc = len(d), int(d.protein_coding.sum())
det_all, det_pc = int(d.protein_detected_HipSci.sum()), int(d[d.protein_coding].protein_detected_HipSci.sum())
say(f"protein detected: {det_all} of {n_all} expressed genes ({100*det_all/n_all:.1f}%); "
    f"{det_pc} of {n_pc} protein-coding ({100*det_pc/n_pc:.1f}%)")
 
# S13a: by decile (same deciles as v3), both denominators
dec = d.groupby("expression_decile").agg(genes=("gene", "size"), protein_coding=("protein_coding", "sum"),
                                         protein_detected=("protein_detected_HipSci", "sum"),
                                         median_count=("mean_norm_count", "median")).reset_index()
dec["pct_detected"] = 100 * dec.protein_detected / dec.genes
dec["pct_detected_of_protein_coding"] = 100 * dec.protein_detected / dec.protein_coding
say(dec.round(1).to_string(index=False))
 
# S13c: protein-supported candidates (same rule as v3)
cand = d[(d.p < .05) & d.protein_detected_HipSci & (d.hipsci_rna_protein_corr >= .3)].sort_values("p")
v3c = pd.read_csv(f"{R}/protein_crossref/S13c_protein_supported_candidates.tsv", sep="\t")
say(f"candidates p<0.05 & protein & r>=0.3: {len(cand)} (v3: {len(v3c)}; overlap {len(set(cand.gene) & set(v3c.gene))})")
say(f"nominal p<0.05 genes with protein evidence: {int(((d.p<.05)&d.protein_detected_HipSci).sum())} of {int((d.p<.05).sum())}")
 
# S13b: PCOS-locus genes
v3b = pd.read_csv(f"{R}/protein_crossref/S13b_protein_PCOS_locus_genes.tsv", sep="\t")
loc = v3b[["gene", "source", "expressed_iPSC"]].copy()
mm = [map_gene(g) for g in loc.gene]; loc["ensembl_gene_id"] = [x[1] for x in mm]
loc = loc.merge(prot[["uniprot_id", "hipsci_mean_protein", "hipsci_protein_nlines", "hipsci_protein_npeptides",
                      "hipsci_rna_protein_corr", "pqtl", "eqtl"]], left_on="ensembl_gene_id", right_index=True, how="left")
loc["protein_detected_HipSci"] = loc.hipsci_mean_protein.notna()
loc["protein_detected_HipSci_v3"] = loc.gene.map(v3b.set_index("gene").protein_detected_HipSci)
e = loc[loc.source.str.startswith("consensus") & loc.expressed_iPSC]
say(f"PCOS consensus genes expressed: {len(e)}; protein-detected: {int(e.protein_detected_HipSci.sum())} "
    f"(v3: {int(e.protein_detected_HipSci_v3.sum())}): {list(e.gene[e.protein_detected_HipSci])}")
say("changed calls:", loc[loc.protein_detected_HipSci != loc.protein_detected_HipSci_v3][["gene", "source", "protein_detected_HipSci_v3", "protein_detected_HipSci", "hipsci_protein_nlines", "hipsci_protein_npeptides"]].to_string(index=False))
 
# S14 literature screen: refresh its protein column
s14 = pd.read_csv(f"{R}/transcriptome/S14_literature_candidate_genes.tsv", sep="\t")
ens = [map_gene(g)[1] for g in s14.gene]
s14_old = s14.protein_detected_HipSci.copy()
s14["protein_detected_HipSci"] = [bool(x) and x in prot.index for x in ens]
say(f"S14 literature genes protein-detected: v3 {int(s14_old.sum())} -> v4 {int(s14.protein_detected_HipSci.sum())}")
 
d.to_csv(f"{OUT}/S13d_protein_crossref_all_genes.tsv.gz", sep="\t", index=False, float_format="%.6g")
dec.to_csv(f"{OUT}/S13a_protein_by_expression_decile.tsv", sep="\t", index=False)
cand.to_csv(f"{OUT}/S13c_protein_supported_candidates.tsv", sep="\t", index=False, float_format="%.6g")
loc.drop(columns="protein_detected_HipSci_v3").to_csv(f"{OUT}/S13b_protein_PCOS_locus_genes.tsv", sep="\t", index=False, float_format="%.6g")
s14.to_csv(f"{OUT}/S14_literature_candidate_genes.tsv", sep="\t", index=False, float_format="%.6g")
log.close()
