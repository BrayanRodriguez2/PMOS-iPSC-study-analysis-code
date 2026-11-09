#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/09_deseq2_limma.py   (original step script: 67b_methods_compare.py)
# Produces: DESeq2 and limma-voom on raw counts and comparison with the OLS.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 67b - compare DESeq2 and limma-voom (raw counts, step 67a) with the OLS on log2 normalised counts (step 65, ours_s2) and with the original Table S4:
genes with p < 0.05 and FDR < 0.1, Spearman correlation of the rankings (signed statistic: OLS t-like = sign(log2FC) x -log10 p), overlap of the nominal sets, and a figure.
Usage: python 06_RNAseq_counts_and_DE_validation/09_deseq2_limma.py <methods_dir> <ols table de_ols.tsv.gz> <the original S4 tsv> <out_dir>
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

md, olsf, s4f, out = sys.argv[1:5]; os.makedirs(os.path.join(out, "figures"), exist_ok=True)
ols = pd.read_csv(olsf, sep="\t").set_index("gene"); s4 = pd.read_csv(s4f, sep="\t").set_index("gene")
rows, rep = [], []
tabs = {}
for lab, olp, olq, olb in [("group", "p", "q_BH", "log2FC"), ("group_xist", "p_adjXIST", "q_BH_adjXIST", "log2FC_adjXIST")]:
    base = pd.DataFrame({"ols_log2FC": ols[olb], "ols_p": ols[olp], "ols_q": ols[olq]}); base["ols_signed"] = np.sign(base.ols_log2FC) * -np.log10(base.ols_p)
    for meth in ["deseq2", "limma_voom"]:
        d = pd.read_csv(os.path.join(md, f"{meth}_{lab}.tsv"), sep="\t").set_index("gene")
        d["signed"] = np.sign(d.log2FC) * -np.log10(d.p)
        j = base.join(d, how="inner").dropna(subset=["p"])
        nom_m, nom_o = set(j.index[j.p < .05]), set(j.index[j.ols_p < .05])
        rho = stats.spearmanr(j.ols_signed, j.signed)[0]; r_fc = np.corrcoef(j.ols_log2FC, j.log2FC)[0, 1]
        rows.append(dict(model=lab, method=meth, genes_tested=len(j), genes_p05=int((j.p < .05).sum()), genes_FDR10=int((j.padj < .1).sum()), min_adj_p=float(j.padj.min()), ols_genes_p05=len(nom_o),
                         overlap_nominal=len(nom_m & nom_o), spearman_signed_p_vs_ols=rho, pearson_log2FC_vs_ols=r_fc))
        tabs[(lab, meth)] = j
    # the original S4
    if lab == "group":
        sj = ols.join(s4, rsuffix="_S4", how="inner")
        rows.append(dict(model=lab, method="the original S4 (the original counts)", genes_tested=len(sj), genes_p05=int((sj.p_S4 < .05).sum()), genes_FDR10=int((sj.q_BH_S4 < .1).sum()), ols_genes_p05=int((sj.p < .05).sum()),
                         overlap_nominal=int(((sj.p_S4 < .05) & (sj.p < .05)).sum()), spearman_signed_p_vs_ols=stats.spearmanr(np.sign(sj.log2FC) * -np.log10(sj.p), np.sign(sj.log2FC_PMOS_vs_control) * -np.log10(sj.p_S4))[0],
                         pearson_log2FC_vs_ols=np.corrcoef(sj.log2FC, sj.log2FC_PMOS_vs_control)[0, 1]))
R = pd.DataFrame(rows); R.to_csv(os.path.join(out, "methods_comparison.tsv"), sep="\t", index=False, float_format="%.4g")
print(R.round(3).to_string(index=False)); open(os.path.join(out, "methods_comparison.txt"), "w").write(R.round(4).to_string(index=False) + "\n")
fig, ax = plt.subplots(2, 2, figsize=(10, 9))
for k, (meth, ttl) in enumerate([("deseq2", "DESeq2 (raw counts, NB Wald)"), ("limma_voom", "limma-voom (raw counts)")]):
    j = tabs[("group", meth)]
    ax[0, k].scatter(j.ols_signed, j.signed, s=3, alpha=.4); ax[0, k].set_xlabel("OLS on log2(norm + 1): sign x -log10 p"); ax[0, k].set_ylabel(f"{ttl}: sign x -log10 p")
    r = R[(R.model == "group") & (R.method == meth)].iloc[0]; ax[0, k].set_title(f"Spearman {r.spearman_signed_p_vs_ols:.2f}; p<0.05: {r.genes_p05} (OLS {r.ols_genes_p05}); FDR<0.1: {r.genes_FDR10}", fontsize=9)
    ax[1, k].scatter(j.ols_log2FC, j.log2FC, s=3, alpha=.4); ax[1, k].set_xlabel("OLS log2FC"); ax[1, k].set_ylabel("log2FC " + meth); ax[1, k].set_title(f"log2FC r = {r.pearson_log2FC_vs_ols:.2f}", fontsize=9)
    lim = [-3, 3]; ax[1, k].plot(lim, lim, "r:", lw=.8); ax[1, k].set_xlim(lim); ax[1, k].set_ylim(lim)
fig.tight_layout(); fig.savefig(os.path.join(out, "figures", "V6_deseq2_limma_vs_ols.png"), dpi=160)
