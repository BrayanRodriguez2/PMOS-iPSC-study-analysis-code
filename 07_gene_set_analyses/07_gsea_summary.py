#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/07_gsea_summary.py   (original step script: 76_gsea_summary.py)
# Produces: Summary tables and calibration figure.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 76 - summary tables and calibration figure of the gene-set analyses.
  gsea_summary_custom_sets.tsv : one row per (ranking, custom set): set size, genes present, NES, p, padj (fgsea, gene permutation) and leading edge; plus, for the two OLS rankings, the exact
                                 sample-label permutation p (step 75).
  gsea_summary_hallmark.tsv    : per ranking, number of Hallmark sets with p < 0.05 / padj < 0.05 (fgsea) and top 5 sets (NES, padj).
  hallmark_labelperm_vs_fgsea.tsv + figures/G1_hallmark_fgsea_vs_labelperm.png : fgsea p (gene permutation) against the exact label-permutation p for the 50 Hallmark sets.
Usage: python 07_gene_set_analyses/07_gsea_summary.py <gsea_dir>
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

gd = sys.argv[1]
f = pd.read_csv(f"{gd}/tables/fgsea_results_all.tsv", sep="\t"); lp = pd.read_csv(f"{gd}/tables/fgsea_labelperm_calibration.tsv", sep="\t")
c = f[f.collection != "Hallmark"].copy()
lpm = {("ols_signed", r.pathway): r for r in lp[(lp.lines == "all14") & (lp.collection != "Hallmark")].itertuples()}
lpm.update({("ols_signed_no350", r.pathway): r for r in lp[(lp.lines == "no350") & (lp.collection != "Hallmark")].itertuples()})
c["p_exact_label_permutation"] = [lpm[(r, p)].p_labelperm if (r, p) in lpm else np.nan for r, p in zip(c.ranking, c.pathway)]
c["NES_label_permutation"] = [lpm[(r, p)].NES_labelperm if (r, p) in lpm else np.nan for r, p in zip(c.ranking, c.pathway)]
c = c.rename(columns={"pval": "p_fgsea", "padj_BH_collection": "padj_fgsea_BH"})
c[["ranking", "collection", "pathway", "n_in_set", "n_present_in_ranking", "NES", "p_fgsea", "padj_fgsea_BH", "NES_label_permutation", "p_exact_label_permutation", "leadingEdge"]].to_csv(f"{gd}/tables/gsea_summary_custom_sets.tsv", sep="\t", index=False, float_format="%.4g")
rows = []
for r, d in f[f.collection == "Hallmark"].groupby("ranking"):
    d = d.sort_values("pval"); top = "; ".join(f"{p.replace('HALLMARK_', '')} ({n:+.2f}, padj {q:.3g})" for p, n, q in zip(d.pathway.head(5), d.NES.head(5), d.padj_BH_collection.head(5)))
    rows.append(dict(ranking=r, sets=len(d), p_lt_0_05=int((d.pval < .05).sum()), padj_lt_0_05=int((d.padj_BH_collection < .05).sum()), positive_NES_padj_lt_0_05=int(((d.NES > 0) & (d.padj_BH_collection < .05)).sum()),
                     negative_NES_padj_lt_0_05=int(((d.NES < 0) & (d.padj_BH_collection < .05)).sum()), top5=top))
pd.DataFrame(rows).to_csv(f"{gd}/tables/gsea_summary_hallmark.tsv", sep="\t", index=False)
h = lp[lp.collection == "Hallmark"].copy(); h.to_csv(f"{gd}/tables/hallmark_labelperm_vs_fgsea.tsv", sep="\t", index=False, float_format="%.5g")
os.makedirs(f"{gd}/figures", exist_ok=True)
fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
for a, (tag, ttl) in zip(ax, [("all14", "14 lines (OLS ranking)"), ("no350", "without HEL350.3 (13 lines)")]):
    d = h[h.lines == tag]; a.scatter(-np.log10(d.p_fgsea.clip(lower=1e-4)), -np.log10(d.p_labelperm.clip(lower=1e-4)), s=18, c="#2c7fb8")
    a.axhline(-np.log10(.05), color="k", ls=":", lw=.7); a.axvline(-np.log10(.05), color="k", ls=":", lw=.7); a.plot([0, 4], [0, 4], "r:", lw=.8)
    a.set_xlabel("fgsea p (gene permutation), -log10"); a.set_ylabel("exact sample-label permutation p, -log10"); a.set_title(f"Hallmark, {ttl}: {int((d.p_fgsea < .05).sum())} sets p < 0.05 by fgsea vs {int((d.p_labelperm < .05).sum())} by label permutation", fontsize=8)
    for r in d.sort_values("p_fgsea").head(4).itertuples(): a.annotate(r.pathway.replace("HALLMARK_", "")[:18], (-np.log10(max(r.p_fgsea, 1e-4)), -np.log10(max(r.p_labelperm, 1e-4))), fontsize=6, xytext=(3, 3), textcoords="offset points")
fig.tight_layout(); fig.savefig(f"{gd}/figures/G1_hallmark_fgsea_vs_labelperm.png", dpi=160)
print(pd.read_csv(f"{gd}/tables/gsea_summary_hallmark.tsv", sep="\t").drop(columns="top5").to_string(index=False))
print(pd.read_csv(f"{gd}/tables/gsea_summary_custom_sets.tsv", sep="\t").drop(columns="leadingEdge").round(3).to_string(index=False))
