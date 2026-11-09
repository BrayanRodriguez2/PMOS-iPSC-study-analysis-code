#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/06_gsea_labelperm.py   (original step script: 75_gsea_labelperm.py)
# Produces: Exact sample-label permutation calibration of the GSEA.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 75 - calibration of the preranked GSEA by EXACT sample-label permutation (own code, numpy).
Why: fgsea permutes GENES, which assumes independent genes; with 14 samples and strongly co-expressed sets (e.g. oxidative phosphorylation) a single coordinated difference in one
donor split gives a tiny gene-permutation p although the same set would also be extreme in many random splits of the 14 lines. Here the null is built from the data themselves:
every possible assignment of the lines to the two groups (all C(14, 7) = 3,432 splits; 7 PMOS vs 7 CTRL; or C(13, 5) = 1,287 without HEL350.3), re-computing for each split the
per-gene OLS statistic (sign x -log10 p, i.e. the 'ols_signed' ranking of step 73) and the weighted enrichment score (ES, GSEA weight 1) of every set. As in GSEA, NES = ES / mean of the
null ES of the same sign; p = fraction of the splits (the observed one included) with |NES| >= the observed |NES|; BH within each collection.
Check: the observed ES must agree with fgsea's ES for the same ranking (step 74; small differences only from tied statistics).
Inputs: expressed-gene normalised counts (step 64, symbols), donor_meta.tsv, custom gene sets and Hallmark (written by 74 as gene_sets/hallmark.tsv), fgsea results (74).
Usage: python 07_gene_set_analyses/06_gsea_labelperm.py <gsea_dir> <norm_symbol.tsv.gz> <donor_meta.tsv>
"""
import itertools, os, re, sys
import numpy as np, pandas as pd
from scipy import stats

gd, normf, metaf = sys.argv[1:4]
meta = pd.read_csv(metaf, sep="\t", index_col=0)
def hel(x):
    m = re.search(r"HEL\s*_?(\d+)[._](\d+)", x); n = f"HEL{m.group(1)}.{m.group(2)}"; return "HEL362.1" if n == "HEL362.3" else n

N = pd.read_csv(normf, sep="\t", index_col=0); N.columns = [hel(c) for c in N.columns]
genes = pd.read_csv(f"{gd}/rankings/ols_signed.tsv", sep="\t").gene.values
cs = pd.read_csv(f"{gd}/gene_sets/custom_gene_sets.tsv", sep="\t"); hm = pd.read_csv(f"{gd}/gene_sets/hallmark.tsv", sep="\t")
sets = {}
for (c, s), g in cs.groupby(["collection", "set"]): sets[(c, s)] = set(g.gene)
for s, g in hm.groupby("set"): sets[("Hallmark", s)] = set(g.gene)
fg = pd.read_csv(f"{gd}/tables/fgsea_results_all.tsv", sep="\t")


def stat_matrix(Y, W, n1, n):
    """Y genes x n (log2), W n x P 0/1 indicator of 'PMOS' per split -> signed -log10 p of the per-gene OLS (two groups), genes x P."""
    n0 = n - n1; s1 = Y @ W; sT = Y.sum(axis=1, keepdims=True); s0 = sT - s1
    q1 = (Y ** 2) @ W; qT = (Y ** 2).sum(axis=1, keepdims=True); q0 = qT - q1
    m1, m0 = s1 / n1, s0 / n0; sse = (q1 - n1 * m1 ** 2) + (q0 - n0 * m0 ** 2)
    se = np.sqrt(np.maximum(sse, 1e-300) / (n - 2) * (1 / n1 + 1 / n0)); t = (m1 - m0) / se
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    return np.sign(t) * -np.log10(np.maximum(p, 1e-300))


def es_for(order_pos, absw, m_idx, Ng):
    """ES of one set for one ranking. order_pos: rank position (0 = highest) of every gene; absw: |stat| per gene; m_idx: gene indices of the set."""
    pos = order_pos[m_idx]; o = np.argsort(pos); pos = pos[o]; w = absw[m_idx][o]; m = len(pos); k = np.arange(1, m + 1)
    NR = w.sum()
    if NR == 0: return 0.0
    phit = np.cumsum(w) / NR; pm_after = (pos + 1 - k) / (Ng - m)
    after = phit - pm_after
    before = np.concatenate([[0.0], phit[:-1]]) - (pos - (k - 1)) / (Ng - m)
    up, dn = max(after.max(), 0.0), min(before.min(), 0.0)
    return up if up >= -dn else dn


out_rows = []
for tag, drop in [("all14", None), ("no350", "HEL350.3")]:
    lines = [l for l in meta.sort_values(["group", "bmi"]).index if l != drop]; n = len(lines)
    Y = np.log2(N.loc[genes, lines].values + 1); g = (meta.loc[lines, "group"] == "PMOS").values; n1 = int(g.sum())
    obsW = g.astype(float)[:, None]; st_obs = stat_matrix(Y, obsW, n1, n)[:, 0]
    combos = list(itertools.combinations(range(n), n1)); P = len(combos)
    W = np.zeros((n, P));
    for j, c in enumerate(combos): W[list(c), j] = 1
    gidx = {gname: i for i, gname in enumerate(genes)}; Ng = len(genes)
    set_idx = {k: np.array([gidx[x] for x in v if x in gidx]) for k, v in sets.items()}
    set_idx = {k: v for k, v in set_idx.items() if 5 <= len(v) <= 500}
    def all_es(st):
        order = np.argsort(-st, kind="stable"); pos = np.empty(Ng, int); pos[order] = np.arange(Ng); absw = np.abs(st)
        return np.array([es_for(pos, absw, v, Ng) for v in set_idx.values()])
    es_obs = all_es(st_obs); es_null = np.empty((P, len(set_idx)))
    for a0 in range(0, P, 400):
        S = stat_matrix(Y, W[:, a0:a0 + 400], n1, n)
        for j in range(S.shape[1]): es_null[a0 + j] = all_es(S[:, j])
    keys = list(set_idx.keys())
    for i, k in enumerate(keys):
        e = es_null[:, i]; pos_m = e[e > 0].mean() if (e > 0).any() else np.nan; neg_m = -e[e < 0].mean() if (e < 0).any() else np.nan
        nes_null = np.where(e > 0, e / pos_m, e / neg_m); nes_obs = es_obs[i] / pos_m if es_obs[i] > 0 else es_obs[i] / neg_m
        pval = np.sum(np.abs(nes_null) >= abs(nes_obs) - 1e-12) / P
        rk = "ols_signed" if tag == "all14" else "ols_signed_no350"
        f = fg[(fg.ranking == rk) & (fg.pathway == k[1]) & (fg.collection == k[0])]
        out_rows.append(dict(lines=tag, n_splits=P, collection=k[0], pathway=k[1], n_present=len(set_idx[k]), ES_obs=es_obs[i], ES_fgsea=(f.ES.iloc[0] if len(f) else np.nan), NES_labelperm=nes_obs,
                             p_labelperm=pval, NES_fgsea=(f.NES.iloc[0] if len(f) else np.nan), p_fgsea=(f.pval.iloc[0] if len(f) else np.nan)))
    print(tag, "done", P, "splits", flush=True)
R = pd.DataFrame(out_rows)
R["padj_labelperm_BH"] = R.groupby(["lines", "collection"]).p_labelperm.transform(lambda p: stats.false_discovery_control(p, method="bh"))
R["padj_fgsea_BH"] = R.groupby(["lines", "collection"]).p_fgsea.transform(lambda p: stats.false_discovery_control(p.fillna(1), method="bh"))
R.to_csv(f"{gd}/tables/fgsea_labelperm_calibration.tsv", sep="\t", index=False, float_format="%.5g")
chk = R.dropna(subset=["ES_fgsea"]); print(f"ES check vs fgsea: max |diff| = {(chk.ES_obs - chk.ES_fgsea).abs().max():.4f}, correlation {np.corrcoef(chk.ES_obs, chk.ES_fgsea)[0, 1]:.5f}")
for tag in ["all14", "no350"]:
    for coll in ["Hallmark"]:
        x = R[(R.lines == tag) & (R.collection == coll)]
        print(f"{tag} {coll}: fgsea padj<0.05 {int((x.padj_fgsea_BH < .05).sum())} | label-permutation padj<0.05 {int((x.padj_labelperm_BH < .05).sum())} (p<0.05: {int((x.p_labelperm < .05).sum())}); smallest label-permutation p {x.p_labelperm.min():.4f}")
    print(R[(R.lines == tag) & (R.collection != "Hallmark")][["collection", "pathway", "n_present", "NES_labelperm", "p_labelperm", "NES_fgsea", "p_fgsea"]].round(3).to_string(index=False))
print(R[(R.lines == "all14") & (R.collection == "Hallmark")].sort_values("p_labelperm").head(12)[["pathway", "NES_labelperm", "p_labelperm", "padj_labelperm_BH", "NES_fgsea", "p_fgsea"]].round(4).to_string(index=False))
