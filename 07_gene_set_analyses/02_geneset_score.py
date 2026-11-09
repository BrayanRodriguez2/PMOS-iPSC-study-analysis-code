#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/02_geneset_score.py   (original step script: 72_geneset_score.py)
# Produces: Replication of the GWAS-gene-set score test with expression-matched random sets.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 72 - independent replication of Fig. 5B / Table S7b of the final manuscript package (PMOS_iPSC_SCRT_061126): score of the PCOS GWAS consensus genes
(33 genes at the 29 loci of Moolhuijsen et al. 2026; Table S7a) compared with random gene sets matched on expression. Own code, the method of Additional file 3 (an3.py):
  * expressed genes: >= 10 normalised counts in >= 7 of 14 lines; y = log2(norm + 1); the consensus genes that are expressed form the gene set (24 in the original data);
  * z-score of every gene across the 14 lines (SD with n - 1); score of a line = mean z of the gene set (Fig. 5B);
  * null: 10,000 random gene sets with the same number of genes in each of 20 expression bins (quantile bins of the mean log2 expression of ALL expressed genes; sampled
    without replacement, the set's own genes are not excluded), numpy default_rng(7), drawn bin by bin in the order of pandas value_counts (as the original);
  * test: Spearman rho between the line score and PMOS status (0/1), AMH, FAI, testosterone, BMI and HOMA-IR; two-sided permutation p = (n(|null rho| >= |rho|) + 1) / (10,000 + 1).
Option --exclude-genes: genes removed from the set (sensitivity: e.g. AMH, expressed in our counts but not in the original).
Extra (descriptive): Mann-Whitney of the line scores PMOS vs control; which consensus genes are expressed in the input.
Inputs: --norm (normalised counts, gene symbols x 14 lines), --meta (donor_meta.tsv of step 65a, corrected clinical file), --s7a / --s7b (the original tables as TSV).
Outputs (--out/<tag>): geneset_score_per_line.tsv, geneset_tests.tsv, consensus_genes_expression.tsv, report.txt
"""
import argparse, os, re
import numpy as np, pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--norm", required=True); ap.add_argument("--meta", required=True); ap.add_argument("--s7a", required=True); ap.add_argument("--s7b", required=True)
ap.add_argument("--tag", required=True); ap.add_argument("--exclude-genes", default=""); ap.add_argument("--out", required=True); ap.add_argument("--n-null", type=int, default=10000)
a = ap.parse_args()
OUT = os.path.join(a.out, a.tag); os.makedirs(OUT, exist_ok=True)
rep = []
log = lambda m="": (print(m, flush=True), rep.append(str(m)))


def hel(x):
    m = re.search(r"HEL\s*_?(\d+)[._](\d+)", x); n = f"HEL{m.group(1)}.{m.group(2)}"; return "HEL362.1" if n == "HEL362.3" else n


meta = pd.read_csv(a.meta, sep="\t", index_col=0)
N = pd.read_csv(a.norm, sep="\t", index_col=0); N.columns = [hel(c) for c in N.columns]
lines = list(meta.sort_values(["group", "bmi"]).index); N = N[lines]
expr = N[(N >= 10).sum(axis=1) >= 7]; lg = np.log2(expr + 1)
s7a = pd.read_csv(a.s7a, sep="\t"); cons = list(s7a.gene)
excl = [x for x in a.exclude_genes.split(",") if x]
gs = [g for g in cons if g in lg.index and g not in excl]
log(f"[{a.tag}] expressed genes {len(lg):,}; consensus genes {len(cons)} (S7a); expressed here: {len(gs)} (the original S7a: {int(s7a.expressed.sum())})")
exp_theirs = set(s7a.gene[s7a.expressed == 1])
log(f"    expressed in this input but not in the original S7a: {sorted(set(gs) - exp_theirs)}; expressed in the original S7a but not here: {sorted(exp_theirs - set(gs))}")
cg = pd.DataFrame({"gene": cons, "mean_norm_count": [N.loc[g].mean() if g in N.index else np.nan for g in cons], "expressed_here": [g in lg.index for g in cons],
                   "expressed_S7a": [g in exp_theirs for g in cons], "S7a_mean_norm_count": s7a.mean_norm_count.values})
cg.to_csv(os.path.join(OUT, "consensus_genes_expression.tsv"), sep="\t", index=False, float_format="%.4g")

Z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1), axis=0)
score = Z.loc[gs].mean()
bins = pd.qcut(lg.mean(axis=1), 20, labels=False); need = bins.loc[gs].value_counts()
pos = {b: np.where(bins.values == b)[0] for b in need.index}
rng = np.random.default_rng(7); Zv = Z.values; NS = a.n_null
null_scores = np.empty((NS, len(lines)))
for i in range(NS):
    pick = np.concatenate([rng.choice(pos[b], n, replace=False) for b, n in need.items()])
    null_scores[i] = Zv[pick].mean(axis=0)
rk = lambda M: np.apply_along_axis(stats.rankdata, -1, M)
NR = rk(null_scores); NR = NR - NR.mean(axis=1, keepdims=True); NR /= np.sqrt((NR ** 2).sum(axis=1, keepdims=True))
traits = {"PMOS status": (meta.loc[lines, "group"] == "PMOS").astype(float).values, "AMH": meta.loc[lines, "amh"].values, "FAI": meta.loc[lines, "fai"].values,
          "Testosterone": meta.loc[lines, "testosterone"].values, "BMI": meta.loc[lines, "bmi"].values, "HOMA-IR": meta.loc[lines, "homa_ir"].values}
rows = []
s7b = pd.read_csv(a.s7b, sep="\t").set_index("test")
for k, x in traits.items():
    xr = stats.rankdata(x); xr = xr - xr.mean(); xr /= np.sqrt((xr ** 2).sum())
    nl = NR @ xr; r = stats.spearmanr(x, score.values)[0]
    p = (np.sum(np.abs(nl) >= abs(r) - 1e-12) + 1) / (NS + 1)
    h = s7b.loc[f"PCOS GWAS gene-set score vs {k}"]
    rows.append(dict(trait=k, n_genes=len(gs), rho=r, p_perm_matched=p, null_rho_SD=float(nl.std()), her_rho=h.rho, her_p=h.p_perm_matched))
T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, "geneset_tests.tsv"), sep="\t", index=False, float_format="%.4g")
log(T.round(3).to_string(index=False))
grp = traits["PMOS status"]; mwp = stats.mannwhitneyu(score.values[grp == 1], score.values[grp == 0], alternative="two-sided").pvalue
sc = pd.DataFrame({"line": lines, "group": meta.loc[lines, "group"].values, "mean_z_score": score.values}); sc.to_csv(os.path.join(OUT, "geneset_score_per_line.tsv"), sep="\t", index=False, float_format="%.4f")
log(f"    mean score PMOS {score.values[grp == 1].mean():+.3f} vs CTRL {score.values[grp == 0].mean():+.3f}; Mann-Whitney p = {mwp:.3f}")
open(os.path.join(OUT, "report.txt"), "w").write("\n".join(rep) + "\n")
