#!/usr/bin/env python3
"""
STEP 8b  Targeted eQTL replication in the 14 iPSC lines
 
Only the short-listed SNP-gene pairs (from 6b, optionally filtered by coloc/AlphaGenome)
are tested, so the multiple-testing burden stays small.
 
Model per pair:  log2(norm count + 1) ~ risk-allele dosage (+ log2 XIST)
  - XIST is included by default because X-inactivation erosion varies strongly
    between the lines (exclude chrX genes from the pair list)
Reports per pair: beta (log2 change per risk allele), p, BH q,
  direction agreement with i2QTL (beta_risk) and with AlphaGenome (if given),
  and overall a sign test: do more pairs agree in direction than expected by chance?
 
Usage:
  python3 08b_targeted_eqtl.py --pairs results/06_i2qtl/i2qtl_best_per_locus.tsv \
     --dosages results/08_targeted/dosages.raw --counts rna/Normalized_counts.tabular \
     --key meta/sample_key.csv --out results/08_targeted
Pair file needs columns: snp, eGene, tag_allele (risk-tagging allele), beta_risk (i2QTL);
optional: alphagenome_dir (+1 / -1 predicted effect of the risk allele).
"""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import read_key, read_counts, ols, bh
 
ap = argparse.ArgumentParser()
ap.add_argument("--pairs", required=True); ap.add_argument("--dosages", required=True)
ap.add_argument("--counts", required=True); ap.add_argument("--key", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--covar", choices=["xist", "none"], default="xist")
ap.add_argument("--min-minor", type=int, default=3, help="min minor-allele copies among the lines")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
 
key = read_key(a.key)
id2line = dict(zip(key.array_iid, key.hel_line))
expr = np.log2(read_counts(a.counts) + 1)
 
raw = pd.read_csv(a.dosages, sep=r"\s+")
raw["line"] = raw.IID.astype(str).map(id2line)
raw = raw.dropna(subset=["line"]).set_index("line")
snpcols = {c.rsplit("_", 1)[0]: c for c in raw.columns[6:]}          # rs123_A -> counted A
 
lines = [l for l in expr.columns if l in raw.index]
missing = set(expr.columns) - set(lines)
if missing: print(f"WARNING: no genotypes for {sorted(missing)}", file=sys.stderr)
print(f"{len(lines)} lines with expression + genotypes")
cov = expr.loc["XIST", lines].values if a.covar == "xist" and "XIST" in expr.index else None
 
P = pd.read_csv(a.pairs, sep="\t")
out = []
for _, r in P.iterrows():
    base = dict(snp=r.snp, eGene=r.eGene, locus=r.get("locus", ""), risk_allele=r.tag_allele,
                i2qtl_beta_risk=r.get("beta_risk", np.nan), alphagenome_dir=r.get("alphagenome_dir", np.nan))
    if r.snp not in snpcols: out.append({**base, "status": "SNP not genotyped"}); continue
    if r.eGene not in expr.index: out.append({**base, "status": "gene not in counts"}); continue
    col = snpcols[r.snp]; counted = col.rsplit("_", 1)[1]
    d = raw.loc[lines, col].astype(float)
    risk = d if counted.upper() == str(r.tag_allele).upper() else 2 - d
    ok = risk.notna().values
    y = expr.loc[r.eGene, lines].values[ok]; x = risk.values[ok]
    hc = np.rint(x); minor = min(hc.sum(), 2 * len(hc) - hc.sum())
    base.update(n=int(ok.sum()), n_0=int((hc == 0).sum()), n_1=int((hc == 1).sum()), n_2=int((hc == 2).sum()),
                mean_log2expr=float(y.mean()))
    if minor < a.min_minor or np.var(x) == 0:
        out.append({**base, "status": "too few risk/non-risk alleles"}); continue
    X = x[:, None] if cov is None else np.column_stack([x, cov[ok]])
    b, se, t, p = ols(y, X)
    out.append({**base, "status": "tested", "beta_log2_per_risk_allele": b[0], "se": se[0], "p": p[0]})
 
R = pd.DataFrame(out)
t = R.status == "tested"
if t.any():
    R.loc[t, "q_BH"] = bh(R.loc[t, "p"])
    R.loc[t, "agrees_i2qtl"] = np.sign(R.loc[t, "beta_log2_per_risk_allele"]) == np.sign(R.loc[t, "i2qtl_beta_risk"])
    if R.alphagenome_dir.notna().any():
        R.loc[t, "agrees_alphagenome"] = np.sign(R.loc[t, "beta_log2_per_risk_allele"]) == np.sign(R.loc[t, "alphagenome_dir"])
R.sort_values(["status", "p"]).to_csv(f"{a.out}/targeted_eqtl_results.tsv", sep="\t", index=False)
 
# ---- overall direction (sign) test
summ = []
for c in ["agrees_i2qtl", "agrees_alphagenome"]:
    if c in R and R.loc[t, c].notna().any():
        v = R.loc[t, c].dropna().astype(bool)
        bt = stats.binomtest(int(v.sum()), len(v), 0.5, alternative="greater")
        summ.append(dict(comparison=c, pairs=len(v), agree=int(v.sum()), fraction=v.mean(), p_sign_test=bt.pvalue))
pd.DataFrame(summ).to_csv(f"{a.out}/direction_agreement.tsv", sep="\t", index=False)
print(R.status.value_counts().to_string()); print(pd.DataFrame(summ).to_string(index=False))
 
# ---- figure: expression by risk-allele count for tested pairs (max 24)
T = R[t].sort_values("p").head(24)
if len(T):
    n = len(T); nc = min(4, n); nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(3 * nc, 2.6 * nr), squeeze=False)
    for ax, (_, r) in zip(axs.flat, T.iterrows()):
        col = snpcols[r.snp]; counted = col.rsplit("_", 1)[1]
        d = raw.loc[lines, col].astype(float)
        x = d if counted.upper() == str(r.risk_allele).upper() else 2 - d
        y = expr.loc[r.eGene, lines]
        ax.scatter(x + np.random.uniform(-.08, .08, len(x)), y, s=22, color="#2471a3", alpha=.8)
        for g in (0, 1, 2):
            m = y[np.rint(x) == g]
            if len(m): ax.hlines(m.median(), g - .25, g + .25, color="#c0392b", lw=2)
        ax.set_xticks([0, 1, 2]); ax.set_xlabel(f"{r.snp} risk alleles", fontsize=8)
        ax.set_title(f"{r.eGene}  p={r.p:.2g}", fontsize=9); ax.set_ylabel("log2 expr", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axs.flat[n:]: ax.axis("off")
    fig.tight_layout(); fig.savefig(f"{a.out}/targeted_eqtl_plots.png", dpi=200)
