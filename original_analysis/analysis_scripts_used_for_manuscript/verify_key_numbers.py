#!/usr/bin/env python3
"""
Independent re-computation of the key numbers of the PMOS iPSC manuscript (v5.1, Stem Cell Research & Therapy version, 4 Nov 2026).
 
Written from scratch with numpy/scipy only; it does NOT import or call any of the manuscript scripts,
so a PASS is an independent confirmation. Runs in about 1-2 minutes on a laptop.
 
  Transcriptome  C01-C17   Genotypes / targeted eQTL  B01-B07   Protein (HipSci result table)  P01-P04
 
Usage:
  python3 verify_key_numbers.py --counts Normalized_counts.tabular --package <package folder>
Output: a PASS/FAIL table on screen and verification_report.tsv in the current folder.
Requirements: Python 3.9+, numpy, scipy, pandas.
"""
import argparse, itertools, re, sys
import numpy as np, pandas as pd
from scipy import stats
 
ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True); ap.add_argument("--package", required=True)
a = ap.parse_args()
R = a.package.rstrip("/") + "/4_Analysis_results"
rows = []
def check(cid, what, expected, observed, tol=0.0):
    if isinstance(expected, str) or isinstance(observed, str): ok = str(expected) == str(observed)
    else: ok = abs(float(observed) - float(expected)) <= tol
    rows.append(dict(id=cid, check=what, expected=expected, observed=observed, result="PASS" if ok else "FAIL"))
    print(f"{'PASS' if ok else 'FAIL'}  {cid:4s} {what:62s} expected {expected!s:>14}  observed {observed!s:>14}", flush=True)
 
# ---------------- data ----------------
cnt = pd.read_csv(a.counts, sep="\t", index_col=0)
def line(c):
    m = re.search(r"HEL\s*_?(\d+)[._](\d+)", c); n = f"HEL{m.group(1)}.{m.group(2)}"
    return {"HEL362.3": "HEL362.1"}.get(n, n)          # RNA-seq label HEL362_3 = clone HEL362.1
cnt.columns = [line(c) for c in cnt.columns]
key = pd.read_csv(f"{R}/metadata/sample_key_14_ipsc_donors.csv").set_index("hel_line")
L = list(key.index); cnt = cnt[L]; n = len(L)
grp = (key.group == "PMOS").astype(float).values
bmi, age = key.bmi.values.astype(float), key.age.values.astype(float)
xist = np.log2(cnt.loc["XIST"].values + 1)
expd = cnt[(cnt >= 10).sum(axis=1) >= 7]; Y = np.log2(expd.values + 1); G = len(Y)
 
def ols(X):
    """per-gene OLS of Y on [1, X]; returns coef, p (two-sided t), SSE"""
    X1 = np.column_stack([np.ones(n), X]); B, *_ = np.linalg.lstsq(X1, Y.T, rcond=None); B = B.T
    E = Y - B @ X1.T; df = n - X1.shape[1]; sse = (E ** 2).sum(1)
    se = np.sqrt(np.outer(sse / df, np.diag(np.linalg.inv(X1.T @ X1))))
    return B[:, 1:], 2 * stats.t.sf(np.abs(B / se), df)[:, 1:], sse
def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
 
if True:
    print("\n=== Transcriptome ===")
    check("C01", "expressed genes (>=10 counts in >=7 of 14 lines)", 16489, G)
    b, p, _ = ols(grp[:, None])
    check("C02", "OLS PMOS vs control: genes p<0.05", 443, int((p[:, 0] < .05).sum()))
    check("C03", "OLS PMOS vs control: genes BH-FDR<0.1", 0, int((bh(p[:, 0]) < .1).sum()))
    _, px, _ = ols(np.column_stack([grp, xist]))
    check("C04", "OLS PMOS + log2 XIST: genes p<0.05", 627, int((px[:, 0] < .05).sum()))
    # Spearman count + exact null over all 7-vs-7 splits
    Rk = np.apply_along_axis(stats.rankdata, 1, Y); Rk -= Rk.mean(1, keepdims=True); Rk /= np.sqrt((Rk ** 2).sum(1, keepdims=True))
    def nsig(Xc):
        Xr = np.apply_along_axis(stats.rankdata, 0, Xc); Xr -= Xr.mean(0); Xr /= np.sqrt((Xr ** 2).sum(0))
        r = np.clip(Rk @ Xr, -0.999999, 0.999999); t = r * np.sqrt((n - 2) / (1 - r ** 2))
        return (2 * stats.t.sf(np.abs(t), n - 2) < .05).sum(0)
    obs = int(nsig(grp[:, None])[0])
    splits = np.array([[1.0 if i in c else 0.0 for i in range(n)] for c in itertools.combinations(range(n), 7)]).T
    null = np.concatenate([nsig(splits[:, k:k + 500]) for k in range(0, splits.shape[1], 500)])
    check("C05", "Spearman PMOS: genes p<0.05", 505, obs)
    check("C06", "exact null (3,432 splits): median", 708, int(np.median(null)))
    check("C07", "exact permutation p (PMOS status)", 0.749, round((np.sum(null >= obs) + 1) / (len(null) + 1), 3), tol=0.002)
    # PCA
    v = np.argsort(Y.var(1))[::-1][:1000]; Z = (Y[v] - Y[v].mean(1, keepdims=True)).T
    s = np.linalg.svd(Z, compute_uv=False); ve = s ** 2 / (s ** 2).sum()
    check("C08", "PCA top-1,000 genes: PC1 % variance", 43.1, round(100 * ve[0], 1), tol=0.1)
    Za = (Y - Y.mean(1, keepdims=True)).T; sa = np.linalg.svd(Za, compute_uv=False)
    check("C09", "PCA all genes: PC1 % variance", 34.6, round(100 * sa[0] ** 2 / (sa ** 2).sum(), 1), tol=0.1)
    # sample correlation
    C = np.corrcoef(Y.T); iu = np.triu_indices(n, 1)
    check("C10", "sample Pearson r range (min-max)", "0.936-0.982", f"{C[iu].min():.3f}-{C[iu].max():.3f}")
    def wb(lab):
        same = lab[:, None] == lab[None, :]; return C[iu][same[iu]].mean() - C[iu][~same[iu]].mean()
    ero = (xist < 8).astype(float); o = wb(ero)
    nx = [wb(np.isin(np.arange(n), c).astype(float)) for c in itertools.combinations(range(n), 3)]
    check("C11", "XIST-eroded lines cluster: exact p (364 splits)", 0.033, round(np.mean(np.array(nx) >= o), 3), tol=0.001)
    check("C12", "XIST-eroded lines (log2 XIST < 8)", "HEL349.4,HEL361.2,HEL363.2", ",".join(sorted(np.array(L)[xist < 8])))
    # variance partition (Type II SS / total SS), observed medians
    SST = ((Y - Y.mean(1, keepdims=True)) ** 2).sum(1); Xd = np.column_stack([grp, bmi, age, xist]); _, _, sse_f = ols(Xd)
    med = [np.median((ols(np.delete(Xd, j, 1))[2] - sse_f) / SST) * 100 for j in range(4)]
    check("C13", "variance partition medians % (PMOS/BMI/age/XIST)", "3.2/5.2/7.5/6.2", "/".join(f"{m:.1f}" for m in med))
    kd = stats.spearmanr(xist, np.log2(cnt.loc["KDM5C"].values + 1))
    check("C14", "XIST vs KDM5C Spearman rho", -0.78, round(kd.statistic, 2), tol=0.01)
    ob = (bmi > 30).astype(float); _, pa, _ = ols(np.column_stack([grp, ob]))
    check("C15", "PMOS + BMI category: genes p<0.05", 465, int((pa[:, 0] < .05).sum()))
    check("C16", "  ... of which shared with unadjusted 443", 362, int(((pa[:, 0] < .05) & (p[:, 0] < .05)).sum()))
    gw = pd.read_csv(f"{R}/transcriptome/S7_PCOS_GWAS_genes_expression.tsv", sep="\t")
    check("C17", "PCOS consensus genes expressed (of 33)", 24, int(gw.gene.isin(expd.index).sum()))
 
if True:
    print("\n=== Genotypes and targeted cis-eQTL ===")
    dos = pd.read_csv(f"{R}/genotype_QC/donor_genotypes_6snps.tsv", sep="\t", index_col=0).loc[L]
    eq = pd.read_csv(f"{R}/genotype_QC/targeted_eqtl_6snps.tsv", sep="\t"); eq = eq[eq.status == "tested"]
    g = dos["rs10876864"]
    check("B01", "rs10876864 risk-allele counts 0/1/2 (13 called)", "3/9/1", "/".join(str(int((g == k).sum())) for k in (0, 1, 2)))
    check("B02", "rs10876864 missing call", "HEL362.1", ",".join(g.index[g.isna()]))
    res = []
    for _, r in eq.iterrows():
        d = dos[r.snp]; ok = d.notna().values
        y = np.log2(cnt.loc[r.gene].values[ok] + 1); X = np.column_stack([np.ones(ok.sum()), d.values[ok], xist[ok]])
        fit = np.linalg.lstsq(X, y, rcond=None)[0]; e = y - X @ fit; df = ok.sum() - 3
        se = np.sqrt((e @ e) / df * np.linalg.inv(X.T @ X)[1, 1]); res.append((r.snp, r.gene, fit[1], 2 * stats.t.sf(abs(fit[1] / se), df)))
    res = pd.DataFrame(res, columns=["snp", "gene", "beta", "p"]); res["q"] = bh(res.p)
    e1 = res[res.gene == "ERBB3"].iloc[0]; e2 = res[res.gene == "RPS26"].iloc[0]
    check("B03", "ERBB3 ~ rs10876864 + XIST: beta per risk allele", 0.309, round(e1.beta, 3), tol=0.002)
    check("B04", "ERBB3: p / q (13 pairs)", "0.024/0.31", f"{e1.p:.3f}/{e1.q:.2f}")
    check("B05", "RPS26 ~ rs10876864 + XIST: beta / p", "0.35/0.093", f"{e2.beta:.2f}/{e2.p:.3f}")
    d = dos["rs10876864"]; ok = d.notna().values; y = np.log2(cnt.loc["ERBB3"].values[ok] + 1)
    def dfit(cols):
        X = np.column_stack([np.ones(ok.sum())] + cols); f = np.linalg.lstsq(X, y, rcond=None)[0]; e = y - X @ f
        df = ok.sum() - X.shape[1]; se = np.sqrt((e @ e) / df * np.linalg.inv(X.T @ X)[1, 1]); return f[1], 2 * stats.t.sf(abs(f[1] / se), df)
    b6, p6 = dfit([d.values[ok], grp[ok]])
    check("B06", "ERBB3 ~ rs10876864 + diagnosis (no XIST): beta / p", "0.27/0.03", f"{b6:.2f}/{p6:.2f}")
    b6, p6 = dfit([d.values[ok], xist[ok], grp[ok]])
    check("B06b", "ERBB3 ~ rs10876864 + XIST + diagnosis: beta / p (v5.1)", "0.17/0.19", f"{b6:.2f}/{p6:.2f}")
    c = (grp[ok] == 0); rc = stats.spearmanr(d.values[ok][c], y[c])
    check("B06c", "controls only: Spearman rho / p (n = 6)", "0.62/0.19", f"{rc.statistic:.2f}/{rc.pvalue:.2f}")
    check("B07", "number of tested variant-gene pairs", 13, len(res))
 
if True:
    print("\n=== Protein cross-reference (from result table) ===")
    d = pd.read_csv(f"{R}/protein_crossref/S13d_protein_crossref_all_genes.tsv.gz", sep="\t"); pc = d[d.protein_coding]
    check("P01", "expressed protein-coding genes (HGNC)", 13241, len(pc))
    check("P02", "  ... with protein detected in HipSci", 8726, int(pc.protein_detected_HipSci.sum()))
    c = d[(d.p < .05) & d.protein_detected_HipSci & (d.hipsci_rna_protein_corr >= .3)]
    check("P03", "protein-supported candidates (p<0.05, r>=0.3)", 49, len(c))
    b = pd.read_csv(f"{R}/protein_crossref/S13b_protein_PCOS_locus_genes.tsv", sep="\t")
    e = b[b.source.str.startswith("consensus") & b.expressed_iPSC]
    check("P04", "expressed PCOS consensus genes with protein", 16, int(e.protein_detected_HipSci.sum()))
 
out = pd.DataFrame(rows); out.to_csv("verification_report.tsv", sep="\t", index=False)
print(f"\n{(out.result == 'PASS').sum()} PASS, {(out.result == 'FAIL').sum()} FAIL -> verification_report.tsv")
sys.exit(0 if (out.result == "PASS").all() else 1)
