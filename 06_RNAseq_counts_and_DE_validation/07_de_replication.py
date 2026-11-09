#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/07_de_replication.py   (original step script: 65_deg_analysis.py)
# Produces: Independent replication of the expression analyses (OLS, PCA, permutations, variance partition, chr20).
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 65 - INDEPENDENT re-implementation (own code, numpy/scipy) of the expression analyses of the manuscript, with the parameters stated in Additional_file_3 and the
figure book. The same code is run on three inputs (--tag): our counts (stranded, the correct orientation), our counts counted unstranded, and the original Normalized_counts.tabular.
Input : normalised counts, rows = gene symbols, columns = 14 lines (hel names); donor_meta.tsv (65a); S4 and S15 of the original Additional_file_1; gene positions (GENCODE v50).
Rules (the original): expressed = >= 10 normalised counts in >= 7 of 14 lines; y = log2(norm + 1); per gene OLS y ~ PMOS (0/1) and y ~ PMOS + log2(XIST + 1), BH over all genes;
 PCA = SVD of centred (not scaled) y on the 1,000 most variable genes and on all expressed genes; Pearson correlation between lines and exact tests of (mean within-group r - mean
 between-group r) over all splits (3,432 for diagnosis and BMI category, C(14, k) for XIST erosion, eroded = log2(XIST + 1) < 8);
 'number of genes with Spearman p < 0.05' (rank-standardised y, t approximation) against the exact null of all C(14,7) splits (PMOS) or 2,000 permutations of each trait (seed 42,
 same call order as the original); variance partition PMOS + BMI + age + log2 XIST (Type II SS / total SS, median over genes) against 1,000 joint permutations of the covariates (seed 2026);
 BMI-category model y ~ PMOS + obese (> 30); sensitivity (Table S8): PMOS test without HEL350.3 + HEL351.3 (792 exact splits) or without the XIST-eroded lines (330 splits).
Extra: OLS without HEL350.3; chromosome-20 expression of HEL350.3 (per-gene z across the 14 lines, median per arm).
Outputs (--out/<tag>): numbers.tsv (every number, key/value), de_ols.tsv.gz, pcs*.tsv, line_correlation.tsv, vp_summary.tsv, nsig_*.tsv, s8.tsv, chr20_*.tsv, report.txt
"""
import argparse, itertools, os, re
import numpy as np, pandas as pd
from scipy import stats
from scipy.cluster import hierarchy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--norm", required=True); ap.add_argument("--tag", required=True); ap.add_argument("--meta", required=True); ap.add_argument("--s4", required=True)
ap.add_argument("--s15", required=True); ap.add_argument("--gene-pos", required=True); ap.add_argument("--out", required=True); ap.add_argument("--figs", action="store_true")
a = ap.parse_args()
OUT = os.path.join(a.out, a.tag); os.makedirs(OUT, exist_ok=True); FIG = os.path.join(OUT, "figures"); os.makedirs(FIG, exist_ok=True)
NUM, rep = [], []


def num(key, val, note=""):
    NUM.append(dict(key=key, value=val, note=note))


def log(m=""):
    print(m, flush=True); rep.append(str(m))


# ------------------------------------------------------------------ helpers (own implementation)
def hel(x):
    m = re.search(r"HEL\s*_?(\d+)[._](\d+)", x); n = f"HEL{m.group(1)}.{m.group(2)}"; return "HEL362.1" if n == "HEL362.3" else n


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(q, 1); return out


def ols(Y, X):
    """Y genes x n, X n x k (no intercept column). Returns beta, p (genes x k), sse (genes) for the k covariates."""
    n = Y.shape[0] if Y.ndim == 1 else Y.shape[1]
    X1 = np.column_stack([np.ones(n), X]); XtXi = np.linalg.inv(X1.T @ X1)
    B = Y @ X1 @ XtXi.T; E = Y - B @ X1.T; df = n - X1.shape[1]; sse = (E ** 2).sum(axis=1)
    se = np.sqrt((sse / df)[:, None] * np.diag(XtXi)[None, :]); t = B / se
    return B[:, 1:], 2 * stats.t.sf(np.abs(t[:, 1:]), df), sse


def sse_only(Y, X):
    n = Y.shape[1]; X1 = np.column_stack([np.ones(n), X]); B = Y @ X1 @ np.linalg.pinv(X1.T @ X1).T; return ((Y - B @ X1.T) ** 2).sum(axis=1)


def rstd(v):
    r = stats.rankdata(v); r = r - r.mean(); return r / np.sqrt((r ** 2).sum())


def rcrit(n):
    t = stats.t.isf(0.025, n - 2); return t / np.sqrt(n - 2 + t ** 2)


# ------------------------------------------------------------------ data
meta = pd.read_csv(a.meta, sep="\t", index_col=0)
N = pd.read_csv(a.norm, sep="\t", index_col=0)
N.columns = [hel(c) for c in N.columns]
lines = list(meta.sort_values(["group", "bmi"]).index); N = N[lines]
if N.index.duplicated().any():
    N = N.groupby(level=0).sum()
grp = (meta.loc[lines, "group"] == "PMOS").astype(float).values
bmi, age = meta.loc[lines, "bmi"].values, meta.loc[lines, "age"].values
obese = (bmi > 30).astype(float); n = 14
expr = N[(N >= 10).sum(axis=1) >= 7]; L = np.log2(expr + 1); Y = L.values; G = len(L)
xist = np.log2(N.loc["XIST", lines].values + 1); eroded = (xist < 8)
log(f"[{a.tag}] genes in the matrix {len(N):,} | expressed (>= 10 in >= 7 of 14 lines): {G:,} | XIST log2 per line: " + ", ".join(f"{l} {x:.2f}" for l, x in zip(lines, xist)))
log(f"    XIST-eroded (log2 < 8): {[l for l, e in zip(lines, eroded) if e]}")
num("n_expressed", G, "expected 16,489"); num("xist_eroded_lines", ";".join(l for l, e in zip(lines, eroded) if e), "expected HEL361.2, HEL349.4, HEL363.2")
s4 = pd.read_csv(a.s4, sep="\t").set_index("gene")
inter = L.index.intersection(s4.index)
num("n_expressed_in_her_S4_list", len(inter), f"of the original {len(s4):,}"); num("expressed_only_ours", G - len(inter)); num("expressed_only_hers", len(s4) - len(inter))
log(f"    overlap of the expressed-gene list with the original Table S4: {len(inter):,} (ours only {G - len(inter):,}, the original only {len(s4) - len(inter):,})")

# ------------------------------------------------------------------ Step 2: gene-wise OLS
b0, p0, _ = ols(Y, grp[:, None]); b1, p1, _ = ols(Y, np.column_stack([grp, xist]))
q0, q1 = bh(p0[:, 0]), bh(p1[:, 0])
de = pd.DataFrame({"gene": L.index, "log2FC": b0[:, 0], "p": p0[:, 0], "q_BH": q0, "log2FC_adjXIST": b1[:, 0], "p_adjXIST": p1[:, 0], "q_BH_adjXIST": q1, "mean_norm_count": expr.mean(axis=1).values})
de.to_csv(os.path.join(OUT, "de_ols.tsv.gz"), sep="\t", index=False, float_format="%.6g")
num("ols_PMOS_p<0.05", int((p0[:, 0] < .05).sum()), "expected 443"); num("ols_PMOS_FDR<0.1", int((q0 < .1).sum()), "expected 0"); num("ols_PMOS_min_FDR", float(q0.min()))
num("ols_PMOS+XIST_p<0.05", int((p1[:, 0] < .05).sum()), "expected 627"); num("ols_PMOS+XIST_FDR<0.1", int((q1 < .1).sum()))
log(f"    OLS ~ PMOS: p < 0.05 {int((p0[:, 0] < .05).sum())} (expected 443), FDR < 0.1 {int((q0 < .1).sum())}, min FDR {q0.min():.3f}; ~ PMOS + log2 XIST: p < 0.05 {int((p1[:, 0] < .05).sum())} (expected 627), FDR < 0.1 {int((q1 < .1).sum())}")
d = de.set_index("gene").loc[inter]; t4 = s4.loc[inter]
for lab, cf, cp, tf, tp in [("PMOS", "log2FC", "p", "log2FC_PMOS_vs_control", "p"), ("PMOS+XIST", "log2FC_adjXIST", "p_adjXIST", "log2FC_PMOS_vs_control_adjXIST", "p_adjXIST")]:
    r_fc = np.corrcoef(d[cf], t4[tf])[0, 1]; r_p = np.corrcoef(-np.log10(d[cp]), -np.log10(t4[tp]))[0, 1]
    nom_o = set(d.index[d[cp] < .05]); nom_t = set(t4.index[t4[tp] < .05])
    num(f"vsS4_{lab}_pearson_log2FC", r_fc); num(f"vsS4_{lab}_pearson_-log10p", r_p); num(f"vsS4_{lab}_spearman_log2FC", stats.spearmanr(d[cf], t4[tf])[0])
    num(f"vsS4_{lab}_nominal_ours", len(nom_o)); num(f"vsS4_{lab}_nominal_hers", len(nom_t)); num(f"vsS4_{lab}_nominal_overlap", len(nom_o & nom_t))
    log(f"    vs the original S4 ({lab}, {len(inter):,} common genes): r(log2FC) {r_fc:.3f}, r(-log10 p) {r_p:.3f}; nominal genes ours {len(nom_o)}, the original {len(nom_t)}, overlap {len(nom_o & nom_t)}")
s15 = pd.read_csv(a.s15, sep="\t", skiprows=1) if False else None
s15_raw = [l.rstrip("\n").split("\t") for l in open(a.s15)]; hdr = next(i for i, r in enumerate(s15_raw) if r and r[0] == "gene"); top40_theirs = [r[0] for r in s15_raw[hdr + 1:] if r and r[0]][:40]
top40_ours = list(de.sort_values("p").gene.head(40)); ov = len(set(top40_ours) & set(top40_theirs))
num("top40_overlap_with_S15", ov); log(f"    top-40 genes by p: overlap with the original Table S15 = {ov}/40")
de.sort_values("p").head(40).to_csv(os.path.join(OUT, "top40_ours.tsv"), sep="\t", index=False, float_format="%.5g")

# ------------------------------------------------------------------ Step 3: PCA
def pca(Z):
    Zc = Z - Z.mean(axis=1, keepdims=True); U, S, _ = np.linalg.svd(Zc.T, full_matrices=False); return U * S, S ** 2 / (S ** 2).sum()
top = L.var(axis=1).sort_values(ascending=False).index[:1000]
sc1, ve1 = pca(L.loc[top].values); sca, vea = pca(Y)
num("pca_top1000_PC1_%", 100 * ve1[0], "expected 43.1"); num("pca_top1000_PC2_%", 100 * ve1[1], "expected 8.5"); num("pca_all_PC1_%", 100 * vea[0], "expected 34.6"); num("pca_all_PC2_%", 100 * vea[1], "expected 9.9")
log(f"    PCA top-1000: PC1 {100 * ve1[0]:.1f} %, PC2 {100 * ve1[1]:.1f} % (expected 43.1 / 8.5); all genes: {100 * vea[0]:.1f} / {100 * vea[1]:.1f} (expected 34.6 / 9.9)")
pd.DataFrame(sc1[:, :5], index=lines, columns=[f"PC{i}" for i in range(1, 6)]).assign(group=meta.loc[lines, "group"].values).to_csv(os.path.join(OUT, "pcs_top1000.tsv"), sep="\t", float_format="%.4f")
traits = {"AMH": meta.loc[lines, "amh"].values, "Testosterone": meta.loc[lines, "testosterone"].values, "FAI": meta.loc[lines, "fai"].values, "FG score": meta.loc[lines, "fg_score"].values,
          "BMI": bmi, "HOMA-IR": meta.loc[lines, "homa_ir"].values, "Age": age, "XIST (log2)": xist, "PMOS status": grp}
rho = pd.DataFrame({k: [stats.spearmanr(v, sc1[:, i])[0] for i in range(5)] for k, v in traits.items()}, index=[f"PC{i}" for i in range(1, 6)])
pv = pd.DataFrame({k: [stats.spearmanr(v, sc1[:, i])[1] for i in range(5)] for k, v in traits.items()}, index=rho.index)
rho.to_csv(os.path.join(OUT, "pc_trait_rho.tsv"), sep="\t", float_format="%.3f"); pv.to_csv(os.path.join(OUT, "pc_trait_p.tsv"), sep="\t", float_format="%.4f")
sig = [(i, k, rho.loc[i, k], pv.loc[i, k]) for i in rho.index for k in rho.columns if pv.loc[i, k] < .05 and k != "PMOS status"]
num("pc_trait_significant_pairs", "; ".join(f"{i}-{k} rho {r:.2f} p {p:.3f}" for i, k, r, p in sig), "expected only PC1-age rho -0.54 p 0.044")
num("PC1_age_rho", rho.loc["PC1", "Age"]); num("PC1_age_p", pv.loc["PC1", "Age"])
log(f"    PC1-5 vs traits (Spearman, p < 0.05, excluding PMOS): {[(i, k, round(r, 2), round(p, 3)) for i, k, r, p in sig]} (expected only PC1-age -0.54, p 0.044)")

# ------------------------------------------------------------------ Step 3: sample correlation and exact tests
C = np.corrcoef(Y.T); iu = np.triu_indices(n, 1)
pd.DataFrame(C, index=lines, columns=lines).to_csv(os.path.join(OUT, "line_correlation.tsv"), sep="\t", float_format="%.4f")
num("line_corr_min", C[iu].min(), "expected 0.936"); num("line_corr_max", C[iu].max(), "expected 0.982")


def blockstat(lab):
    same = lab[:, None] == lab[None, :]; return C[iu][same[iu]].mean() - C[iu][~same[iu]].mean()


def exact_p(lab, k):
    obs = blockstat(lab); cnt = tot = 0
    for c in itertools.combinations(range(n), k):
        z = np.zeros(n); z[list(c)] = 1; tot += 1; cnt += blockstat(z) >= obs - 1e-15
    return obs, cnt / tot, tot


for nm, lab, k, exp in [("diagnosis", grp, 7, "0.79"), ("BMI_category", obese, 7, "0.08"), ("XIST_erosion", eroded.astype(float), int(eroded.sum()), "0.033")]:
    obs, p, tot = exact_p(lab, k); num(f"corr_exact_p_{nm}", p, f"expected {exp}; splits {tot}; stat {obs:.5f}")
    log(f"    Pearson between lines: within - between r for {nm}: {obs:.5f}, exact p = {p:.4f} over {tot} splits (expected {exp})")

# ------------------------------------------------------------------ Step 3: genes with Spearman p < 0.05 vs permutation null (Fig 3D, Table S8)
R = np.apply_along_axis(stats.rankdata, 1, Y); R = R - R.mean(axis=1, keepdims=True); R /= np.sqrt((R ** 2).sum(axis=1, keepdims=True))


def nsig_cols(Rm, X):
    """X: m x n (rank-standardised trait vectors); returns number of genes with p < 0.05 for each row."""
    rc = rcrit(Rm.shape[1]); return (np.abs(Rm @ X.T) > rc).sum(axis=0)


def exact_split_null(Rm, k):
    nn = Rm.shape[1]; rows = []
    for c in itertools.combinations(range(nn), k):
        z = np.zeros(nn); z[list(c)] = 1; rows.append((z - z.mean()) / np.sqrt(((z - z.mean()) ** 2).sum()))
    return nsig_cols(Rm, np.array(rows))


obs_pmos = int(nsig_cols(R, rstd(grp)[None, :])[0]); null_p = exact_split_null(R, 7)
pp = (np.sum(null_p >= obs_pmos) + 1) / (len(null_p) + 1)
num("nsig_PMOS_obs", obs_pmos, "expected 505"); num("nsig_PMOS_null_median", float(np.median(null_p)), "expected 708"); num("nsig_PMOS_perm_p", pp, "expected 0.75"); num("nsig_PMOS_n_splits", len(null_p))
log(f"    genes with Spearman p < 0.05 for PMOS: {obs_pmos} vs exact null median {np.median(null_p):.0f} (3,432 splits), p = {pp:.3f} (expected 505 vs 708, p 0.75)")
rng = np.random.default_rng(42); nrows = [dict(trait="PMOS status", genes_p05=obs_pmos, null_median=float(np.median(null_p)), null_95pct=float(np.percentile(null_p, 95)), perm_p=pp, n_perm=len(null_p))]
for k_, x in traits.items():
    if k_ in ("PMOS status", "XIST (log2)"): continue
    obs = int(nsig_cols(R, rstd(x)[None, :])[0])
    nl = nsig_cols(R, np.array([rstd(rng.permutation(x)) for _ in range(2000)]))
    nrows.append(dict(trait=k_, genes_p05=obs, null_median=float(np.median(nl)), null_95pct=float(np.percentile(nl, 95)), perm_p=(np.sum(nl >= obs) + 1) / (len(nl) + 1), n_perm=len(nl)))
nsd = pd.DataFrame(nrows); nsd.to_csv(os.path.join(OUT, "nsig_permutation.tsv"), sep="\t", index=False, float_format="%.4g")
for r in nsd.itertuples():
    num(f"nsig_{r.trait}_obs", r.genes_p05); num(f"nsig_{r.trait}_perm_p", r.perm_p)
log("    trait associations (genes with Spearman p < 0.05, permutation p):\n" + nsd.round(3).to_string(index=False))

# ------------------------------------------------------------------ Step 3: variance partition (Type II SS) with 1,000 joint permutations
SST = ((Y - Y.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
Xd = np.column_stack([grp, bmi, age, xist])


def vp(X, Yv=Y, sst=SST):
    full = sse_only(Yv, X); out = [(sse_only(Yv, np.delete(X, j, 1)) - full) / sst for j in range(X.shape[1])]; out.append(full / sst); return out


V = vp(Xd); rng = np.random.default_rng(2026); nm_ = np.zeros((1000, 4))
for k in range(1000):
    Vk = vp(Xd[rng.permutation(n)]); nm_[k] = [np.median(Vk[j]) for j in range(4)]
rows = []
for j, f in enumerate(["PMOS status", "BMI", "Age", "XIST (log2)"]):
    m = np.median(V[j]); rows.append(dict(factor=f, median_pct=100 * m, null_median_pct=100 * np.median(nm_[:, j]), perm_p=(np.sum(nm_[:, j] >= m) + 1) / 1001))
rows.append(dict(factor="Residual", median_pct=100 * np.median(V[4])))
vps = pd.DataFrame(rows); vps.to_csv(os.path.join(OUT, "vp_summary.tsv"), sep="\t", index=False, float_format="%.4g")
for r in vps.itertuples():
    num(f"vp_median_{r.factor}", r.median_pct, "expected PMOS 3.2, BMI 5.2, age 7.5, XIST 6.2, residual 64.1");
    if r.factor != "Residual": num(f"vp_perm_p_{r.factor}", r.perm_p, "expected 0.72, 0.22, 0.04, 0.11")
log("    variance partition (median % of total SS; permutation p):\n" + vps.round(3).to_string(index=False))
np.save(os.path.join(OUT, "vp_nullmed.npy"), nm_); pd.DataFrame({f: V[j] for j, f in enumerate(["PMOS", "BMI", "Age", "XIST", "Residual"])}, index=L.index).to_csv(os.path.join(OUT, "vp_per_gene.tsv.gz"), sep="\t", float_format="%.5f")

# ------------------------------------------------------------------ Step 3: BMI-category model
bA, pA, _ = ols(Y, np.column_stack([grp, obese])); u, v_ = p0[:, 0] < .05, pA[:, 0] < .05
num("bmi_cat_PMOS_p<0.05", int(v_.sum()), "expected 465"); num("bmi_cat_overlap_with_main", int((u & v_).sum()), "expected 362 shared with 443"); num("bmi_cat_FDR<0.1", int((bh(pA[:, 0]) < .1).sum()), "expected 0")
log(f"    PMOS + BMI category: p < 0.05 {int(v_.sum())} (expected 465), shared with the main model {int((u & v_).sum())} (expected 362 of 443), FDR < 0.1 {int((bh(pA[:, 0]) < .1).sum())}")

# ------------------------------------------------------------------ Step 3: XIST and KDM5C
if "KDM5C" in N.index:
    rr, pr_ = stats.spearmanr(N.loc["XIST", lines].values, N.loc["KDM5C", lines].values); num("XIST_KDM5C_rho", rr, "expected -0.78"); num("XIST_KDM5C_p", pr_, "expected 0.001")
    log(f"    XIST vs KDM5C: Spearman rho {rr:.3f}, p {pr_:.4f} (expected -0.78, 0.001)")
er_ctrl = int(eroded[grp == 0].sum()); er_pm = int(eroded[grp == 1].sum())
fp = stats.fisher_exact([[er_ctrl, 7 - er_ctrl], [er_pm, 7 - er_pm]])[1]; num("xist_erosion_fisher_p", fp, f"eroded: {er_ctrl} CTRL / {er_pm} PMOS; expected 0.19")
log(f"    XIST-eroded lines: {er_ctrl} control / {er_pm} PMOS; Fisher exact p = {fp:.3f} (expected 0.19 for 3/0)")

# ------------------------------------------------------------------ Step 4: Table S8 sensitivity (exact splits)
s8 = []


def subset_test(keep, label):
    idx = [i for i, l in enumerate(lines) if l in keep]; Rs = np.apply_along_axis(stats.rankdata, 1, Y[:, idx]); Rs = Rs - Rs.mean(axis=1, keepdims=True); Rs /= np.sqrt((Rs ** 2).sum(axis=1, keepdims=True))
    g = grp[idx]; k = int(g.sum()); nn = len(idx); xs = (g - g.mean()) / np.sqrt(((g - g.mean()) ** 2).sum()); rc = rcrit(nn)
    obs = int((np.abs(Rs @ xs) > rc).sum()); nl = exact_split_null_n(Rs, k); p_ = (np.sum(nl >= obs) + 1) / (len(nl) + 1)
    _, pols, _ = ols(Y[:, idx], g[:, None]); fdr = int((bh(pols[:, 0]) < .1).sum())
    return dict(analysis=label, n_control=nn - k, n_PMOS=k, genes_p05=obs, null_median=float(np.median(nl)), perm_p=p_, n_splits=len(nl), genes_FDR10=fdr, ols_genes_p05=int((pols[:, 0] < .05).sum()))


def exact_split_null_n(Rm, k):
    nn = Rm.shape[1]; rows = []
    for c in itertools.combinations(range(nn), k):
        z = np.zeros(nn); z[list(c)] = 1; rows.append((z - z.mean()) / np.sqrt(((z - z.mean()) ** 2).sum()))
    rc = rcrit(nn); return (np.abs(Rm @ np.array(rows).T) > rc).sum(axis=0)


s8.append(subset_test(set(lines), "all 14"))
s8.append(subset_test(set(lines) - {"HEL350.3", "HEL351.3"}, "excl. HEL350.3, HEL351.3"))
er_lines = {l for l, e in zip(lines, eroded) if e}
s8.append(subset_test(set(lines) - er_lines, f"excl. XIST-eroded lines ({len(er_lines)})"))
s8.append(subset_test(set(lines) - {"HEL361.2", "HEL349.4", "HEL363.2"}, "excl. HEL361.2, HEL349.4, HEL363.2 (the original three)"))
s8d = pd.DataFrame(s8); s8d.to_csv(os.path.join(OUT, "s8.tsv"), sep="\t", index=False, float_format="%.4g")
for r in s8d.itertuples():
    num(f"S8_{r.analysis}_genes_p05", r.genes_p05); num(f"S8_{r.analysis}_null_median", r.null_median); num(f"S8_{r.analysis}_perm_p", r.perm_p); num(f"S8_{r.analysis}_n_splits", r.n_splits)
log("    Table S8 (expected: all 14 505/708/0.749; excl. 350.3 + 351.3 491/625/0.663 (792); excl. 3 eroded 642/580/0.447 (330)):\n" + s8d.round(3).to_string(index=False))

# ------------------------------------------------------------------ Step 5: OLS without HEL350.3 and chr20 of HEL350.3
keep = [i for i, l in enumerate(lines) if l != "HEL350.3"]; b_ns, p_ns, _ = ols(Y[:, keep], grp[keep][:, None]); b_ns1, p_ns1, _ = ols(Y[:, keep], np.column_stack([grp[keep], xist[keep]]))
num("noHEL350.3_ols_p<0.05", int((p_ns[:, 0] < .05).sum())); num("noHEL350.3_ols_FDR<0.1", int((bh(p_ns[:, 0]) < .1).sum())); num("noHEL350.3_ols_log2FC_corr_with_full", np.corrcoef(b_ns[:, 0], b0[:, 0])[0, 1])
num("noHEL350.3_ols+XIST_p<0.05", int((p_ns1[:, 0] < .05).sum())); num("noHEL350.3_ols+XIST_FDR<0.1", int((bh(p_ns1[:, 0]) < .1).sum()))
log(f"    without HEL350.3: OLS p < 0.05 {int((p_ns[:, 0] < .05).sum())} (full {int((p0[:, 0] < .05).sum())}), FDR < 0.1 {int((bh(p_ns[:, 0]) < .1).sum())}, r(log2FC, full) {np.corrcoef(b_ns[:, 0], b0[:, 0])[0, 1]:.3f}; + XIST: p < 0.05 {int((p_ns1[:, 0] < .05).sum())}, FDR < 0.1 {int((bh(p_ns1[:, 0]) < .1).sum())}")
pos = pd.read_csv(a.gene_pos, sep="\t").drop_duplicates("gene_name").set_index("gene_name")
Z = (Y - Y.mean(axis=1, keepdims=True)) / Y.std(axis=1, ddof=1, keepdims=True); Zd = pd.DataFrame(Z, index=L.index, columns=lines)
pz = pos.reindex(L.index); pz["mid"] = (pz.start + pz.end) / 2
chr_of = pz.chr.fillna("NA")
c20 = pz[(pz.chr == "chr20")]
arm = pd.Series("NA", index=L.index)
arm[(pz.chr == "chr20") & (pz.end <= 25_600_000)] = "20p"; arm[(pz.chr == "chr20") & (pz.start >= 30_400_000)] = "20q"; arm[(pz.chr == "chr20") & (arm == "NA")] = "20 pericentromeric"
rows = []
for l in lines:
    for ar in ["20p", "20q"]:
        rows.append(dict(line=l, arm=ar, n_genes=int((arm == ar).sum()), median_z=float(Zd.loc[arm == ar, l].median()), mean_log2_diff_vs_others=float((L.loc[arm == ar, l] - L.loc[arm == ar, [x for x in lines if x != l]].mean(axis=1)).mean())))
ca = pd.DataFrame(rows); ca.to_csv(os.path.join(OUT, "chr20_arm_median_z.tsv"), sep="\t", index=False, float_format="%.4f")
# null: random gene sets of the same size taken from the other chromosomes (1,000 draws) for HEL350.3
rngc = np.random.default_rng(20); other = np.where(~chr_of.isin(["chr20", "NA"]).values)[0]; nullp = {}
for ar in ["20p", "20q"]:
    k = int((arm == ar).sum()); nullp[ar] = np.array([np.median(Zd["HEL350.3"].values[rngc.choice(other, k, replace=False)]) for _ in range(1000)])
    o = ca[(ca.line == "HEL350.3") & (ca.arm == ar)].median_z.iloc[0]
    pe = (np.sum(np.abs(nullp[ar]) >= abs(o)) + 1) / 1001; rank = int((ca[ca.arm == ar].median_z >= o).sum())
    num(f"chr20_{ar}_HEL350.3_median_z", o, f"{k} genes; null SD {nullp[ar].std():.3f}; empirical p {pe:.3f}; rank {rank} of 14 (1 = highest)")
    num(f"chr20_{ar}_HEL350.3_empirical_p", pe); num(f"chr20_{ar}_HEL350.3_rank_of_14", rank)
    log(f"    chr{ar}: HEL350.3 median z = {o:+.3f} over {k} genes (null from random gene sets: SD {nullp[ar].std():.3f}, empirical p = {pe:.3f}; rank {rank}/14 among the lines; the other lines' arm medians range {ca[(ca.arm == ar) & (ca.line != 'HEL350.3')].median_z.min():+.3f} to {ca[(ca.arm == ar) & (ca.line != 'HEL350.3')].median_z.max():+.3f})")
chm = pd.DataFrame({"chr": chr_of.values, "z": Zd["HEL350.3"].values}).query("chr != 'NA'").groupby("chr").z.agg(["median", "size"]); chm.to_csv(os.path.join(OUT, "HEL350.3_median_z_by_chromosome.tsv"), sep="\t", float_format="%.4f")
bins = pd.cut(c20.loc[c20.index.intersection(L.index), "mid"] / 1e6, bins=np.arange(0, 66, 3))
prof = pd.DataFrame({l: Zd.loc[c20.index.intersection(L.index), l].groupby(bins, observed=True).median() for l in lines}); prof.to_csv(os.path.join(OUT, "chr20_profile_3Mb.tsv"), sep="\t", float_format="%.4f")

# ------------------------------------------------------------------ figures (only for the main tag)
if a.figs:
    col = {"PMOS": "#c0392b", "CTRL": "#2471a3"}; lc = [col[meta.loc[l, "group"]] for l in lines]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for k_, (lab, bb, pp_, qq) in enumerate([("PMOS", b0[:, 0], p0[:, 0], q0), ("PMOS + log2 XIST", b1[:, 0], p1[:, 0], q1)]):
        ax[k_].scatter(bb, -np.log10(pp_), s=4, c=np.where(pp_ < .05, "#c0392b", "#aaaaaa"), alpha=.6); ax[k_].axhline(-np.log10(.05), color="k", ls=":", lw=.7)
        ax[k_].set_xlabel("log2 fold change (PMOS vs control)"); ax[k_].set_ylabel("-log10 p"); ax[k_].set_title(f"{lab}: p < 0.05 = {int((pp_ < .05).sum())}, FDR < 0.1 = {int((qq < .1).sum())}", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "V1_volcano.png"), dpi=160); plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for k_, (sc_, ve_, ttl) in enumerate([(sc1, ve1, "1,000 most variable genes"), (sca, vea, f"all {G:,} expressed genes")]):
        ax[k_].scatter(sc_[:, 0], sc_[:, 1], c=lc, s=45)
        for i, l in enumerate(lines): ax[k_].annotate(l.replace("HEL", ""), (sc_[i, 0], sc_[i, 1]), fontsize=7, xytext=(3, 3), textcoords="offset points")
        ax[k_].set_xlabel(f"PC1 ({100 * ve_[0]:.1f} %)"); ax[k_].set_ylabel(f"PC2 ({100 * ve_[1]:.1f} %)"); ax[k_].set_title(ttl, fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "V2_pca.png"), dpi=160); plt.close(fig)
    lk = hierarchy.linkage(C[iu] * 0 + (1 - C[iu]), "average"); order = hierarchy.leaves_list(lk)
    fig, ax = plt.subplots(figsize=(6.5, 5.5)); im = ax.imshow(C[np.ix_(order, order)], cmap="viridis", vmin=C[iu].min()); ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels([lines[i] for i in order], rotation=90, fontsize=7); ax.set_yticklabels([lines[i] for i in order], fontsize=7); fig.colorbar(im, fraction=.046); ax.set_title("Pearson r between lines (log2, expressed genes)", fontsize=9)
    for t_, i in zip(ax.get_xticklabels(), order): t_.set_color(lc[i])
    for t_, i in zip(ax.get_yticklabels(), order): t_.set_color(lc[i])
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "V3_line_correlation.png"), dpi=160); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.5, 4.3)); fs = ["PMOS status", "BMI", "Age", "XIST (log2)", "Residual"]
    ax.bar(range(5), vps.median_pct, color="#2c7fb8"); ax.scatter(range(4), vps.null_median_pct[:4], marker="_", s=400, color="k", label="permutation null median"); ax.set_xticks(range(5)); ax.set_xticklabels(fs, fontsize=8)
    ax.set_ylabel("median % of total variance (Type II SS)"); ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(FIG, "V4_variance_partition.png"), dpi=160); plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.3)); lo_ = np.log2(1)
    ax[0].scatter(t4.log2FC_PMOS_vs_control, d.log2FC, s=3, alpha=.4); lim = [min(t4.log2FC_PMOS_vs_control.min(), d.log2FC.min()), max(t4.log2FC_PMOS_vs_control.max(), d.log2FC.max())]
    ax[0].plot(lim, lim, "r:", lw=.8); ax[0].set_xlabel("the original Table S4 log2FC"); ax[0].set_ylabel(f"our log2FC ({a.tag})"); ax[0].set_title(f"log2FC, r = {np.corrcoef(d.log2FC, t4.log2FC_PMOS_vs_control)[0, 1]:.3f} ({len(inter):,} genes)", fontsize=9)
    ax[1].scatter(-np.log10(t4.p), -np.log10(d.p), s=3, alpha=.4); ax[1].set_xlabel("original -log10 p"); ax[1].set_ylabel("our -log10 p"); ax[1].set_title("p values", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "V5_log2FC_vs_S4.png"), dpi=160); plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    xs_ = np.array([i.mid / 1 for i in prof.index]) if False else np.arange(len(prof)) * 3 + 1.5
    for l in lines:
        ax[0].plot(xs_, prof[l].values, color="#bbbbbb", lw=.8)
    ax[0].plot(xs_, prof["HEL350.3"].values, color="#c0392b", lw=2, label="HEL350.3"); ax[0].axvspan(25.6, 30.4, color="#f4d03f", alpha=.3); ax[0].axhline(0, color="k", lw=.5)
    ax[0].set_xlabel("chr20 position (Mb, 3-Mb bins)"); ax[0].set_ylabel("median per-gene z across the 14 lines"); ax[0].legend(frameon=False); ax[0].set_title("chromosome 20 profile (grey = other lines; yellow = centromere region)", fontsize=9)
    ar_ = ["20p", "20q"]; w = .35
    for k_, ar in enumerate(ar_):
        v = ca[ca.arm == ar].set_index("line").median_z.reindex(lines)
        ax[1].scatter([k_] * 13, v.drop("HEL350.3"), color="#999999", s=25); ax[1].scatter([k_], v["HEL350.3"], color="#c0392b", s=70, zorder=3)
    ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(ar_); ax[1].axhline(0, color="k", lw=.5); ax[1].set_ylabel("median z of the arm's genes"); ax[1].set_title("arm medians per line (red = HEL350.3)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "V7_chr20_HEL350.3.png"), dpi=160); plt.close(fig)

pd.DataFrame(NUM).to_csv(os.path.join(OUT, "numbers.tsv"), sep="\t", index=False)
open(os.path.join(OUT, "report.txt"), "w").write("\n".join(rep) + "\n")
