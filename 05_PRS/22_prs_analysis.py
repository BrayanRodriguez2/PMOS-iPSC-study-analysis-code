#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/22_prs_analysis.py   (original step script: 52_prs_analysis.py)
# Produces: Case-control statistics, stability, comparison with an external PRS, figures.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 52 - analysis and figures of the PRS computed on the IMPUTED clean dataset (steps 49-51). DESCRIPTIVE only (few samples, no covariates).

Scores (SCORE1_SUM of plink2 --score on imputed dosages, 135 clean samples; see 05_PRS/18_score_imputed.sbatch):
  V1 (main; PRS-CS, UKBB panel), V2 (1000G panel), V3 (UKBB panel, whole-panel bim), V1_r2ge0.8 (V1 weights restricted to SNPs with imputation R2 >= 0.8), GRS29 (29 ST2 leads).
Reference PRS of an external group (clumping + p-value thresholds 1e-3, 1e-4, 1e-5, 1e-6; Oulu bank; 200 kb windows; column PRS.3..PRS.6), scored on all 140 samples.
Sample sets (z-scores are always made WITHIN the set that is analysed: z = (x - mean) / SD, SD with n-1):
  main_129 : the 129 QC-pass samples (keep_129.txt)      s133 : keep_135 without the contaminated samples      main_minus_lowlevel : main_129 without the low-level-contamination sample
Case-control labels: DNA-extraction sheet (step 48b); PMOS = case, CTRL = control; ALL donors are labelled, including 22-30 (the asterisk / Nanodrop note of the sheet refers to the DNA reading, not to the diagnosis; correction 2026-11-05).
Statistics per score and set: difference of mean z (PMOS - CTRL) with Welch 95 % CI; Wilcoxon rank-sum (Mann-Whitney U, two-sided); logistic regression
  (PMOS ~ z, OR per SD, Wald 95 % CI and p); AUC (PMOS = positive, higher score = more risk) with DeLong 95 % CI.
Also: Pearson and Spearman correlations between scores and with the reference PRS of an external group (common samples), the 14 iPSC donors (percentile and z within the
  129; two donors are not in the clean file and are NOT scored), and figures.
Outputs (--outdir): tables (TSV), summary_report.txt, figures/*.png
"""
import argparse, os, re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--scoredir", required=True); ap.add_argument("--weightdir", required=True); ap.add_argument("--clean-dir", required=True)
ap.add_argument("--labels", required=True); ap.add_argument("--donor-key", required=True); ap.add_argument("--extprs", required=True)
ap.add_argument("--bim-counts", required=True); ap.add_argument("--grs-orientation", required=True); ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(os.path.join(a.outdir, "figures"), exist_ok=True)
FIG = os.path.join(a.outdir, "figures")
SCORES = ["V1", "V2", "V3", "V1_r2ge0.8", "GRS29"]
EXTPRS = {"PRS.3": "ExtPRS_p1e-3", "PRS.4": "ExtPRS_p1e-4", "PRS.5": "ExtPRS_p1e-5", "PRS.6": "ExtPRS_p1e-6"}
COL = {"PMOS": "#c0392b", "CTRL": "#2471a3", "unlabelled": "#999999"}
rep = []


def log(m=""):
    print(m, flush=True); rep.append(str(m))


# ------------------------------------------------------------------------------------------------ inputs
def read_keep(fn):
    d = pd.read_csv(fn, sep=r"\s+"); return set(d.iloc[:, 1].astype(int))   # header "#FID IID" is space-separated, rows tab-separated


k129 = read_keep(os.path.join(a.clean_dir, "keep_129.txt")); k135 = read_keep(os.path.join(a.clean_dir, "keep_135.txt"))
_ids = lambda n: {int(x) for x in __import__("os").environ.get(n, "").split(",") if x}   # sample IDs come from the environment (config/config.sh)
SETS = {"main_129": k129, "s133": k135 - _ids("SENSITIVITY_EXCLUDE_IDS"), "main_minus_lowlevel": k129 - _ids("LOWLEVEL_SAMPLE_IDS")}
assert len(k129) == 129 and len(k135) == 135 and k129 <= k135
sc = pd.DataFrame(index=sorted(k135))
allele_ct = {}
for s in SCORES:
    d = pd.read_csv(os.path.join(a.scoredir, f"{s}.sscore"), sep="\t"); d.columns = [c.lstrip("#") for c in d.columns]
    d = d.set_index("IID")
    assert set(d.index) == k135, f"{s}: scored samples differ from keep_135"
    assert d.SCORE1_SUM.notna().all()
    allele_ct[s] = sorted(d.ALLELE_CT.unique().tolist())
    sc[s] = d.SCORE1_SUM
log(f"samples scored: {len(sc)} (all of keep_135; none missing). ALLELE_CT per sample (one value = no missing dosages): " + "; ".join(f"{k}: {v}" for k, v in allele_ct.items()))
lab = pd.read_csv(a.labels, sep="\t").set_index("esf_id")
sc["group_cc"] = [lab.group_cc.get(i) if pd.notna(lab.group_cc.get(i)) and lab.group_cc.get(i) != "" else "unlabelled" for i in sc.index]
sc["nanodrop_flag_22_30"] = [lab.nanodrop_flag_22_30.get(i, "") for i in sc.index]
log("case-control labels in main_129: " + str(sc.loc[sorted(k129), "group_cc"].value_counts().to_dict()))   # expected: no 'unlabelled'
li = pd.read_csv(a.extprs, sep="\t").set_index("Sample.ID")
li.index = li.index.astype(int)
for c, n in EXTPRS.items():
    sc[n] = li[c].reindex(sc.index)
log(f"an external group PRS: {len(li)} samples; present among the 135: {int(sc[list(EXTPRS.values())[0]].notna().sum())}")


def zs(x, ref):
    return (x - x[ref].mean()) / x[ref].std(ddof=1)


# ------------------------------------------------------------------------------------------------ statistics helpers
def delong(y, s):
    pos, neg = s[y == 1], s[y == 0]; m, n = len(pos), len(neg)
    v10 = np.array([(np.sum(neg < p) + 0.5 * np.sum(neg == p)) / n for p in pos])
    v01 = np.array([(np.sum(pos > q) + 0.5 * np.sum(pos == q)) / m for q in neg])
    auc = v10.mean(); var = v10.var(ddof=1) / m + v01.var(ddof=1) / n
    se = np.sqrt(var)
    return auc, max(0.0, auc - 1.96 * se), min(1.0, auc + 1.96 * se)


def logit(y, z):
    X = np.column_stack([np.ones(len(z)), z]); b = np.zeros(2)
    for _ in range(100):
        p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
        H = X.T @ (X * W[:, None]); g = X.T @ (y - p)
        step = np.linalg.solve(H, g); b = b + step
        if np.abs(step).max() < 1e-10:
            break
    cov = np.linalg.inv(X.T @ (X * (p * (1 - p))[:, None])); se = np.sqrt(np.diag(cov))
    zval = b[1] / se[1]; pval = 2 * (1 - stats.norm.cdf(abs(zval)))
    return np.exp(b[1]), np.exp(b[1] - 1.96 * se[1]), np.exp(b[1] + 1.96 * se[1]), pval


def casectrl(setname, name):
    S = sorted(SETS[setname])
    v = sc[name].loc[S].astype(float); ok = v.notna()
    v = v[ok]; S = list(v.index)
    z = (v - v.mean()) / v.std(ddof=1)
    g = sc.group_cc.loc[S]; keep = g.isin(["PMOS", "CTRL"]); zz = z[keep]; y = (g[keep] == "PMOS").astype(int).values
    zp, zc = zz[y == 1].values, zz[y == 0].values
    diff = zp.mean() - zc.mean(); sed = np.sqrt(zp.var(ddof=1) / len(zp) + zc.var(ddof=1) / len(zc))
    dfw = sed ** 4 / ((zp.var(ddof=1) / len(zp)) ** 2 / (len(zp) - 1) + (zc.var(ddof=1) / len(zc)) ** 2 / (len(zc) - 1)); tcrit = stats.t.ppf(0.975, dfw)
    u = stats.mannwhitneyu(zp, zc, alternative="two-sided")
    orr, olo, ohi, op = logit(y.astype(float), zz.values)
    auc, alo, ahi = delong(y, zz.values)
    return dict(set=setname, score=name, n_set=len(S), n_PMOS=len(zp), n_CTRL=len(zc), mean_z_PMOS=zp.mean(), mean_z_CTRL=zc.mean(), diff_z=diff, diff_lo=diff - tcrit * sed,
                diff_hi=diff + tcrit * sed, wilcoxon_p=u.pvalue, OR_per_SD=orr, OR_lo=olo, OR_hi=ohi, logit_p=op, AUC=auc, AUC_lo=alo, AUC_hi=ahi), zz, y


def corr(x, y):
    ok = x.notna() & y.notna(); x, y = x[ok], y[ok]; n = len(x)
    r, p = stats.pearsonr(x, y); rho, ps = stats.spearmanr(x, y)
    fz = np.arctanh(r); lo, hi = np.tanh(fz - 1.96 / np.sqrt(n - 3)), np.tanh(fz + 1.96 / np.sqrt(n - 3))
    return dict(n=n, pearson_r=r, pearson_lo=lo, pearson_hi=hi, pearson_p=p, spearman_rho=rho, spearman_p=ps)


# ------------------------------------------------------------------------------------------------ (a) stability
M = sorted(k129)
rows = []
for x, y in [("V1", "V2"), ("V1", "V3"), ("V2", "V3"), ("V1", "V1_r2ge0.8"), ("V1", "GRS29"), ("V2", "GRS29"), ("V3", "GRS29")]:
    rows.append(dict(a=x, b=y, **corr(sc.loc[M, x], sc.loc[M, y])))
cor_a = pd.DataFrame(rows); cor_a.to_csv(os.path.join(a.outdir, "correlations_between_scores.tsv"), sep="\t", index=False, float_format="%.4g")
log("\n(a) correlations between scores, 129 QC samples (before imputation V1-V2 0.874, V1-V3 0.730, V2-V3 0.714):")
log(cor_a.round(4).to_string(index=False))

# ------------------------------------------------------------------------------------------------ (b), (e) case vs control
cc, zstore = [], {}
for setname in SETS:
    for name in SCORES + list(EXTPRS.values()):
        r, zz, y = casectrl(setname, name); cc.append(r); zstore[(setname, name)] = (zz, y)
ccdf = pd.DataFrame(cc); ccdf.to_csv(os.path.join(a.outdir, "case_control_stats.tsv"), sep="\t", index=False, float_format="%.4g")
log("\n(b) case (PMOS) vs control (CTRL), main_129 (z within the 129; all donors labelled, including 22-30):")
log(ccdf[ccdf.set == "main_129"][["score", "n_PMOS", "n_CTRL", "diff_z", "diff_lo", "diff_hi", "wilcoxon_p", "OR_per_SD", "OR_lo", "OR_hi", "logit_p", "AUC", "AUC_lo", "AUC_hi"]].round(3).to_string(index=False))
log("\n(e) sensitivity: V1 and GRS29 in the other sets, and V1 restricted to R2 >= 0.8:")
log(ccdf[ccdf.score.isin(["V1", "V1_r2ge0.8", "GRS29"])][["set", "score", "n_set", "n_PMOS", "n_CTRL", "diff_z", "wilcoxon_p", "OR_per_SD", "OR_lo", "OR_hi", "logit_p", "AUC", "AUC_lo", "AUC_hi"]].round(3).to_string(index=False))

# ------------------------------------------------------------------------------------------------ (c) ExtPRS
rows = []
for name in SCORES:
    for n in EXTPRS.values():
        rows.append(dict(score=name, extprs=n, **corr(sc.loc[M, name], sc.loc[M, n])))
cl = pd.DataFrame(rows); cl.to_csv(os.path.join(a.outdir, "correlation_with_extprs_prs.tsv"), sep="\t", index=False, float_format="%.4g")
log("\n(c) correlation with the PRS of an external group, 129 samples:")
log(cl.round(4).to_string(index=False))
liis = pd.DataFrame([dict(a=n1, b=n2, **corr(sc.loc[M, n1], sc.loc[M, n2])) for i, n1 in enumerate(EXTPRS.values()) for n2 in list(EXTPRS.values())[i + 1:]])
log("\n  (for reference) ExtPRS's thresholds among themselves:\n" + liis[["a", "b", "pearson_r", "spearman_rho"]].round(3).to_string(index=False))

# ------------------------------------------------------------------------------------------------ (d) 14 donors
key = pd.read_csv(a.donor_key, sep="\t")
pct = lambda name, i: stats.percentileofscore(sc.loc[M, name], sc.loc[i, name], kind="mean")
drows = []
for r in key.itertuples():
    i = int(r.esf_id); d = dict(ipsc_line=r.ipsc_line, esf_id=i, group=("PMOS" if r.group in ("PCOS", "PMOS") else "CTRL"), group_dna_sheet=lab.group_sheet.get(i, ""),
                                 in_clean_file=i in k135, in_main_129=i in k129)
    for name in SCORES:
        if i in k135:
            ref = sc.loc[M, name]; d[f"{name}_raw"] = sc.loc[i, name]; d[f"{name}_z"] = (sc.loc[i, name] - ref.mean()) / ref.std(ddof=1); d[f"{name}_pct"] = pct(name, i)
        else:
            d[f"{name}_raw"] = d[f"{name}_z"] = d[f"{name}_pct"] = np.nan
    note = []
    if i not in k135:
        note.append("NOT in the clean file (excluded at QC: contamination / call rate): no score computed")
    if i in _ids("LOWLEVEL_SAMPLE_IDS"): note.append("low-level contamination estimate: kept in main_129, see the sensitivity set without it")
    d["note"] = "; ".join(note)
    drows.append(d)
dn = pd.DataFrame(drows); dn.to_csv(os.path.join(a.outdir, "donors_14_prs_imputed.tsv"), sep="\t", index=False, float_format="%.4g")
log("\n(d) the 14 donors of the iPSC lines (z and percentile within the 129; donors absent from the clean file are not scored):")
log(dn[["ipsc_line", "esf_id", "group", "in_main_129", "V1_z", "V1_pct", "V2_pct", "V3_pct", "V1_r2ge0.8_pct", "GRS29_z", "GRS29_pct"]].round(2).to_string(index=False))
d12 = dn[dn.in_main_129]
for name in ["V1", "GRS29"]:
    p = d12[d12.group == "PMOS"][f"{name}_pct"]; c = d12[d12.group == "CTRL"][f"{name}_pct"]
    log(f"  {name}: mean percentile of the {len(p)} PMOS donors scored = {p.mean():.1f}, of the {len(c)} CTRL = {c.mean():.1f}")

# ------------------------------------------------------------------------------------------------ per-sample table
out = sc.copy(); out.index.name = "esf_id"
out["in_main_129"] = [i in k129 for i in out.index]; out["in_s133"] = [i in SETS["s133"] for i in out.index]
for name in SCORES:
    ref = sc.loc[M, name]; out[f"{name}_z129"] = (sc[name] - ref.mean()) / ref.std(ddof=1)
out.to_csv(os.path.join(a.outdir, "prs_scores_per_sample.tsv"), sep="\t", float_format="%.6g")

# ------------------------------------------------------------------------------------------------ SNP counts
bc = pd.read_csv(a.bim_counts, sep="\t", header=None, names=["step", "UKBB_EUR", "1000G_EUR"])
bc.to_csv(os.path.join(a.outdir, "snp_counts_bim_steps.tsv"), sep="\t", index=False)
srows = []
for s in SCORES:
    w = pd.read_csv(os.path.join(a.weightdir, f"{s}_weights.txt"), sep="\t")
    n_scored = sum(1 for _ in open(os.path.join(a.scoredir, f"{s}.sscore.vars")))
    srows.append(dict(set=s, n_weights_file=len(w), n_variants_scored_by_plink2=n_scored, n_weights_nonzero=int((w.BETA != 0).sum()), sum_abs_beta=w.BETA.abs().sum(), sum_beta_squared=(w.BETA ** 2).sum()))
sn = pd.DataFrame(srows); sn.to_csv(os.path.join(a.outdir, "snp_counts_scored.tsv"), sep="\t", index=False, float_format="%.5g")
log("\nSNPs per score:\n" + sn.to_string(index=False))
log("\nSNPs of the LD panel at each step of the bim construction:\n" + bc.to_string(index=False))

# ------------------------------------------------------------------------------------------------ figures
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def donor_marks(ax, name, dtab):
    xs = []
    for r in dtab[dtab.in_main_129].sort_values(f"{name}_z").itertuples():
        xs.append((getattr(r, f"{name}_z"), r.ipsc_line, r.group))
    ymax = ax.get_ylim()[1]
    for k, (x, ln, g) in enumerate(xs):
        yb = ymax * (0.05 + 0.12 * (k % 6))          # six staggered heights so that neighbouring donors do not overprint their labels
        ax.plot([x, x], [0, yb], color=COL[g], lw=1.3); ax.plot(x, yb, "o", color=COL[g], ms=5)
        ax.text(x + 0.02, yb + ymax * 0.01, ln, rotation=90, ha="left", va="bottom", fontsize=6.5, color=COL[g])
    ax.set_ylim(0, ymax * 1.9)


fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
for ax, name, ttl in zip(axs, ["V1", "GRS29"], ["PRS-CS V1 (imputed, UKBB panel)", "GRS of the 29 lead SNPs (ST2)"]):
    z = (sc.loc[M, name] - sc.loc[M, name].mean()) / sc.loc[M, name].std(ddof=1)
    ax.hist(z, bins=24, color="#d5d8dc", edgecolor="white"); ax.set_xlabel("z-score within the 129 QC samples"); ax.set_ylabel("donors"); ax.set_title(ttl)
    donor_marks(ax, name, dn)
fig.text(0.5, 0.005, "iPSC donors marked (red = PMOS, blue = CTRL); donors that are not in the clean file are not scored", ha="center", fontsize=8)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(FIG, "F1_distribution_14_donors.png"), dpi=160); plt.close(fig)


def roc_points(y, s):
    o = np.argsort(-s); ys = y[o]; tpr = np.concatenate([[0], np.cumsum(ys) / ys.sum()]); fpr = np.concatenate([[0], np.cumsum(1 - ys) / (1 - ys).sum()])
    return fpr, tpr


def fig_cc(name, ttl, fn):
    zz, y = zstore[("main_129", name)]; r = ccdf[(ccdf.set == "main_129") & (ccdf.score == name)].iloc[0]
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    data = [zz[y == 0].values, zz[y == 1].values]
    vp = ax[0].violinplot(data, positions=[0, 1], showextrema=False, widths=0.8)
    for b, c in zip(vp["bodies"], [COL["CTRL"], COL["PMOS"]]): b.set_facecolor(c); b.set_alpha(0.25)
    ax[0].boxplot(data, positions=[0, 1], widths=0.25, showfliers=False, medianprops=dict(color="k"))
    rng = np.random.default_rng(1)
    for k, c in enumerate([COL["CTRL"], COL["PMOS"]]):
        ax[0].scatter(k + rng.uniform(-0.12, 0.12, len(data[k])), data[k], s=12, color=c, alpha=0.8)
    ax[0].set_xticks([0, 1]); ax[0].set_xticklabels([f"CTRL (n={r.n_CTRL})", f"PMOS (n={r.n_PMOS})"]); ax[0].set_ylabel("z-score (129 QC samples)")
    ax[0].set_title(f"{ttl}\nWilcoxon p = {r.wilcoxon_p:.3g}; OR/SD = {r.OR_per_SD:.2f} ({r.OR_lo:.2f}-{r.OR_hi:.2f})", fontsize=9)
    fpr, tpr = roc_points(y, zz.values)
    ax[1].plot(fpr, tpr, color="k", lw=1.6); ax[1].plot([0, 1], [0, 1], ":", color="grey")
    ax[1].set_xlabel("1 - specificity"); ax[1].set_ylabel("sensitivity"); ax[1].set_title(f"ROC: AUC = {r.AUC:.3f} (95 % CI {r.AUC_lo:.3f}-{r.AUC_hi:.3f})", fontsize=9)
    ax[1].set_aspect("equal")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, fn), dpi=160); plt.close(fig)


fig_cc("V1", "PRS-CS V1 (imputed)", "F2_V1_cases_vs_controls.png")
fig_cc("GRS29", "GRS 29 SNPs (ST2)", "F2b_GRS29_cases_vs_controls.png")


def scatter_grid(pairs, fn, ncol, ttl=None):
    n = len(pairs); nrow = int(np.ceil(n / ncol)); fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4 * nrow), squeeze=False)
    for ax, (xn, yn, xl, yl) in zip(axs.ravel(), pairs):
        d = sc.loc[M]; x, y = d[xn], d[yn]
        ok = x.notna() & y.notna(); xz = (x - x[ok].mean()) / x[ok].std(); yz = (y - y[ok].mean()) / y[ok].std()
        for g in ["unlabelled", "CTRL", "PMOS"]:
            m = (d.group_cc == g) & ok
            if not m.any():
                continue                 # no empty legend entry (there is no unlabelled group any more)
            ax.scatter(xz[m], yz[m], s=14, color=COL[g], alpha=0.75, label=g if g != "unlabelled" else "unlabelled (22-30)")
        r = corr(x, y); lim = [min(xz.min(), yz.min()) - 0.3, max(xz.max(), yz.max()) + 0.3]
        ax.plot(lim, lim, ":", color="grey", lw=0.8); ax.set_xlabel(xl + " (z)"); ax.set_ylabel(yl + " (z)")
        ax.set_title(f"Pearson r = {r['pearson_r']:.2f}, Spearman = {r['spearman_rho']:.2f} (n = {r['n']})", fontsize=8.5)
    axs.ravel()[0].legend(fontsize=7, frameon=False)
    for ax in axs.ravel()[n:]: ax.axis("off")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, fn), dpi=160); plt.close(fig)


scatter_grid([("V1", "V2", "PRS V1", "PRS V2 (1000G panel)"), ("V1", "V3", "PRS V1", "PRS V3 (whole-panel bim)"), ("V1", "V1_r2ge0.8", "PRS V1", "PRS V1, SNPs R2 >= 0.8")], "F3_V1_vs_V2_V3_R2.png", 3)
scatter_grid([("V1", "GRS29", "PRS V1", "GRS 29 SNPs")], "F4_V1_vs_GRS29.png", 1)
lab_l = lambda n: f"an external group PRS (p < {n.split('_p')[1]})"
scatter_grid([("V1", n, "PRS V1", lab_l(n)) for n in EXTPRS.values()], "F5_V1_vs_ExtPRS.png", 2)
scatter_grid([("GRS29", n, "GRS 29 SNPs", lab_l(n)) for n in EXTPRS.values()], "F5b_GRS29_vs_ExtPRS.png", 2)

# sensitivity forest
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
lbl = [("main_129", "V1"), ("s133", "V1"), ("main_minus_lowlevel", "V1"), ("main_129", "V1_r2ge0.8")]
for ax, names, ttl in [(axs[0], lbl, "PRS V1: AUC (95 % CI)"), (axs[1], [(s, "GRS29") for s in SETS], "GRS 29 SNPs: AUC (95 % CI)")]:
    for k, (s, n) in enumerate(names):
        r = ccdf[(ccdf.set == s) & (ccdf.score == n)].iloc[0]
        ax.errorbar(r.AUC, k, xerr=[[r.AUC - r.AUC_lo], [r.AUC_hi - r.AUC]], fmt="o", color="k", capsize=3)
    ax.axvline(0.5, ls=":", color="grey"); ax.set_yticks(range(len(names))); ax.set_yticklabels([f"{s}\n{n}" for s, n in names], fontsize=7.5); ax.invert_yaxis(); ax.set_xlabel("AUC"); ax.set_title(ttl, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "F6_sensitivity_AUC.png"), dpi=160); plt.close(fig)

open(os.path.join(a.outdir, "summary_report.txt"), "w").write("\n".join(rep) + "\n")
os.system(f"cp '{a.grs_orientation}' '{os.path.join(a.outdir, 'grs29_orientation.tsv')}'")
log("\nfigures written: " + ", ".join(sorted(os.listdir(FIG))))
