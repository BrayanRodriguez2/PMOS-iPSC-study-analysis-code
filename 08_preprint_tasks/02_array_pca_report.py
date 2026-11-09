#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/02_array_pca_report.py   (original step script: 80b_array_pca_report.py)
# Produces: Variance explained, outlier check, comparison with the external table, figures.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 80b - report and figures of the array genotype PCA (step 80): variance explained, position of the 12 iPSC donors that pass the QC (two donors failed it and are not
in the 129), whether any donor lies outside the cloud of the cohort, whether PC1 is still dominated by a related pair / outlier (one related sample is already excluded), comparison with the PCs of the original analysis' Table S2,
and the check of the array ID <-> iPSC line key against the original Table S1.
  * % of variance of PC k = eigenvalue k / sum of all eigenvalues (plink2 --pca 128 gives all components of the 129-sample matrix).
  * 'outside the cloud': a donor is flagged if |robust z| > 3.5 (median / 1.4826 MAD of the 129) on any of PC1-PC5, or if the squared Mahalanobis distance in PC1-PC5 (cohort covariance) has chi2(5) p < 0.001.
  * PC1 domination: share of PC1 sum of squares carried by the top 1 and top 2 samples, and the top samples; kinship pairs >= 0.0884 among the samples of the PCA (king.kin0 of step 09).
Usage: python 08_preprint_tasks/02_array_pca_report.py <array_pca_dir> <donor_key.tsv> <rna_samples.tsv> <sample_qc_table.tsv> <king.kin0> <the original S2 tsv> <the original S1 tsv> <out_dir>
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pdir, keyf, rnaf, qcf, kingf, s2f, s1f, out = sys.argv[1:9]
os.makedirs(f"{out}/figures", exist_ok=True)
rep = []
log = lambda m="": (print(m), rep.append(str(m)))
COL = {"PMOS": "#c0392b", "CTRL": "#2471a3"}


def load(tag):
    v = pd.read_csv(f"{pdir}/pca{tag}.eigenvec", sep="\t"); v = v.rename(columns={"#FID": "FID"}).set_index("IID")
    e = pd.read_csv(f"{pdir}/pca{tag}.eigenval", header=None)[0].values
    return v, 100 * e / e.sum()


v129, ve129 = load("129"); v140, ve140 = load("140")
key = pd.read_csv(keyf, sep="\t"); key["group"] = key.group.replace({"PCOS": "PMOS"}); key = key.set_index("esf_id")
pcs = [c for c in v129.columns if c.startswith("PC")]
log(f"PCA of the array data (autosomes, MAF >= 0.05, LD-pruned 200/50/0.2, long-range LD regions excluded): {len(v129)} QC samples; variance explained (%) PC1-PC10: " + ", ".join(f"{x:.2f}" for x in ve129[:10]))
log(f"   same SNPs, all 140 samples: PC1-PC10 (%): " + ", ".join(f"{x:.2f}" for x in ve140[:10]))
pd.DataFrame({"PC": range(1, 11), "pct_var_129_QC": ve129[:10], "pct_var_140": ve140[:10]}).to_csv(f"{out}/array_pca_variance.tsv", sep="\t", index=False, float_format="%.3f")
X = v129[pcs[:5]].values; med = np.median(X, axis=0); mad = 1.4826 * np.median(np.abs(X - med), axis=0)
cov = np.cov(X.T); ci = np.linalg.inv(cov)
rows = []
for esf, r in key.iterrows():
    if esf in v129.index:
        x = v129.loc[esf, pcs[:5]].values.astype(float); z = (x - med) / mad; d2 = float((x - X.mean(0)) @ ci @ (x - X.mean(0))); p = stats.chi2.sf(d2, 5)
        flag = bool((np.abs(z) > 3.5).any() or p < 0.001)
        rows.append(dict(ipsc_line=r.ipsc_line, esf_id=esf, group=r.group, in_QC_129=True, **{f"PC{i + 1}": x[i] for i in range(5)}, **{f"z_PC{i + 1}": z[i] for i in range(5)}, mahalanobis_d2=d2, chi2_p=p, outside_cloud=flag))
    else:
        rows.append(dict(ipsc_line=r.ipsc_line, esf_id=esf, group=r.group, in_QC_129=False))
D = pd.DataFrame(rows); D.to_csv(f"{out}/array_pca_donors.tsv", sep="\t", index=False, float_format="%.4g")
d = D[D.in_QC_129]
log(f"iPSC donors in the 129 QC samples: {len(d)} of 14 (not in the QC set: {list(D[~D.in_QC_129].ipsc_line)} = donors that failed the QC)")
log(f"   donors outside the cloud (|robust z| > 3.5 on PC1-5 or Mahalanobis p < 0.001): {list(d[d.outside_cloud].ipsc_line) or 'none'}; max |z| over PC1-5 per donor: " + ", ".join(f"{l} {np.abs(d[[f'z_PC{i}' for i in range(1, 6)]].loc[i_]).max():.1f}" for i_, l in zip(d.index, d.ipsc_line)))
# PC1 domination
for tag, v in [("129 QC", v129), ("140", v140)]:
    out_s = []
    for k in (1, 2):
        s = v[f"PC{k}"].values; ss = s ** 2 / (s ** 2).sum(); o = np.argsort(-ss)
        out_s.append(f"PC{k}: top sample {v.index[o[0]]} = {100 * ss[o[0]]:.1f} % of its sum of squares, top 2 ({v.index[o[0]]}, {v.index[o[1]]}) = {100 * (ss[o[0]] + ss[o[1]]):.1f} %; kurtosis {stats.kurtosis(s):.1f}")
    log(f"   [{tag}] " + " | ".join(out_s))
kin = pd.read_csv(kingf, sep=r"\s+"); kc = [c for c in kin.columns if c.upper().lstrip("#") in ("IID1", "ID1")][0]; k2 = [c for c in kin.columns if c.upper().lstrip("#") in ("IID2", "ID2")][0]
kv = [c for c in kin.columns if c.upper() in ("KINSHIP", "PI_HAT")][0]
kin = kin[kin[kv] >= 0.0884]; qs = set(v129.index)
log("   kinship >= 0.0884 pairs among the 140 (king.kin0): " + "; ".join(f"{a}-{b} {k:.3f}{' (both in the 129)' if a in qs and b in qs else ''}" for a, b, k in zip(kin[kc], kin[k2], kin[kv])))
# comparison with the original S2 PCs
s2 = pd.read_csv(s2f, sep="\t").set_index("Genotype_array_ID"); s2.index = s2.index.astype(int)
common = v140.index.intersection(s2.index); log(f"   the original Table S2 PCs vs ours (140-sample PCA, {len(common)} samples): |Pearson r| between the original PC i (rows) and our PC j (columns, PC1-5):")
cm = pd.DataFrame({f"our_PC{j}": [abs(np.corrcoef(s2.loc[common, f"PC{i}"].astype(float), v140.loc[common, f"PC{j}"])[0, 1]) for i in range(1, 6)] for j in range(1, 6)}, index=[f"her_PC{i}" for i in range(1, 6)])
cm.to_csv(f"{out}/array_pca_vs_her_S2.tsv", sep="\t", float_format="%.3f"); log(cm.round(2).to_string())
# key check
s1 = pd.read_csv(s1f, sep="\t"); s1 = s1[s1.iPSC_line_RNAseq.astype(str).str.startswith("HEL")].set_index("iPSC_line_RNAseq")
rna = pd.read_csv(rnaf, sep="\t").set_index("ipsc_line")
kr = []
for esf, r in key.iterrows():
    kr.append(dict(ipsc_line=r.ipsc_line, ours_donor_key=esf, her_S1_array_ID=int(s1.loc[r.ipsc_line, "Genotype_array_ID"]), rna_identity_key=int(rna.loc[r.ipsc_line, "esf_id"]),
                   group_ours=r.group, group_her_S1="PMOS" if s1.loc[r.ipsc_line, "Group"] == "PMOS" else "CTRL"))
K = pd.DataFrame(kr); K["match"] = (K.ours_donor_key == K.her_S1_array_ID) & (K.ours_donor_key == K.rna_identity_key) & (K.group_ours == K.group_her_S1); K.to_csv(f"{out}/array_id_line_key_check.tsv", sep="\t", index=False)
log(f"   array ID <-> line key: ours = the original S1 = RNA-seq identity key for {int(K.match.sum())}/14 lines ("
    + ", ".join(f"{e} -> {K[K.ours_donor_key == e].ipsc_line.iloc[0]}" for e in [int(x) for x in __import__("os").environ.get("KEY_CHECK_SAMPLE_IDS", "").split(",") if x]) + ")")
v129.assign(group=[key.group.get(i, "") for i in v129.index], ipsc_line=[key.ipsc_line.get(i, "") for i in v129.index]).iloc[:, :].to_csv(f"{out}/array_pca_scores_129.tsv", sep="\t", float_format="%.5g")


def panel(ax, v, ve, xk, yk, donors, extra=None):
    ax.scatter(v[f"PC{xk}"], v[f"PC{yk}"], s=16, c="#c8c8c8", label=f"cohort (n = {len(v)})")
    for esf, r in donors.iterrows():
        if esf in v.index:
            ax.scatter(v.loc[esf, f"PC{xk}"], v.loc[esf, f"PC{yk}"], s=55, c=COL[r.group], edgecolor="k", lw=.6, zorder=3)
            ax.annotate(r.ipsc_line, (v.loc[esf, f"PC{xk}"], v.loc[esf, f"PC{yk}"]), fontsize=6.5, xytext=(4, 3), textcoords="offset points", color=COL[r.group])
    if extra:
        for esf, lab in extra.items():
            if esf in v.index: ax.scatter(v.loc[esf, f"PC{xk}"], v.loc[esf, f"PC{yk}"], s=40, marker="X", c="k", zorder=4); ax.annotate(lab, (v.loc[esf, f"PC{xk}"], v.loc[esf, f"PC{yk}"]), fontsize=6.5, xytext=(4, -9), textcoords="offset points")
    ax.set_xlabel(f"PC{xk} ({ve[xk - 1]:.1f} %)"); ax.set_ylabel(f"PC{yk} ({ve[yk - 1]:.1f} %)")


fig, ax = plt.subplots(1, 2, figsize=(11, 5))
panel(ax[0], v129, ve129, 2, 3, key); panel(ax[1], v129, ve129, 1, 2, key); ax[0].set_title("PC2 vs PC3 (as in the original Fig. 2G)", fontsize=9); ax[1].set_title("PC1 vs PC2", fontsize=9)
ax[0].scatter([], [], c=COL["PMOS"], edgecolor="k", s=40, label="iPSC donor, PMOS"); ax[0].scatter([], [], c=COL["CTRL"], edgecolor="k", s=40, label="iPSC donor, control"); ax[0].legend(frameon=False, fontsize=7)
fig.text(0.5, 0.005, "Array data (not imputed); 129 QC samples; 78,749 LD-pruned autosomal SNPs, MAF >= 0.05, long-range LD regions excluded. Donors that failed the QC are not shown.", ha="center", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(f"{out}/figures/P1_array_PCA_129_QC.png", dpi=170); plt.close(fig)
fig, ax = plt.subplots(1, 2, figsize=(11, 5))
ex = {int(x): x for x in __import__("os").environ.get("EXTRA_PLOT_SAMPLE_IDS", "").split(",") if x}
panel(ax[0], v140, ve140, 2, 3, key, ex); panel(ax[1], v140, ve140, 1, 2, key, ex); ax[0].set_title("all 140 samples: PC2 vs PC3", fontsize=9); ax[1].set_title("all 140 samples: PC1 vs PC2 (marked samples: first-degree pair; X = samples not in the 129)", fontsize=8)
fig.tight_layout(); fig.savefig(f"{out}/figures/P2_array_PCA_140_all.png", dpi=170); plt.close(fig)
open(f"{out}/array_pca_report.txt", "w").write("\n".join(rep) + "\n")
