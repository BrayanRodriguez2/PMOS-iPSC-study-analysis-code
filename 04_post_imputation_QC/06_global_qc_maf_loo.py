#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/06_global_qc_maf_loo.py   (original step script: 41_qc_global_maf_loo.py)
# Produces: Global R2, R2 by MAF interval and leave-one-out ER2 of genotyped SNPs.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 41 - Post-imputation QC, parts (a), (b), (c): global summary, R2 by MAF interval and leave-one-out accuracy of the genotyped variants.
Diagnostic only; reads the tables of step 40 and writes new tables / figures.

(a) variants per chromosome; Typed / Imputed (and Typed_only if present) counts; distribution of R2.
(b) R2 by MAF interval [<0.001, 0.001-0.005, 0.005-0.01, 0.01-0.05, 0.05-0.2, 0.2-0.5]: n, mean R2, median R2, % with R2 >= 0.3 and >= 0.8,
    for all variants and for imputed variants only. MAF is the server's estimate over the 140 uploaded samples.
    IMPORTANT: the job was run with the server filter "R2 >= 0.3", so every delivered variant has R2 >= 0.3: "% R2 >= 0.3" is 100 % by
    construction, the R2 distribution is truncated at 0.3 and the mean R2 per MAF interval is biased upwards (variants with R2 < 0.3
    were dropped, mostly the rarest). Typed variants are reported with the same filter; their R2 minimum is printed.
(c) ER2 = empirical leave-one-out R2 (the server provides it only for genotyped variants; there is no EmpR column in this server version):
    distribution, by MAF interval, and for the subsets of the QC-only review: the 75 kept frequency-discordant variants
    (keep_fin_variants_75_v2.tsv) and the 358 kept discordant ones (compare_1000g_446.tsv, |diff vs TOPMed| <= 0.2) versus the other genotyped SNPs.
Inputs : info_tables/chr*.tsv.gz (step 40); results/topmed_qc_only/{keep_fin_variants_88.tsv, compare_1000g_446.tsv} (steps 34-35)
Outputs: results/imputation_qc/{summary_by_chrom.tsv, r2_by_maf.tsv, r2_by_maf_imputed_only.tsv, r2_histogram.tsv, er2_summary.tsv, er2_by_maf.tsv,
         er2_subsets.tsv, er2_low_outliers.tsv, qc_global_report.txt} and figures/{F1_r2_distribution, F2_r2_by_maf, F3_er2_genotyped}.png
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tables-dir", required=True); ap.add_argument("--keep75", required=True); ap.add_argument("--compare446", required=True)
ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(os.path.join(a.outdir, "figures"), exist_ok=True)
F = os.path.join(a.outdir, "figures")
parts = []
for c in range(1, 23):
    d = pd.read_csv(os.path.join(a.tables_dir, f"chr{c}.tsv.gz"), sep="\t", usecols=["CHROM", "POS", "REF", "ALT", "MAF", "R2", "ER2", "TYPE"],
                    dtype={"CHROM": "category", "REF": "category", "ALT": "category", "TYPE": "category"})
    parts.append(d)
t = pd.concat(parts, ignore_index=True)
t["chrn"] = t.CHROM.astype(str).str.replace("chr", "").astype(int)
L = [f"variants delivered by the server (R2 >= 0.3): {len(t):,} on 22 autosomes"]
L.append("TYPE counts: " + str(t.TYPE.value_counts().to_dict()))

# ------------------------------------------------------------------ (a)
by = t.groupby("chrn").agg(n_variants=("POS", "size"), n_typed=("TYPE", lambda s: int((s == "Typed").sum())),
                           n_imputed=("TYPE", lambda s: int((s == "Imputed").sum())), mean_R2=("R2", "mean"), median_R2=("R2", "median"),
                           pct_R2_ge_0_8=("R2", lambda s: 100 * float((s >= 0.8).mean())))
by.loc["total"] = [len(t), int((t.TYPE == "Typed").sum()), int((t.TYPE == "Imputed").sum()), t.R2.mean(), t.R2.median(), 100 * float((t.R2 >= 0.8).mean())]
by.to_csv(os.path.join(a.outdir, "summary_by_chrom.tsv"), sep="\t", float_format="%.4f")
q = t.groupby("TYPE", observed=True).R2.describe(percentiles=[.05, .25, .5, .75, .95])
L.append("R2 quantiles by TYPE:\n" + q.round(3).to_string())
L.append(f"R2 min overall {t.R2.min():.3f}; Typed R2 min {t[t.TYPE == 'Typed'].R2.min():.3f} | share R2 >= 0.5: {100 * (t.R2 >= 0.5).mean():.1f} % | >= 0.8: {100 * (t.R2 >= 0.8).mean():.1f} % | >= 0.9: {100 * (t.R2 >= 0.9).mean():.1f} %")
bins = np.round(np.arange(0.3, 1.0001, 0.05), 2)
hist = pd.DataFrame({"R2_from": bins[:-1], "R2_to": bins[1:]})
for nm, sub in [("all", t), ("imputed", t[t.TYPE == "Imputed"]), ("typed", t[t.TYPE == "Typed"])]:
    hist[nm] = np.histogram(sub.R2, bins=bins)[0]
hist.to_csv(os.path.join(a.outdir, "r2_histogram.tsv"), sep="\t", index=False)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for nm, col in [("all", "#444444"), ("imputed", "#2c7fb8"), ("typed", "#c0392b")]:
    ax[0].step(bins[:-1], hist[nm], where="post", label=nm, color=col)
    ax[1].step(bins[:-1], 100 * hist[nm] / hist[nm].sum(), where="post", label=nm, color=col)
ax[0].set_yscale("log"); ax[0].set_ylabel("variants"); ax[1].set_ylabel("% of the group")
for x in ax: x.set_xlabel("R2 (server estimate), truncated at 0.3 by the R2 >= 0.3 filter"); x.legend()
fig.suptitle("Distribution of imputation R2 (140 samples, TOPMed r3)", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(F, "F1_r2_distribution.png"), dpi=160); plt.close(fig)

# ------------------------------------------------------------------ (b)
edges = [0, 0.001, 0.005, 0.01, 0.05, 0.2, 0.5000001]
labels = ["<0.001", "0.001-0.005", "0.005-0.01", "0.01-0.05", "0.05-0.2", "0.2-0.5"]
t["maf_bin"] = pd.cut(t.MAF, edges, labels=labels, right=False)


def by_maf(d):
    g = d.groupby("maf_bin", observed=False).R2
    return pd.DataFrame({"n": g.size(), "mean_R2": g.mean(), "median_R2": g.median(),
                         "pct_R2_ge_0.3": g.apply(lambda s: 100 * float((s >= 0.3).mean()) if len(s) else np.nan),
                         "pct_R2_ge_0.8": g.apply(lambda s: 100 * float((s >= 0.8).mean()) if len(s) else np.nan)})


bm_all, bm_imp = by_maf(t), by_maf(t[t.TYPE == "Imputed"])
bm_all.to_csv(os.path.join(a.outdir, "r2_by_maf.tsv"), sep="\t", float_format="%.4f"); bm_imp.to_csv(os.path.join(a.outdir, "r2_by_maf_imputed_only.tsv"), sep="\t", float_format="%.4f")
L.append("R2 by MAF interval, ALL delivered variants:\n" + bm_all.round(3).to_string()); L.append("R2 by MAF interval, IMPUTED only:\n" + bm_imp.round(3).to_string())
fig, ax = plt.subplots(1, 2, figsize=(11, 4)); x = np.arange(len(labels))
ax[0].plot(x, bm_all.mean_R2, "o-", color="#444444", label="all"); ax[0].plot(x, bm_imp.mean_R2, "s--", color="#2c7fb8", label="imputed only")
ax[0].set_ylabel("mean R2"); ax[0].set_ylim(0.3, 1.0)
ax[1].plot(x, bm_all["pct_R2_ge_0.8"], "o-", color="#444444", label="% R2 >= 0.8 (all)"); ax[1].plot(x, bm_imp["pct_R2_ge_0.8"], "s--", color="#2c7fb8", label="% R2 >= 0.8 (imputed)")
ax[1].plot(x, bm_all["pct_R2_ge_0.3"], ":", color="#999999", label="% R2 >= 0.3 (100 % by construction)"); ax[1].set_ylabel("% of variants that pass"); ax[1].set_ylim(0, 105)
for k, x_ in enumerate(ax):
    x_.set_xticks(x); x_.set_xticklabels(labels, rotation=30); x_.set_xlabel("MAF interval (server estimate, 140 samples)"); x_.legend(fontsize=8)
    for i in range(len(labels)): x_.annotate(f"{int(bm_all.n.iloc[i]):,}", (i, (bm_all.mean_R2 if k == 0 else bm_all["pct_R2_ge_0.8"]).iloc[i]), fontsize=6, xytext=(0, 6), textcoords="offset points", ha="center")
fig.suptitle("Imputation R2 by MAF (numbers = variants in the interval; R2 < 0.3 variants were already removed by the server)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(F, "F2_r2_by_maf.png"), dpi=160); plt.close(fig)

# ------------------------------------------------------------------ (c)
ty = t[(t.TYPE == "Typed") & t.ER2.notna()].copy()
key = lambda d: d.CHROM.astype(str) + "|" + d.POS.astype(str) + "|" + d.REF.astype(str) + "|" + d.ALT.astype(str)
ty["key"] = key(ty)
k75 = pd.read_csv(a.keep75, sep="\t", dtype={"CHROM": str}); k75["key"] = k75.CHROM + "|" + k75.POS.astype(str) + "|" + k75.REF + "|" + k75.ALT
c446 = pd.read_csv(a.compare446, sep="\t", dtype={"CHROM": str}); c446["key"] = c446.CHROM + "|" + c446.POS.astype(str) + "|" + c446.REF + "|" + c446.ALT
k358 = c446[~c446.abs_diff_topmed_gt_0_2.astype(bool)]
ty["subset"] = "other genotyped"
ty.loc[ty.key.isin(set(k358.key)), "subset"] = "358 discordant kept"
ty.loc[ty.key.isin(set(k75.key)), "subset"] = "75 kept (FIN-supported)"
L.append(f"genotyped variants with ER2: {len(ty):,} of {int((t.TYPE == 'Typed').sum()):,} Typed | of the 75: {int((ty.subset == '75 kept (FIN-supported)').sum())} found, of the 358: {int((ty.subset == '358 discordant kept').sum())} found")
qs = [.01, .05, .25, .5, .75, .95]


def er2_row(d):
    return {"n": len(d), "mean": d.ER2.mean(), **{f"q{int(q * 100):02d}": d.ER2.quantile(q) for q in qs},
            "pct_ge_0.3": 100 * float((d.ER2 >= 0.3).mean()), "pct_ge_0.8": 100 * float((d.ER2 >= 0.8).mean()), "pct_ge_0.9": 100 * float((d.ER2 >= 0.9).mean()),
            "n_below_0.3": int((d.ER2 < 0.3).sum()), "n_below_0.5": int((d.ER2 < 0.5).sum())}


es = pd.DataFrame({"all genotyped": er2_row(ty), **{s: er2_row(d) for s, d in ty.groupby("subset")}}).T
es.to_csv(os.path.join(a.outdir, "er2_summary.tsv"), sep="\t", float_format="%.4f")
L.append("ER2 (leave-one-out R2) of the genotyped variants:\n" + es.round(3).to_string())
ty["maf_bin"] = pd.cut(ty.MAF, edges, labels=labels, right=False)
em = ty.groupby("maf_bin", observed=False).ER2.agg(n="size", mean="mean", median="median", pct_ge_0_8=lambda s: 100 * float((s >= 0.8).mean()) if len(s) else np.nan)
em.to_csv(os.path.join(a.outdir, "er2_by_maf.tsv"), sep="\t", float_format="%.4f")
L.append("ER2 by MAF interval (genotyped):\n" + em.round(3).to_string())
ty[["subset"]].value_counts().to_csv(os.path.join(a.outdir, "er2_subsets.tsv"), sep="\t")
low = ty[ty.ER2 < 0.5][["CHROM", "POS", "REF", "ALT", "MAF", "R2", "ER2", "subset"]].sort_values("ER2")
low.to_csv(os.path.join(a.outdir, "er2_low_outliers.tsv"), sep="\t", index=False, float_format="%.4g")
L.append(f"genotyped variants with ER2 < 0.5: {len(low):,} (of which in the 75: {int((low.subset == '75 kept (FIN-supported)').sum())}, in the 358: {int((low.subset == '358 discordant kept').sum())})")
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
b2 = np.linspace(0, 1, 41)
for s, col in [("other genotyped", "#7f8c8d"), ("358 discordant kept", "#2c7fb8"), ("75 kept (FIN-supported)", "#c0392b")]:
    d = ty[ty.subset == s].ER2
    if len(d):
        ax[0].hist(d, bins=b2, density=True, histtype="step", color=col, label=f"{s} (n={len(d):,})", lw=1.6)
        ax[1].plot(np.sort(d), np.arange(1, len(d) + 1) / len(d), color=col, label=s)
ax[0].set_xlabel("ER2 (leave-one-out R2)"); ax[0].set_ylabel("density"); ax[0].set_yscale("log"); ax[0].legend(fontsize=7)
ax[1].set_xlabel("ER2"); ax[1].set_ylabel("cumulative fraction"); ax[1].legend(fontsize=7)
ax[2].plot(np.arange(len(labels)), em["mean"], "o-", color="#444444"); ax[2].set_xticks(np.arange(len(labels))); ax[2].set_xticklabels(labels, rotation=30)
ax[2].set_xlabel("MAF interval"); ax[2].set_ylabel("mean ER2 (genotyped)"); ax[2].set_ylim(0, 1)
fig.suptitle("Leave-one-out accuracy (ER2) of the genotyped variants", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(F, "F3_er2_genotyped.png"), dpi=160); plt.close(fig)
open(os.path.join(a.outdir, "qc_global_report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
