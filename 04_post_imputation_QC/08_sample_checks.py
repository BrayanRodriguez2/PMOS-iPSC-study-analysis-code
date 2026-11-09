#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/08_sample_checks.py   (original step script: 43_sample_checks.py)
# Produces: Concordance, heterozygosity and outlier checks of the samples after imputation.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 43 - Post-imputation QC, part (e): sample behaviour in the imputed data (all 140 samples). Diagnostic only.

Sub-commands
  concord : one chromosome. Compares, at the GENOTYPED (Typed) sites, the original genotype (v2 upload VCF, GT) with the imputed dosage
            (dose VCF, DS rounded to 0/1/2). Per sample: sites compared, concordant, original het -> imputed hom, original hom -> imputed
            het, original hom -> the other hom, and sum |DS - original dosage|.
            Note: at typed sites the imputed dose comes from the phased haplotypes (Eagle) + reference panel, so it is mostly the observed
            genotype smoothed by the panel; a contaminated or noisy array sample shows MORE discordance (excess het calls the panel does not
            support) - this is exactly what is being checked. The server's per-variant leave-one-out ER2 is used in step 41.
  het     : heterozygosity of the best-guess imputed genotypes (GT of the dose VCF) at variants with R2 >= 0.8 and MAF >= 0.05 in 66 windows of
            300 kb (same windows as step 42): per sample het rate and inbreeding coefficient F = (E[het] - O[het]) / E[het] with E[het] = sum 2p(1-p)
            using the allele frequency p of the 129 QC samples; robust z-scores (median/MAD) over the 129. The windows are contiguous (LD), so F is
            a relative measure, comparable between samples but not to the F of the LD-pruned array data (step 09) in absolute terms.
  merge   : joins everything with the sample QC table (step 09), writes the table, the figures and the report: do the 11 excluded samples
            (especially the 5 contaminated samples) stand out, and does any other sample become an outlier?
Outputs: results/imputation_qc/{sample_checks_imputed.tsv, sample_checks_report.txt}, figures/{F5_sample_concordance, F6_het_imputed_vs_array}.png
"""
import argparse, glob, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)
c = sub.add_parser("concord"); c.add_argument("--orig", required=True); c.add_argument("--dose", required=True); c.add_argument("--samples", required=True); c.add_argument("--out", required=True)
h = sub.add_parser("het"); h.add_argument("--gt", required=True); h.add_argument("--samples", required=True); h.add_argument("--keep129", required=True); h.add_argument("--out", required=True)
m = sub.add_parser("merge"); m.add_argument("--counts-dir", required=True); m.add_argument("--het", required=True); m.add_argument("--qc-table", required=True); m.add_argument("--outdir", required=True)
a = ap.parse_args()
GT2D = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/0": 1, "1/1": 2, "1|1": 2}

if a.cmd == "concord":
    samples = [l.strip() for l in open(a.samples) if l.strip()]
    o = pd.read_csv(a.orig, sep="\t", header=None, dtype=str); o.columns = ["POS", "REF", "ALT"] + samples
    d = pd.read_csv(a.dose, sep="\t", header=None, dtype=str); d.columns = ["POS", "REF", "ALT"] + samples
    o["k"] = o.POS + "|" + o.REF + "|" + o.ALT; d["k"] = d.POS + "|" + d.REF + "|" + d.ALT
    o = o.drop_duplicates("k").set_index("k"); d = d.drop_duplicates("k").set_index("k")
    common = o.index.intersection(d.index)
    O = np.vstack([o.loc[common, s].map(GT2D).fillna(-1).to_numpy(dtype=int) for s in samples]).T
    D = np.vstack([pd.to_numeric(d.loc[common, s], errors="coerce").to_numpy() for s in samples]).T
    Dr = np.clip(np.rint(D), 0, 2); ok = (O >= 0) & ~np.isnan(D)
    rows = []
    for j, s in enumerate(samples):
        v = ok[:, j]; Oj, Dj, Drj = O[v, j], D[v, j], Dr[v, j]
        rows.append(dict(IID=s, n_typed_in_both=len(common), n_compared=int(v.sum()), n_orig_missing=int((O[:, j] < 0).sum()), n_concordant=int((Oj == Drj).sum()),
                         n_orig_het_imp_hom=int(((Oj == 1) & (Drj != 1)).sum()), n_orig_hom_imp_het=int(((Oj != 1) & (Drj == 1)).sum()),
                         n_orig_hom_imp_otherhom=int((((Oj == 0) & (Drj == 2)) | ((Oj == 2) & (Drj == 0))).sum()),
                         n_orig_het=int((Oj == 1).sum()), sum_abs_err=float(np.abs(Dj - Oj).sum())))
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    print(f"{os.path.basename(a.out)}: {len(common):,} typed sites in both files")

elif a.cmd == "het":
    samples = [l.strip() for l in open(a.samples) if l.strip()]
    keep = set(l.strip() for l in open(a.keep129) if l.strip())
    g = pd.read_csv(a.gt, sep="\t", header=None, dtype=str)
    G = np.vstack([g[j].map(GT2D).fillna(-1).to_numpy(dtype=int) for j in range(g.shape[1])]).T      # variants x samples
    m129 = np.array([s in keep for s in samples])
    called = G >= 0
    p = np.where(called[:, m129].sum(1) > 0, np.where(called[:, m129], G[:, m129], 0).sum(1) / (2 * np.maximum(called[:, m129].sum(1), 1)), np.nan)
    use = (np.minimum(p, 1 - p) >= 0.05) & ~np.isnan(p)
    G, called, p = G[use], called[use], p[use]
    eh = 2 * p * (1 - p)
    rows = []
    for j, s in enumerate(samples):
        v = called[:, j]; ohet = int((G[v, j] == 1).sum()); ehet = float(eh[v].sum())
        rows.append(dict(IID=s, n_variants=int(v.sum()), het_rate_imputed=ohet / max(v.sum(), 1), F_imputed=(ehet - ohet) / ehet if ehet else np.nan))
    t = pd.DataFrame(rows)
    ref = t[t.IID.isin(keep)]
    for col in ["het_rate_imputed", "F_imputed"]:
        med = ref[col].median(); mad = 1.4826 * (ref[col] - med).abs().median(); t[col + "_robust_z"] = (t[col] - med) / mad
    t.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    print(f"het: {len(t)} samples, {int(use.sum()):,} variants (R2 >= 0.8, imputed MAF >= 0.05 in the 129)")

else:
    cnt = pd.concat([pd.read_csv(f, sep="\t", dtype={"IID": str}) for f in sorted(glob.glob(os.path.join(a.counts_dir, "concord_chr*.tsv")))])
    s = cnt.groupby("IID").sum(numeric_only=True)
    s["concordance"] = s.n_concordant / s.n_compared
    s["het_to_hom_rate"] = s.n_orig_het_imp_hom / s.n_orig_het
    s["hom_to_het_rate"] = s.n_orig_hom_imp_het / (s.n_compared - s.n_orig_het)
    s["mean_abs_err"] = s.sum_abs_err / s.n_compared
    het = pd.read_csv(a.het, sep="\t", dtype={"IID": str}).set_index("IID")
    qc = pd.read_csv(a.qc_table, sep="\t", dtype={"IID": str}).set_index("IID")
    t = s[["n_compared", "concordance", "het_to_hom_rate", "hom_to_het_rate", "mean_abs_err"]].join(het[["het_rate_imputed", "F_imputed", "F_imputed_robust_z"]]).join(
        qc[["call_rate", "het_F", "excluded", "reasons", "is_donor_of_ipsc_line"]])
    ref = t[t.excluded == "no"]
    for col in ["concordance", "het_to_hom_rate", "hom_to_het_rate"]:
        med = ref[col].median(); mad = 1.4826 * (ref[col] - med).abs().median(); t[col + "_robust_z"] = (t[col] - med) / mad
    t = t.reset_index().rename(columns={"index": "IID"})
    t.to_csv(os.path.join(a.outdir, "sample_checks_imputed.tsv"), sep="\t", index=False, float_format="%.6g")
    cont = [x for x in __import__("os").environ.get("CONTAMINATED_SAMPLE_IDS", "").split(",") if x]
    L = [f"samples: {len(t)} | typed sites compared per sample: median {int(t.n_compared.median()):,}",
         f"concordance (rounded imputed dose = original genotype at typed sites): 129 QC-pass median {ref.concordance.median():.5f} (min {ref.concordance.min():.5f}); "
         f"11 excluded median {t[t.excluded == 'yes'].concordance.median():.5f} (min {t[t.excluded == 'yes'].concordance.min():.5f})",
         f"original het -> imputed hom rate: 129 QC-pass median {ref.het_to_hom_rate.median():.5f}; 5 contaminated {[round(t[t.IID == i].het_to_hom_rate.iloc[0], 5) for i in cont]}",
         f"F in the imputed windows: 129 QC-pass median {ref.F_imputed.median():.4f}; excluded: " + ", ".join(f"{r.IID}={r.F_imputed:.3f} (z {r.F_imputed_robust_z:.1f})" for r in t[t.excluded == 'yes'].itertuples())]
    for col, thr, name in [("F_imputed_robust_z", 3, "F (imputed)"), ("concordance_robust_z", -3, "concordance"), ("het_to_hom_rate_robust_z", 3, "het->hom rate")]:
        out = t[(t.excluded == "no") & ((t[col] < thr) if thr < 0 else (t[col].abs() > thr))]
        L.append(f"QC-pass samples that are outliers by {name} (robust z {'<' if thr < 0 else '|z| >'} {abs(thr)}): " + (", ".join(f"{r.IID} ({getattr(r, col):.1f})" for r in out.itertuples()) or "none"))
    for col, thr, name in [("F_imputed_robust_z", 3, "F (imputed)"), ("concordance_robust_z", -3, "concordance")]:
        ex = t[t.excluded == "yes"]; flagged = ex[(ex[col] < thr) if thr < 0 else (ex[col].abs() > thr)]
        L.append(f"excluded samples flagged by {name}: {', '.join(flagged.IID) or 'none'} ({len(flagged)}/11)")
    L.append("lowest 8 by concordance:\n" + t.sort_values("concordance").head(8)[["IID", "concordance", "het_to_hom_rate", "F_imputed", "het_F", "call_rate", "excluded", "reasons"]].round(5).to_string(index=False))
    open(os.path.join(a.outdir, "sample_checks_report.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    os.makedirs(os.path.join(a.outdir, "figures"), exist_ok=True)
    col = np.where(t.IID.isin(cont), "#c0392b", np.where(t.excluded == "yes", "#e67e22", "#2c7fb8"))
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    o = t.sort_values("concordance").reset_index(drop=True)
    cc = np.where(o.IID.isin(cont), "#c0392b", np.where(o.excluded == "yes", "#e67e22", "#2c7fb8"))
    ax[0].scatter(range(len(o)), o.concordance, c=cc, s=14); ax[0].set_ylabel("concordance at genotyped sites (original GT vs rounded dose)"); ax[0].set_xlabel("samples sorted by concordance")
    for i, r in o.head(8).iterrows(): ax[0].annotate(r.IID, (i, r.concordance), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax[1].scatter(t.het_to_hom_rate * 100, t.concordance, c=col, s=14); ax[1].set_xlabel("original het -> imputed hom (% of original hets)"); ax[1].set_ylabel("concordance")
    for r in t[t.excluded == "yes"].itertuples(): ax[1].annotate(r.IID, (r.het_to_hom_rate * 100, r.concordance), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax[2].scatter(t.het_F, t.F_imputed, c=col, s=14); ax[2].axhline(0, color="#999999", lw=.5); ax[2].axvline(0, color="#999999", lw=.5)
    ax[2].set_xlabel("F in the original array data (step 09)"); ax[2].set_ylabel("F in the imputed best-guess genotypes (windows, R2 >= 0.8)")
    for r in t[t.excluded == "yes"].itertuples(): ax[2].annotate(r.IID, (r.het_F, r.F_imputed), fontsize=7, xytext=(3, 3), textcoords="offset points")
    fig.suptitle("Samples after imputation (red = 5 contaminated samples; orange = other excluded; blue = QC-pass)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "figures", "F5_sample_concordance.png"), dpi=160); plt.close(fig)
