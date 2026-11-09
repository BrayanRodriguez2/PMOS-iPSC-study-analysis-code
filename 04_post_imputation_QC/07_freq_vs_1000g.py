#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/07_freq_vs_1000g.py   (original step script: 42_freq_vs_1000g.py)
# Produces: Imputed allele frequencies against 1000 Genomes EUR/FIN in sampled windows.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 42 - Post-imputation QC, part (d): imputed allele frequency (the 129 QC-pass samples) versus 1000 Genomes EUR (503, includes FIN) and FIN
(99), in common variants. Diagnostic only.

Design (documented limitation): reading all ~11.7 M imputed variants for 129 samples and the 1000G genotypes of all of them is not needed to
see whether the imputed frequencies are sound, so a SAMPLE of the genome is used: 3 windows of 300 kb per autosome (centred at the 20 %, 50 %
and 80 % quantiles of the variant positions of that chromosome) = 66 windows, ~20 Mb. In each window: imputed allele frequency =
mean(DS)/2 over the 129 samples (DS = imputed dosage of the dose VCF), and the 1000G allele frequency of the SAME allele computed from
the genotypes of the 1000G samples (EUR: 503; FIN: 99) with allele matching (same pair, or swapped -> 1 - AF; otherwise no comparison).
"Common" = imputed MAF >= 0.05 in the 129. Reported: n, Pearson r, mean difference, number with |diff| > 0.2, for all common variants and
for R2 >= 0.8; separately for genotyped (Typed) and imputed variants. The 1000G frequencies come from a different, larger and
multi-country European sample, so small differences are expected (FIN differs most because it is only 99 individuals).
Inputs : windows (made here from info_tables), dose_windows.tsv (bcftools query of the 129 samples), g1000_windows_gt.tsv (+ g1000_samples.txt,
         1kg.panel)
Outputs: results/imputation_qc/{freq_vs_1000g_variants.tsv, freq_vs_1000g_summary.tsv, freq_vs_1000g_report.txt}, figures/F4_freq_vs_1000g.png
Sub-commands: windows (writes the window list) | compare (needs the two extracted tables)
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)
w = sub.add_parser("windows"); w.add_argument("--tables-dir", required=True); w.add_argument("--out", required=True); w.add_argument("--width", type=int, default=300_000)
c = sub.add_parser("compare")
c.add_argument("--dose", required=True); c.add_argument("--g1000", required=True); c.add_argument("--g1000-samples", required=True)
c.add_argument("--panel", required=True); c.add_argument("--outdir", required=True)
a = ap.parse_args()

if a.cmd == "windows":
    rows = []
    for ch in range(1, 23):
        pos = pd.read_csv(os.path.join(a.tables_dir, f"chr{ch}.tsv.gz"), sep="\t", usecols=["POS"]).POS.to_numpy()
        for q in (0.2, 0.5, 0.8):
            mid = int(np.quantile(pos, q)); rows.append((f"chr{ch}", max(1, mid - a.width // 2), mid + a.width // 2))
    pd.DataFrame(rows).to_csv(a.out, sep="\t", header=False, index=False)
    print(f"{len(rows)} windows written")
    raise SystemExit(0)

GT2D = {"0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2, "0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}
d = pd.read_csv(a.dose, sep="\t", header=None, dtype=str)
nd = d.shape[1] - 9
d.columns = ["CHROM", "POS", "ID", "REF", "ALT", "R2", "MAF", "IMPUTED", "TYPED"] + list(range(nd))
for col in ["POS", "R2", "MAF"]: d[col] = pd.to_numeric(d[col])
ds = d[list(range(nd))].apply(pd.to_numeric, errors="coerce").to_numpy()
d["af_imputed_129"] = np.nanmean(ds, axis=1) / 2
d["type"] = np.where(d.TYPED == "1", "Typed", "Imputed")
d = d[["CHROM", "POS", "ID", "REF", "ALT", "R2", "MAF", "type", "af_imputed_129"]].copy()
d["maf_129"] = np.minimum(d.af_imputed_129, 1 - d.af_imputed_129)
print(f"imputed variants in the windows: {len(d):,} (129 samples: {nd})", flush=True)
samples = [l.strip() for l in open(a.g1000_samples) if l.strip()]
panel = pd.read_csv(a.panel, sep="\t").set_index("sample")
fin = [i for i, s in enumerate(samples) if panel.loc[s, "pop"] == "FIN"]; eur = [i for i, s in enumerate(samples) if panel.loc[s, "super_pop"] == "EUR"]
g = pd.read_csv(a.g1000, sep="\t", header=None, dtype=str); g.columns = ["CHROM", "POS", "REF", "ALT"] + samples
g["CHROM"] = "chr" + g.CHROM.str.replace("chr", "", regex=False); g["POS"] = g.POS.astype(int)
G = np.vstack([g[s].map(GT2D).fillna(np.nan).to_numpy(dtype=float) for s in samples]).T


def af(rows, cols):
    x = G[rows][:, cols]; return np.nanmean(x, axis=1) / 2


g["af_eur_alt"] = af(slice(None), eur); g["af_fin_alt"] = af(slice(None), fin)
g["k"] = g.CHROM + "|" + g.POS.astype(str)
idx = {}
for i, k in enumerate(g.k): idx.setdefault(k, []).append(i)
e_af, f_af, st = [], [], []
for r in d.itertuples():
    hit = None
    for i in idx.get(f"{r.CHROM}|{r.POS}", []):
        r0, a0 = g.REF.iat[i], g.ALT.iat[i]
        if (r0, a0) == (r.REF, r.ALT): hit = (g.af_eur_alt.iat[i], g.af_fin_alt.iat[i]); break
        if (r0, a0) == (r.ALT, r.REF): hit = (1 - g.af_eur_alt.iat[i], 1 - g.af_fin_alt.iat[i]); break
    e_af.append(hit[0] if hit else np.nan); f_af.append(hit[1] if hit else np.nan); st.append("matched" if hit else "no_1000G_match")
d["af_eur"], d["af_fin"], d["g1000"] = e_af, f_af, st
d["diff_eur"] = d.af_imputed_129 - d.af_eur; d["diff_fin"] = d.af_imputed_129 - d.af_fin
d.to_csv(os.path.join(a.outdir, "freq_vs_1000g_variants.tsv"), sep="\t", index=False, float_format="%.5g")
L = [f"variants in the 66 windows: {len(d):,} | matched to a 1000G biallelic SNP with the same alleles: {int((d.g1000 == 'matched').sum()):,}"]
rows = []
for nm, m in [("common, all", d.maf_129 >= 0.05), ("common, R2 >= 0.8", (d.maf_129 >= 0.05) & (d.R2 >= 0.8)),
              ("common, Typed", (d.maf_129 >= 0.05) & (d.type == "Typed")), ("common, Imputed", (d.maf_129 >= 0.05) & (d.type == "Imputed")),
              ("common, Imputed R2 0.3-0.8", (d.maf_129 >= 0.05) & (d.type == "Imputed") & (d.R2 < 0.8))]:
    s = d[m & (d.g1000 == "matched")]
    for ref, col, dc in [("EUR", "af_eur", "diff_eur"), ("FIN", "af_fin", "diff_fin")]:
        rows.append(dict(set=nm, reference=ref, n=len(s), pearson_r=float(np.corrcoef(s.af_imputed_129, s[col])[0, 1]) if len(s) > 2 else np.nan,
                         mean_diff=s[dc].mean(), mean_abs_diff=s[dc].abs().mean(), n_abs_diff_gt_0_2=int((s[dc].abs() > 0.2).sum()),
                         pct_abs_diff_gt_0_2=100 * float((s[dc].abs() > 0.2).mean()) if len(s) else np.nan))
sm = pd.DataFrame(rows); sm.to_csv(os.path.join(a.outdir, "freq_vs_1000g_summary.tsv"), sep="\t", index=False, float_format="%.4f")
L.append(sm.round(4).to_string(index=False))
big = d[(d.maf_129 >= 0.05) & (d.g1000 == "matched") & ((d.diff_eur.abs() > 0.2) | (d.diff_fin.abs() > 0.2))]
L.append(f"common variants with |diff| > 0.2 vs EUR or FIN: {len(big)} (Typed {int((big.type == 'Typed').sum())}, Imputed {int((big.type == 'Imputed').sum())}); R2 median {big.R2.median() if len(big) else float('nan'):.3f}")
fig, ax = plt.subplots(1, 2, figsize=(11, 5))
cm = (d.maf_129 >= 0.05) & (d.g1000 == "matched")
for x, (col, nm) in zip(ax, [("af_eur", "1000G EUR (503)"), ("af_fin", "1000G FIN (99)")]):
    s = d[cm]; r = np.corrcoef(s.af_imputed_129, s[col])[0, 1]
    x.scatter(s[col], s.af_imputed_129, s=2, alpha=.25, color="#2c7fb8"); x.plot([0, 1], [0, 1], "k-", lw=.6)
    x.plot([0, 0.8], [0.2, 1], "r:", lw=.6); x.plot([0.2, 1], [0, 0.8], "r:", lw=.6)
    x.set_xlabel(f"allele frequency in {nm}"); x.set_ylabel("imputed allele frequency (129 QC samples)")
    x.set_title(f"n = {len(s):,} common variants | r = {r:.4f} | |diff| > 0.2: {int(((s.af_imputed_129 - s[col]).abs() > 0.2).sum())}", fontsize=9)
fig.suptitle("Imputed frequency vs 1000 Genomes in 66 windows of 300 kb (red dotted = |diff| 0.2)", fontsize=10)
fig.tight_layout(); os.makedirs(os.path.join(a.outdir, "figures"), exist_ok=True); fig.savefig(os.path.join(a.outdir, "figures", "F4_freq_vs_1000g.png"), dpi=160)
open(os.path.join(a.outdir, "freq_vs_1000g_report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
