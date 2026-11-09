#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/12_summarize_prs_array.py   (original step script: 16_summarize_prs.py)
# Produces: Standardised scores, correlations between versions and figures (array PRS).
# Run: python 05_PRS/12_summarize_prs_array.py <args>   (or from a SLURM script; see the usage line below)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 16c - standardise the PRS, compare versions and draw the descriptive figures.

Inputs (from `plink2 --score`, step 16 .sbatch): <scoredir>/V1.sscore, V2.sscore, V3.sscore
  columns used: IID, ALLELE_CT, SCORE1_SUM  (SUM of weight x allele dosage; missing genotypes were
  mean-imputed with allele frequencies of the QC cohort, see the .sbatch)
Standardisation: z = (score - mean) / SD, where mean and SD are computed ONLY over the samples that
pass sample QC (step 09), so that the reference distribution is the QC cohort. Samples that failed QC
are still scored (same weights, same z formula) but they are always labelled "QC_fail" and never
enter the mean/SD, the correlations or the histogram bars.
Comparisons: Pearson and Spearman correlations between V1, V2 and V3 over the QC-pass samples.
The 14 donors of the iPSC lines (donor_key.tsv; the IID of the .fam is the eSF ID) are marked by group
(PCOS / CTRL) on the histograms. THIS IS DESCRIPTIVE ONLY: no association test is done (n = 14).
SNP counts per version: number of SNPs PRS-CS saw in the reference & sumstats & bim (from its logs),
number of weights written, and number of variants that were scored.
Outputs (--outdir): prs_scores.tsv, donors_14_prs.tsv, snps_per_version.tsv, prs_correlations.tsv,
prs_distribution.png/.pdf, prs_version_scatter.png/.pdf, prs_summary_report.txt
Needs numpy, scipy and matplotlib.
"""
import argparse
import glob
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

VERSIONS = ["V1", "V2", "V3"]
LABEL = {"V1": "V1 main: UKBB panel + QC array", "V2": "V2 sensitivity: 1000G panel + QC array",
         "V3": "V3 article-like: UKBB panel, all HapMap3 SNPs; weights applied to array"}
COLOR = {"PCOS": "#D55E00", "CTRL": "#0072B2"}  # colour-blind-safe pair


def read_table(path):
    with open(path) as fh:
        h = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in fh]


def count_lines(path):
    with open(path) as fh:
        return sum(1 for _ in fh)


def prscs_common_snps(logdir):
    """Sum over chromosomes of 'N common SNPs in the reference, sumstats, and validation set' per version
    (latest log per (version, chromosome))."""
    best = {}
    for fn in glob.glob(os.path.join(logdir, "prscs_*.out")):
        txt = open(fn, errors="replace").read()
        m = re.search(r"### (V\d) chr(\d+)", txt)
        c = re.search(r"\.\.\. (\d+) common SNPs in the reference, sumstats, and validation set", txt)
        if m and c and "Traceback" not in txt:
            key = (m.group(1), int(m.group(2)))
            if key not in best or os.path.getmtime(fn) > best[key][0]:
                best[key] = (os.path.getmtime(fn), int(c.group(1)))
    out = {}
    for (v, ch), (_, n) in best.items():
        out.setdefault(v, {})[ch] = n
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scoredir", required=True)
    ap.add_argument("--qc-table", required=True)
    ap.add_argument("--donor-key", required=True)
    ap.add_argument("--logdir", required=True)
    ap.add_argument("--weightdir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rep = []

    def log(m=""):
        print(m, flush=True)
        rep.append(m)

    qc = {r["IID"]: r for r in read_table(args.qc_table)}
    donors = {r["esf_id"]: r for r in read_table(args.donor_key)}
    raw = {}
    for v in VERSIONS:
        for r in read_table(os.path.join(args.scoredir, f"{v}.sscore")):
            raw.setdefault(r["IID"], {})[v] = float(r["SCORE1_SUM"])
    iids = sorted(raw, key=int)
    passing = [i for i in iids if qc[i]["excluded"] == "no"]
    log(f"[samples] scored: {len(iids)} | pass sample QC: {len(passing)} | fail: {len(iids) - len(passing)}")

    z = {v: {} for v in VERSIONS}
    for v in VERSIONS:
        ref = np.array([raw[i][v] for i in passing])
        mu, sd = ref.mean(), ref.std(ddof=1)
        for i in iids:
            z[v][i] = (raw[i][v] - mu) / sd
        log(f"[{v}] QC cohort: mean {mu:.6g} | SD {sd:.6g} | min {ref.min():.6g} | max {ref.max():.6g}")

    with open(os.path.join(args.outdir, "prs_scores.tsv"), "w") as out:
        out.write("IID\tqc_status\tqc_reasons\tdonor_ipsc_line\tgroup\t" + "\t".join(f"{v}_sum" for v in VERSIONS) + "\t"
                  + "\t".join(f"{v}_z" for v in VERSIONS) + "\n")
        for i in iids:
            d = donors.get(i)
            out.write(f"{i}\t{'pass' if qc[i]['excluded'] == 'no' else 'QC_fail'}\t{qc[i]['reasons']}\t"
                      f"{d['ipsc_line'] if d else ''}\t{d['group'] if d else ''}\t"
                      + "\t".join(f"{raw[i][v]:.6g}" for v in VERSIONS) + "\t" + "\t".join(f"{z[v][i]:.4f}" for v in VERSIONS) + "\n")
    with open(os.path.join(args.outdir, "donors_14_prs.tsv"), "w") as out:
        out.write("ipsc_line\tesf_id\tgroup\tqc_status\t" + "\t".join(f"{v}_z" for v in VERSIONS) + "\t"
                  + "\t".join(f"{v}_percentile_in_QC_cohort" for v in VERSIONS) + "\n")
        for esf, d in sorted(donors.items(), key=lambda x: x[1]["ipsc_line"]):
            if esf not in raw:
                out.write(f"{d['ipsc_line']}\t{esf}\t{d['group']}\tnot_genotyped\n")
                continue
            pct = []
            for v in VERSIONS:
                ref = np.array([z[v][i] for i in passing])
                pct.append(100.0 * (ref < z[v][esf]).mean())
            out.write(f"{d['ipsc_line']}\t{esf}\t{d['group']}\t{'pass' if qc[esf]['excluded'] == 'no' else 'QC_fail'}\t"
                      + "\t".join(f"{z[v][esf]:.3f}" for v in VERSIONS) + "\t" + "\t".join(f"{p:.1f}" for p in pct) + "\n")

    # ---- SNP counts --------------------------------------------------------------------
    common = prscs_common_snps(args.logdir)
    with open(os.path.join(args.outdir, "snps_per_version.tsv"), "w") as out:
        out.write("version\tchromosomes_with_prscs_log\tsnps_seen_by_prscs\tweights_written_for_array\tvariants_scored_plink2\n")
        for v in VERSIONS:
            cm = common.get(v, {})
            wfile = os.path.join(args.weightdir, f"{v}_weights_array.txt")
            nw = count_lines(wfile) - 1 if os.path.exists(wfile) else -1
            vs = os.path.join(args.scoredir, f"{v}.sscore.vars")
            ns = count_lines(vs) if os.path.exists(vs) else -1
            out.write(f"{v}\t{len(cm)}\t{sum(cm.values())}\t{nw}\t{ns}\n")
            log(f"[SNPs {v}] PRS-CS saw {sum(cm.values()):,} SNPs over {len(cm)} chromosomes | weights written for the array: {nw:,} "
                f"| variants scored by plink2: {ns:,}")

    # ---- correlations -----------------------------------------------------------------------
    pairs = [("V1", "V2"), ("V1", "V3"), ("V2", "V3")]
    with open(os.path.join(args.outdir, "prs_correlations.tsv"), "w") as out:
        out.write("pair\tn_samples\tpearson_r\tspearman_rho\n")
        for a, b in pairs:
            x = np.array([raw[i][a] for i in passing])
            y = np.array([raw[i][b] for i in passing])
            r, rho = stats.pearsonr(x, y)[0], stats.spearmanr(x, y)[0]
            out.write(f"{a}-{b}\t{len(passing)}\t{r:.4f}\t{rho:.4f}\n")
            log(f"[corr {a}-{b}] n={len(passing)} | Pearson r = {r:.3f} | Spearman rho = {rho:.3f}")

    # ---- figures ---------------------------------------------------------------------------------
    don_ids = [i for i in donors if i in raw]
    # One block per version: histogram of the QC cohort on top and, below it on the same x axis, a
    # dot plot with one row per donor (sorted by PRS) so that the 14 donors never overlap.
    fig = plt.figure(figsize=(9, 15.5))
    outer = fig.add_gridspec(3, 1, hspace=0.30, top=0.94, bottom=0.07)
    for k_v, v in enumerate(VERSIONS):
        inner = outer[k_v].subgridspec(2, 1, height_ratios=[1.6, 2.4], hspace=0.08)
        ax = fig.add_subplot(inner[0])
        ax2 = fig.add_subplot(inner[1], sharex=ax)
        vals = np.array([z[v][i] for i in passing])
        ax.hist(vals, bins=20, color="#BBBBBB", edgecolor="white")
        ax.set_ylabel("donors (QC cohort)")
        ax.set_title(f"{LABEL[v]}\n(grey bars: {len(passing)} donors passing sample QC)", fontsize=9.5, loc="left")
        plt.setp(ax.get_xticklabels(), visible=False)
        order = sorted(don_ids, key=lambda i: z[v][i])
        for row, i in enumerate(order):
            d, ok = donors[i], qc[i]["excluded"] == "no"
            ax2.plot(z[v][i], row, marker="o" if ok else "D", ms=8, color=COLOR[d["group"]],
                     mfc=COLOR[d["group"]] if ok else "white", mew=1.6, ls="none")
        ax2.set_yticks(range(len(order)))
        ax2.set_yticklabels([donors[i]["ipsc_line"] + ("" if qc[i]["excluded"] == "no" else " *") for i in order], fontsize=7.5)
        for tick, i in zip(ax2.get_yticklabels(), order):
            tick.set_color(COLOR[donors[i]["group"]])
        ax2.axvline(0, color="#999999", lw=0.8, ls=":")
        ax2.set_ylim(-0.7, len(order) - 0.3)
        ax2.grid(axis="y", color="#EEEEEE", lw=0.6)
        ax2.set_xlabel("PRS (z-score, mean 0 / SD 1 in the QC cohort)")
    handles = [plt.Line2D([], [], marker="o", ls="none", color=COLOR["PCOS"], label="PCOS donor (passes QC)"),
               plt.Line2D([], [], marker="o", ls="none", color=COLOR["CTRL"], label="CTRL donor (passes QC)"),
               plt.Line2D([], [], marker="D", ls="none", mfc="white", mew=1.6, color="#444444",
                          label="* donor failing sample QC (hollow diamond; colour = group)")]
    fig.legend(handles=handles, loc="lower center", ncol=1, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("PCOS PRS (PRS-CS, phi = 1e-4) in the donor cohort - descriptive, no association test", fontsize=11, y=0.992)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.outdir, f"prs_distribution.{ext}"), dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
    for ax, (a, b) in zip(axes, pairs):
        ax.scatter([z[a][i] for i in passing], [z[b][i] for i in passing], s=12, color="#999999")
        for i in don_ids:
            d, ok = donors[i], qc[i]["excluded"] == "no"
            ax.scatter(z[a][i], z[b][i], s=40, marker="o" if ok else "D", facecolors=COLOR[d["group"]] if ok else "white",
                       edgecolors=COLOR[d["group"]], linewidths=1.5)
        r = stats.pearsonr([raw[i][a] for i in passing], [raw[i][b] for i in passing])[0]
        ax.set_xlabel(f"{a} (z)")
        ax.set_ylabel(f"{b} (z)")
        ax.set_title(f"{a} vs {b}: Pearson r = {r:.3f} (n={len(passing)})", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.outdir, f"prs_version_scatter.{ext}"), dpi=200)
    plt.close(fig)
    with open(os.path.join(args.outdir, "prs_summary_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
