#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/06_qc_by_chip.py   (original step script: 09b_qc_by_chip.py)
# Produces: Call rate and heterozygosity summarised by array chip.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 09b - sample call rate and autosomal heterozygosity (F) summarised by genotyping CHIP
(SentrixBarcode_A of the original sample sheet AKurg_06.10.22_b38.csv), with one figure.

Nothing is recomputed: the per-sample metrics are those of step 09 (results/qc/samples/sample_qc_table.tsv);
this script only joins them to the sample sheet. Join key: Sample_ID of the sheet == IID of the .fam
(the sheet has 140 rows with Sample_ID 1-141 without 37, the same IDs as the .fam).

Outputs (--outdir):
  chip_sample_map.tsv    every sample: chip, plate well, chip position, call rate, het F, excluded, reasons, donor
  chip_summary.tsv       one row per chip (n, Sample_ID range, excluded, call-rate and F summaries)
  chip_tests.txt         exploratory tests (Kruskal-Wallis across chips; Fisher exact and Mann-Whitney for the chip
                         with the most exclusions vs all other chips)
  qc_by_chip.png/.pdf    call rate and heterozygosity F by chip, thresholds drawn, excluded samples and the 14
                         iPSC donors labelled (same symbols as the other QC figures)
CAVEAT (also in chip_tests.txt): chips hold CONTIGUOUS blocks of Sample_ID (e.g. consecutive IDs on one chip), so a chip
effect cannot be separated from a plate/collection/extraction-batch effect; with 140 samples and 7 chips (one of
only 4 samples) the tests are exploratory, not confirmatory.
Needs numpy, scipy, matplotlib.
"""
import argparse
import csv
import os
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

MIN_CALL_RATE, HET_SD = 0.98, 3.0
COLOR = {"PCOS": "#D55E00", "CTRL": "#0072B2"}
GREY = "#9A9A9A"


def read_tsv(path):
    with open(path) as fh:
        h = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in fh]


def read_sheet(path):
    """Rows of the [Data] section: Sample_ID -> (plate well, chip barcode, chip position)."""
    out, header = {}, None
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if header is None:
                if row and row[0] == "Sample_ID":
                    header = row
                continue
            if not row or not row[0].strip():
                continue
            i = {k: header.index(k) for k in ("Sample_ID", "Sample_Well", "SentrixBarcode_A", "SentrixPosition_A")}
            out[row[i["Sample_ID"]].strip()] = (row[i["Sample_Well"]], row[i["SentrixBarcode_A"]], row[i["SentrixPosition_A"]])
    return out


def short_reasons(reasons):
    out = []
    for t in reasons.split(";"):
        if t.startswith("low_call_rate"):
            out.append("call rate")
        elif t.startswith("heterozygosity"):
            out.append("het")
        elif t.startswith("related_to"):
            out.append("kinship")
        elif "sex" in t or "female" in t:
            out.append("sex")
    return " + ".join(out)


def spread(ys, gap):
    """Push label heights apart (ascending order) so that consecutive labels are >= gap apart."""
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    out = list(ys)
    prev = None
    for i in order:
        if prev is not None and out[i] < out[prev] + gap:
            out[i] = out[prev] + gap
        prev = i
    shift = (sum(out) - sum(ys)) / max(len(ys), 1) * 0.5  # re-centre a little so labels stay near their points
    return [v - shift for v in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-sheet", required=True)
    ap.add_argument("--qc-table", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sheet = read_sheet(args.sample_sheet)
    S = []
    for r in read_tsv(args.qc_table):
        well, chip, pos = sheet.get(r["IID"], ("NA", "NA", "NA"))
        donor, group = ("", "")
        if r["is_donor_of_ipsc_line"]:
            donor, group = r["is_donor_of_ipsc_line"].rstrip(")").split("(")
        S.append(dict(iid=r["IID"], chip=chip, well=well, pos=pos, call=float(r["call_rate"]), het=float(r["het_F"]),
                      z=float(r["het_z"]), excl=r["excluded"] == "yes", reasons=r["reasons"], donor=donor, group=group))
    missing = [s["iid"] for s in S if s["chip"] == "NA"]
    chips = sorted({s["chip"] for s in S if s["chip"] != "NA"}, key=lambda c: min(int(s["iid"]) for s in S if s["chip"] == c))
    rep = [f"samples joined to the sheet: {len(S) - len(missing)} of {len(S)} | not found in the sheet: {missing}",
           f"chips: {len(chips)} (ordered by lowest Sample_ID)"]
    # het limits as used by step 09 (recovered from het_z)
    sd, mu = np.polyfit([s["z"] for s in S], [s["het"] for s in S], 1)
    lo, hi = mu - HET_SD * sd, mu + HET_SD * sd

    with open(os.path.join(args.outdir, "chip_sample_map.tsv"), "w") as out:
        out.write("IID\tchip\tsample_well\tchip_position\tcall_rate\thet_F\texcluded\treasons\tdonor_of_ipsc_line\n")
        for s in sorted(S, key=lambda s: (chips.index(s["chip"]) if s["chip"] in chips else 99, int(s["iid"]))):
            out.write(f"{s['iid']}\t{s['chip']}\t{s['well']}\t{s['pos']}\t{s['call']:.4f}\t{s['het']:.4f}\t"
                      f"{'yes' if s['excl'] else 'no'}\t{s['reasons']}\t{s['donor']}\n")
    head = ["chip", "n_samples", "sample_id_range", "n_excluded", "excluded_IIDs(reason)", "call_rate_mean", "call_rate_median",
            "call_rate_min", "n_call_rate_below_0.98", "hetF_mean", "hetF_SD", "hetF_min", "hetF_max", "n_het_outliers",
            "n_ipsc_donors"]
    summ = []
    for c in chips:
        g = [s for s in S if s["chip"] == c]
        ids = sorted(int(s["iid"]) for s in g)
        cr, hf = [s["call"] for s in g], [s["het"] for s in g]
        ex = [s for s in g if s["excl"]]
        summ.append([c, len(g), f"{ids[0]}-{ids[-1]}", len(ex), "; ".join(f"{s['iid']}({short_reasons(s['reasons'])})" for s in ex) or "-",
                     f"{statistics.mean(cr):.4f}", f"{statistics.median(cr):.4f}", f"{min(cr):.4f}",
                     sum(1 for v in cr if v < MIN_CALL_RATE), f"{statistics.mean(hf):.4f}",
                     f"{statistics.stdev(hf):.4f}" if len(hf) > 1 else "NA", f"{min(hf):.4f}", f"{max(hf):.4f}",
                     sum(1 for s in g if abs(s["z"]) > HET_SD), sum(1 for s in g if s["donor"])])
    with open(os.path.join(args.outdir, "chip_summary.tsv"), "w") as out:
        out.write("\t".join(head) + "\n")
        for r in summ:
            out.write("\t".join(map(str, r)) + "\n")
    for r in summ:
        rep.append("  " + " | ".join(f"{h}={v}" for h, v in zip(head, r)))

    # ---- exploratory tests ---------------------------------------------------------------------
    tests = []
    groups_cr = [[s["call"] for s in S if s["chip"] == c] for c in chips]
    groups_hf = [[s["het"] for s in S if s["chip"] == c] for c in chips]
    h1, p1 = stats.kruskal(*groups_cr)
    h2, p2 = stats.kruskal(*groups_hf)
    tests.append(f"Kruskal-Wallis across {len(chips)} chips: call rate H = {h1:.2f}, p = {p1:.3g} | heterozygosity F H = {h2:.2f}, p = {p2:.3g}")
    top = max(chips, key=lambda c: sum(1 for s in S if s["chip"] == c and s["excl"]))
    a = [s for s in S if s["chip"] == top]
    b = [s for s in S if s["chip"] != top]
    ea, eb = sum(s["excl"] for s in a), sum(s["excl"] for s in b)
    odds, pf = stats.fisher_exact([[ea, len(a) - ea], [eb, len(b) - eb]])
    tests.append(f"Chip with most exclusions: {top}: {ea}/{len(a)} excluded vs {eb}/{len(b)} on the other chips | "
                 f"Fisher exact (two-sided) odds ratio {odds:.2f}, p = {pf:.3g}")
    for lab, key in (("call rate", "call"), ("heterozygosity F", "het")):
        u, pu = stats.mannwhitneyu([s[key] for s in a], [s[key] for s in b], alternative="two-sided")
        tests.append(f"Mann-Whitney {top} vs the rest, {lab}: median {statistics.median([s[key] for s in a]):.4f} vs "
                     f"{statistics.median([s[key] for s in b]):.4f}, p = {pu:.3g}")
    pos_ex = [f"{s['iid']}: well {s['well']}, chip position {s['pos']}" for s in a if s["excl"]]
    tests.append(f"Excluded samples on {top}: " + "; ".join(pos_ex))
    tests.append("CAVEAT: chips are contiguous Sample_ID blocks, so chip is confounded with plate / extraction batch; one chip has only "
                 "4 samples; tests are exploratory.")
    rep += tests
    with open(os.path.join(args.outdir, "chip_tests.txt"), "w") as out:
        out.write("\n".join(tests) + "\n")

    # ---- figure ---------------------------------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(15, 11))
    rng = np.random.default_rng(1)  # fixed seed: reproducible jitter
    for ax, key, ylab, (ymin, ymax) in ((axes[0], "call", "sample call rate", (0.958, 0.9935)),
                                        (axes[1], "het", "autosomal heterozygosity F", (-0.15, 0.10))):
        for xi, c in enumerate(chips):
            g = [s for s in S if s["chip"] == c]
            for s in g:
                s["x"] = xi + rng.uniform(-0.22, 0.22)
            vals = [s[key] for s in g]
            ax.hlines(statistics.median(vals), xi - 0.32, xi + 0.32, color="#444444", lw=1.5, zorder=1)
            ax.text(xi, ymax - (ymax - ymin) * 0.01, f"{sum(s['excl'] for s in g)}/{len(g)} excluded", ha="center", va="top", fontsize=8.5)
            lab = [s for s in g if s["excl"] or s["donor"]]
            gap = (ymax - ymin) * 0.034
            ys = spread([s[key] for s in lab], gap)
            for s, ly in zip(lab, ys):
                col = (COLOR[s["group"]] if s["donor"] else "black")
                txt = (s["donor"] + f" ({s['iid']})" if s["donor"] else f"IID {s['iid']}") + (f" *{short_reasons(s['reasons'])}" if s["excl"] else "")
                ax.annotate(txt, xy=(s["x"], s[key]), xytext=(xi + 0.36, ly), textcoords="data", fontsize=6.5, color=col, va="center",
                            arrowprops=dict(arrowstyle="-", color=col, lw=0.4, alpha=0.6, shrinkA=0, shrinkB=2))
            for s in sorted(g, key=lambda s: (s["donor"] != "", s["excl"])):
                if s["donor"]:
                    c_ = COLOR[s["group"]]
                    ax.scatter(s["x"], s[key], marker="D" if s["excl"] else "o", s=55, facecolors="white" if s["excl"] else c_,
                               edgecolors=c_, linewidths=1.5, zorder=5)
                elif s["excl"]:
                    ax.scatter(s["x"], s[key], marker="X", s=55, color="black", zorder=4)
                else:
                    ax.scatter(s["x"], s[key], s=10, color=GREY, zorder=2)
        ax.set_xlim(-0.6, len(chips) - 0.1)
        ax.set_ylim(ymin, ymax)
        ax.set_xticks(range(len(chips)))
        ax.set_xticklabels([f"{c}\nSample_ID {r[2]}, n = {r[1]}" for c, r in zip(chips, summ)], fontsize=8)
        ax.set_ylabel(ylab)
        ax.grid(axis="y", color="#EEEEEE", lw=0.6)
    axes[0].axhline(MIN_CALL_RATE, color="k", ls="--", lw=1)
    axes[0].text(-0.58, MIN_CALL_RATE, f" threshold {MIN_CALL_RATE}", fontsize=8, va="bottom")
    for v, t in ((lo, f"mean - {HET_SD:g} SD = {lo:.3f}"), (hi, f"mean + {HET_SD:g} SD = {hi:.3f}")):
        axes[1].axhline(v, color="k", ls="--", lw=1)
        axes[1].text(-0.58, v, " " + t, fontsize=8, va="bottom" if v > 0 else "top")
    axes[0].set_title("Sample call rate by chip (SentrixBarcode_A); grey line = chip median", fontsize=10, loc="left")
    axes[1].set_title("Autosomal heterozygosity F by chip", fontsize=10, loc="left")
    handles = [plt.Line2D([], [], marker="o", ls="none", color=COLOR["PCOS"], label="PCOS donor of an iPSC line (passes QC)"),
               plt.Line2D([], [], marker="o", ls="none", color=COLOR["CTRL"], label="CTRL donor of an iPSC line (passes QC)"),
               plt.Line2D([], [], marker="D", ls="none", mfc="white", mew=1.5, color="#555555", label="donor failing sample QC (hollow)"),
               plt.Line2D([], [], marker="X", ls="none", color="black", label="excluded sample, not an iPSC donor"),
               plt.Line2D([], [], marker="o", ls="none", color=GREY, ms=3, label="other samples (pass QC)")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.outdir, f"qc_by_chip.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    with open(os.path.join(args.outdir, "qc_by_chip_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
