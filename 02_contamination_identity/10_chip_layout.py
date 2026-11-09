#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/10_chip_layout.py   (original step script: 29_chip_layout.py)
# Produces: Tests and layout figure of the chip with the contamination signal.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 29 - Chip 206149980019: is the whole chip worse, or only the 5 failed samples?  (diagnostic, nothing is modified)

(1) Mann-Whitney U test (two-sided) of call rate and heterozygosity F: the 18 samples of chip 206149980019 that PASS the sample QC
    (step 09) versus all samples that pass the QC on the other chips (n = 111). Also run for alpha_hat (contamination estimate of
    step 24), the share of intermediate-BAF homozygous calls and the LRR SD as additional (non-requested) descriptors.
    Effect size = rank-biserial correlation r = 2U/(n1*n2) - 1 (U of the chip-019 group; r < 0: chip-019 values are lower).
    Medians and IQR are reported per group. A per-chip table of the passing samples is written for context.
(2) Chip layout figure (SentrixPosition R01..R12 = 12 rows, C01..C02 = 2 columns; R01C01 top-left) for chip 206149980019 and for a
    "normal" chip (206482920107: 24 sample positions, 22 with alpha_hat < 1 %): left panel coloured by call rate, right panel by
    alpha_hat. Each cell shows the eSF ID (or "NA" for an IDAT position that has no row in the sample sheet); QC-failed samples get a
    thick red frame, donors of an iPSC line get a star and the line name. Call rate = PLINK call rate of step 09 (all SNPs); for
    the unassigned positions (not in the QC table) the autosomal GTC call rate of step 24 is used (the two agree to < 0.001 in the
    assigned samples).

Input : validation/results/contamination/sample_intensity_metrics.tsv (step 24; includes the QC-table columns)
Output: validation/results/contamination/chip_layout/{mannwhitney_chip019_vs_other_chips.tsv, chip_summary_passing.tsv,
        figures/C1_chip_layout_019_vs_normal.png/.pdf}
"""
import argparse, os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ap = argparse.ArgumentParser()
ap.add_argument("--metrics", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--chip", default="206149980019")
ap.add_argument("--normal-chip", default="206482920107")
args = ap.parse_args()
os.makedirs(os.path.join(args.outdir, "figures"), exist_ok=True)

d = pd.read_csv(args.metrics, sep="\t", dtype={"IID": str, "chip": str})
d["row"] = d.chip_position.str.extract(r"_R(\d+)C")[0].astype(int)
d["col"] = d.chip_position.str.extract(r"C(\d+)$")[0].astype(int)
d["assigned"] = d.IID != "NA_unassigned"
d["cr_plot"] = d.call_rate.fillna(d.call_rate_auto)

# ---------------------------------------------------------------- (1) Mann-Whitney
passing = d[d.excluded == "no"]
a = passing[passing.chip == args.chip]; b = passing[passing.chip != args.chip]
print(f"chip {args.chip}: {len(a)} passing samples | other chips: {len(b)} passing samples")
rows = []
for var, label in [("call_rate", "call rate (PLINK, all SNPs)"), ("het_F", "heterozygosity F"),
                   ("alpha_hat", "alpha_hat (contamination estimate)"),
                   ("frac_hom_intermediate_baf", "share of homozygous calls with intermediate BAF"), ("lrr_sd", "LRR SD")]:
    x, y = a[var].dropna(), b[var].dropna()
    U, p = mannwhitneyu(x, y, alternative="two-sided")
    rows.append(dict(variable=label, n_chip019=len(x), n_other=len(y),
                     median_chip019=x.median(), IQR_chip019=f"{x.quantile(.25):.4f}-{x.quantile(.75):.4f}",
                     median_other=y.median(), IQR_other=f"{y.quantile(.25):.4f}-{y.quantile(.75):.4f}",
                     U=U, p_value=p, rank_biserial_r=2 * U / (len(x) * len(y)) - 1))
mw = pd.DataFrame(rows)
mw.to_csv(os.path.join(args.outdir, "mannwhitney_chip019_vs_other_chips.tsv"), sep="\t", index=False, float_format="%.5g")
print(mw.to_string(index=False), flush=True)
cs = passing.groupby("chip").agg(n_pass=("IID", "size"), median_call_rate=("call_rate", "median"), median_het_F=("het_F", "median"),
                                 median_alpha_hat=("alpha_hat", "median"), max_alpha_hat=("alpha_hat", "max"))
cs.to_csv(os.path.join(args.outdir, "chip_summary_passing.tsv"), sep="\t", float_format="%.5g")
print(cs.round(4).to_string(), flush=True)

# ---------------------------------------------------------------- (2) layout figure
def layout(ax, chip, value, vmin, vmax, cmap, title, fmt):
    c = d[d.chip == chip]
    ax.set_xlim(0, 2); ax.set_ylim(12, 0); ax.set_aspect(0.55); ax.axis("off")
    ax.set_title(title, fontsize=9, pad=18)
    cm = plt.get_cmap(cmap)
    for r in c.itertuples():
        v = getattr(r, value)
        colour = cm(np.clip((v - vmin) / (vmax - vmin), 0, 1)) if np.isfinite(v) else "#dddddd"
        x0, y0 = r.col - 1, r.row - 1
        ax.add_patch(Rectangle((x0 + .03, y0 + .03), .94, .94, fc=colour, ec="#888888", lw=.6))
        fail = r.excluded == "yes"
        if fail:
            ax.add_patch(Rectangle((x0 + .05, y0 + .05), .90, .90, fc="none", ec="#c0392b", lw=3))
        if not r.assigned:
            ax.add_patch(Rectangle((x0 + .03, y0 + .03), .94, .94, fc="none", ec="#333333", lw=.8, hatch="///"))
        txt = ("NA" if not r.assigned else r.IID) + ("  ★" if isinstance(r.is_donor_of_ipsc_line, str) else "")
        lum = 0.299 * colour[0] + 0.587 * colour[1] + 0.114 * colour[2] if np.isfinite(v) else 1   # white text on dark cells
        tc = "white" if lum < 0.45 else "black"
        ax.text(x0 + .5, y0 + .38, txt, ha="center", va="center", fontsize=8.5, fontweight="bold", color=tc)
        sub = fmt(v) + (" | " + r.is_donor_of_ipsc_line.split("(")[0] if isinstance(r.is_donor_of_ipsc_line, str) else "")
        ax.text(x0 + .5, y0 + .72, sub, ha="center", va="center", fontsize=6.5, color=tc)
    for rr in range(12): ax.text(-.04, rr + .5, f"R{rr + 1:02d}", ha="right", va="center", fontsize=6)
    ax.text(.5, -.12, "C01", ha="center", fontsize=7); ax.text(1.5, -.12, "C02", ha="center", fontsize=7)


fig, axes = plt.subplots(1, 4, figsize=(17, 9.5))
chips = [args.chip, args.normal_chip]
for k, chip in enumerate(chips):
    nm = "chip" if chip == args.chip else "normal chip"
    layout(axes[2 * k], chip, "cr_plot", 0.95, 1.0, "viridis", f"{nm} {chip}\ncall rate", lambda v: f"{v:.4f}")
    layout(axes[2 * k + 1], chip, "alpha_hat", 0, 0.12, "Reds", f"{nm} {chip}\ncontamination estimate alpha_hat", lambda v: f"{v * 100:.1f} %")
for ax, (vmin, vmax, cmap, lab) in zip([axes[0], axes[1]], [(0.95, 1.0, "viridis", "call rate (<=0.95 saturates)"), (0, 12, "Reds", "alpha_hat (%)")]):
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax)); sm.set_array([])
    fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=.03, pad=.02, label=lab)
fig.suptitle("Chip layout (SentrixPosition, 12 rows x 2 columns). Red frame = failed sample QC | ★ = donor of an iPSC line (line name shown) | "
             "hatched, 'NA' = IDAT position with no sample in the sheet", fontsize=9)
fig.tight_layout(rect=(0, 0, 1, .96))
fig.savefig(os.path.join(args.outdir, "figures", "C1_chip_layout_019_vs_normal.png"), dpi=170)
fig.savefig(os.path.join(args.outdir, "figures", "C1_chip_layout_019_vs_normal.pdf")); plt.close(fig)
print("figure written", flush=True)
