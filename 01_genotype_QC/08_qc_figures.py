#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/08_qc_figures.py   (original step script: 10b_qc_figures.py)
# Produces: Standard QC figures drawn from the QC tables.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 10b - standard genotype-QC figures, drawn ONLY from tables that already exist
(nothing is recomputed: no plink2, no filtering, no new statistics except the exact threshold
lines that are recovered from the tables themselves).
It covers the outputs of step 09 (sample QC) and step 10 (variant QC).

Figures (PNG + PDF, in --outdir):
  01_sample_call_rate_hist   histogram of the sample call rate + one row per labelled sample
  02_het_vs_call_rate        autosomal heterozygosity F against sample call rate
  03_chrX_F_by_sample        chrX F statistic (sex check): histogram + one row per labelled sample
  04_kinship_pairs           KING kinship of the related pairs that were saved
  05_snp_call_rate           per-variant call rate of the variants that passed the filter (see note)
  06_snp_maf_hist            per-variant MAF before QC (140 samples) and in the final dataset (129 samples)
Labelled in every sample-level figure: the samples EXCLUDED by step 09 (with the reason: call rate, het,
kinship, sex) and the 14 donors of the iPSC lines (donor_of_ipsc_line column of sample_qc_table.tsv).
Symbols (same code as the PRS figures):
  filled circle, colour = group (orange PCOS / blue CTRL) : donor of an iPSC line that PASSES QC
  hollow diamond, colour = group                          : donor of an iPSC line that FAILS sample QC (*)
  black cross                                             : excluded sample that is not an iPSC donor
  small grey dot                                          : any other sample (passes QC)

Thresholds drawn (the ones used in steps 09/10; see the READMEs):
  sample call rate >= 0.98 | heterozygosity F within mean +/- 3 SD | chrX F: female <= 0.2, male >= 0.8
  KING: 0.0884 second degree (reported only), 0.177 first degree (sample removed), 0.354 duplicate/MZ
  variant call rate >= 0.98 (= at most 2 missing calls among 129 samples) | MAF >= 0.01
The mean and SD of the heterozygosity F used by step 09 are recovered EXACTLY from the table by fitting
het_F = mean + SD * het_z (het_z is stored in the table), so the lines are the ones the filter used.

Known limits of the existing tables (also stated on the figures):
  * Only variants that PASSED the call-rate filter are in variants/b38_qc.vmiss. If variants/raw_129.vmiss
    exists (step 10c: per-variant missingness of the RAW array on the 129 QC-pass samples) figure 05 uses
    it and shows the removed tail too (30,849 variants, cross-checked with qc_counts.tsv; chrY is left out
    because in an all-female cohort its call rate is undefined). Otherwise it falls back to b38_qc.vmiss.
    With 129 samples the call rate is discrete, so it is also drawn per number of missing calls.
  * The KING table keeps only pairs with kinship >= 0.0884, so figure 04 shows those pairs, not the full
    distribution.
  * MAF "before QC" comes from the plink2 --freq table of the raw 140 samples (results/bim_check/freq.afreq).

Inputs : --qc-dir (results/qc), --afreq-raw, --afreq-final, --outdir.   Needs numpy and matplotlib.
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# ---- thresholds used by steps 09/10 ------------------------------------------------------------
MIN_CALL_RATE = 0.98
HET_SD = 3.0
XF_FEMALE_MAX, XF_MALE_MIN = 0.2, 0.8
KING_2ND, KING_1ST, KING_DUP = 0.0884, 0.177, 0.354
SNP_CALL_RATE, SNP_MAF = 0.98, 0.01

COLOR = {"PCOS": "#D55E00", "CTRL": "#0072B2"}  # colour-blind-safe pair, as in the PRS figures
GREY = "#9A9A9A"


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


def short_reasons(reasons):
    """'low_call_rate;heterozygosity_outlier_low' -> 'call rate + het'."""
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


def load_samples(qc_dir):
    rows = read_tsv(os.path.join(qc_dir, "samples", "sample_qc_table.tsv"))
    out = []
    for r in rows:
        donor, group = "", ""
        if r["is_donor_of_ipsc_line"]:
            donor, group = r["is_donor_of_ipsc_line"].rstrip(")").split("(")
        out.append(dict(iid=r["IID"], pedsex=r["pedsex"], xf=float(r["xF"]), call=float(r["call_rate"]),
                        het=float(r["het_F"]), z=float(r["het_z"]), donor=donor, group=group,
                        excluded=r["excluded"] == "yes", reasons=r["reasons"]))
    return out


def label_of(s):
    txt = f"{s['donor']} (IID {s['iid']})" if s["donor"] else f"IID {s['iid']}"
    if s["excluded"]:
        txt += (" *" if s["donor"] else " excl.") + f" ({short_reasons(s['reasons'])})"
    if s["pedsex"] == "0":
        txt += " [.fam sex=0]"
    return txt


def style(s):
    """(marker, facecolor, edgecolor, size, zorder, text colour) for a sample."""
    if s["donor"]:
        c = COLOR[s["group"]]
        if s["excluded"]:
            return "D", "white", c, 70, 5, c
        return "o", c, c, 60, 4, c
    if s["excluded"]:
        return "X", "black", "black", 60, 4, "black"
    return "o", GREY, GREY, 10, 2, GREY


def draw_points(ax, samples, xkey, ykey):
    for s in sorted(samples, key=lambda s: style(s)[4]):
        m, fc, ec, sz, zo, _ = style(s)
        ax.scatter(s[xkey], s[ykey], marker=m, s=sz, facecolors=fc, edgecolors=ec, linewidths=1.5 if m == "D" else 0.6, zorder=zo)


def labelled(samples):
    return [s for s in samples if s["donor"] or s["excluded"]]


def fan_labels(ax, samples, xkey, ykey, fontsize=7):
    """Labels in a column right of the axes (ordered by y), joined to their points by thin lines:
    avoids the overlaps that appear when many donors sit in the same cluster."""
    pts = sorted(labelled(samples), key=lambda s: s[ykey])
    n = len(pts)
    ys = np.linspace(0.02, 0.98, n) if n > 1 else [0.5]
    for s, yy in zip(pts, ys):
        col = style(s)[5]
        ax.annotate(label_of(s), xy=(s[xkey], s[ykey]), xycoords="data", xytext=(1.03, yy), textcoords="axes fraction",
                    fontsize=fontsize, color=col, va="center", ha="left", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.5, alpha=0.6, shrinkA=0, shrinkB=3))


def legend_handles():
    return [
        plt.Line2D([], [], marker="o", ls="none", color=COLOR["PCOS"], label="PCOS donor of an iPSC line (passes QC)"),
        plt.Line2D([], [], marker="o", ls="none", color=COLOR["CTRL"], label="CTRL donor of an iPSC line (passes QC)"),
        plt.Line2D([], [], marker="D", ls="none", mfc="white", mew=1.5, color="#555555",
                   label="donor of an iPSC line that FAILS sample QC (hollow, *)"),
        plt.Line2D([], [], marker="X", ls="none", color="black", label="excluded sample that is not an iPSC donor"),
        plt.Line2D([], [], marker="o", ls="none", color=GREY, ms=3, label="other samples (pass QC)"),
    ]


def save(fig, outdir, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---- figures -------------------------------------------------------------------------------------
def hist_strip(samples, key, bins, xlabel, title, outdir, name, vlines, spans=(), xlim=None):
    """Histogram of all samples (top) + one row per labelled sample sorted by `key` (bottom), same x axis.
    vlines: [(x, text)] threshold lines; spans: [(x0, x1, colour)] shaded regions."""
    lab = sorted(labelled(samples), key=lambda s: s[key])
    fig = plt.figure(figsize=(9.5, 4 + 0.27 * len(lab)))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.5, 0.27 * len(lab) / 1.6 + 1.0], hspace=0.08)
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax)
    ax.hist(np.array([s[key] for s in samples]), bins=bins, color="#BBBBBB", edgecolor="white")
    for a in (ax, ax2):
        for x0, x1, c in spans:
            a.axvspan(x0, x1, color=c, zorder=0, alpha=0.8)
        for x, _ in vlines:
            a.axvline(x, color="k", ls="--", lw=1)
    top = ax.get_ylim()[1]
    for x, txt in vlines:
        ax.text(x, top * 0.95, " " + txt, fontsize=8, va="top")
    if xlim:
        ax.set_xlim(*xlim)
    ax.set_ylabel("samples")
    ax.set_title(title, fontsize=10, loc="left")
    plt.setp(ax.get_xticklabels(), visible=False)
    for row, s in enumerate(lab):
        m, fc, ec, sz, zo, _ = style(s)
        ax2.scatter(s[key], row, marker=m, s=max(sz, 45), facecolors=fc, edgecolors=ec, linewidths=1.5 if m == "D" else 0.6, zorder=3)
    ax2.set_yticks(range(len(lab)))
    ax2.set_yticklabels([label_of(s) for s in lab], fontsize=7.5)
    for tick, s in zip(ax2.get_yticklabels(), lab):
        tick.set_color(style(s)[5])
    ax2.set_ylim(-0.7, len(lab) - 0.3)
    ax2.grid(axis="y", color="#EEEEEE", lw=0.6)
    ax2.set_xlabel(xlabel)
    fig.legend(handles=legend_handles(), loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=7.5, frameon=False)
    save(fig, outdir, name)


def fig_call_rate(samples, outdir):
    calls = np.array([s["call"] for s in samples])
    n_below = int((calls < MIN_CALL_RATE).sum())
    hist_strip(samples, "call", np.arange(0.960, calls.max() + 0.0015, 0.001), "sample call rate",
               f"Sample call rate (autosomal variants with call rate >= 95 %), n = {len(samples)}",
               outdir, "01_sample_call_rate_hist", [(MIN_CALL_RATE, f"threshold {MIN_CALL_RATE}\n {n_below} samples below")])


def het_fit(samples):
    """Recover mean and SD used by step 09 from het_F = mean + SD * het_z."""
    z = np.array([s["z"] for s in samples])
    f = np.array([s["het"] for s in samples])
    sd, mu = np.polyfit(z, f, 1)
    return mu, sd


def fig_het(samples, outdir, rep):
    mu, sd = het_fit(samples)
    rep.append(f"heterozygosity F: mean {mu:.5f}, SD {sd:.5f} (recovered from het_z); limits {mu - HET_SD * sd:.4f} .. {mu + HET_SD * sd:.4f}")
    fig, ax = plt.subplots(figsize=(11, 6.5))
    draw_points(ax, samples, "call", "het")
    ax.axvline(MIN_CALL_RATE, color="k", ls="--", lw=1)
    ax.axhline(mu, color="#777777", ls=":", lw=0.8)
    for k in (-HET_SD, HET_SD):
        ax.axhline(mu + k * sd, color="k", ls="--", lw=1)
    ax.text(ax.get_xlim()[0], mu + HET_SD * sd, f" mean + {HET_SD:g} SD = {mu + HET_SD * sd:.3f}", fontsize=8, va="bottom")
    ax.text(ax.get_xlim()[0], mu - HET_SD * sd, f" mean - {HET_SD:g} SD = {mu - HET_SD * sd:.3f}", fontsize=8, va="top")
    ax.text(MIN_CALL_RATE, ax.get_ylim()[1], f" call rate {MIN_CALL_RATE}", fontsize=8, va="top")
    ax.set_xlabel("sample call rate")
    ax.set_ylabel("autosomal heterozygosity F (LD-pruned, MAF >= 5 %)")
    ax.set_title("Heterozygosity vs call rate: a sample is excluded if it is left of the vertical line or outside the horizontal lines",
                 fontsize=9.5, loc="left")
    fan_labels(ax, samples, "call", "het")
    fig.subplots_adjust(right=0.70)
    fig.legend(handles=legend_handles(), loc="upper center", bbox_to_anchor=(0.45, 0.02), ncol=2, fontsize=7.5, frameon=False)
    save(fig, outdir, "02_het_vs_call_rate")


def fig_xf(samples, outdir):
    n_f = sum(1 for s in samples if s["xf"] <= XF_FEMALE_MAX)
    hist_strip(samples, "xf", np.arange(-0.2, 1.0001, 0.02), "chrX inbreeding coefficient F",
               f"Sex check (chrX F): {n_f} of {len(samples)} samples are called female, including the 2 with sex=0 in the .fam",
               outdir, "03_chrX_F_by_sample",
               [(XF_FEMALE_MAX, f"max F for female = {XF_FEMALE_MAX}"), (XF_MALE_MIN, f"min F for male = {XF_MALE_MIN}")],
               spans=[(XF_FEMALE_MAX, XF_MALE_MIN, "#FFF3C4"), (XF_MALE_MIN, 1.0, "#E1ECF7")], xlim=(-0.2, 1.0))


def fig_kinship(qc_dir, samples, outdir):
    rows = read_tsv(os.path.join(qc_dir, "samples", "related_pairs.tsv"))
    info = {s["iid"]: s for s in samples}
    fig, ax = plt.subplots(figsize=(9, 2.6 + 0.5 * len(rows)))
    for thr, txt in ((KING_2ND, "second degree\n(reported only)"), (KING_1ST, "first degree\n(one sample removed)"),
                     (KING_DUP, "duplicate / MZ")):
        ax.axvline(thr, color="k", ls="--", lw=1)
        ax.text(thr, len(rows) - 0.35, f" {thr}\n {txt}", fontsize=7.5, va="bottom")
    names = []
    for i, r in enumerate(sorted(rows, key=lambda r: float(r["kinship"]))):
        k = float(r["kinship"])
        ax.scatter(k, i, s=70, color="black" if k >= KING_1ST else GREY, zorder=3)
        nm = []
        for iid in (r["IID1"], r["IID2"]):
            s = info[iid]
            nm.append(f"{iid}" + (f" ({s['donor']}{' *' if s['excluded'] else ''})" if s["donor"] else "")
                      + ("" if not s["excluded"] or s["donor"] else " excl."))
        names.append(f"{nm[0]}  -  {nm[1]}")
        ax.text(k + 0.006, i, f"{k:.3f}", va="center", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_ylim(-0.7, len(rows) + 0.5)
    ax.set_xlim(0, 0.42)
    ax.set_xlabel("KING-robust kinship")
    ax.set_title("Related pairs saved by step 09 (kinship >= 0.0884). Pairs below the table filter are not stored.\n"
                 "IID = sample ID; * = fails sample QC; excl. = excluded", fontsize=9, loc="left")
    save(fig, outdir, "04_kinship_pairs")


def fig_snp_call_rate_raw(raw_vmiss, qc_dir, outdir, rep):
    """Figure 05 with the REMOVED tail: per-variant call rate of the raw array on the 129 QC-pass samples
    (results/qc/variants/raw_129.vmiss, step 10c). chrY is left out: in an all-female cohort plink2 has no
    observations there (OBS_CT = 0, call rate undefined) and --geno does not act on it (they leave at the
    autosome filter instead). With chrY excluded, the number of variants above 2 % missing equals the
    variants removed by filter S2 of step 10 (qc_counts.tsv): this is checked and printed."""
    rows = [r for r in read_tsv(raw_vmiss) if r["CHROM"] != "Y"]
    miss = np.array([int(r["MISSING_CT"]) for r in rows])
    fmiss = np.array([float(r["F_MISS"]) for r in rows])
    obs = int(rows[0]["OBS_CT"])
    removed_tab = [int(r["removed_in_step"]) for r in read_tsv(os.path.join(qc_dir, "variants", "qc_counts.tsv"))
                   if r["step"].startswith("S2")][0]
    n_fail = int((fmiss > 1 - SNP_CALL_RATE + 1e-12).sum())
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5))
    ax = axes[0]
    top_lvl = 15
    levels = np.arange(0, top_lvl + 1)
    n = np.array([(miss == m).sum() for m in levels])
    over = int((miss > top_lvl).sum())
    max_ok = int(np.floor(obs * (1 - SNP_CALL_RATE)))
    colors = ["#BBBBBB" if m <= max_ok else "#C0504D" for m in levels] + ["#C0504D"]
    ax.bar(list(levels) + [top_lvl + 1.2], list(n) + [over], color=colors, edgecolor="white", width=0.8)
    ax.set_yscale("log")
    ax.set_ylim(top=n.max() * 6)
    for x, c in zip(list(levels) + [top_lvl + 1.2], list(n) + [over]):
        ax.text(x, c * 1.15, f"{c:,}", ha="center", va="bottom", fontsize=6.5, rotation=90)
    ax.set_xticks(list(levels) + [top_lvl + 1.2])
    ax.set_xticklabels([str(m) for m in levels] + [f">{top_lvl}"], fontsize=8)
    ax.axvline(max_ok + 0.5, color="k", ls="--", lw=1)
    ax.text(max_ok + 0.6, n.max() * 3, f" threshold {SNP_CALL_RATE}: <= {max_ok} missing calls\n of {obs} (grey = kept, red = removed)", fontsize=8, va="top")
    ax.set_xlabel(f"missing calls per variant among the {obs} samples that pass QC")
    ax.set_ylabel("variants (log scale)")
    ax.set_title("Number of missing calls per variant", fontsize=9.5, loc="left")
    ax = axes[1]
    ax.hist(1 - fmiss, bins=np.arange(0, 1.0001, 0.02), color="#BBBBBB", edgecolor="white")
    ax.set_yscale("log")
    ax.axvline(SNP_CALL_RATE, color="k", ls="--", lw=1)
    ax.text(SNP_CALL_RATE, ax.get_ylim()[1], f" threshold {SNP_CALL_RATE}", fontsize=8, va="top", ha="right")
    ax.set_xlabel("variant call rate")
    ax.set_ylabel("variants (log scale)")
    ax.set_title("Whole distribution (2 % bins)", fontsize=9.5, loc="left")
    fig.suptitle(f"Per-variant call rate of the raw array ({len(rows):,} variants, chrY excluded), 129 samples passing QC: "
                 f"{n_fail:,} variants fall below {SNP_CALL_RATE} (filter S2 of step 10 removed {removed_tab:,}: "
                 f"{'consistent' if n_fail == removed_tab else 'NOT consistent'})", fontsize=9.5)
    fig.tight_layout()
    rep.append(f"RAW SNP call rate (chrY excluded, n={len(rows):,}): below {SNP_CALL_RATE}: {n_fail:,} | S2 removed (qc_counts): "
               f"{removed_tab:,} | consistent: {n_fail == removed_tab} | variants per missing-call count 0..5: "
               + ", ".join(f"{m}->{int((miss == m).sum()):,}" for m in range(6)) + f" | > {top_lvl}: {over:,}")
    save(fig, outdir, "05_snp_call_rate")


def fig_snp_call_rate(qc_dir, outdir, rep):
    rows = read_tsv(os.path.join(qc_dir, "variants", "b38_qc.vmiss"))
    miss = np.array([int(r["MISSING_CT"]) for r in rows])
    obs = int(rows[0]["OBS_CT"])
    counts = {r["step"]: int(r["removed_in_step"]) for r in read_tsv(os.path.join(qc_dir, "variants", "qc_counts.tsv"))
              if r["removed_in_step"] != "NA"}
    removed = [v for k, v in counts.items() if k.startswith("S2")][0]
    levels = np.arange(0, miss.max() + 1)
    n = np.array([(miss == m).sum() for m in levels])
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.bar(levels, n, color="#BBBBBB", edgecolor="white", width=0.7)
    for m, c in zip(levels, n):
        ax.text(m, c, f"{c:,}\n({100 * c / len(miss):.1f} %)", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(levels)
    ax.set_xticklabels([f"{m} missing\n(call rate {100 * (1 - m / obs):.2f} %)" for m in levels])
    ax.set_ylim(0, n.max() * 1.25)
    ax.set_ylabel("variants")
    ax.set_xlabel(f"missing calls per variant among the {obs} samples that pass QC")
    ax.set_title(f"Per-variant call rate of the {len(miss):,} variants in b38_qc (autosomal, polymorphic, HWE and MAF filters applied)",
                 fontsize=9, loc="left")
    max_ok = int(np.floor(obs * (1 - SNP_CALL_RATE)))
    ax.text(0.98, 0.96, f"threshold {SNP_CALL_RATE}: at most {max_ok} missing calls of {obs}\n(3 missing = {100 * (1 - 3 / obs):.2f} % would be removed)",
            transform=ax.transAxes, fontsize=8.5, va="top", ha="right")
    ax.text(0.98, 0.78, f"NOTE: the {removed:,} variants removed by the call-rate filter (S2) are NOT\nin the saved table, so the removed tail cannot be shown.",
            transform=ax.transAxes, fontsize=8, va="top", ha="right", color="#8A1C1C")
    rep.append(f"SNP call rate: {len(miss):,} variants in b38_qc.vmiss; missing calls per variant: "
               + ", ".join(f"{m} -> {c:,}" for m, c in zip(levels, n)) + f" (of {obs} samples); S2 removed {removed:,}")
    save(fig, outdir, "05_snp_call_rate")


def read_maf(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        ia = header.index("ALT_FREQS")
        out = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            p = float(f[ia]) if f[ia] not in ("NA", "nan") else float("nan")
            out.append(min(p, 1 - p))
    return np.array(out)


def fig_snp_maf(afreq_raw, afreq_final, outdir, rep):
    raw, fin = read_maf(afreq_raw), read_maf(afreq_final)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for ax, maf, title in ((axes[0], raw, f"Before QC: raw array, 140 samples, {len(raw):,} variants"),
                           (axes[1], fin, f"Final dataset: 129 samples, {len(fin):,} variants")):
        ax.hist(maf[~np.isnan(maf)], bins=np.arange(0, 0.5001, 0.01), color="#BBBBBB", edgecolor="white")
        ax.set_yscale("log")
        ax.set_ylim(bottom=1e3)
        ax.axvline(SNP_MAF, color="k", ls="--", lw=1)
        ax.set_xlabel("minor allele frequency")
        ax.set_ylabel("variants (log scale)")
        ax.set_title(title, fontsize=9, loc="left")
        mono = int(np.sum(maf == 0))
        low = int(np.sum((maf > 0) & (maf < SNP_MAF)))
        ok = int(np.sum(maf >= SNP_MAF))
        ax.text(0.98, 0.97, f"monomorphic (MAF = 0): {mono:,}\n0 < MAF < {SNP_MAF}: {low:,}\nMAF >= {SNP_MAF}: {ok:,}",
                transform=ax.transAxes, fontsize=8, va="top", ha="right")
        ax.text(SNP_MAF, ax.get_ylim()[1], f" threshold {SNP_MAF}", fontsize=8, va="top")
        rep.append(f"MAF [{title}]: monomorphic {mono:,}; 0<MAF<{SNP_MAF} {low:,}; >= {SNP_MAF} {ok:,}")
    fig.suptitle("Minor allele frequency per variant (the first bin, MAF 0-0.01, contains the monomorphic variants)", fontsize=10)
    fig.tight_layout()
    save(fig, outdir, "06_snp_maf_hist")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--afreq-raw", required=True)
    ap.add_argument("--afreq-final", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    samples = load_samples(args.qc_dir)
    rep = [f"samples in table: {len(samples)} | excluded: {sum(s['excluded'] for s in samples)} | "
           f"iPSC donors: {sum(1 for s in samples if s['donor'])} (of which excluded: {sum(1 for s in samples if s['donor'] and s['excluded'])})"]
    fig_call_rate(samples, args.outdir)
    fig_het(samples, args.outdir, rep)
    fig_xf(samples, args.outdir)
    fig_kinship(args.qc_dir, samples, args.outdir)
    raw_vmiss = os.path.join(args.qc_dir, "variants", "raw_129.vmiss")
    if os.path.exists(raw_vmiss):
        fig_snp_call_rate_raw(raw_vmiss, args.qc_dir, args.outdir, rep)  # with the removed tail (step 10c)
    else:
        fig_snp_call_rate(args.qc_dir, args.outdir, rep)  # kept variants only
    fig_snp_maf(args.afreq_raw, args.afreq_final, args.outdir, rep)
    with open(os.path.join(args.outdir, "qc_figures_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
