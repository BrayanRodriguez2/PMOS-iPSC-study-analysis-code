#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/04_baf_contamination.py   (original step script: 24_baf_contamination.py)
# Produces: BAF/LRR summaries and the DNA-contamination estimate per sample.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 24 - Part A of the contamination diagnostic: is the excess heterozygosity caused by DNA contamination?

Reads the per-sample GTC files written by step 23 (genotype calls, B allele frequency = BAF, Log R Ratio = LRR),
computes per-sample intensity metrics on the autosomes and a contamination estimate, and draws the BAF distribution
of the samples of interest next to clean control samples.

THIS IS A DIAGNOSTIC ONLY: it reads GTC / sample-sheet / QC tables and writes new files under
validation/results/contamination/. No QC dataset is modified.

Libraries: numpy, pandas, matplotlib (conda env "idat"). The GTC files are read with the small parser defined below
(Illumina's BeadArrayFiles is not on PyPI and was not downloaded from GitHub); the locus table (chromosome per locus,
in manifest order) comes from the BPM through `bcftools +idat2gtc --snp-map` (step 23b).

Method (all on autosomal loci; chromosome from the BPM manifest):
  * q  = frequency of the B allele among the 140 donors (GenCall calls: AA=1, AB=2, BB=3); loci with q in [0.05, 0.95]
         and call rate >= 0.9 are "informative".
  * Contamination estimate ("BAF regression"): if a fraction alpha of the DNA comes from another person drawn from
    the same population, then at loci where the host is AA the expected BAF is  ~ alpha * q  (the contaminant brings
    B alleles with frequency q) and at BB loci  1-BAF ~ alpha * (1-q). We bin the AA-called loci by q (bins of 0.05),
    take the mean BAF per bin and fit a weighted straight line: the SLOPE estimates alpha (the intercept absorbs the
    cluster noise that every clean sample also has). The same is done for BB calls with (1-q, 1-BAF).
    alpha_hat = mean of the two slopes. REVISION: the first run used all q in 0.05-0.95 (column alpha_v1_allq, kept):
    it gave 3-8 % in samples whose other metrics are clean (e.g. IID 20: 7 % with 0.3 % intermediate BAF) because the
    extreme-q bins hold few loci and mis-called hets dominate the fit. The final alpha_hat uses only calls homozygous for
    the MAJOR allele (x = minor allele frequency in 0.05-0.5, 9 bins) and drops off-band loci (BAF > 0.3 / < 0.7).
    It is an approximation (BAF is a non-linear function of the B fraction and the
    contaminant is assumed to come from the same population), so it is reported next to the same statistic of the
    clean samples, which is the real reference. Rule of thumb in the literature: >= ~3-5 % is clearly contaminated.
  * Intermediate-BAF fraction among homozygous calls: share of AA/BB calls with BAF in [0.10, 0.35] U [0.65, 0.90]
    (contamination moves hom calls towards the middle; clean samples have almost none).
  * Het BAF spread: SD of BAF among AB calls (a mixture of two genotypes widens/splits the central band).
  * LRR SD (autosomes): general intensity noise; bad DNA gives high SD without BAF bands, contamination gives both.
  * call rate and het rate from the same calls (sanity check against the PLINK-based QC table of step 09).
Samples: all 144 IDAT positions; the sample sheet (AKurg_06.10.22_b38.csv) maps chip+position -> eSF ID. IDAT positions
that have no row in the sheet are labelled "NA_unassigned".
"""
import argparse, glob, os, struct, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


class GenotypeCalls:
    """Minimal reader of an Illumina GTC file (format v3-v5): only the per-locus arrays this analysis needs.
    Layout (checked on the files written by step 23): b"gtc" + version byte, int32 number of TOC entries, then
    TOC entries (uint16 id, int32 offset). Per-locus arrays sit at `offset` as int32 count + values:
    1002 genotype calls (uint8: 0 = no call, 1 = AA, 2 = AB, 3 = BB), 1012 B allele frequency (float32),
    1013 Log R Ratio (float32). Id 1 holds the number of loci directly in the offset field.
    (Illumina's own reader, BeadArrayFiles, is not on PyPI; ids as in https://github.com/Illumina/BeadArrayFiles.)"""
    def __init__(self, path):
        with open(path, "rb") as fh:
            self.b = fh.read()
        if self.b[:3] != b"gtc" or self.b[3] not in (3, 4, 5):
            raise ValueError(f"{path}: not a GTC v3-v5 file")
        n_toc = struct.unpack_from("<i", self.b, 4)[0]
        self.toc = dict(struct.unpack_from("<Hi", self.b, 8 + 6 * k) for k in range(n_toc))

    def _array(self, i, dtype):
        off = self.toc[i]; n = struct.unpack_from("<i", self.b, off)[0]
        return np.frombuffer(self.b, dtype=dtype, count=n, offset=off + 4)

    def get_num_snps(self): return self.toc[1]
    def get_genotypes(self): return self._array(1002, "u1")
    def get_ballele_freqs(self): return self._array(1012, "<f4")
    def get_logr_ratios(self): return self._array(1013, "<f4")


ap = argparse.ArgumentParser()
ap.add_argument("--gtc-dir", required=True, help="PRS/intensity/gtc (one sub-folder per chip)")
ap.add_argument("--snp-map", required=True, help="PRS/intensity/bpm_snp_map.csv (step 23b: locus table in manifest order)")
ap.add_argument("--sample-sheet", required=True)
ap.add_argument("--qc-table", required=True, help="results/qc/samples/sample_qc_table.tsv (step 09)")
ap.add_argument("--outdir", required=True)
ap.add_argument("--focal", default=__import__("os").environ.get("FOCAL_SAMPLE_IDS", ""), help="IIDs to plot (comma separated)")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)
fig_dir = os.path.join(args.outdir, "figures"); os.makedirs(fig_dir, exist_ok=True)

# ---------------------------------------------------------------- manifest (locus table from the BPM) + sample sheet
snpmap = pd.read_csv(args.snp_map, dtype={"Chromosome": str})
chroms = snpmap["Chromosome"].to_numpy()
autosome = np.isin(chroms, [str(i) for i in range(1, 23)])
print(f"BPM loci: {len(chroms):,} | autosomal: {autosome.sum():,}", flush=True)

sheet = pd.read_csv(args.sample_sheet, skiprows=8, dtype=str).dropna(subset=["SentrixBarcode_A"])
sheet["key"] = sheet["SentrixBarcode_A"] + "_" + sheet["SentrixPosition_A"]
key2iid = dict(zip(sheet["key"], sheet["Sample_ID"]))
qc = pd.read_csv(args.qc_table, sep="\t", dtype={"IID": str}).set_index("IID")

# ---------------------------------------------------------------- load all GTC files
files = sorted(glob.glob(os.path.join(args.gtc_dir, "*", "*.gtc")))
print(f"GTC files: {len(files)}", flush=True)
keys, GT, BAF, LRRSD, CALLRATE_ALL = [], [], [], [], []
focal_set = set(x for x in args.focal.split(",") if x)
focal_store = {}
for i, f in enumerate(files):
    key = os.path.basename(f)[:-4]
    g = GenotypeCalls(f)
    if g.get_num_snps() != len(chroms):
        sys.exit(f"{key}: {g.get_num_snps()} loci in GTC vs {len(chroms)} in BPM")
    gt = g.get_genotypes()
    baf = g.get_ballele_freqs().astype(np.float32)
    lrr = g.get_logr_ratios().astype(np.float32)
    ok = autosome & np.isfinite(lrr)
    keys.append(key); GT.append(gt[autosome].astype(np.uint8)); BAF.append(baf[autosome])
    LRRSD.append(float(np.nanstd(lrr[ok])))
    iid = key2iid.get(key)
    if iid in focal_set:
        focal_store[iid] = (baf[autosome], lrr[autosome], gt[autosome])
    if (i + 1) % 20 == 0: print(f"  read {i + 1}/{len(files)}", flush=True)
GT = np.vstack(GT); BAF = np.vstack(BAF); n = len(keys)

# ---------------------------------------------------------------- allele frequency of the B allele (assigned donors only)
assigned = np.array([k in key2iid for k in keys])
called = (GT[assigned] > 0)
nB = ((GT[assigned] == 2) * 1 + (GT[assigned] == 3) * 2).sum(0)
nc = called.sum(0)
q = np.where(nc > 0, nB / (2 * np.maximum(nc, 1)), np.nan)
cr_snp = nc / assigned.sum()
informative = (q >= 0.05) & (q <= 0.95) & (cr_snp >= 0.9)
print(f"informative autosomal loci: {informative.sum():,}", flush=True)


def slope(x, y, bins=np.linspace(0, 1, 21)):
    """weighted LS slope of the per-bin mean of y on the bin mean of x (see module docstring)"""
    idx = np.digitize(x, bins) - 1
    mx, my, w = [], [], []
    for b in range(len(bins) - 1):
        m = idx == b
        if m.sum() >= 50:
            mx.append(x[m].mean()); my.append(y[m].mean()); w.append(m.sum())
    if len(mx) < 4: return np.nan, np.nan
    s, ic = np.polyfit(mx, my, 1, w=np.sqrt(w))
    return float(s), float(ic)


rows = []
for j, key in enumerate(keys):
    gt, baf = GT[j], BAF[j]
    inf = informative & np.isfinite(baf)
    aa = inf & (gt == 1); bb = inf & (gt == 3); ab = (gt == 2) & np.isfinite(baf)
    # v1 (first version, kept for transparency): all informative q, 20 bins on [0, 1] -> noisy and biased upwards in
    # clean samples because the extreme-q bins hold few loci and a few mis-called hets dominate the weighted fit
    s_aa1, _ = slope(q[aa], baf[aa]); s_bb1, _ = slope(1 - q[bb], 1 - baf[bb])
    # final: only calls homozygous for the MAJOR allele (AA with q<=0.5, BB with q>=0.5; x = minor allele frequency),
    # x in [0.05, 0.5], loci whose "off-band" BAF exceeds 0.3 dropped (mis-called hets; real contamination <= ~0.15),
    # 9 bins of 0.05; separately for AA and BB calls
    aa2 = aa & (q <= 0.5) & (baf <= 0.3); bb2 = bb & (q >= 0.5) & (baf >= 0.7)
    mb = np.linspace(0.05, 0.5, 10)
    s_aa, i_aa = slope(q[aa2], baf[aa2], mb)
    s_bb, i_bb = slope(1 - q[bb2], 1 - baf[bb2], mb)
    hom = aa | bb
    hb = np.where(aa, baf, 1 - baf)[hom]
    mid = np.mean((hb >= 0.10) & (hb <= 0.35))
    nc_j = (gt > 0).sum()
    rows.append(dict(
        chip_position=key, IID=key2iid.get(key, "NA_unassigned"),
        n_called=int(nc_j), call_rate_auto=nc_j / len(gt),
        het_rate=float((gt == 2).sum() / max(nc_j, 1)),
        lrr_sd=LRRSD[j],
        alpha_AA=s_aa, alpha_BB=s_bb, alpha_hat=np.nanmean([s_aa, s_bb]),
        alpha_v1_allq=np.nanmean([s_aa1, s_bb1]),
        intercept_AA=i_aa, intercept_BB=i_bb,
        mean_baf_AA=float(baf[aa].mean()), mean_1mbaf_BB=float((1 - baf[bb]).mean()),
        frac_hom_intermediate_baf=float(mid), het_baf_sd=float(np.std(baf[ab])),
    ))
df = pd.DataFrame(rows)
df["chip"] = df["chip_position"].str.slice(0, 12)
# robust z-score of alpha_hat among the assigned samples (median/MAD)
a = df.loc[df.IID != "NA_unassigned", "alpha_hat"]
med, mad = a.median(), 1.4826 * (a - a.median()).abs().median()
df["alpha_robust_z"] = (df["alpha_hat"] - med) / mad
# join with the PLINK-based QC table of step 09 (call rate over all SNPs, het F, exclusion reasons, iPSC line)
df = df.merge(qc[["call_rate", "het_F", "het_z", "excluded", "reasons", "is_donor_of_ipsc_line"]],
              left_on="IID", right_index=True, how="left")
df = df.sort_values(["chip", "chip_position"])
cols = ["chip_position", "chip", "IID", "is_donor_of_ipsc_line", "excluded", "reasons", "call_rate", "het_F",
        "n_called", "call_rate_auto", "het_rate", "lrr_sd", "alpha_AA", "alpha_BB", "alpha_hat", "alpha_v1_allq", "alpha_robust_z",
        "intercept_AA", "intercept_BB", "mean_baf_AA", "mean_1mbaf_BB", "frac_hom_intermediate_baf", "het_baf_sd"]
out_tsv = os.path.join(args.outdir, "sample_intensity_metrics.tsv")
df[cols].to_csv(out_tsv, sep="\t", index=False, float_format="%.5g")
print("wrote", out_tsv, flush=True)
print(df[cols].sort_values("alpha_hat", ascending=False).head(12).to_string(index=False), flush=True)

# ---------------------------------------------------------------- controls (clean samples) and figures
clean = df[(df.excluded == "no") & (df.het_F.abs() < 0.01) & (df.call_rate >= 0.99)
           & (df.frac_hom_intermediate_baf < 0.002)].copy()
clean["IID_n"] = clean["IID"].astype(int)
same_chip = clean[clean.chip == "206149980019"].sort_values("IID_n")
other = clean[clean.chip != "206149980019"].sort_values("IID_n")
controls = [same_chip.IID.iloc[0], other.IID.iloc[0]]
print("clean controls:", controls, flush=True)
pos = {r.IID: r.chip_position for r in df.itertuples()}
need = [i for i in controls if i not in focal_store]
for i in need:
    g = GenotypeCalls(glob.glob(os.path.join(args.gtc_dir, "*", pos[i] + ".gtc"))[0])
    focal_store[i] = (g.get_ballele_freqs().astype(np.float32)[autosome],
                      g.get_logr_ratios().astype(np.float32)[autosome], g.get_genotypes()[autosome])
order = [i for i in args.focal.split(",") if i] + controls
np.savez_compressed(os.path.join(args.outdir, "focal_baf_lrr.npz"),
                    **{f"baf_{i}": focal_store[i][0] for i in order}, **{f"gt_{i}": focal_store[i][2] for i in order})
meta = df.set_index("IID")

# Figure 1: BAF density (all autosomal loci) of focal samples vs clean controls
fig, axes = plt.subplots(len(order), 1, figsize=(8, 2.0 * len(order)), sharex=True)
for ax, i in zip(axes, order):
    b = focal_store[i][0]; b = b[np.isfinite(b)]
    ax.hist(b, bins=200, range=(0, 1), color="#c0392b" if i in focal_set else "#2c7fb8")
    ax.set_yscale("log"); ax.set_ylabel("loci")
    m = meta.loc[i]
    lab = "CONTROL (clean)" if i in controls else "focal"
    ax.set_title(f"IID {i} ({lab}) | chip {m.chip} | het F={m.het_F:.3f} | alpha_hat={m.alpha_hat * 100:.1f} % | "
                 f"intermediate hom BAF={m.frac_hom_intermediate_baf * 100:.2f} %", fontsize=8, loc="left")
axes[-1].set_xlabel("B allele frequency (BAF), all autosomal loci")
fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "A1_baf_distribution_focal_vs_clean.png"), dpi=160)
fig.savefig(os.path.join(fig_dir, "A1_baf_distribution_focal_vs_clean.pdf")); plt.close(fig)

# Figure 2: zoom on the homozygous bands (BAF of AA calls and 1-BAF of BB calls) - where contamination shows up
fig, axes = plt.subplots(len(order), 1, figsize=(8, 2.0 * len(order)), sharex=True)
for ax, i in zip(axes, order):
    b, _, g = focal_store[i]
    hb = np.concatenate([b[(g == 1) & informative], 1 - b[(g == 3) & informative]]); hb = hb[np.isfinite(hb)]
    ax.hist(hb, bins=100, range=(0, 0.5), color="#c0392b" if i in focal_set else "#2c7fb8")
    ax.set_yscale("log"); ax.set_ylabel("loci")
    ax.set_title(f"IID {i}: distance of homozygous calls from their pure band", fontsize=8, loc="left")
axes[-1].set_xlabel("BAF of AA calls  /  1-BAF of BB calls (0 = pure; shifted bands = contamination)")
fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "A2_hom_band_zoom.png"), dpi=160)
fig.savefig(os.path.join(fig_dir, "A2_hom_band_zoom.pdf")); plt.close(fig)

# Figure 3: estimate per sample (all 144 positions), coloured by chip; focal samples labelled
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
chips = sorted(df.chip.unique()); cmap = plt.get_cmap("tab10")
for k, c in enumerate(chips):
    d = df[df.chip == c]
    axes[0].scatter(d.alpha_AA * 100, d.alpha_BB * 100, s=18, color=cmap(k), label=c, alpha=.8)
    axes[1].scatter(d.het_F, d.alpha_hat * 100, s=18, color=cmap(k), alpha=.8)
    axes[2].scatter(d.lrr_sd, d.alpha_hat * 100, s=18, color=cmap(k), alpha=.8)
for ax, (xl) in zip(axes, ["alpha from AA calls (%)", "heterozygosity F (PLINK, step 09)", "LRR SD (autosomes)"]):
    ax.set_xlabel(xl); ax.set_ylabel("alpha_hat (%)" if ax is not axes[0] else "alpha from BB calls (%)")
    for i in order[:3]:
        r = meta.loc[i]
        x = {0: r.alpha_AA * 100, 1: r.het_F, 2: r.lrr_sd}[list(axes).index(ax)]
        y = r.alpha_BB * 100 if ax is axes[0] else r.alpha_hat * 100
        ax.annotate(i, (x, y), fontsize=9, color="#c0392b", xytext=(4, 4), textcoords="offset points")
axes[0].legend(fontsize=6, title="chip", title_fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "A3_contamination_estimate_all_samples.png"), dpi=160)
fig.savefig(os.path.join(fig_dir, "A3_contamination_estimate_all_samples.pdf")); plt.close(fig)
print("figures written to", fig_dir, flush=True)
