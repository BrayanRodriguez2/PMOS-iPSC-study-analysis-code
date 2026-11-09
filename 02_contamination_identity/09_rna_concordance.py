#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/09_rna_concordance.py   (original step script: 28_rna_concordance.py)
# Produces: Line x donor genotype concordance (identity check).
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 28 - Part B of the contamination diagnostic: is each iPSC line really from its donor?  (RNA-seq genotypes vs array genotypes)

DIAGNOSTIC ONLY: reads the RNA genotype calls of step 27b, the donor array genotypes extracted from the TOPMed-prepared VCF
(step 20: 140 donors, no sample QC) and writes new tables/figures under validation/results/rna_identity/. No QC dataset is modified.

For each iPSC line L (14) and each donor d (140), over the array SNPs inside exons where L has >= MIN_DP reads and d has a call:
  concordance[L,d] = fraction of SNPs with the same genotype (0/0, 0/1, 1/1)        -> 14 x 140 matrix (headline, as requested)
  ibs0[L,d]        = fraction of SNPs where one is 0/0 and the other 1/1              (the cleanest discordance: no het ambiguity)
  concordance_common[L,d] = concordance restricted to common SNPs (donor-cohort MAF >= 0.05) where hom/hom matches are less trivial
The line's own donor should stand out from the 139 others (same person -> ~0.9+ concordance; unrelated -> clearly lower).

Stratified view for the donors of interest (split by the DONOR's array genotype, for the line vs its own donor):
  * donor-HOMOZYGOUS sites: concordance should be high whatever happens (a contaminant mostly turns hom calls into hets, not into
    the other hom), and "RNA het at donor-hom" is the reverse error;
  * donor-HETEROZYGOUS sites: the RNA must also be 0/1. If the array DNA is contaminated, part of the array hets are FALSE hets
    (true genotype hom + contaminant allele) that the RNA does not confirm -> lower het concordance with intact hom concordance.
  RNA-only effects also lower het concordance (allele-specific expression, low depth, mapping bias): therefore the other 12 lines
  (own donors) are the baseline, and the lines of interest are compared with their mean +- SD.
Also: array het rate of the donor vs RNA het rate of the line on the SAME sites.

Inputs : rna_samples.tsv (step 26), genotypes/<fastq_name>.calls.tsv.gz (27b), donor_gt.tsv + donor_samples.txt (made by the sbatch
         wrapper with bcftools from the array VCF), sample_qc_table.tsv (step 09; flags only)
Outputs: concordance_matrix.tsv, ibs0_matrix.tsv, nsites_matrix.tsv, concordance_common_matrix.tsv, best_donor_per_line.tsv,
         het_hom_split.tsv, report_B.txt, figures/B1_concordance_heatmap.*, B2_het_hom_split.*
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--rna-samples", required=True)
ap.add_argument("--calls-dir", required=True)
ap.add_argument("--donor-gt", required=True, help="CHROM POS REF ALT then one GT column per donor (no header)")
ap.add_argument("--donor-samples", required=True, help="donor IIDs in the column order of --donor-gt")
ap.add_argument("--qc-table", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--min-dp", type=int, default=10)
ap.add_argument("--focal-donors", default=__import__("os").environ.get("FOCAL_SAMPLE_IDS", ""))
args = ap.parse_args()
os.makedirs(os.path.join(args.outdir, "figures"), exist_ok=True)
fig_dir = os.path.join(args.outdir, "figures")
focal = [x for x in args.focal_donors.split(",") if x]
GT2D = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/0": 1, "1/1": 2, "1|1": 2}

# ------------------------------------------------------------------ donors
donors = [l.strip() for l in open(args.donor_samples) if l.strip()]
dg = pd.read_csv(args.donor_gt, sep="\t", header=None, dtype=str)
assert dg.shape[1] == 4 + len(donors), (dg.shape, len(donors))
site_key = dg[0] + ":" + dg[1]
assert site_key.is_unique
D = np.vstack([dg[4 + j].map(GT2D).fillna(-1).astype(np.int8).to_numpy() for j in range(len(donors))]).T  # sites x donors
site_idx = {k: i for i, k in enumerate(site_key)}
called = D >= 0
alt_af = np.where(called.sum(1) > 0, np.where(called, D, 0).sum(1) / (2 * np.maximum(called.sum(1), 1)), np.nan)
maf = np.minimum(alt_af, 1 - alt_af)
common = maf >= 0.05
print(f"donor genotype matrix: {D.shape[0]:,} exonic array SNPs x {len(donors)} donors; common (MAF>=0.05): {common.sum():,}", flush=True)

# ------------------------------------------------------------------ lines
rs = pd.read_csv(args.rna_samples, sep="\t", dtype=str)
qc = pd.read_csv(args.qc_table, sep="\t", dtype={"IID": str}).set_index("IID")
R, depth_ok, alt_frac = {}, {}, {}
for r in rs.itertuples():
    c = pd.read_csv(os.path.join(args.calls_dir, f"{r.fastq_name}.calls.tsv.gz"), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "alt", "gt", "dp", "ad"], dtype={"chrom": str, "gt": str, "ad": str})
    c = c[c.dp >= args.min_dp]
    key = (c.chrom + ":" + c.pos.astype(str)).map(site_idx)
    c = c[key.notna()]; key = key[key.notna()].astype(int)
    dos = c["gt"].map(GT2D)      # (c.gt would be the DataFrame method "greater than")
    ok = dos.notna().to_numpy()
    v = np.full(D.shape[0], -1, dtype=np.int8); v[key.to_numpy()[ok]] = dos.to_numpy()[ok].astype(np.int8)
    af = np.full(D.shape[0], np.nan)
    ad = c.ad.str.split(",", expand=True).reindex(columns=[0, 1]).apply(pd.to_numeric, errors="coerce")  # ref, alt read counts
    af[key.to_numpy()] = (ad[1] / (ad[0] + ad[1])).to_numpy()
    R[r.fastq_name] = v; alt_frac[r.fastq_name] = af
    print(f"{r.fastq_name} ({r.ipsc_line}, expected donor {r.esf_id}): {int((v >= 0).sum()):,} SNPs genotyped with >= {args.min_dp} reads", flush=True)

lines = list(rs.fastq_name)
conc = pd.DataFrame(np.nan, index=lines, columns=donors); ibs0 = conc.copy(); nsite = conc.copy(); conc_c = conc.copy()
for L in lines:
    v = R[L]
    for j, d in enumerate(donors):
        m = (v >= 0) & (D[:, j] >= 0)
        n = int(m.sum()); nsite.loc[L, d] = n
        if n:
            conc.loc[L, d] = np.mean(v[m] == D[m, j]); ibs0.loc[L, d] = np.mean(np.abs(v[m] - D[m, j]) == 2)
        mc = m & common
        if mc.sum(): conc_c.loc[L, d] = np.mean(v[mc] == D[mc, j])
label = {r.fastq_name: f"{r.ipsc_line} (donor {r.esf_id}, {r.group})" for r in rs.itertuples()}
for name, mat in [("concordance_matrix", conc), ("ibs0_matrix", ibs0), ("nsites_matrix", nsite), ("concordance_common_matrix", conc_c)]:
    out = mat.copy(); out.index = [label[i] for i in out.index]
    out.to_csv(os.path.join(args.outdir, name + ".tsv"), sep="\t", float_format="%.5f")

# ------------------------------------------------------------------ best donor per line
rows = []
for r in rs.itertuples():
    s = conc.loc[r.fastq_name].astype(float); si = ibs0.loc[r.fastq_name].astype(float)
    order = s.sort_values(ascending=False)
    exp = r.esf_id
    others = s.drop(exp)
    rows.append(dict(fastq_name=r.fastq_name, ipsc_line=r.ipsc_line, group=r.group, expected_donor=exp,
                     n_sites_expected=int(nsite.loc[r.fastq_name, exp]),
                     conc_expected=s[exp], ibs0_expected=si[exp],
                     best_donor=order.index[0], conc_best=order.iloc[0], second_best_donor=order.index[1], conc_second=order.iloc[1],
                     margin_best_vs_second=order.iloc[0] - order.iloc[1],
                     rank_of_expected_by_concordance=int(list(order.index).index(exp)) + 1,
                     best_donor_by_min_ibs0=si.idxmin(), ibs0_best=si.min(),
                     z_expected_vs_139_others=(s[exp] - others.mean()) / others.std(),
                     max_conc_among_139_others=others.max(),
                     identity="CORRECT (own donor is the best match)" if order.index[0] == exp else "MISMATCH (best match is another donor)"))
best = pd.DataFrame(rows)
best.to_csv(os.path.join(args.outdir, "best_donor_per_line.tsv"), sep="\t", index=False, float_format="%.5f")
print(best[["ipsc_line", "expected_donor", "conc_expected", "best_donor", "conc_best", "conc_second", "rank_of_expected_by_concordance",
            "n_sites_expected", "identity"]].to_string(index=False), flush=True)

# ------------------------------------------------------------------ stratified view vs own donor
srows = []
for r in rs.itertuples():
    j = donors.index(r.esf_id); v = R[r.fastq_name]; d = D[:, j]; m = (v >= 0) & (d >= 0)
    hom = m & (d != 1); het = m & (d == 1)
    srows.append(dict(
        fastq_name=r.fastq_name, ipsc_line=r.ipsc_line, donor=r.esf_id, n_sites=int(m.sum()),
        n_donor_hom=int(hom.sum()), conc_at_donor_hom=np.mean(v[hom] == d[hom]),
        rna_het_at_donor_hom=np.mean(v[hom] == 1),
        n_donor_het=int(het.sum()), conc_at_donor_het=np.mean(v[het] == 1),
        rna_hom_at_donor_het=np.mean(v[het] != 1),
        rna_homref_at_donor_het=np.mean(v[het] == 0), rna_homalt_at_donor_het=np.mean(v[het] == 2),
        mean_alt_frac_rna_at_donor_het=np.nanmean(alt_frac[r.fastq_name][het]),
        donor_het_rate_same_sites=np.mean(d[m] == 1), rna_het_rate_same_sites=np.mean(v[m] == 1)))
sp = pd.DataFrame(srows)
sp["het_rate_gap_array_minus_rna"] = sp.donor_het_rate_same_sites - sp.rna_het_rate_same_sites
sp["donor_excluded_by_qc"] = [qc.loc[d, "excluded"] if d in qc.index else "NA" for d in sp.donor]
sp["donor_het_F"] = [qc.loc[d, "het_F"] if d in qc.index else np.nan for d in sp.donor]
# baseline = lines whose donors are NOT the focal ones
base = sp[~sp.donor.isin(focal)]
for col in ["conc_at_donor_hom", "conc_at_donor_het", "rna_hom_at_donor_het", "rna_het_at_donor_hom", "het_rate_gap_array_minus_rna"]:
    sp[col + "_z_vs_other_lines"] = (sp[col] - base[col].mean()) / base[col].std()
sp.to_csv(os.path.join(args.outdir, "het_hom_split.tsv"), sep="\t", index=False, float_format="%.5f")
print(sp[["ipsc_line", "donor", "n_donor_hom", "conc_at_donor_hom", "n_donor_het", "conc_at_donor_het", "rna_het_at_donor_hom",
          "donor_het_rate_same_sites", "rna_het_rate_same_sites", "donor_het_F"]].round(4).to_string(index=False), flush=True)

# ------------------------------------------------------------------ figures
# B1 heatmaps
order_d = sorted(donors, key=int)
fig, axes = plt.subplots(2, 1, figsize=(16, 9), gridspec_kw={"height_ratios": [1, 1]})
ylab = [label[i] for i in lines]
for ax, mat, title, cmap in [(axes[0], conc, "Genotype concordance (all SNPs with >= %d RNA reads)" % args.min_dp, "viridis"),
                              (axes[1], ibs0, "IBS0 rate (one genotype 0/0, the other 1/1; lower = more similar)", "magma_r")]:
    M = mat[order_d].astype(float).to_numpy()
    im = ax.imshow(M, aspect="auto", cmap=cmap, interpolation="nearest")
    ax.set_yticks(range(len(lines))); ax.set_yticklabels(ylab, fontsize=7)
    ax.set_xticks(range(0, len(order_d), 5)); ax.set_xticklabels(order_d[::5], fontsize=6, rotation=90)
    for i, r in enumerate(rs.itertuples()):
        x = order_d.index(r.esf_id)
        ax.plot(x, i, "o", mfc="none", mec="white" if cmap == "viridis" else "black", ms=7)
    for i, r in enumerate(rs.itertuples()):
        if r.esf_id in focal:
            ax.get_yticklabels()[i].set_color("#c0392b"); ax.get_yticklabels()[i].set_fontweight("bold")
    ax.set_title(title + "   (circle = expected donor)", fontsize=9); plt.colorbar(im, ax=ax, pad=0.01)
axes[1].set_xlabel("donor IID (eSF ID)")
fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "B1_concordance_heatmap.png"), dpi=170)
fig.savefig(os.path.join(fig_dir, "B1_concordance_heatmap.pdf")); plt.close(fig)

# B2 stratified view
fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
x = np.arange(len(sp)); col = ["#c0392b" if d in focal else "#2c7fb8" for d in sp.donor]
axes[0].bar(x, sp.conc_at_donor_hom, color=col); axes[0].set_ylim(sp.conc_at_donor_hom.min() - 0.02, 1.0)
axes[0].set_title("concordance at donor-HOMOZYGOUS sites", fontsize=9)
axes[1].bar(x, sp.conc_at_donor_het, color=col); axes[1].set_ylim(max(0, sp.conc_at_donor_het.min() - 0.1), 1.0)
axes[1].set_title("concordance at donor-HETEROZYGOUS sites (RNA must be 0/1)", fontsize=9)
axes[2].bar(x - 0.2, sp.donor_het_rate_same_sites, width=0.4, color="#7f8c8d", label="donor array het rate")
axes[2].bar(x + 0.2, sp.rna_het_rate_same_sites, width=0.4, color=col, label="RNA het rate (same sites)")
axes[2].set_title("het rate on the same SNPs: array vs RNA", fontsize=9); axes[2].legend(fontsize=7)
for ax in axes:
    ax.set_xticks(x); ax.set_xticklabels([f"{l}\n({d})" for l, d in zip(sp.ipsc_line, sp.donor)], rotation=90, fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "B2_het_hom_split.png"), dpi=170)
fig.savefig(os.path.join(fig_dir, "B2_het_hom_split.pdf")); plt.close(fig)

with open(os.path.join(args.outdir, "report_B.txt"), "w") as fh:
    fh.write("Part B - identity of the iPSC lines against the 140 array donors (RNA-seq genotypes, >= %d reads)\n\n" % args.min_dp)
    fh.write(best.to_string(index=False) + "\n\n")
    fh.write(sp.round(4).to_string(index=False) + "\n")
print("done", flush=True)
