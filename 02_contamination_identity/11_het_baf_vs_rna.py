#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/11_het_baf_vs_rna.py   (original step script: 30_het_baf_vs_rna.py)
# Produces: BAF of array-heterozygous sites confirmed or not by RNA-seq.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 30 - Joins Part A (array intensities) and Part B (RNA genotypes) for the donors of the 14 iPSC lines (diagnostic only).

Question: the array heterozygous calls that the RNA does NOT confirm (array 0/1, RNA homozygous at >= 10 reads) - are they
"false hets" caused by DNA contamination?  If the array DNA is a mixture, a false het comes from a truly homozygous site pushed
towards the middle by the contaminant's alleles, so its BAF is OFF-CENTRE (typically 0.2-0.35 or 0.65-0.8), whereas a het that the RNA
misses for biological/technical reasons (allele-specific expression, mapping bias, low depth) is a true het with BAF ~ 0.5.
For each of the 14 donors we compare |BAF - 0.5| of the array-het sites split into RNA-confirmed (RNA 0/1) and RNA-discordant
(RNA homozygous), and the share of sites with |BAF - 0.5| > 0.15. The 12 donors that pass the QC with F ~ 0 are the baseline.

Site -> BAF link: array VCF site (GRCh38 chr:pos, step 20) -> locus name through the PLINK b38.bim of the project export (same
names as the manifest) -> row of the BPM locus table (step 23b, manifest order) -> BAF in the donor's GTC (step 23). Positions with
more than one bim name, or names not unique in the locus table, are skipped (reported).
Inputs : results/rna_identity/work/{donor_gt.tsv, donor_samples.txt}, rnaseq/{rna_samples.tsv, genotypes/*.calls.tsv.gz},
         PRS/genotypes/plink_bin/b38.bim, PRS/intensity/{bpm_snp_map.csv, gtc/}, results/contamination/sample_intensity_metrics.tsv
Outputs: results/rna_identity/het_baf_vs_rna.tsv, figures/B3_het_baf_rna_confirmed_vs_discordant.{png,pdf}
"""
import argparse, glob, os, struct
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--rna-dir", required=True); ap.add_argument("--rna-samples", required=True)
ap.add_argument("--identity-dir", required=True); ap.add_argument("--bim", required=True)
ap.add_argument("--snp-map", required=True); ap.add_argument("--gtc-dir", required=True)
ap.add_argument("--metrics", required=True); ap.add_argument("--min-dp", type=int, default=10)
args = ap.parse_args()
GT2D = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/0": 1, "1/1": 2, "1|1": 2}


def read_baf(path):
    """BAF array (float32, manifest order) of a GTC v3-v5 file: TOC (uint16 id, int32 offset) -> id 1012 = int32 n + float32[n]"""
    b = open(path, "rb").read()
    n_toc = struct.unpack_from("<i", b, 4)[0]
    toc = dict(struct.unpack_from("<Hi", b, 8 + 6 * k) for k in range(n_toc))
    off = toc[1012]; n = struct.unpack_from("<i", b, off)[0]
    return np.frombuffer(b, dtype="<f4", count=n, offset=off + 4)


# ---- donor genotypes (exonic array SNPs) and site -> manifest row
donors = [l.strip() for l in open(os.path.join(args.identity_dir, "work", "donor_samples.txt")) if l.strip()]
dg = pd.read_csv(os.path.join(args.identity_dir, "work", "donor_gt.tsv"), sep="\t", header=None, dtype=str)
D = np.vstack([dg[4 + j].map(GT2D).fillna(-1).astype(np.int8).to_numpy() for j in range(len(donors))]).T
site_key = (dg[0].str.replace("chr", "", regex=False) + ":" + dg[1]).to_numpy()
bim = pd.read_csv(args.bim, sep=r"\s+", header=None, names=["chr", "name", "cm", "pos", "a1", "a2"], dtype={"chr": str, "pos": str})
bim["key"] = bim.chr + ":" + bim.pos
snpmap = pd.read_csv(args.snp_map, dtype={"Chromosome": str})
row_of = snpmap.reset_index().groupby("Name")["index"].agg(["first", "count"])
cnt = bim.groupby("key").name.agg(["first", "count"])
site_name = pd.Series(site_key).map(cnt["first"]); site_n = pd.Series(site_key).map(cnt["count"])
site_row = site_name.map(row_of["first"]); site_rn = site_name.map(row_of["count"])
ok_site = (site_n == 1) & (site_rn == 1) & site_row.notna()
print(f"exonic array SNPs: {len(site_key):,} | linked to exactly one manifest locus: {int(ok_site.sum()):,}", flush=True)
site_row = site_row.fillna(-1).astype(int).to_numpy(); ok_site = ok_site.to_numpy()

met = pd.read_csv(args.metrics, sep="\t", dtype={"IID": str, "chip": str}).set_index("IID")
rs = pd.read_csv(args.rna_samples, sep="\t", dtype=str)
focal = {x for x in __import__("os").environ.get("FOCAL_SAMPLE_IDS", "").split(",") if x}
rows, store = [], {}
for r in rs.itertuples():
    c = pd.read_csv(os.path.join(args.rna_dir, f"{r.fastq_name}.calls.tsv.gz"), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "alt", "gt", "dp", "ad"], dtype={"chrom": str, "gt": str, "ad": str})
    c = c[c.dp >= args.min_dp]
    key = (c.chrom.str.replace("chr", "", regex=False) + ":" + c.pos.astype(str))
    idx = pd.Series(range(len(site_key)), index=site_key)
    ix = key.map(idx); good = ix.notna() & c["gt"].map(GT2D).notna()
    v = np.full(len(site_key), -1, dtype=np.int8); v[ix[good].astype(int).to_numpy()] = c.loc[good, "gt"].map(GT2D).to_numpy().astype(np.int8)
    j = donors.index(r.esf_id); d = D[:, j]
    baf = read_baf(glob.glob(os.path.join(args.gtc_dir, "*", met.loc[r.esf_id, "chip_position"] + ".gtc"))[0])
    het = (d == 1) & (v >= 0) & ok_site
    conf = het & (v == 1); disc = het & (v != 1)
    dev_c = np.abs(baf[site_row[conf]] - 0.5); dev_d = np.abs(baf[site_row[disc]] - 0.5)
    U, p = mannwhitneyu(dev_d, dev_c, alternative="two-sided")
    rows.append(dict(ipsc_line=r.ipsc_line, donor=r.esf_id, donor_F=met.loc[r.esf_id, "het_F"], donor_alpha_hat=met.loc[r.esf_id, "alpha_hat"],
                     n_het_rna_confirmed=int(conf.sum()), n_het_rna_discordant=int(disc.sum()),
                     median_absdev_confirmed=np.median(dev_c), median_absdev_discordant=np.median(dev_d),
                     frac_off_centre_confirmed=np.mean(dev_c > 0.15), frac_off_centre_discordant=np.mean(dev_d > 0.15),
                     p_mannwhitney=p))
    store[r.esf_id] = (baf[site_row[conf]], baf[site_row[disc]])
t = pd.DataFrame(rows)
base = t[~t.donor.isin(focal)]
t["frac_off_centre_discordant_z_vs_baseline"] = (t.frac_off_centre_discordant - base.frac_off_centre_discordant.mean()) / base.frac_off_centre_discordant.std()
t.to_csv(os.path.join(args.identity_dir, "het_baf_vs_rna.tsv"), sep="\t", index=False, float_format="%.5g")
print(t.round(4).to_string(index=False), flush=True)

show = sorted(focal) + [x for x in __import__("os").environ.get("COMPARISON_SAMPLE_IDS", "").split(",") if x]
fig, axes = plt.subplots(1, 4, figsize=(16, 3.8), sharey=True)
for ax, d in zip(axes, show):
    c, x = store[d]
    ax.hist(c, bins=40, range=(0, 1), density=True, alpha=.6, color="#2c7fb8", label=f"RNA-confirmed het (n={len(c)})")
    ax.hist(x, bins=40, range=(0, 1), density=True, alpha=.6, color="#c0392b", label=f"RNA homozygous (n={len(x)})")
    ax.set_title(f"donor {d} | F={met.loc[d, 'het_F']:.3f} | alpha_hat={met.loc[d, 'alpha_hat'] * 100:.1f} %", fontsize=9)
    ax.set_xlabel("array BAF at array-het sites"); ax.legend(fontsize=7)
axes[0].set_ylabel("density")
fig.tight_layout(); os.makedirs(os.path.join(args.identity_dir, "figures"), exist_ok=True)
fig.savefig(os.path.join(args.identity_dir, "figures", "B3_het_baf_rna_confirmed_vs_discordant.png"), dpi=170)
fig.savefig(os.path.join(args.identity_dir, "figures", "B3_het_baf_rna_confirmed_vs_discordant.pdf")); plt.close(fig)
print("done", flush=True)
