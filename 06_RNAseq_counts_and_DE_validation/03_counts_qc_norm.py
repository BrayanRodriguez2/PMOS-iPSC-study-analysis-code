#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/03_counts_qc_norm.py   (original step script: 64_counts_qc_norm.py)
# Produces: Library orientation, per-line QC, size-factor normalisation, comparison with the external file.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 64 - own gene counts (featureCounts, step 63) -> library orientation, per-line QC table, DESeq2-style size-factor normalisation (median of ratios, as Galaxy's
DESeq2) and comparison with the original analysis' Normalized_counts.tabular.

1. Strandedness: fraction of reads assigned to genes with featureCounts -s 0 / 1 / 2 (single-end). The setting with by far the highest fraction is the library
   orientation (unstranded libraries give ~ equal fractions for 0, and about half of it for -s 1 and -s 2).
2. QC table per line (STAR Log.final.out + featureCounts summary of the chosen setting): input reads, uniquely mapped %, multi-mapped %, unmapped %,
   assigned % (of input reads and of the uniquely mapped), counts assigned, genes detected (>= 1 and >= 10 raw counts), % of counts on mitochondrial genes, rRNA genes and
   protein-coding genes; a line is flagged when a metric lies > 3 robust SD (median/MAD) from the other lines or an absolute red line is crossed.
3. Size factors: DESeq2's default "ratio" estimator (geometric mean over genes with non-zero counts in all lines; factor = median of count / geometric mean), normalised
   counts = counts / size factor. Rows are Ensembl gene IDs (no version) with the GENCODE symbol attached.
4. Comparison with the original file (rows = gene SYMBOLS, 62,813 genes): genes in common by symbol (our IDs with the same symbol are summed after normalisation), per-line correlation of
   log2(norm + 1), the 14 x 14 cross-correlation (our line vs the original line: identity check), per-line median scaling ratio, genes with the largest mean difference, and the
   genes least correlated across lines.
Outputs (--out): strandedness.tsv, qc_table.tsv, size_factors.tsv, counts_raw.tsv.gz, norm_counts_gene_id.tsv.gz, norm_counts_symbol.tsv.gz, compare_*.tsv, compare_report.txt, figures/F1_*.png
"""
import argparse, os, re
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--fc-prefix", required=True); ap.add_argument("--strand", type=int, default=None); ap.add_argument("--star-dir", required=True)
ap.add_argument("--rna-samples", required=True); ap.add_argument("--theirs", required=True); ap.add_argument("--their-s1", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(os.path.join(a.out, "figures"), exist_ok=True)
FIG = os.path.join(a.out, "figures")
rep = []


def log(m=""):
    print(m, flush=True); rep.append(str(m))


rs = pd.read_csv(a.rna_samples, sep="\t").set_index("fastq_name")          # fastq_name -> ipsc_line (HEL362_3 is the clone HEL362.1)


def read_fc(s):
    d = pd.read_csv(f"{a.fc_prefix}{s}.txt", sep="\t", comment="#")
    d.columns = [re.sub(r".*/([^/]+)\.bam$", r"\1", c) if c.endswith(".bam") else c for c in d.columns]
    summ = pd.read_csv(f"{a.fc_prefix}{s}.txt.summary", sep="\t")
    summ.columns = [re.sub(r".*/([^/]+)\.bam$", r"\1", c) if c.endswith(".bam") else c for c in summ.columns]
    return d, summ.set_index("Status")


# ------------------------------------------------------------------ 1. strandedness
srows = []
fcs = {}
for s in (0, 1, 2):
    d, summ = read_fc(s); fcs[s] = (d, summ)
    tot = summ.sum(axis=0); asg = summ.loc["Assigned"]
    for l in summ.columns:
        srows.append(dict(strand_setting=s, fastq_name=l, line=rs.loc[l, "ipsc_line"], assigned=int(asg[l]), alignments_counted=int(tot[l]), frac_assigned=asg[l] / tot[l]))
st = pd.DataFrame(srows)
stw = st.pivot(index="line", columns="strand_setting", values="frac_assigned"); stw.columns = [f"s{c}" for c in stw.columns]
stw.to_csv(os.path.join(a.out, "strandedness.tsv"), sep="\t", float_format="%.4f")
m = stw.mean()
log("(1) library orientation: mean fraction of alignments assigned to a gene, featureCounts -s 0 / 1 / 2 (single-end, GENCODE v50, exon, gene_id):")
log("    " + " | ".join(f"-s {c[1]}: {v:.3f}" for c, v in m.items()))
best = int(m.idxmax()[1]); S = a.strand if a.strand is not None else best
log(f"    highest: -s {best}; used for everything below: -s {S}. (unstranded would give similar s1 and s2 ~ half of s0; here s{best} / s0 = {m.max() / m['s0']:.2f}, "
    f"s{'1' if best == 2 else '2'} / s0 = {m[f's{1 if best == 2 else 2}'] / m['s0']:.2f})")

# ------------------------------------------------------------------ 2. QC
fc, summ = fcs[S]
meta = fc[["Geneid", "Chr", "gene_name", "gene_type"]].copy()
samples = [c for c in fc.columns if c in rs.index]
counts = fc[samples].copy(); counts.index = fc.Geneid.str.replace(r"\.\d+(_PAR_Y)?$", lambda m_: "_PAR_Y" if m_.group(1) else "", regex=True)
meta.index = counts.index
counts.columns = [rs.loc[c, "ipsc_line"] for c in samples]
counts.insert(0, "gene_name", meta.gene_name); counts.insert(1, "gene_type", meta.gene_type)
counts.to_csv(os.path.join(a.out, "counts_raw.tsv.gz"), sep="\t")
C = counts.drop(columns=["gene_name", "gene_type"]); lines = list(C.columns)


def star_stats(fq):
    t = open(os.path.join(a.star_dir, fq, "Log.final.out")).read()
    g = lambda pat: float(re.search(pat + r"\s*\|\s*([\d.]+)", t).group(1))
    return dict(input_reads=int(g("Number of input reads")), unique_pct=g("Uniquely mapped reads %"), multi_pct=g("% of reads mapped to multiple loci"),
                too_many_loci_pct=g("% of reads mapped to too many loci"), unmapped_short_pct=g("% of reads unmapped: too short"), unmapped_other_pct=g("% of reads unmapped: other"),
                mismatch_rate_pct=g("Mismatch rate per base, %"))


q = []
isMT = (meta.Chr == "chrM") | meta.gene_type.str.startswith("Mt_"); isrRNA = meta.gene_type.str.contains("rRNA") & ~meta.gene_type.str.startswith("Mt_")
isPC = meta.gene_type == "protein_coding"
for fq in samples:
    ln = rs.loc[fq, "ipsc_line"]; sS = star_stats(fq); col = C[ln]
    asg = summ.loc["Assigned", fq]; noF = summ.loc["Unassigned_NoFeatures", fq]; amb = summ.loc["Unassigned_Ambiguity", fq]
    q.append(dict(line=ln, fastq=fq, **sS, assigned_reads=int(asg), assigned_pct_of_input=100 * asg / sS["input_reads"], assigned_pct_of_unique=100 * asg / (asg + noF + amb),
                  ambiguous_pct_of_unique=100 * amb / (asg + noF + amb), no_feature_pct_of_unique=100 * noF / (asg + noF + amb),
                  genes_ge1=int((col >= 1).sum()), genes_ge10=int((col >= 10).sum()), mito_pct_of_counts=100 * col[isMT.values].sum() / col.sum(),
                  rRNA_pct_of_counts=100 * col[isrRNA.values].sum() / col.sum(), protein_coding_pct_of_counts=100 * col[isPC.values].sum() / col.sum()))
qc = pd.DataFrame(q).set_index("line")
flags = {}
for c in ["input_reads", "unique_pct", "assigned_pct_of_input", "assigned_pct_of_unique", "genes_ge1", "genes_ge10", "mito_pct_of_counts", "rRNA_pct_of_counts", "mismatch_rate_pct"]:
    x = qc[c]; mad = 1.4826 * np.median(np.abs(x - x.median()))
    z = (x - x.median()) / (mad if mad > 0 else x.std())
    for l in x.index[np.abs(z) > 3]:
        flags.setdefault(l, []).append(f"{c} {x[l]:.4g} (robust z {z[l]:+.1f})")
for l in qc.index:
    if qc.loc[l, "unique_pct"] < 70: flags.setdefault(l, []).append("unique % < 70")
    if qc.loc[l, "assigned_pct_of_input"] < 50: flags.setdefault(l, []).append("assigned % of input < 50")
    if qc.loc[l, "input_reads"] < 0.5 * qc.input_reads.median(): flags.setdefault(l, []).append("input reads < 50 % of the median")
qc["flags"] = [("; ".join(flags[l]) if l in flags else "") for l in qc.index]
qc.to_csv(os.path.join(a.out, "qc_table.tsv"), sep="\t", float_format="%.4g")
log("\n(2) QC per line:")
log(qc[["fastq", "input_reads", "unique_pct", "multi_pct", "assigned_pct_of_input", "assigned_pct_of_unique", "genes_ge1", "genes_ge10", "mito_pct_of_counts", "rRNA_pct_of_counts", "flags"]].round(2).to_string())

# ------------------------------------------------------------------ 3. size factors and normalisation
nz = (C > 0).all(axis=1); lg = np.log(C[nz].values.astype(float)); geo = lg.mean(axis=1)
sf = np.exp(np.median(lg - geo[:, None], axis=0)); sf = pd.Series(sf, index=lines)
sf.rename("size_factor").to_csv(os.path.join(a.out, "size_factors.tsv"), sep="\t", float_format="%.6g")
N = C / sf
log(f"\n(3) size factors (median of ratios on {int(nz.sum()):,} genes with counts in all 14 lines): " + ", ".join(f"{l} {v:.3f}" for l, v in sf.items()))
Nout = N.copy(); Nout.insert(0, "gene_name", meta.gene_name); Nout.to_csv(os.path.join(a.out, "norm_counts_gene_id.tsv.gz"), sep="\t", float_format="%.5f")
Ns = N.groupby(meta.gene_name.values).sum()                      # symbol level (sum of Ensembl IDs sharing a symbol)
nsym_dup = int((meta.gene_name.value_counts() > 1).sum())
Ns.to_csv(os.path.join(a.out, "norm_counts_symbol.tsv.gz"), sep="\t", float_format="%.5f")
log(f"    gene IDs {len(N):,} -> unique symbols {len(Ns):,} ({nsym_dup:,} symbols shared by > 1 gene ID, summed)")

# ------------------------------------------------------------------ 4. comparison with the original file
T = pd.read_csv(a.theirs, sep="\t", index_col=0)
def hel(x):
    m_ = re.search(r"HEL\s*_?(\d+)[._](\d+)", x); n = f"HEL{m_.group(1)}.{m_.group(2)}"; return "HEL362.1" if n == "HEL362.3" else n
T.columns = [hel(c) for c in T.columns]
dupT = int(T.index.duplicated().sum()); T = T.groupby(level=0).sum() if dupT else T
common = Ns.index.intersection(T.index)
log(f"\n(4) the original file: {len(T):,} rows (duplicated symbols summed: {dupT}) | our symbols {len(Ns):,} | in common {len(common):,} ({100 * len(common) / len(T):.1f} % of the original, {100 * len(common) / len(Ns):.1f} % of ours)")
onlyT = T.index.difference(Ns.index); onlyO = Ns.index.difference(T.index)
log(f"    only in the original file: {len(onlyT):,} (mean norm count > 10 in the original: {int((T.loc[onlyT].mean(axis=1) > 10).sum())}) | only in ours: {len(onlyO):,} (mean norm count > 10 in ours: {int((Ns.loc[onlyO].mean(axis=1) > 10).sum())})")
top_only = pd.concat([T.loc[onlyT].mean(axis=1).sort_values(ascending=False).head(25).rename("mean_norm_count").to_frame().assign(only_in="theirs"),
                      Ns.loc[onlyO].mean(axis=1).sort_values(ascending=False).head(25).rename("mean_norm_count").to_frame().assign(only_in="ours")])
top_only.to_csv(os.path.join(a.out, "compare_genes_only_in_one.tsv"), sep="\t", float_format="%.3f")
lo, lt = np.log2(Ns.loc[common, lines] + 1), np.log2(T.loc[common, lines] + 1)
expr = ((Ns.loc[common, lines] >= 10).sum(axis=1) >= 7) | ((T.loc[common, lines] >= 10).sum(axis=1) >= 7)
both = ((Ns.loc[common, lines] >= 10).sum(axis=1) >= 7) & ((T.loc[common, lines] >= 10).sum(axis=1) >= 7)
rows = []
for l in lines:
    r_all = np.corrcoef(lo[l], lt[l])[0, 1]; rho_all = stats.spearmanr(lo[l], lt[l])[0]
    r_e = np.corrcoef(lo.loc[both, l], lt.loc[both, l])[0, 1]; rho_e = stats.spearmanr(lo.loc[both, l], lt.loc[both, l])[0]
    ratio = np.median((Ns.loc[common, l][both] + 1) / (T.loc[common, l][both] + 1))
    rows.append(dict(line=l, pearson_log2_all_common=r_all, spearman_all_common=rho_all, pearson_log2_expressed_in_both=r_e, spearman_expressed_in_both=rho_e, median_ratio_ours_over_theirs=ratio))
cl = pd.DataFrame(rows).set_index("line"); cl.to_csv(os.path.join(a.out, "compare_per_line.tsv"), sep="\t", float_format="%.5g")
log(f"    genes expressed (>= 10 in >= 7 lines) in ours or the original: {int(expr.sum()):,}; in both: {int(both.sum()):,}")
log(cl.round(4).to_string())
X = np.corrcoef(lo.loc[both, lines].values.T, lt.loc[both, lines].values.T)[:14, 14:]
Xd = pd.DataFrame(X, index=[f"ours {l}" for l in lines], columns=[f"the original {l}" for l in lines]); Xd.to_csv(os.path.join(a.out, "compare_cross_correlation.tsv"), sep="\t", float_format="%.5f")
diag_best = int(sum(np.argmax(X[i]) == i for i in range(14)))
log(f"    identity check (our line i vs the original line j, genes expressed in both): the best match is the same line for {diag_best}/14 lines; mean diagonal r {np.mean(np.diag(X)):.4f}, mean off-diagonal {np.mean(X[~np.eye(14, dtype=bool)]):.4f}")
# genes: mean difference and per-gene correlation across the 14 lines
md = (lo.mean(axis=1) - lt.mean(axis=1))
cg = pd.Series([np.corrcoef(lo.loc[g].values, lt.loc[g].values)[0, 1] if lo.loc[g].std() > 0 and lt.loc[g].std() > 0 else np.nan for g in common[both.values]], index=common[both.values])
gtab = pd.DataFrame({"mean_log2_ours": lo.mean(axis=1), "mean_log2_theirs": lt.mean(axis=1), "mean_diff": md, "mean_norm_ours": Ns.loc[common, lines].mean(axis=1), "mean_norm_theirs": T.loc[common, lines].mean(axis=1)})
gtab["corr_across_lines"] = cg; gtab["expressed_in_both"] = both
sel = gtab[(gtab.mean_norm_ours >= 10) | (gtab.mean_norm_theirs >= 10)]
sel.reindex(sel.mean_diff.abs().sort_values(ascending=False).index).head(60).to_csv(os.path.join(a.out, "compare_genes_largest_mean_diff.tsv"), sep="\t", float_format="%.4g")
gtab[gtab.expressed_in_both].sort_values("corr_across_lines").head(60).to_csv(os.path.join(a.out, "compare_genes_lowest_corr.tsv"), sep="\t", float_format="%.4g")
gtab.to_csv(os.path.join(a.out, "compare_per_gene.tsv.gz"), sep="\t", float_format="%.4g")
log(f"    expressed in both ({int(both.sum()):,} genes): median |mean log2 diff| {md[both].abs().median():.3f}; share with |diff| > 0.5: {100 * (md[both].abs() > 0.5).mean():.1f} %, > 1: {100 * (md[both].abs() > 1).mean():.1f} %; "
    f"per-gene correlation across the 14 lines: median {cg.median():.3f}, share < 0.9: {100 * (cg < 0.9).mean():.1f} %, < 0.5: {100 * (cg < 0.5).mean():.1f} %")
log("    genes with the largest mean difference (mean norm >= 10 in either):\n" + sel.reindex(sel.mean_diff.abs().sort_values(ascending=False).index).head(15)[["mean_norm_ours", "mean_norm_theirs", "mean_diff"]].round(2).to_string())
# XIST vs the original S1
s1 = pd.read_csv(a.their_s1, sep="\t").set_index("iPSC_line_RNAseq")
xo = Ns.loc["XIST", lines]; xs = s1.loc[lines, "XIST_norm_count"].astype(float)
log(f"    XIST: ours vs the original Table S1 norm count, r (log2) = {np.corrcoef(np.log2(xo + 1), np.log2(xs + 1))[0, 1]:.4f}; eroded (log2 < 8) ours: {[l for l in lines if np.log2(xo[l] + 1) < 8]}; the original: {[l for l in lines if np.log2(xs[l] + 1) < 8]}")
pd.DataFrame({"XIST_ours": xo, "XIST_S1": xs, "XIST_theirs_file": T.loc["XIST", lines]}).to_csv(os.path.join(a.out, "compare_XIST.tsv"), sep="\t", float_format="%.2f")

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
ax[0].bar(range(14), cl.pearson_log2_expressed_in_both, color="#2c7fb8"); ax[0].set_xticks(range(14)); ax[0].set_xticklabels(lines, rotation=90, fontsize=7)
ax[0].set_ylim(min(0.9, cl.pearson_log2_expressed_in_both.min() - 0.01), 1); ax[0].set_ylabel("Pearson r, log2(norm + 1)"); ax[0].set_title(f"per line, {int(both.sum()):,} genes expressed in both", fontsize=9)
l0 = lines[0]; ax[1].hexbin(lt.loc[both, l0], lo.loc[both, l0], gridsize=70, bins="log", cmap="viridis", mincnt=1); ax[1].plot([0, 18], [0, 18], "r:", lw=.8)
ax[1].set_xlabel(f"the original file: log2(norm + 1), {l0}"); ax[1].set_ylabel("our counts: log2(norm + 1)"); ax[1].set_title(f"{l0}: r = {cl.loc[l0, 'pearson_log2_expressed_in_both']:.3f}", fontsize=9)
im = ax[2].imshow(X, vmin=X.min(), vmax=1, cmap="magma"); ax[2].set_xticks(range(14)); ax[2].set_xticklabels(lines, rotation=90, fontsize=6); ax[2].set_yticks(range(14)); ax[2].set_yticklabels(lines, fontsize=6)
ax[2].set_xlabel("the original file"); ax[2].set_ylabel("our counts"); ax[2].set_title("cross-correlation of lines (identity check)", fontsize=9); fig.colorbar(im, ax=ax[2], fraction=0.046)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "F1_counts_vs_theirs.png"), dpi=160); plt.close(fig)
open(os.path.join(a.out, "compare_report.txt"), "w").write("\n".join(rep) + "\n")
