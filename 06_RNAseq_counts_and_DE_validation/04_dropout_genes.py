#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/04_dropout_genes.py   (original step script: 64b_dropout_genes.py)
# Produces: Genes lost in one of the two count tables.
# Run: python 06_RNAseq_counts_and_DE_validation/04_dropout_genes.py <args>   (or from a SLURM script; see the usage line below)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 64b - which genes are (nearly) absent from the original counts although well expressed in ours (e.g. TUBB 52,537 vs 0.18) and vice versa? Uses compare_per_gene.tsv.gz (step 64) and the
GENCODE v50 gene positions (65_deg_run: gene_pos.tsv). Lists the genes with mean normalised count >= 100 in one table and < 5 in the other (common symbols), their chromosomes
(relative to the share of all common genes) and gene types, to test the explanation "reads lost as multi-mappers/ambiguous in the original pipeline" (e.g. genes with copies on ALT/patch
contigs or overlapping antisense genes) instead of a normalisation problem.
Usage: python 06_RNAseq_counts_and_DE_validation/04_dropout_genes.py <compare_per_gene.tsv.gz> <gene_pos.tsv> <out_dir>
"""
import os, sys
import pandas as pd

cp, gp, out = sys.argv[1:4]
g = pd.read_csv(cp, sep="\t", index_col=0); pos = pd.read_csv(gp, sep="\t").drop_duplicates("gene_name").set_index("gene_name")
g = g.join(pos[["chr", "gene_type"]], how="left")
lost = g[(g.mean_norm_ours >= 100) & (g.mean_norm_theirs < 5)].sort_values("mean_norm_ours", ascending=False)
gain = g[(g.mean_norm_theirs >= 100) & (g.mean_norm_ours < 5)].sort_values("mean_norm_theirs", ascending=False)
lost.to_csv(os.path.join(out, "dropout_genes_lost_in_hers.tsv"), sep="\t", float_format="%.4g"); gain.to_csv(os.path.join(out, "dropout_genes_lost_in_ours.tsv"), sep="\t", float_format="%.4g")
base = g.chr.value_counts(normalize=True)
rows = []
for nm, d in [("expressed in ours, ~0 in the original", lost), ("expressed in the original, ~0 in ours", gain)]:
    sh = d.chr.value_counts(normalize=True)
    rows.append(f"{nm}: {len(d)} genes; top chromosomes {dict((d.chr.value_counts().head(6)).items())}; chr6 share {100 * sh.get('chr6', 0):.1f} % (all common genes {100 * base.get('chr6', 0):.1f} %); "
                f"gene types {dict(d.gene_type.value_counts().head(4).items())}")
rows.append("examples lost in the original: " + ", ".join(f"{i} ({r.mean_norm_ours:.0f} vs {r.mean_norm_theirs:.1f}, {r.chr})" for i, r in lost.head(25).iterrows()))
rows.append("examples lost in ours: " + ", ".join(f"{i} ({r.mean_norm_ours:.1f} vs {r.mean_norm_theirs:.0f}, {r.chr})" for i, r in gain.head(15).iterrows()))
rows.append(f"genes with mean norm >= 100 in both tables: share with |mean log2 diff| > 1: {100 * ((g.mean_norm_ours >= 100) & (g.mean_norm_theirs >= 100) & (g.mean_diff.abs() > 1)).sum() / ((g.mean_norm_ours >= 100) & (g.mean_norm_theirs >= 100)).sum():.2f} %")
open(os.path.join(out, "dropout_genes_report.txt"), "w").write("\n".join(rows) + "\n"); print("\n".join(rows))
