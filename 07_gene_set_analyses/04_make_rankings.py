#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/04_make_rankings.py   (original step script: 73_make_rankings.py)
# Produces: Ranking statistics and custom gene sets for the GSEA.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 73 - ranking statistics and custom gene sets for the GSEA (step 74).
Rankings (one TSV each: gene symbol, stat; positive = higher in PMOS than control; only genes of our expressed set = the 16,495 genes tested in 06_RNAseq_counts_and_DE_validation):
  deseq2_wald        main: Wald statistic of DESeq2 ~ group on RAW stranded counts (06_RNAseq_counts_and_DE_validation, methods/deseq2_group.tsv)
  deseq2_wald_xist   DESeq2 ~ log2 XIST + group (methods/deseq2_group_xist.tsv)
  ols_signed         the original method: sign(log2FC) x -log10 p of the per-gene OLS on log2(norm + 1) (analysis/ours_s2/de_ols.tsv.gz)
  and the same three WITHOUT HEL350.3 (suffix _no350): DESeq2 re-run on 13 lines (73a), OLS re-fitted here on 13 lines (own OLS, numpy; XIST version not needed for the OLS ranking).
Gene sets (long TSV: collection, set, gene):
  GWAS_consensus_33   the 33 consensus genes of the 29 loci (Table S7a of the final package; all 33, whether expressed or not)
  Literature_60       the 60 literature genes (Table S14; one set; the expressed ones enter the ranking)
  Min_2018_up / Min_2018_down / Min_2019_up / Min_2019_down   ONLY if the files gene_lists/min2018_up.txt etc. exist (one symbol per line; NOT invented: see the report)
Usage: python 07_gene_set_analyses/04_make_rankings.py <deg_results_dir> <package_scrt_dir> <gene_lists_dir> <out_dir>
"""
import os, re, sys
import numpy as np, pandas as pd
from scipy import stats

R, PK, GL, out = sys.argv[1:5]
os.makedirs(f"{out}/rankings", exist_ok=True); os.makedirs(f"{out}/gene_sets", exist_ok=True)
rep = []
log = lambda m="": (print(m), rep.append(str(m)))


def hel(x):
    m = re.search(r"HEL\s*_?(\d+)[._](\d+)", x); n = f"HEL{m.group(1)}.{m.group(2)}"; return "HEL362.1" if n == "HEL362.3" else n


def save(name, d):
    d = d.dropna().drop_duplicates("gene")
    d.to_csv(f"{out}/rankings/{name}.tsv", sep="\t", index=False, float_format="%.6g"); log(f"ranking {name}: {len(d):,} genes, stat range {d.stat.min():.2f} .. {d.stat.max():.2f}, ties {int(d.stat.duplicated().sum())}")


for name, f in [("deseq2_wald", "methods/deseq2_group.tsv"), ("deseq2_wald_xist", "methods/deseq2_group_xist.tsv"), ("deseq2_wald_no350", "methods_no350/deseq2_group.tsv"), ("deseq2_wald_xist_no350", "methods_no350/deseq2_group_xist.tsv")]:
    if os.path.exists(f"{R}/{f}"):
        d = pd.read_csv(f"{R}/{f}", sep="\t"); save(name, d[["gene", "stat"]])
    else:
        log(f"(missing {f})")
ols = pd.read_csv(f"{R}/analysis/ours_s2/de_ols.tsv.gz", sep="\t")
save("ols_signed", pd.DataFrame({"gene": ols.gene, "stat": np.sign(ols.log2FC) * -np.log10(ols.p)}))
# OLS without HEL350.3 (own fit on the 13 remaining lines; same expressed genes)
meta = pd.read_csv(f"{R}/clinical/donor_meta.tsv", sep="\t", index_col=0)
N = pd.read_csv(f"{R}/step1/norm_counts_symbol.tsv.gz", sep="\t", index_col=0); N.columns = [hel(c) for c in N.columns]
lines = [l for l in meta.sort_values(["group", "bmi"]).index if l != "HEL350.3"]
Y = np.log2(N.loc[ols.gene, lines].values + 1); g = (meta.loc[lines, "group"] == "PMOS").astype(float).values; n = len(lines)
X1 = np.column_stack([np.ones(n), g]); XtXi = np.linalg.inv(X1.T @ X1); B = Y @ X1 @ XtXi.T; E = Y - B @ X1.T; df = n - 2
se = np.sqrt(((E ** 2).sum(axis=1) / df) * XtXi[1, 1]); t = B[:, 1] / se; p = 2 * stats.t.sf(np.abs(t), df)
save("ols_signed_no350", pd.DataFrame({"gene": ols.gene, "stat": np.sign(B[:, 1]) * -np.log10(p)}))
# gene sets
rows = []
s7a = pd.read_csv(f"{PK}/sheet10_S7a_GWAS_genes.tsv", sep="\t"); rows += [("GWAS_consensus_33", "GWAS_consensus_33", x) for x in s7a.gene.unique()]
s14 = pd.read_csv(f"{PK}/sheet18_S14_literature_genes.tsv", sep="\t", skiprows=1); rows += [("Literature_60", "Literature_60", x) for x in s14.gene.unique()]
log(f"GWAS consensus genes: {s7a.gene.nunique()} | literature genes: {s14.gene.nunique()}")
for nm in ["min2018_up", "min2018_down", "min2019_up", "min2019_down"]:
    f = f"{GL}/{nm}.txt"
    if os.path.exists(f):
        gl = [x.strip() for x in open(f) if x.strip()]; rows += [("Min_et_al", nm, x) for x in gl]; log(f"{nm}: {len(gl)} genes")
    else:
        log(f"Min et al. list {nm}: NOT FOUND ({f}); not used")
pd.DataFrame(rows, columns=["collection", "set", "gene"]).to_csv(f"{out}/gene_sets/custom_gene_sets.tsv", sep="\t", index=False)
open(f"{out}/rankings_and_sets_report.txt", "w").write("\n".join(rep) + "\n")
