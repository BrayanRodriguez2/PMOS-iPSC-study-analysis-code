#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/03_min_direction.py   (original step script: 81_min_direction.py)
# Produces: Direction check of 20 genes named in a published PCOS iPSC study.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 81 - direction check of the 20 genes named in Tables 2 and 3 of Min et al. 2018 (Stem Cell Res Ther 9:210; top 10 genes per direction, PCOS vs control iPSCs, 3 vs 3, microarray).
The COMPLETE lists of differential genes of Min et al. 2018 and 2019 are NOT public (supplements: primers, figures, methods; no GEO deposit), so no GSEA is possible and this check is limited
to these 20 genes. It is DESCRIPTIVE: 20 pre-selected extreme genes cannot test the signature.
  down in PCOS (Min): FN1, NTS, CER1, SPP1, SLC7A3, SLC2A3, ZFP42, HAS2, PTPRZ1, CDH1
  up in PCOS (Min):   IFI16, CAPN6, LAMA4, IL18, FOLH1B, FOLH1, TBX5, FBP1, AGL, KIAA1324 (= ELAPOR1 in current HGNC/GENCODE)
For every gene: present in GENCODE v50?; expressed in our lines (section-28 filter: >= 10 normalised counts in >= 7 of 14 lines; stranded counts); log2FC (PMOS - control) and p from (a) our OLS on
log2(norm + 1), (b) DESeq2 on raw counts, (c) the original analysis' Table S4 (the original counts), each with the 14 lines and WITHOUT HEL350.3 (OLS refitted; DESeq2 of step 73a; S4 only with 14 lines); concordance of the
direction with Min. Summary: concordant / expressed genes with a two-sided exact binomial test against 50 % (descriptive), genes with nominal p < 0.05 and FDR < 0.1 (BH over all tested genes of each table).
Usage: python 08_preprint_tasks/03_min_direction.py <deg_results_dir> <gene_pos.tsv> <S4 tsv> <out_dir>
"""
import os, re, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R, GP, S4F, out = sys.argv[1:5]; os.makedirs(f"{out}/figures", exist_ok=True)
DOWN = ["FN1", "NTS", "CER1", "SPP1", "SLC7A3", "SLC2A3", "ZFP42", "HAS2", "PTPRZ1", "CDH1"]
UP = ["IFI16", "CAPN6", "LAMA4", "IL18", "FOLH1B", "FOLH1", "TBX5", "FBP1", "AGL", "KIAA1324"]
ALIAS = {"KIAA1324": "ELAPOR1"}
rep = []
log = lambda m="": (print(m), rep.append(str(m)))


def hel(x):
    m = re.search(r"HEL\s*_?(\d+)[._](\d+)", x); n = f"HEL{m.group(1)}.{m.group(2)}"; return "HEL362.1" if n == "HEL362.3" else n


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1); q = np.minimum.accumulate(q[::-1])[::-1]; r = np.empty(n); r[o] = np.minimum(q, 1); return r


def ols_group(Y, g):
    n = len(g); X1 = np.column_stack([np.ones(n), g]); XtXi = np.linalg.inv(X1.T @ X1); B = Y @ X1 @ XtXi.T; E = Y - B @ X1.T; df = n - 2
    se = np.sqrt(((E ** 2).sum(axis=1) / df) * XtXi[1, 1]); t = B[:, 1] / se
    return B[:, 1], 2 * stats.t.sf(np.abs(t), df)


meta = pd.read_csv(f"{R}/clinical/donor_meta.tsv", sep="\t", index_col=0)
N = pd.read_csv(f"{R}/step1/norm_counts_symbol.tsv.gz", sep="\t", index_col=0); N.columns = [hel(c) for c in N.columns]
gp = pd.read_csv(GP, sep="\t"); gencode = set(gp.gene_name)
lines = list(meta.sort_values(["group", "bmi"]).index); N = N[lines]
exp_mask = (N >= 10).sum(axis=1) >= 7; expr = N[exp_mask]; L = np.log2(expr + 1)
g14 = (meta.loc[lines, "group"] == "PMOS").astype(float).values
b14, p14 = ols_group(L.values, g14); keep13 = [i for i, l in enumerate(lines) if l != "HEL350.3"]; b13, p13 = ols_group(L.values[:, keep13], g14[keep13])
ols14 = pd.DataFrame({"log2FC": b14, "p": p14, "q": bh(p14)}, index=L.index); ols13 = pd.DataFrame({"log2FC": b13, "p": p13, "q": bh(p13)}, index=L.index)
d14 = pd.read_csv(f"{R}/methods/deseq2_group.tsv", sep="\t").set_index("gene"); d13 = pd.read_csv(f"{R}/methods_no350/deseq2_group.tsv", sep="\t").set_index("gene")
for d in (d14, d13): d["q"] = bh(d.p.fillna(1))
s4 = pd.read_csv(S4F, sep="\t").set_index("gene")
rows = []
for direction, genes in [("down", DOWN), ("up", UP)]:
    for gname in genes:
        sym = ALIAS.get(gname, gname); r = dict(min_gene=gname, symbol_used=sym, min_direction_in_PCOS=direction, in_GENCODE_v50=sym in gencode)
        if sym in N.index:
            r.update(mean_norm_count=float(N.loc[sym].mean()), lines_ge10=int((N.loc[sym] >= 10).sum()), expressed=bool(exp_mask[sym]))
        else:
            r.update(mean_norm_count=np.nan, lines_ge10=np.nan, expressed=False)
        # descriptive OLS for any gene present in the counts (also not expressed)
        if sym in N.index:
            y = np.log2(N.loc[sym, lines].values + 1)[None, :]; b, p = ols_group(y, g14); r.update(OLS14_log2FC_any=b[0], OLS14_p_any=p[0])
        for tag, tab in [("OLS14", ols14), ("OLS13", ols13), ("DESeq2_14", d14), ("DESeq2_13", d13)]:
            if sym in tab.index and r["expressed"]:
                r[f"{tag}_log2FC"] = tab.loc[sym, "log2FC"]; r[f"{tag}_p"] = tab.loc[sym, "p"]; r[f"{tag}_q"] = tab.loc[sym, "q"]
        if sym in s4.index:
            r.update(S4_log2FC=s4.loc[sym, "log2FC_PMOS_vs_control"], S4_p=s4.loc[sym, "p"], S4_q=s4.loc[sym, "q_BH"])
        rows.append(r)
T = pd.DataFrame(rows)
for tag in ["OLS14", "OLS13", "DESeq2_14", "DESeq2_13", "S4"]:
    col = f"{tag}_log2FC"; sgn = np.sign(T[col]); want = np.where(T.min_direction_in_PCOS == "up", 1, -1)
    T[f"{tag}_agrees_with_Min"] = np.where(T[col].isna(), "", np.where(sgn == want, "yes", "no"))
T.to_csv(f"{out}/min_20_genes_direction.tsv", sep="\t", index=False, float_format="%.4g")
log(f"Genes present in GENCODE v50 (after alias KIAA1324 -> ELAPOR1): {int(T.in_GENCODE_v50.sum())} of 20; missing: {list(T.min_gene[~T.in_GENCODE_v50])}")
log(f"Expressed in our lines (section-28 filter): {int(T.expressed.sum())} of 20: up {int(T[T.min_direction_in_PCOS == 'up'].expressed.sum())}/10, down {int(T[T.min_direction_in_PCOS == 'down'].expressed.sum())}/10; "
    f"not expressed: {list(T.min_gene[~T.expressed])}")
srows = []
for tag, lab in [("OLS14", "our OLS, 14 lines"), ("OLS13", "our OLS, without HEL350.3"), ("DESeq2_14", "our DESeq2, 14 lines"), ("DESeq2_13", "our DESeq2, without HEL350.3"), ("S4", "the manuscript authors Table S4 (the original counts)")]:
    d = T[T[f"{tag}_agrees_with_Min"] != ""]; k = int((d[f"{tag}_agrees_with_Min"] == "yes").sum()); n = len(d)
    bt = stats.binomtest(k, n, 0.5).pvalue if n else np.nan
    nom = d[d[f"{tag}_p"] < 0.05]; fdr = d[d[f"{tag}_q"] < 0.1] if f"{tag}_q" in d else d.iloc[0:0]
    srows.append(dict(analysis=lab, genes_tested=n, same_direction_as_Min=k, pct=100 * k / n if n else np.nan, binomial_p_vs_50pct=bt, nominal_p_lt_0_05=len(nom),
                      nominal_genes="; ".join(f"{g} ({'concordant' if a == 'yes' else 'discordant'})" for g, a in zip(nom.min_gene, nom[f"{tag}_agrees_with_Min"])), FDR_lt_0_1=len(fdr)))
S = pd.DataFrame(srows); S.to_csv(f"{out}/min_20_genes_summary.tsv", sep="\t", index=False, float_format="%.4g")
log(S.round(3).to_string(index=False))
for direction in ["up", "down"]:
    d = T[(T.min_direction_in_PCOS == direction) & T.expressed]
    k = int((d.OLS14_agrees_with_Min == "yes").sum()); log(f"  Min '{direction}' genes (expressed, n = {len(d)}): same direction in our OLS {k}/{len(d)}")
open(f"{out}/min_20_genes_report.txt", "w").write("\n".join(rep) + "\n")
# figure: log2FC of the expressed genes (OLS 14, OLS 13, DESeq2 14, S4) coloured by Min direction
E = T[T.expressed].copy(); E["order"] = E.min_direction_in_PCOS.map({"down": 0, "up": 1}); E = E.sort_values(["order", "OLS14_log2FC"])
fig, ax = plt.subplots(figsize=(8, 0.35 * len(E) + 1.5)); ypos = np.arange(len(E))
for tag, mk, off, col in [("OLS14", "o", -0.2, "#1f77b4"), ("OLS13", "s", 0.0, "#6baed6"), ("DESeq2_14", "^", 0.2, "#d62728"), ("S4", "x", 0.4, "#555555")]:
    ax.scatter(E[f"{tag}_log2FC"], ypos + off, marker=mk, c=col, s=26, label={"OLS14": "our OLS (14 lines)", "OLS13": "our OLS without HEL350.3", "DESeq2_14": "our DESeq2 (14)", "S4": "the original Table S4"}[tag])
ax.axvline(0, color="k", lw=.8); ax.set_yticks(ypos); ax.set_yticklabels([f"{g} (Min: {d})" for g, d in zip(E.min_gene, E.min_direction_in_PCOS)], fontsize=7)
ax.set_xlabel("log2 fold change, PMOS vs control (positive = higher in PMOS)"); ax.legend(frameon=False, fontsize=7, loc="lower right"); ax.set_title("Min et al. 2018, Tables 2-3: direction in our iPSC lines (expressed genes only)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{out}/figures/M1_min_20_genes_direction.png", dpi=160)
