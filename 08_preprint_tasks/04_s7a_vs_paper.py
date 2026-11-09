#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/04_s7a_vs_paper.py   (original step script: 82_s7a_vs_paper.py)
# Produces: Consensus genes of the manuscript table against the GWAS supplementary tables.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 82 - check of the 33 consensus genes of Table S7a (final package of the manuscript) against the supplementary material of Moolhuijsen et al. 2026 (Nat Genet 58:1040):
  ST2  (29 signals; columns SNP, Variant, 'Nearest Gene', 'Consensus Gene (based on biological and literature evidence)')
  ST9  (harmonisation of the bio-informatic gene (ST8) and the literature genes: 'Bio-informatics top gene' and 'Best consensus gene'; part A novel signals, part B replicated signals).
One row per gene of S7a: rsID of the locus, gene in S7a and its 'consensus' label, nearest / consensus gene(s) in ST2, bioinformatic top gene and best consensus gene in ST9, whether the S7a gene
is named in the paper's consensus gene field (token match, '/' separated, e.g. 'IRF1 / RAD50'), and whether it is named in the ST9 free-text description of the locus's candidate genes
('[GENE]:' headers). Locus-level checks: the 29 rsIDs of ST2 against the 29 loci of S7a; number of genes (33) and loci with two genes.
Usage: python 08_preprint_tasks/04_s7a_vs_paper.py <S7a tsv> <ST2 tsv> <ST9 tsv> <out_dir>
"""
import os, re, sys
import pandas as pd

s7af, st2f, st9f, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True)
rep = []
log = lambda m="": (print(m), rep.append(str(m)))
s7a = pd.read_csv(s7af, sep="\t"); s7a["locus_rsid"] = s7a.locus_rsid.str.strip()
raw = [l.rstrip("\n").split("\t") for l in open(st2f, encoding="utf-8")]; hi = next(i for i, r in enumerate(raw) if r and r[0].strip() == "SNP"); hdr = [h.strip() for h in raw[hi]]
st2 = pd.DataFrame([r + [""] * (len(hdr) - len(r)) for r in raw[hi + 1:] if r and r[0].strip().startswith("rs")], columns=hdr)
st2["rsid"] = st2.SNP.str.strip(); ccol = [c for c in st2.columns if c.startswith("Consensus Gene")][0]; ncol = [c for c in st2.columns if c.startswith("Nearest Gene")][0]
st2 = st2.set_index("rsid")
rows9 = {}
for l in open(st9f, encoding="utf-8"):
    f = l.rstrip("\n").split("\t")
    if f and re.match(r"^rs\d+", f[0].strip()) and len(f) >= 7:
        rows9[f[0].strip()] = dict(nearest_ST9=f[2].strip(), bio_top=f[4].strip(), function_text=f[5], best_consensus=f[6].strip())
st9 = pd.DataFrame(rows9).T
log(f"ST2: {len(st2)} signals | ST9: {len(st9)} loci | S7a: {len(s7a)} genes in {s7a.locus_rsid.nunique()} loci ({int((s7a.groupby('locus_rsid').size() > 1).sum())} loci with two genes)")
log(f"loci in S7a not in ST2: {sorted(set(s7a.locus_rsid) - set(st2.index))} | in ST2 not in S7a: {sorted(set(st2.index) - set(s7a.locus_rsid))}")
rows = []
tok = lambda s: [t for t in re.split(r"[/,;]| and ", s.replace("(", " ").replace(")", " ")) if t.strip()]
for r in s7a.itertuples():
    rs = r.locus_rsid; g = r.gene
    c2 = st2.loc[rs, ccol].strip() if rs in st2.index else ""; n2 = st2.loc[rs, ncol].strip() if rs in st2.index else ""
    c9 = st9.loc[rs, "best_consensus"] if rs in st9.index else ""; bt = st9.loc[rs, "bio_top"] if rs in st9.index else ""
    ft = st9.loc[rs, "function_text"] if rs in st9.index else ""
    names9 = [t.strip() for t in tok(c9)]; names2 = [t.strip() for t in tok(c2)]
    in_cons = g in names9 or g in names2 or any(t.strip() in (g, g.rstrip("0123456789")) and len(t.strip()) >= 3 and t.strip() == g for t in names9)
    in_text = bool(re.search(r"\[" + re.escape(g) + r"\]|\b" + re.escape(g) + r"\b", ft))
    if in_cons:
        verdict = "match"
    elif c9 and g.startswith(c9.split()[0].rstrip("123456789")) and "complex" in c9.lower():
        verdict = "DIFFERENCE: the paper gives a complex ('" + c9 + "'); S7a picks specific gene(s)"
    else:
        verdict = "DIFFERENCE"
    rows.append(dict(locus_rsid=rs, S7a_consensus_label=r.consensus, S7a_gene=g, ST2_nearest_gene=n2, ST2_consensus_gene=c2, ST9_best_consensus_gene=c9, ST9_bioinformatics_top_gene=bt,
                     S7a_gene_in_paper_consensus_field=in_cons, S7a_gene_named_in_ST9_text=in_text, verdict=verdict))
T = pd.DataFrame(rows); T.to_csv(f"{out}/s7a_vs_moolhuijsen_genes.tsv", sep="\t", index=False)
log(f"S7a genes whose symbol appears in the paper's consensus-gene field (ST2 or ST9): {int(T.S7a_gene_in_paper_consensus_field.sum())} of {len(T)}")
log("differences:\n" + (T[T.verdict != "match"][["locus_rsid", "S7a_consensus_label", "S7a_gene", "ST2_consensus_gene", "ST9_best_consensus_gene", "S7a_gene_named_in_ST9_text", "verdict"]].to_string(index=False) or "none"))
# loci where the paper's consensus names more genes than S7a includes
miss = []
for rs in st9.index:
    pg = [t.strip() for t in tok(st9.loc[rs, "best_consensus"])]; sg = list(s7a[s7a.locus_rsid == rs].gene)
    nm = [x for x in pg if x not in sg and "complex" not in x.lower()]
    if nm or (not sg): miss.append((rs, st9.loc[rs, "best_consensus"], sg))
log(f"loci where the paper's consensus gene field names a gene that is not in S7a: {miss or 'none'}")
log("consensus label vs paper consensus (locus level, S7a 'consensus' column vs ST9 best consensus gene): "
    + ("identical for all 29 loci" if all(str(a).strip() == str(b).strip() for a, b in zip(s7a.drop_duplicates("locus_rsid").sort_values("locus_rsid").consensus, st9.loc[s7a.drop_duplicates("locus_rsid").sort_values("locus_rsid").locus_rsid].best_consensus)) else
       "different at: " + "; ".join(f"{rs}: S7a '{c}' vs ST9 '{st9.loc[rs, 'best_consensus']}'" for rs, c in s7a.drop_duplicates("locus_rsid")[["locus_rsid", "consensus"]].values if str(c).strip() != st9.loc[rs, "best_consensus"].strip())))
open(f"{out}/s7a_vs_moolhuijsen_report.txt", "w").write("\n".join(rep) + "\n")
