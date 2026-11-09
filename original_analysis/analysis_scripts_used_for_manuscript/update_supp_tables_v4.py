"""v4: rebuild sheet S13 (HipSci protein cross-reference, now S13a-e) and refresh the protein column of S14
in Supplementary_Tables (v3 -> v4). Same layout conventions as v3: bold section title rows, bold header row,
blank row between sections. v4 changes are highlighted in light turquoise.
Usage: python3 update_supp_tables_v4.py <v3.xlsx> <v4_results_dir> <out.xlsx>"""
import sys
import pandas as pd, openpyxl
from openpyxl.styles import Font, PatternFill
SRC, V4, OUT = sys.argv[1:4]
TQ = PatternFill("solid", fgColor="FFB4E5F0")
B = Font(bold=True)
wb = openpyxl.load_workbook(SRC)
idx = wb.sheetnames.index("S13_protein_HipSci")
del wb["S13_protein_HipSci"]
ws = wb.create_sheet("S13_protein_HipSci", idx)
 
secs = [
 ("S13a. Protein detected in the HipSci iPSC proteome (Mirauta et al. 2020, eLife 9:e57390; TMT-MS, 202 lines) by expression decile in our lines. "
  "v4: original HipSci data (eLife Figure 1 source data 2-3), matched via HGNC/Ensembl; both denominators shown",
  pd.read_csv(f"{V4}/S13a_protein_by_expression_decile.tsv", sep="\t")),
 ("S13b. PCOS-locus genes and eQTL-tested neighbours: protein detection, HipSci RNA-protein correlation and pQTL/eQTL flags (v4: MYO10, ID4, NEIL2, FANCC, ERBB3 and RPGRIP1L now detected)",
  pd.read_csv(f"{V4}/S13b_protein_PCOS_locus_genes.tsv", sep="\t")),
 ("S13c. Protein-supported candidates: nominal PMOS vs control p < 0.05, protein detected and HipSci RNA-protein r >= 0.3 (49 genes; includes all 35 of v3). None reaches FDR < 0.1",
  pd.read_csv(f"{V4}/S13c_protein_supported_candidates.tsv", sep="\t")),
 ("S13e. GO / pathway enrichment (g:Profiler, g:SCS < 0.05, custom backgrounds; GO:BP/MF/CC, KEGG, Reactome; 4 Nov 2026)",
  pd.read_csv(f"{V4}/S13e_GO_enrichment.tsv", sep="\t")),
 ("S13d. All 16,489 expressed genes with HipSci reference values (hipsci_mean_rna_SD2 in source-data units; r = Pearson correlation of HipSci RNA and log10 protein across 201 lines)",
  pd.read_csv(f"{V4}/S13d_protein_crossref_all_genes.tsv.gz", sep="\t")),
]
r = 1
for title, df in secs:
    c = ws.cell(r, 1, title); c.font = B; c.fill = TQ
    r += 1
    for j, col in enumerate(df.columns, 1):
        ws.cell(r, j, col).font = B
    r += 1
    for row in df.itertuples(index=False):
        for j, v in enumerate(row, 1):
            if isinstance(v, float) and pd.isna(v): continue
            ws.cell(r, j, v.item() if hasattr(v, "item") else v)
        r += 1
    r += 1
ws.column_dimensions["A"].width = 13
 
# S14: refresh protein column, highlight changed cells
s14 = pd.read_csv(f"{V4}/S14_literature_candidate_genes.tsv", sep="\t").set_index("gene")
w14 = wb["S14_literature_genes"]
hdr = [c.value for c in w14[2]]; gcol = hdr.index("gene") + 1; pcol = hdr.index("protein_detected_HipSci") + 1
changed = []
for rr in range(3, w14.max_row + 1):
    g = w14.cell(rr, gcol).value
    if g in s14.index:
        new = bool(s14.loc[g, "protein_detected_HipSci"])
        if w14.cell(rr, pcol).value != new:
            w14.cell(rr, pcol).value = new; w14.cell(rr, pcol).fill = TQ; changed.append(g)
print("S14 changed:", changed)
 
# README
wr = wb["README"]
for rr in range(1, wr.max_row + 1):
    if wr.cell(rr, 1).value == "S13":
        wr.cell(rr, 2).value = "Protein-level cross-reference with the HipSci iPSC proteome, rebuilt from the original HipSci data, with GO enrichment (v4)"
        wr.cell(rr, 2).fill = TQ
n = wr.max_row + 1
wr.cell(n, 1, "v4, 4 Nov 2026").fill = TQ
wr.cell(n, 2, "Changes highlighted in light turquoise. S13 rebuilt from the original HipSci data (Mirauta et al. 2020, eLife Figure 1 source data 2-3, "
              "Figure 2 source data 1-2) instead of values copied from the Sept 2026 supplement (which covered ~5,500 genes); genes matched through HGNC "
              "symbols/previous symbols to Ensembl IDs. Where both versions have a value they agree (RNA-protein r 0.999; protein abundance 1.000; "
              "pQTL/eQTL flags 100%). New S13e: GO/pathway enrichment. S14 protein column refreshed.").fill = TQ
wb.save(OUT)
print("saved", OUT)
