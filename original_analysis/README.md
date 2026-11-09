# Original analysis of the manuscript

This folder holds the code that produced the results, figures and tables of the manuscript (Donor-associated variation outweighs PMOS status in patient-derived iPSC transcriptomes: a genotyped, BMI-adjusted study of 14 lines, Masuma Khatun, Brayan Stiveen Rodriguez Morales, Vijayachitra Modhukur, Apostol Apostol, Karolina Lundin, Riikka K. Arffman, Timo Tuuri, Terhi T. Piltonen, Juha S. Tapanainen, Andres Salumets). It is the code of Additional file 3 ("Detailed methods and code", version of 6 November 2026), reproduced here as files.

**Scripts are published unmodified, with one exception requested by the author: in `analysis_scripts_used_for_manuscript/figs.py` (line 56) the array IDs of a pair of related samples omitted from the genotype-PC panel were replaced by the environment variable `RELATED_PAIR_SAMPLE_IDS` (comma-separated, empty by default; IDs held by the authors; without it the panel shows all genotyped samples). This is the only modification of the scripts. They otherwise contain array sample IDs used in the original analysis.**

Because the files are otherwise unchanged:

- Several scripts contain **absolute paths of the original analysis session** (for example `/home/claude/ms`, `/home/claude/pheno`, `/home/claude/geno`, `/mnt/user-data/attachments/...`). To reproduce an analysis, edit those paths or recreate the same folder layout. The later scripts (`an_restore.py`, `an_protein_v4.py`, `figs_v3.py`, `figs_v4_protein.py`, `verify_key_numbers.py`, `update_supp_tables_v4.py`, `go_enrichment_v4.py`) take the package folder and the count matrix as arguments; see the usage line at the top of each script.
- `config.sh` is a template (`/scratch/project_XXXXXXX/...`); the genotype scripts in `scripts/` read their settings from it.
- The input data are **not** included (normalised counts, donor key and phenotypes, genotypes, HipSci tables, i2QTL results). See Section 2 of Additional file 3 for their sources.
- The code is provided as text extracted from Additional file 3; line endings and the trailing blank line of each file may differ from the original files.

Software (Additional file 3): Python 3.13 (NumPy, SciPy, pandas, Matplotlib, openpyxl), PLINK 1.9/2 and bcftools for the genotype pipeline, g:Profiler web service (queried on 4 November 2026). Random seeds: trait permutations 42 (`an2.py`); gene-set null 7 (`an3.py`); variance-partition null 2026 (`an_restore.py`); figure jitter 1-5.

## Map: script to figure / table

| Script | What it does | Figures / tables |
|---|---|---|
| `scripts/common.py` | Shared helpers: count-matrix reader, line-name harmonisation (`hel_name`), per-gene OLS (`ols`), Benjamini-Hochberg (`bh`). | used by the transcriptome scripts |
| `analysis_scripts_used_for_manuscript/an1.py` | Expressed-gene filter (>= 10 normalised counts in >= 7 of 14 lines), gene-wise PMOS models with and without XIST, trait associations, principal components. | Fig 3A-C, S1, S4; Tables S4, S5a, S5c |
| `.../an2.py` | Transcriptome-wide permutation tests (exact 3,432 splits for PMOS; 2,000 permutations for continuous traits). | Fig 3D; Table S5b |
| `.../an3.py` | Expression of the 33 GWAS consensus genes and gene-set score against expression-matched random sets. | Fig 5; Tables S7a, S7b |
| `.../an_restore.py` | Sample-to-sample correlation with exact permutation, variance partition, BMI adjustment and strata, literature candidates, top 40 genes. | Fig 3B, 4, S2; Tables S9, S10, S11, S12, S14, S15 (Table 1 together with `figs.py`) |
| `.../an_protein_v4.py` | Cross-reference with the HipSci iPSC proteome (gene matching, RNA-protein correlation, protein-supported candidates). | Fig 5 (protein markers), S3; Tables S13a-d, S14 (protein column) |
| `.../go_enrichment_v4.py`, `.../gprofiler_browser_call_v4.js` | g:Profiler GO / KEGG / Reactome enrichment (Python version and the browser version that was actually used). | Fig S3C; Table S13e |
| `.../bpm.py` | Probe positions from the Illumina BPM manifest (used to match lead variants to the array). | Fig 2G, S2, S3 (with `01_genotype_qc.sh`); Fig 6A-D, Table S6 |
| `scripts/01_genotype_qc.sh` | First-pass genotype QC pipeline of the original analysis (thresholds superseded by the final QC described in Section 3.14 of Additional file 3). | Fig 2G, S2, S3 (as listed in Section 4 of Additional file 3) |
| `scripts/08a_extract_dosages.sh`, `scripts/08b_targeted_eqtl.py` | Targeted cis-eQTL at the genotyped lead variants. | Fig 6A-D; Table S6 |
| `.../fig6_prs_panels.py` | Polygenic-score panels. | Fig 6E-F |
| `.../power_detectable_effect.py` | Smallest detectable log2 fold change at 80% power. | Table S17 |
| `.../figs.py` | Figures 1, 2 and 6. | Fig 1, 2, 6 |
| `.../figs_v3.py` | Figures 3, 4, S1 and S2. | Fig 3, 4, S1, S2 |
| `.../figs_v4_protein.py` | Figure 5 (protein markers) and Figure S3. | Fig 5, S3 |
| `.../update_supp_tables_v4.py` | Rebuilds sheet S13 and the S14 protein column in the supplementary workbook. | Tables S13, S14 |
| `.../verify_key_numbers.py` | Independent re-implementation (NumPy/SciPy only) of 30 key numbers of the manuscript. | verification |
| `.../fix_erbb3_v5_1.py` | Correction of the ERBB3 diagnosis-adjusted estimate (refitted with XIST), found by `verify_key_numbers.py`. | corrected value in the manuscript (v5.1) |

Table 1 and Fig 2A-F are produced by `figs.py` and `an_restore.py`. Table S8 (PMOS test after excluding lines) was produced by the statistic of `an2.py` on line subsets; Additional file 3 states that its exact code was not preserved in the package. Figure 1 and Fig 6E-F, Table 3, Tables S16 and S17 also rely on the genotype QC, imputation and PRS-CS analyses described in Section 3.14, whose SLURM scripts are in the numbered folders of the parent directory.

## Not included

Scripts that assemble manuscript documents rather than analyse data (cover letter, manuscript text edits, verification guide, checklist, figure book, Additional file 2 builder) and the early pipeline scripts that did not produce results cited in the manuscript or were replaced by the final analysis (`02_prepare_imputation.sh`, `03_post_imputation.sh`, `04a`-`04c`, `05_prs_summary.py`, `06a`, `06b`, `07_coloc.R`, `09_prs_expression.py`). They are in Additional file 3.

## Relationship to the independent validation

The numbered folders of the parent directory contain a second, independent analysis written separately from this one (see the main README):

- `06_RNAseq_counts_and_DE_validation/` re-quantifies the RNA-seq reads (STAR, featureCounts) and re-implements the gene-wise analyses of `an1.py`, `an2.py` and `an_restore.py`; `08_external_pipeline_check.*` runs copies of `common.py`, `an1.py`, `an2.py` and sections 1-4 of `an_restore.py` of this folder on the original counts (only the absolute paths of the original session were replaced and `an_restore.py` was cut after section 4); `10_verdicts.py` compares the numbers of the manuscript with both analyses. `09_deseq2_limma.*` adds DESeq2 and limma-voom on raw counts.
- `07_gene_set_analyses/` replicates the gene-set score of `an3.py` (Fig 5B, Table S7b) and adds fgsea analyses.
- `01_` to `05_` repeat the genotype QC, imputation, post-imputation QC and PRS-CS analysis summarised in Section 3.14 of Additional file 3.
- `08_preprint_tasks/` contains the array PCA behind Fig 2G and the updated Tables S2 and S3.

## Licence

Same as the rest of the repository (see `../LICENSE`).
