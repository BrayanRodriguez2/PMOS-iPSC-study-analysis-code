# PMOS iPSC study: analysis code

Code for the genotype, imputation, polygenic-score and RNA-seq analyses of a study of undifferentiated iPSC lines derived from women with polyendocrine metabolic ovarian syndrome (PMOS; formerly polycystic ovary syndrome, PCOS) and controls, with donor genotypes from an Illumina GSA array.

**Related manuscript:** Donor-associated variation outweighs PMOS status in patient-derived iPSC transcriptomes: a genotyped, BMI-adjusted study of 14 lines, Masuma Khatun, Brayan Stiveen Rodriguez Morales, Vijayachitra Modhukur, Apostol Apostol, Karolina Lundin, Riikka K. Arffman, Timo Tuuri, Terhi T. Piltonen, Juha S. Tapanainen, Andres Salumets. Manuscript in preparation. Citation will be added upon publication. Journal: intended for submission to Stem Cell Research & Therapy.

**Version:** 1.0.1 (2026-11-09). Changes since 1.0.0: one line of `original_analysis/` (figs.py) now reads the IDs of a pair of samples from an environment variable instead of listing them (see `original_analysis/README.md`); the sample-ID variable `RELATED_PAIR_SAMPLE_IDS` was added to `config/config.example.sh`; comments mentioning that pair and a log message listing four array-ID-to-line pairs (02_array_pca_report.py) were made generic; version and ORCID added to `CITATION.cff`.

## Two parts: original analysis and independent validation

This repository has two complementary parts.

1. **Original analysis (`original_analysis/`).** The code that produced the results, figures and tables of the manuscript, published as reproduced in Additional file 3 of the manuscript (unmodified except for one change requested by the author, described in its README): gene-wise expression models, permutation tests, variance partition, gene-set score, HipSci protein cross-reference and GO enrichment, targeted cis-eQTL, figure scripts, detectable-effect analysis and an independent check of key numbers. It has its own `README.md` with a map from each script to the figures and tables it produced. Absolute paths in it are those of the original analysis session and must be edited to reproduce it.
2. **Independent validation (folders `01_` to `08_`).** A second analysis, written separately, that re-does the genotype QC, contamination and identity checks, TOPMed imputation and post-imputation QC, PRS-CS polygenic scores with the imputed data, a re-implementation of the expression analyses on re-quantified RNA-seq counts (including DESeq2 and limma-voom), gene-set/GSEA analyses and several preprint checks. `06_RNAseq_counts_and_DE_validation/` and `08_preprint_tasks/` refer to the supplementary tables of the manuscript; `06_.../08_external_pipeline_check.*` runs copies of the original scripts of `original_analysis/` on the original counts.

## What is in this repository

| Folder | Content | Run order |
|---|---|---|
| `original_analysis/` | Original analysis code of the manuscript (from Additional file 3; one change requested by the author, see its README). | see its README |
| `00_common/` | Small utilities used by several steps (xlsx to TSV, PDF to text). | any |
| `01_genotype_QC/` | Array genotype QC: PED/MAP to PLINK, sample QC (call rate, heterozygosity, sex, kinship), SNP QC, rsID renaming, QC by chip, QC figures. | 1 |
| `02_contamination_identity/` | DNA-contamination estimate from raw intensities (IDAT/GTC, BAF), STAR alignment of the RNA-seq data, genotypes called from RNA-seq and line-to-donor identity check. | 2 (can run in parallel with 3) |
| `03_imputation_TOPMed/` | Preparation of the upload files for the TOPMed Imputation Server (GRCh38 alignment, strand resolution, frequency filters), review of the server's QC-only report and the cleaned upload VCFs. | 3 |
| `04_post_imputation_QC/` | Download and verification of the imputation results, post-imputation QC (R2, ER2, frequencies, sample checks, PRS-CS panel coverage) and the clean PLINK2 dosage dataset. | 4 |
| `05_PRS/` | PRS-CS (UK Biobank and 1000 Genomes EUR LD panels) on the array data and on the imputed data, scoring with PLINK2, a 29-SNP genetic risk score and the descriptive analysis. | 5 |
| `06_RNAseq_counts_and_DE_validation/` | Gene counts (featureCounts), library orientation, QC and normalisation, independent re-implementation of the gene-wise expression analyses, DESeq2/limma-voom checks. | 6 |
| `07_gene_set_analyses/` | Gene-set score replication, GSEA (fgsea) with an exact sample-label permutation calibration. | 7 |
| `08_preprint_tasks/` | Array PCA, direction check of published genes, consensus-gene check, updated sample-QC and lead-variant tables. | 8 |
| `config/` | `config.example.sh`: every path and site-specific setting as a variable. | - |
| `envs/` | Conda environment exports (versions pinned, no build strings). | - |
| `results/` | Aggregated tables only (no donor-level data). | - |
| `docs/` | `script_index.tsv` (release script, original script, what it produces, run mode), `results_index.tsv` and `release_changes.tsv`. | - |

Scripts of folders `01_` to `08_` are numbered in execution order inside each folder; a Python script and the SLURM script that calls it share the same number. The original step numbering of the project (used in some comments as "Step NN") is given in `docs/script_index.tsv`.

## How to run

All scripts were written for a SLURM cluster and a conda-based software stack. Nothing is hard-coded to a machine: paths and site-specific settings are variables.

1. Copy `config/config.example.sh` to `config/config.sh`, edit it and `source` it in the shell from which you submit jobs (SLURM exports the environment to the jobs).
2. In the SLURM scripts replace the placeholders `YOUR_SLURM_ACCOUNT` and `YOUR_PARTITION` (for example `sed -i 's/YOUR_SLURM_ACCOUNT/myaccount/; s/YOUR_PARTITION/myqueue/' */*.sbatch`) or override them on the command line (`sbatch --account=... --partition=...`). Time and memory requests are those used for this dataset.
3. Create a log folder in the directory from which you submit (`mkdir -p logs`; the scripts write `logs/%x_%j.out`).
4. Run the scripts of each folder in numerical order (`sbatch 01_genotype_QC/03_sample_qc.sbatch`, ...). Array jobs are documented in the header of each script (`--array=...`, `--export=ALL,VERSION=V1`, `SETS=...`). Python and R scripts print their usage in the docstring/header.
5. Scripts that depend on a decision or on external steps say so in their header: the upload to the TOPMed server is manual; the private share links of the server and the password of the encrypted results are never stored here (variables in the config file); lists of sample IDs used for sensitivity analyses are read from environment variables (see `config.example.sh`).

### Expected working directory (`PROJECT_DIR`)

The scripts keep the layout of the original analysis below `PROJECT_DIR`; the folders are created by the scripts that write to them.

```
PROJECT_DIR/
  gwas/GCST90570579_harmonised/GCST90570579.h.tsv.gz      # GWAS Catalog harmonised summary statistics
  PRS/
    genotypes/{plink_bin,plink_qc,qc_work,topmed}/        # array data and QC datasets (not included)
    reference/{ldblk_ukbb_eur,ldblk_1kg_eur,grch38,star_gencode_v50}/
    metadata/                                             # donor key, sample labels (not included)
    imputation/{topmed_imputation_v2,clean}/              # server output and clean dataset (not included)
    prscs/{inputs,inputs_imputed,out,out_imputed,weights,weights_imputed,scores,scores_imputed}/
    rnaseq/{star,genotypes,sites}/  intensity/  software/PRScs/  env/
    validation/results/                                   # all outputs of the analysis steps
  DEG_validation/{package,counts,results}/                # expression validation
  GSEA_validation/{package_scrt,gene_lists,results}/      # gene-set analyses
```

## Software and versions

Versions are those installed for the analyses (exact pins in `envs/*.yml`; Linux x86-64, conda-forge and bioconda channels).

| Tool | Version | Environment file |
|---|---|---|
| PLINK 2 | 2.0.0-a.7.13 | `environment_prs.yml`, `environment_topmed.yml` |
| bcftools / htslib | 1.22 / 1.23 (contamination); 1.25 (RNA-seq, TOPMed preparation) | `environment_idat.yml`, `environment_rnaseq.yml`, `environment_topmed.yml` |
| gtc2vcf `idat2gtc` plugin (bcftools) | from https://github.com/freeseek/gtc2vcf (see `02_contamination_identity/01_idat_env.sbatch`) | `environment_idat.yml` |
| PRS-CS | commit `fbc9ec017652f0b40a6ab47f63747d97fa7d033b` (https://github.com/getian107/PRScs); run in an environment with NumPy 1.26.4, SciPy 1.13.1, h5py 3.16.0 | `environment_prscs.yml` |
| STAR | 2.7.10b | `environment_rnaseq.yml` |
| samtools | 1.25 | `environment_rnaseq.yml`, `environment_degval.yml` |
| featureCounts (Subread) | 2.1.1 | `environment_degval.yml` |
| RSeQC | 5.0.6 | `environment_degval.yml` |
| R | 4.3.3 | `environment_degval.yml`, `environment_gsea.yml` |
| DESeq2 / limma / edgeR | 1.42.0 / 3.58.1 / 4.0.16 | `environment_degval.yml` |
| fgsea / msigdbr | 1.28.0 / 25.1.1 (MSigDB 2025.1 Hs, Hallmark) | `environment_gsea.yml` |
| Python | 3.11 (NumPy 2.4.6, SciPy 1.17.1, pandas 3.1.0, Matplotlib 3.10-3.11, openpyxl 3.1.5) | `environment_degval.yml`, `environment_idat.yml`, others |
| poppler (`pdftotext`) | 26.10.0 | `environment_pdftools.yml` |

Create an environment with `conda env create -f envs/environment_<name>.yml`, or let the `*setup*`/`*_env.sbatch` scripts create it (they use strict channel priority). The exports were made on the cluster used for the analysis; if a pinned build is not available for your platform, relax the pin of that package.

## External data needed

| Data | Source | Used in |
|---|---|---|
| PCOS GWAS summary statistics, European, age-adjusted, harmonised to GRCh38 (GCST90570579) | GWAS Catalog, https://www.ebi.ac.uk/gwas/studies/GCST90570579 (Moolhuijsen et al., Nat Genet 2026) | `05_PRS/` |
| Supplementary tables of the GWAS paper (lead SNPs, consensus genes) | Supplementary information of the paper | `05_PRS/`, `08_preprint_tasks/` |
| PRS-CS LD reference panels (UK Biobank EUR, 1000 Genomes EUR) | links in the PRS-CS README (https://github.com/getian107/PRScs); download scripts in `05_PRS/` | `05_PRS/` |
| TOPMed Imputation Server (panel TOPMed r3, hg38) | https://imputation.biodatacatalyst.nhlbi.nih.gov (account required; results via private links) | `03_imputation_TOPMed/`, `04_post_imputation_QC/` |
| 1000 Genomes (biallelic SNV/INDEL release 20190312, GRCh38; sample panel) | https://ftp.1000genomes.ebi.ac.uk | `03_...`, `04_...`, `05_PRS/` |
| GRCh38 reference FASTA (hg38) | UCSC, https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips | `03_imputation_TOPMed/` |
| GENCODE v50 (primary-assembly FASTA and comprehensive GTF) | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/ | `02_...`, `06_...` |
| MSigDB Hallmark gene sets | R package `msigdbr` (installed by `07_gene_set_analyses/01_setup_gsea_env.sbatch`) | `07_gene_set_analyses/` |
| Illumina GSA manifest (.bpm) and cluster file (.egt), IDAT files, GenomeStudio export | project data (controlled access) | `01_...`, `02_...` |

## Data availability

The individual-level data of this study (genotypes, imputed dosages, RNA-seq reads and counts per line, clinical data and sample-level QC results) are not included in this repository. They are available under controlled access or from the authors on reasonable request, subject to ethical approval and data-use agreements.

The `results/` folder contains aggregated tables only (counts of variants at each step, summary statistics, per-gene statistics of group comparisons, gene-set results and direction checks of published genes). It does not contain rows per donor or per sample.

## Notes for reviewers of the code

- In folders `00_` to `08_` individual sample IDs are not part of the scripts or their configuration files; the few analyses that annotate or exclude specific samples read their IDs from environment variables (empty by default, see `config/config.example.sh`). Analyses that depend on them need those lists from the data holders.
- The scripts that verify results against the original analysis of the manuscript (`06_RNAseq_counts_and_DE_validation/08_*`, `10_verdicts.py`) refer to its supplementary tables and, for `08_*`, run the scripts of `original_analysis/` (copy `original_analysis/` or point the script to it).
- Changes made when preparing folders `00_` to `08_` for this release (paths to variables, placeholders, sample IDs to environment variables, removal of private links) are listed in `docs/release_changes.tsv`. The scripts in `original_analysis/` were not changed, apart from the single change described in its README, and contain array sample IDs used in the original analysis.

## AI assistance

The code of folders `00_` to `08_` and the analyses they contain were written with the assistance of an AI model (Claude, Anthropic) and were reviewed by the analyst. The scripts in original_analysis/ were written by Masuma Khatun with assistance from Claude (Anthropic) and are included unmodified except for the single change described in its README. The author reviewed and takes responsibility for that code.

## License and citation

Code: MIT License (`LICENSE`). Tables in `results/`: CC BY 4.0 (`LICENSE-results.txt`). 
If you use this code, please cite the archived version of this repository: Khatun M., et al. (2026). Analysis code for: Donor-associated variation outweighs PMOS status in patient-derived iPSC transcriptomes: a genotyped, BMI-adjusted study of 14 lines. Zenodo. https://doi.org/10.5281/zenodo.24025021 (concept DOI, all versions; version 1.0.1). See CITATION.cff. The citation of the associated manuscript will be added upon publication.