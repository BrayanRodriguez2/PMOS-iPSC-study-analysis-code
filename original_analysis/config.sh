#!/usr/bin/env bash
# =====================================================================
#  PMOS iPSC genetics pipeline - configuration
#  Edit the paths below once, then every script reads this file.
#  Designed for CSC Puhti (Linux, Slurm). Keep genotype data on CSC.
# =====================================================================
 
# ---- Project folders ------------------------------------------------
PROJECT=/scratch/project_XXXXXXX/pcos_ipsc          # <- your CSC project scratch
OUT=$PROJECT/results                                  # all outputs go here
 
# ---- Inputs you already have -----------------------------------------
RAW_PLINK=$PROJECT/raw/b38                 # b38.ped + b38.map (GenomeStudio export, GRCh38)
COUNTS=$PROJECT/rna/Normalized_counts.tabular   # DESeq2-normalised counts, genes x 14 lines
SAMPLE_KEY=$PROJECT/meta/sample_key.csv    # template: templates/sample_key_template.csv
 
# ---- Inputs to obtain ------------------------------------------------
# PCOS GWAS summary statistics (Moolhuijsen et al., Nat Genet 2026): GWAS Catalog GCST90570579
# (European, age-adjusted, without 23andMe; public). Full three-strata data are not public.
SUMSTATS=$PROJECT/gwas/pcos_sumstats.txt.gz
GWAS_N=456570          # VERIFY: set to the N of GCST90570579 (European, no 23andMe)
GWAS_CASES=15634       # VERIFY: number of cases in GCST90570579 (for coloc)
# 29 lead variants from Supplementary Table (template: templates/lead_variants_template.tsv)
LEADS=$PROJECT/gwas/lead_variants.tsv
# i2QTL iPSC cis-eQTL nominal results (Bonder et al., Nat Genet 2021), tab-separated
I2QTL=$PROJECT/eqtl/i2qtl_nominal.tsv.gz
I2QTL_N=1367           # adjust to the sample size stated in the file you download
 
# ---- Reference files -------------------------------------------------
FASTA=$PROJECT/ref/GRCh38.fa               # GRCh38 reference (chr-prefixed), indexed (.fai)
DBSNP=$PROJECT/ref/dbsnp_b38.vcf.gz        # dbSNP GRCh38 VCF (chr-prefixed), tabix-indexed
PRSCS_DIR=$PROJECT/tools/PRScs             # git clone https://github.com/getian107/PRScs
PRSCS_LDREF=$PROJECT/ref/ldblk_ukbb_eur    # PRS-CS European LD reference (see PRScs README)
 
# ---- Imputation output (downloaded from TOPMed/Michigan server) ------
IMPUTED_DIR=$PROJECT/imputed               # chr*.dose.vcf.gz after unzipping
 
# ---- Parameters -----------------------------------------------------
THREADS=8
PRSCS_PHI=1e-4         # value used in the paper; set to "auto" to let PRS-CS learn it
IMPUTE_R2=0.8          # post-imputation quality filter
MAF=0.01
 
# ---- Tools (CSC modules or your own installs) ------------------------
PLINK=plink            # PLINK 1.9
PLINK2=plink2          # PLINK 2.0
BCFTOOLS=bcftools
PYTHON=python3
RSCRIPT=Rscript
 
mkdir -p "$OUT"
