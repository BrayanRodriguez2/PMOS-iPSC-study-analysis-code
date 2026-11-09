#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/14_er2_filter_list.sbatch   (original step script: 47a_er2_list.sbatch)
# Produces: List of genotyped variants with ER2 < 0.5 and their PRS-CS panel membership.
# Run: sbatch 04_post_imputation_QC/14_er2_filter_list.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=er2_list
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
# Step 47a - list of the genotyped variants with ER2 < 0.5 to remove from the clean file + PRS-CS panel membership (see 04_post_imputation_QC/14_er2_filter_list.py). Nothing is removed here.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/04_post_imputation_QC/14_er2_filter_list.py" --clean-dir "$PRS/imputation/clean" --qc-only-dir "$PRS/validation/results/topmed_qc_only" \
  --snpinfo "$PRS/reference/ldblk_ukbb_eur/snpinfo_ukbb_hm3" --sumstats-gz "${PROJECT_DIR}/gwas/GCST90570579_harmonised/GCST90570579.h.tsv.gz" \
  --sumstats-prscs "$PRS/prscs/inputs/sumstats_prscs.txt"
