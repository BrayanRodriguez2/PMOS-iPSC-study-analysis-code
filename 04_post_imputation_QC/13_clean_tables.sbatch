#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/13_clean_tables.sbatch   (original step script: 45c_clean_tables.sbatch)
# Produces: Variant table and counts of the clean dataset.
# Run: sbatch 04_post_imputation_QC/13_clean_tables.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=imp_tables
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
# Step 45c - variant table and filter-count table of the clean imputed fileset (python part of 45b; the "prs" env of 45b has plink2 but no pandas,
# so this runs in the "idat" env). See 04_post_imputation_QC/13_clean_tables.py. Input: imputation/clean/{v2_imputed_clean.pvar, stats129.*, counts_chr*.tsv}, info_tables.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/04_post_imputation_QC/13_clean_tables.py" --clean-dir "$PRS/imputation/clean" --tables-dir "$PRS/imputation/topmed_imputation_v2/info_tables" --unzipped-dir "$PRS/imputation/topmed_imputation_v2/unzipped" \
  --qc-only-dir "$PRS/validation/results/topmed_qc_only"
