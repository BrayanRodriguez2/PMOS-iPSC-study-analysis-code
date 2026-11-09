#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/05_info_tables.sbatch   (original step script: 40_info_tables.sbatch)
# Produces: Per-chromosome tables from the server info files.
# Run: sbatch 04_post_imputation_QC/05_info_tables.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=info_tab
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-22
#SBATCH --output=logs/%x_%A_%a.out
# Step 40 - info.gz (sites VCF with INFO R2/ER2/AF/MAF/flags) -> plain table, one chromosome per array task (see 04_post_imputation_QC/05_info_tables.py).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
C=$SLURM_ARRAY_TASK_ID
PRS="${PROJECT_DIR}/PRS"; D="$PRS/imputation/topmed_imputation_v2"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/04_post_imputation_QC/05_info_tables.py" --info "$D/unzipped/chr$C.info.gz" --out "$D/info_tables/chr$C.tsv.gz"
