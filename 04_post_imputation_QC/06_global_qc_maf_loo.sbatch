#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/06_global_qc_maf_loo.sbatch   (original step script: 41_qc_global_maf_loo.sbatch)
# Produces: Global R2, R2 by MAF interval and leave-one-out ER2 of genotyped SNPs.
# Run: sbatch 04_post_imputation_QC/06_global_qc_maf_loo.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=qc_global
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%j.out
# Step 41 - Post-imputation QC (a) global summary, (b) R2 by MAF interval, (c) leave-one-out ER2 of the genotyped variants (see 04_post_imputation_QC/06_global_qc_maf_loo.py).
# Input: info_tables/chr*.tsv.gz (step 40). Output: validation/results/imputation_qc/. Memory: ~12M variants x 8 columns.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; D="$PRS/imputation/topmed_imputation_v2"; Q="$PRS/validation/results/topmed_qc_only"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/04_post_imputation_QC/06_global_qc_maf_loo.py" --tables-dir "$D/info_tables" --keep75 "$Q/keep_fin_variants_88.tsv" \
  --compare446 "$Q/compare_1000g_446.tsv" --outdir "$PRS/validation/results/imputation_qc"
