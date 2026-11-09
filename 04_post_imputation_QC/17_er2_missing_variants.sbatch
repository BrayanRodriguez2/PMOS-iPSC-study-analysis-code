#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/17_er2_missing_variants.sbatch   (original step script: 71_er2_missing7.sbatch)
# Produces: Why a few kept variants have no ER2 (absent from the server output).
# Run: sbatch 04_post_imputation_QC/17_er2_missing_variants.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=er2_missing7
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# Step 71 - Why do 7 variants of the 75 / 358 subsets have no ER2? (follow-up of the post-imputation QC, see 04_post_imputation_QC/17_er2_missing_variants.py).
# READ-ONLY: only reads info tables, variant tables and upload VCFs; writes 3 small files to results/imputation_qc/er2_missing7/.
# Resources: info tables are read per chromosome (<= ~1.2 M rows) and the variant tables in 1 M-row chunks, so 16 GB is ample; 1 CPU (pandas is single-threaded).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; D="$PRS/imputation/topmed_imputation_v2"; Q="$PRS/validation/results/topmed_qc_only"; C="$PRS/imputation/clean"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/04_post_imputation_QC/17_er2_missing_variants.py" --keep75 "$Q/keep_fin_variants_88.tsv" --compare446 "$Q/compare_1000g_446.tsv" \
  --info-dir "$D/info_tables" --clean-table "$C/variant_table.tsv.gz" --pre-table "$C/variant_table_pre_er2filter.tsv.gz" \
  --removed-ids "$C/removed_er2_lt_0.5_ids.txt" --vcf-dir "$PRS/topmed/upload_v2" \
  --outdir "$PRS/validation/results/imputation_qc/er2_missing7"
