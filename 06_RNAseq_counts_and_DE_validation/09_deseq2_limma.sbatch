#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/09_deseq2_limma.sbatch   (original step script: 67_methods.sbatch)
# Produces: DESeq2 and limma-voom on raw counts and comparison with the OLS.
# Run: sbatch 06_RNAseq_counts_and_DE_validation/09_deseq2_limma.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=deg_methods
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# Step 67 - DESeq2 and limma-voom on raw counts (06_RNAseq_counts_and_DE_validation/09_deseq2_limma.R) and comparison with the OLS and the original S4 (06_RNAseq_counts_and_DE_validation/09_deseq2_limma.py).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; D="$H/DEG_validation"; S="${CODE_DIR}"; R="$D/results"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
mkdir -p "$R/methods"
Rscript "$S/06_RNAseq_counts_and_DE_validation/09_deseq2_limma.R" "$R/step1/counts_raw.tsv.gz" "$R/clinical/donor_meta.tsv" "$R/analysis/ours_s2/de_ols.tsv.gz" "$R/step1/size_factors.tsv" "$R/methods"
python -I "$S/06_RNAseq_counts_and_DE_validation/09_deseq2_limma.py" "$R/methods" "$R/analysis/ours_s2/de_ols.tsv.gz" "$D/package/af1_tsv/sheet05_S4_DE.tsv" "$R/methods"
echo done
