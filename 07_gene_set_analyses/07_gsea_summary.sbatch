#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/07_gsea_summary.sbatch   (original step script: 76_gsea_summary.sbatch)
# Produces: Summary tables and calibration figure.
# Run: sbatch 07_gene_set_analyses/07_gsea_summary.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=gsea_summary
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# Step 76 - summary tables and calibration figure (see 07_gene_set_analyses/07_gsea_summary.py).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
python -I "${CODE_DIR}/07_gene_set_analyses/07_gsea_summary.py" "$H/GSEA_validation/results/gsea"
