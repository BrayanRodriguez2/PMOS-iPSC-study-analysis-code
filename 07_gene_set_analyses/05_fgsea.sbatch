#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/05_fgsea.sbatch   (original step script: 74_fgsea.sbatch)
# Produces: fgsea preranked GSEA (custom sets and MSigDB Hallmark).
# Run: sbatch 07_gene_set_analyses/05_fgsea.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=gsea_fgsea
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# Step 74 - fgsea preranked on every ranking of step 73 (see 07_gene_set_analyses/05_fgsea.R). Environment "gsea" (step 72a).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; G="$H/GSEA_validation"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate gsea
Rscript "${CODE_DIR}/07_gene_set_analyses/05_fgsea.R" "$G/results/gsea"
