#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/04_make_rankings.sbatch   (original step script: 73_rankings.sbatch)
# Produces: Ranking statistics and custom gene sets for the GSEA.
# Run: sbatch 07_gene_set_analyses/04_make_rankings.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=gsea_rank
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# Step 73 - rankings and custom gene sets for the GSEA (see 07_gene_set_analyses/04_make_rankings.py). Min et al. lists are read from GSEA_validation/gene_lists/ if present (they were not found in the project).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; G="$H/GSEA_validation"; mkdir -p "$G/gene_lists"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
python -I "${CODE_DIR}/07_gene_set_analyses/04_make_rankings.py" "$H/DEG_validation/results" "$G/package_scrt" "$G/gene_lists" "$G/results/gsea"
