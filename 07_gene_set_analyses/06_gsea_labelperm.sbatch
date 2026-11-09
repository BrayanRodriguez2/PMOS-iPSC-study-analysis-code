#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/06_gsea_labelperm.sbatch   (original step script: 75_gsea_labelperm.sbatch)
# Produces: Exact sample-label permutation calibration of the GSEA.
# Run: sbatch 07_gene_set_analyses/06_gsea_labelperm.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=gsea_labelperm
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
# Step 75 - calibration of the GSEA by exact sample-label permutation (see 07_gene_set_analyses/06_gsea_labelperm.py). Needs the output of 74 (gene_sets/hallmark.tsv, tables/fgsea_results_all.tsv).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; G="$H/GSEA_validation/results/gsea"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
python -I "${CODE_DIR}/07_gene_set_analyses/06_gsea_labelperm.py" "$G" "$H/DEG_validation/results/step1/norm_counts_symbol.tsv.gz" "$H/DEG_validation/results/clinical/donor_meta.tsv"
