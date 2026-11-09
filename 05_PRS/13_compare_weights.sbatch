#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/13_compare_weights.sbatch   (original step script: 17_compare_weights.sbatch)
# Produces: Comparison of the weights between PRS-CS versions.
# Run: sbatch 05_PRS/13_compare_weights.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_wcompare
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# Step 17 - compare V1/V2/V3 weights (see 05_PRS/13_compare_weights.py). Read-only, tiny.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
python -I "${CODE_DIR}/05_PRS/13_compare_weights.py" --weightdir "$PRS/prscs/weights" \
       --out "$PRS/validation/results/prs/weight_comparison.txt"
