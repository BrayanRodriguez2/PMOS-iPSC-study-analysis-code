#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/10_chip_layout.sbatch   (original step script: 29_chip_layout.sbatch)
# Produces: Tests and layout figure of the chip with the contamination signal.
# Run: sbatch 02_contamination_identity/10_chip_layout.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=chip_layout
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 29 - chip 206149980019 vs the other chips: Mann-Whitney tests and chip-layout figure (see 02_contamination_identity/10_chip_layout.py).
# Input : validation/results/contamination/sample_intensity_metrics.tsv (step 24, final version)
# Output: validation/results/contamination/chip_layout/
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/02_contamination_identity/10_chip_layout.py" \
  --metrics "$PRS/validation/results/contamination/sample_intensity_metrics.tsv" \
  --outdir "$PRS/validation/results/contamination/chip_layout" --chip 206149980019 --normal-chip 206482920107
