#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/12_ibs0_donor_profile.sbatch   (original step script: 31_ibs0_donor_profile.sbatch)
# Produces: IBS0 profile of focal donors against the rest of the cohort.
# Run: sbatch 02_contamination_identity/12_ibs0_donor_profile.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=ibs0_prof
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%j.out
# Step 31 - IBS0 profile of each donor against the lines of other donors (see 02_contamination_identity/12_ibs0_donor_profile.py). Reads step 28 tables only.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/02_contamination_identity/12_ibs0_donor_profile.py" --identity-dir "$PRS/validation/results/rna_identity" \
  --qc-table "$PRS/validation/results/qc/samples/sample_qc_table.tsv"
