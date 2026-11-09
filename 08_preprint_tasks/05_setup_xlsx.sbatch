#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/05_setup_xlsx.sbatch   (original step script: 83a_setup_xlsx.sbatch)
# Produces: Adds openpyxl to the 'degval' environment.
# Run: sbatch 08_preprint_tasks/05_setup_xlsx.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=xlsx_env
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# Step 83a - add openpyxl (xlsx writer) to the "degval" environment (conda-forge) to write the updated Tables S2 / S3 as .xlsx.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda install -y -n degval --override-channels --strict-channel-priority -c conda-forge openpyxl
conda activate degval
python -c "import openpyxl; print('openpyxl', openpyxl.__version__)"
