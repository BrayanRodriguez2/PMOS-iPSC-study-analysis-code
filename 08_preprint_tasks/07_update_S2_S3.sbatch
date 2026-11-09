#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/07_update_S2_S3.sbatch   (original step script: 83_update_S2_S3.sbatch)
# Produces: Updated supplementary tables S2 (sample QC) and S3 (lead variants).
# Run: sbatch 08_preprint_tasks/07_update_S2_S3.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=upd_s2s3
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# Step 83 - updated Tables S2 and S3 with our results (see 08_preprint_tasks/07_update_S2_S3.py; needs steps 80, 83a (openpyxl), 83b).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; G="$H/GSEA_validation"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
mkdir -p "$G/results/tareas/S2_S3"
python -I "${CODE_DIR}/08_preprint_tasks/07_update_S2_S3.py" "$H/PRS" "$G/package_scrt/sheet03_S2_genotype_QC.tsv" "$G/package_scrt/sheet04_S3_lead_variants.tsv" "$G/results/tareas/S2_S3"
