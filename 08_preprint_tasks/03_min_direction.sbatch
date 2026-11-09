#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/03_min_direction.sbatch   (original step script: 81_min_direction.sbatch)
# Produces: Direction check of 20 genes named in a published PCOS iPSC study.
# Run: sbatch 08_preprint_tasks/03_min_direction.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=min_dir
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# Step 81 - direction check of the 20 genes of Min et al. 2018 Tables 2-3 (see 08_preprint_tasks/03_min_direction.py). Uses the section-28 results and the original analysis' Table S4 of the final package.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; R="$H/DEG_validation/results"; G="$H/GSEA_validation"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
python -I "${CODE_DIR}/08_preprint_tasks/03_min_direction.py" "$R" "$R/gene_pos.tsv" "$G/package_scrt/sheet05_S4_DE.tsv" "$G/results/min_check"
