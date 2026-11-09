#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/04_s7a_vs_paper.sbatch   (original step script: 82_s7a_vs_paper.sbatch)
# Produces: Consensus genes of the manuscript table against the GWAS supplementary tables.
# Run: sbatch 08_preprint_tasks/04_s7a_vs_paper.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=s7a_check
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# Step 82 - Table S7a (final package) against ST2 / ST9 of Moolhuijsen et al. 2026 (see 08_preprint_tasks/04_s7a_vs_paper.py).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; T="$H/PRS/metadata/supp_tables"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
mkdir -p "$H/GSEA_validation/results/tareas/s7a_check"
python -I "${CODE_DIR}/08_preprint_tasks/04_s7a_vs_paper.py" "$H/GSEA_validation/package_scrt/sheet10_S7a_GWAS_genes.tsv" "$(ls "$T"/sheet03_ST2*.tsv)" "$(ls "$T"/sheet10_ST9*.tsv)" "$H/GSEA_validation/results/tareas/s7a_check"
