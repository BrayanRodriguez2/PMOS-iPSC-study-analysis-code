#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/19_grs29.sbatch   (original step script: 51b_grs29.sbatch)
# Produces: Weights and allele-orientation table of the 29-SNP genetic risk score.
# Run: sbatch 05_PRS/19_grs29.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=grs29
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%j.out
# Step 51b - 29-SNP GRS: weights (ST2 risk allele and effect) mapped to the clean file + orientation table (see 05_PRS/19_grs29.py).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/05_PRS/19_grs29.py" --st2 "$(ls "$PRS"/metadata/supp_tables/sheet03_ST2*.tsv)" \
  --lookup "$PRS/validation/results/loci/lead_snps_lookup.tsv" --variant-table "$PRS/imputation/clean/variant_table.tsv.gz" --outdir "$PRS/prscs/weights_imputed"
