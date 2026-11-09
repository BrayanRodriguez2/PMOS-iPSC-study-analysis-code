#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/06_qc_by_chip.sbatch   (original step script: 09b_qc_by_chip.sbatch)
# Produces: Call rate and heterozygosity summarised by array chip.
# Run: sbatch 01_genotype_QC/06_qc_by_chip.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_qcchip
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 09b - call rate and heterozygosity (F) by genotyping chip (SentrixBarcode_A of the original
# sample sheet) + figure. Joins the step-09 table to the sheet (Sample_ID == IID); nothing is recomputed.
# Input (read-only): ${RAW_GENOTYPE_DIR}/AKurg_06.10.22_b38.csv
#                    results/qc/samples/sample_qc_table.tsv
# Output: results/qc/by_chip/{chip_sample_map.tsv, chip_summary.tsv, chip_tests.txt, qc_by_chip.png/.pdf}
# Resources: 140 rows; 1 CPU / 4 GB on `testing`, seconds.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
python -I "${CODE_DIR}/01_genotype_QC/06_qc_by_chip.py" \
  --sample-sheet ${RAW_GENOTYPE_DIR}/AKurg_06.10.22_b38.csv \
  --qc-table "$PRS/validation/results/qc/samples/sample_qc_table.tsv" \
  --outdir "$PRS/validation/results/qc/by_chip"
