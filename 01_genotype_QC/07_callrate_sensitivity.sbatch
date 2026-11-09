#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/07_callrate_sensitivity.sbatch   (original step script: 09c_callrate_sensitivity.sbatch)
# Produces: Sensitivity of the sample set to the call-rate threshold.
# Run: sbatch 01_genotype_QC/07_callrate_sensitivity.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_crsens
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 09c - sensitivity of the sample exclusions to the call-rate threshold (0.95 ... 0.98),
# keeping the heterozygosity and relatedness filters of step 09. Pure table arithmetic on
# sample_qc_table.tsv (see 01_genotype_QC/07_callrate_sensitivity.py for why no plink2 re-run is needed).
# Output: results/qc/sensitivity/callrate_sensitivity.tsv + report. Resources: trivial.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
python -I "${CODE_DIR}/01_genotype_QC/07_callrate_sensitivity.py" \
  --qc-table "$PRS/validation/results/qc/samples/sample_qc_table.tsv" \
  --outdir "$PRS/validation/results/qc/sensitivity"
