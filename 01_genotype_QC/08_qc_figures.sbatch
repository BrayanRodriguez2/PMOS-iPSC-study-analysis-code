#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/08_qc_figures.sbatch   (original step script: 10b_qc_figures.sbatch)
# Produces: Standard QC figures drawn from the QC tables.
# Run: sbatch 01_genotype_QC/08_qc_figures.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_qcfig
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 10b - standard genotype-QC figures (sample call rate, heterozygosity vs call rate, chrX F,
# KING kinship, SNP call rate, SNP MAF), drawn from the tables that steps 09/10 already wrote.
# NOTHING is recomputed (no plink2 run). See 01_genotype_QC/08_qc_figures.py for the exact tables used and
# for the limits of those tables (e.g. variants removed by the call-rate filter are not stored).
# Inputs (read-only): results/qc/{samples,variants}/* (incl. variants/raw_129.vmiss from step 10c, which
#   lets figure 05 show the variants removed by the call-rate filter), results/bim_check/freq.afreq (raw MAF, step 05),
#   prscs/scores/qc_cohort.afreq (MAF of the final dataset, step 16).
# Output: results/qc/figures/*.png|pdf + qc_figures_report.txt
# Resources: ~777 k MAF values to parse in Python -> 1 CPU / 8 GB, ~1 min, `testing` partition.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
# figure 05 was renamed (05_snp_call_rate_hist -> 05_snp_call_rate); drop the stale files of that first draft
find "$PRS/validation/results/qc/figures" -name '05_snp_call_rate_hist.*' -delete 2>/dev/null || true
python -I "${CODE_DIR}/01_genotype_QC/08_qc_figures.py" \
  --qc-dir "$PRS/validation/results/qc" \
  --afreq-raw "$PRS/validation/results/bim_check/freq.afreq" \
  --afreq-final "$PRS/prscs/scores/qc_cohort.afreq" \
  --outdir "$PRS/validation/results/qc/figures"
