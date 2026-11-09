#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/09_exclusion_rule.sbatch   (original step script: 35_exclusion_rule.sbatch)
# Produces: Rule that removes problematic genotyped SNPs before the final upload.
# Run: sbatch 03_imputation_TOPMed/09_exclusion_rule.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=excl_rule
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# Step 35 - Exclusion list for the v2 VCFs + lead/proxy safety check (see 03_imputation_TOPMed/09_exclusion_rule.py). Writes lists only; removes nothing.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/03_imputation_TOPMed/09_exclusion_rule.py" --outdir "$PRS/validation/results/topmed_qc_only" \
  --leads "$PRS/validation/results/loci/lead_snps_lookup.tsv" --locus-cov "$PRS/validation/results/loci/locus_coverage.tsv"
