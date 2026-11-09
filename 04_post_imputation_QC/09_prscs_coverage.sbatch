#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/09_prscs_coverage.sbatch   (original step script: 44_prscs_coverage.sbatch)
# Produces: Coverage of the PRS-CS LD panel and of the 29 lead SNPs in the imputed data.
# Run: sbatch 04_post_imputation_QC/09_prscs_coverage.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prscs_cov
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
# Step 44 - Post-imputation QC (f): coverage of the PRS-CS panel (UKBB EUR HapMap3) and of the 29 lead SNPs in the imputed data (see 04_post_imputation_QC/09_prscs_coverage.py).
# Memory: 11.5 M imputed variants + 9.4 M sumstats rows.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; D="$PRS/imputation/topmed_imputation_v2"; L="$PRS/validation/results/loci"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/04_post_imputation_QC/09_prscs_coverage.py" --tables-dir "$D/info_tables" --snpinfo "$PRS/reference/ldblk_ukbb_eur/snpinfo_ukbb_hm3" \
  --sumstats-gz "${PROJECT_DIR}/gwas/GCST90570579_harmonised/GCST90570579.h.tsv.gz" --sumstats-prscs "$PRS/prscs/inputs/sumstats_prscs.txt" \
  --leads "$L/lead_snps_lookup.tsv" --locus-cov "$L/locus_coverage.tsv" \
  --outdir "$PRS/validation/results/imputation_qc"
