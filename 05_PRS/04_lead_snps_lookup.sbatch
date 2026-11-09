#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/04_lead_snps_lookup.sbatch   (original step script: 11_lead_snps_lookup.sbatch)
# Produces: The 29 lead SNPs in the summary statistics and genome-wide significant loci.
# Run: sbatch 05_PRS/04_lead_snps_lookup.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_leads
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 11 - Part B.1: the 29 lead SNPs of the PCOS GWAS in the deposited sumstats
# (GCST90570579 harmonised, GRCh38) + distance-based cross-check of significant loci.
# See 05_PRS/04_lead_snps_lookup.py. Streaming pass over 9.1 M rows: 1 CPU / 8 GB, ~1-2 min.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
python -I "${CODE_DIR}/05_PRS/04_lead_snps_lookup.py" \
  --sumstats "${PROJECT_DIR}/gwas/GCST90570579_harmonised/GCST90570579.h.tsv.gz" \
  --outdir "$PRS/validation/results/loci"
