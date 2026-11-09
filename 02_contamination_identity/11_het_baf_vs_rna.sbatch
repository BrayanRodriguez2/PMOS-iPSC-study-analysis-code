#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/11_het_baf_vs_rna.sbatch   (original step script: 30_het_baf_vs_rna.sbatch)
# Produces: BAF of array-heterozygous sites confirmed or not by RNA-seq.
# Run: sbatch 02_contamination_identity/11_het_baf_vs_rna.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=het_baf_rna
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 30 - BAF of the array-het sites that the RNA confirms vs does not confirm (joins Part A and Part B; see 02_contamination_identity/11_het_baf_vs_rna.py)
# Input : steps 23, 23b, 24, 27b, 28 outputs and the PLINK b38.bim of the project export (read-only copy in PRS/genotypes/plink_bin)
# Output: validation/results/rna_identity/{het_baf_vs_rna.tsv, figures/B3_*}
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/02_contamination_identity/11_het_baf_vs_rna.py" \
  --rna-dir "$PRS/rnaseq/genotypes" --rna-samples "$PRS/rnaseq/rna_samples.tsv" \
  --identity-dir "$PRS/validation/results/rna_identity" --bim "$PRS/genotypes/plink_bin/b38.bim" \
  --snp-map "$PRS/intensity/bpm_snp_map.csv" --gtc-dir "$PRS/intensity/gtc" \
  --metrics "$PRS/validation/results/contamination/sample_intensity_metrics.tsv"
