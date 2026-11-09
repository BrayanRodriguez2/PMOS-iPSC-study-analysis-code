#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/03_count_overlap.sbatch   (original step script: 04_count_overlap.sbatch)
# Produces: Overlap of array SNPs with the PRS-CS panel (by rsID, with allele checks).
# Run: sbatch 05_PRS/03_count_overlap.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_overlap
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 04 - Count rsIDs in the sumstats (a), LD panel (b) and array .bim (c),
#           their intersections, and rsID recovery by chr:pos(GRCh38)+alleles.
# See 05_PRS/03_count_overlap.py for the full method. NOTHING is modified: the .bim and
# all inputs are only read; results go to PRS/validation/results/overlap/.
# Resources: single-threaded Python holding ~10 M rsIDs in memory -> 24 GB, 1 CPU.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail

PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prs

SNPINFO=$(find "$PRS/reference/ldblk_1kg_eur" -maxdepth 1 -name 'snpinfo*' | head -n 1)
echo "### snpinfo: $SNPINFO"

python -I "${CODE_DIR}/05_PRS/03_count_overlap.py" \
  --sumstats "${PROJECT_DIR}/gwas/GCST90570579_harmonised/GCST90570579.h.tsv.gz" \
  --snpinfo  "$SNPINFO" \
  --bim      "$PRS/genotypes/plink_bin/b38.bim" \
  --outdir   "$PRS/validation/results/overlap"
