#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/03_deseq2_without_one_line.sbatch   (original step script: 73a_deseq2_no350.sbatch)
# Produces: DESeq2/limma re-run without one line (sensitivity).
# Run: sbatch 07_gene_set_analyses/03_deseq2_without_one_line.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=deseq_no350
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# Step 73a - DESeq2 / limma-voom (06_RNAseq_counts_and_DE_validation/09_deseq2_limma.R) re-run WITHOUT HEL350.3 (13 lines; idic(20) mosaic) for the GSEA sensitivity.
# Size factors are those estimated on all 14 lines (median of ratios over ~19 k genes; one line with a chr20 arm change moves them only marginally); the XIST covariate is the 14-line value of the 13 lines.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; D="$H/DEG_validation/results"; S="${CODE_DIR}"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
mkdir -p "$D/methods_no350"
Rscript "$S/06_RNAseq_counts_and_DE_validation/09_deseq2_limma.R" "$D/step1/counts_raw.tsv.gz" "$D/clinical/donor_meta.tsv" "$D/analysis/ours_s2/de_ols.tsv.gz" "$D/step1/size_factors.tsv" "$D/methods_no350" HEL350.3
