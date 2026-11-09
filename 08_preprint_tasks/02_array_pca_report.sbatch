#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/02_array_pca_report.sbatch   (original step script: 80b_array_pca_report.sbatch)
# Produces: Variance explained, outlier check, comparison with the external table, figures.
# Run: sbatch 08_preprint_tasks/02_array_pca_report.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=pca_report
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# Step 80b - report and figures of the array PCA (see 08_preprint_tasks/02_array_pca_report.py).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; P="$H/PRS"; V="$P/validation/results"; G="$H/GSEA_validation"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
python -I "${CODE_DIR}/08_preprint_tasks/02_array_pca_report.py" "$V/array_pca" "$P/metadata/donor_key.tsv" "$P/rnaseq/rna_samples.tsv" "$V/qc/samples/sample_qc_table.tsv" "$V/qc/samples/king.kin0" \
  "$G/package_scrt/sheet03_S2_genotype_QC.tsv" "$G/package_scrt/sheet02_S1_donors.tsv" "$G/results/tareas/array_pca"
