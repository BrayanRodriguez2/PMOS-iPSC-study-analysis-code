#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/03_counts_qc_norm.sbatch   (original step script: 64_counts_qc_norm.sbatch)
# Produces: Library orientation, per-line QC, size-factor normalisation, comparison with the external file.
# Run: sbatch 06_RNAseq_counts_and_DE_validation/03_counts_qc_norm.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=counts_qc
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# Step 64 - orientation, QC table, normalisation and comparison with the original Normalized_counts.tabular (see 06_RNAseq_counts_and_DE_validation/03_counts_qc_norm.py). STRAND (optional) and OUT (optional) come from the environment: default = the best strand setting and results/step1; STRAND=0 OUT=results/step1_s0 repeats it with unstranded counting.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; D="$H/DEG_validation"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
python -I "${CODE_DIR}/06_RNAseq_counts_and_DE_validation/03_counts_qc_norm.py" --fc-prefix "$D/counts/fc_s" ${STRAND:+--strand $STRAND} --star-dir "$H/PRS/rnaseq/star" \
  --rna-samples "$H/PRS/rnaseq/rna_samples.tsv" --theirs "$D/package/Normalized_counts.tabular" --their-s1 "$D/package/af1_tsv/sheet02_S1_donors.tsv" --out "${OUT:-$D/results/step1}"
