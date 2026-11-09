#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/04_baf_contamination.sbatch   (original step script: 24_baf_contamination.sbatch)
# Produces: BAF/LRR summaries and the DNA-contamination estimate per sample.
# Run: sbatch 02_contamination_identity/04_baf_contamination.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=baf_contam
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 24 - Part A of the contamination diagnostic: BAF/LRR summaries, contamination estimate and BAF figures.
# Input : PRS/intensity/gtc/<chip>/*.gtc (step 23), PRS/intensity/bpm_snp_map.csv (step 23b), the sample sheet,
#         results/qc/samples/sample_qc_table.tsv (step 09)
# Output: validation/results/contamination/{sample_intensity_metrics.tsv, focal_baf_lrr.npz, figures/A1..A3*}
# Memory: genotype calls (uint8) + BAF (float32) of 144 samples x ~730k autosomal loci ~ 0.6 GB; 16 GB is plenty.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
G=${RAW_GENOTYPE_DIR}
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/02_contamination_identity/04_baf_contamination.py" \
  --gtc-dir "$PRS/intensity/gtc" \
  --snp-map "$PRS/intensity/bpm_snp_map.csv" \
  --sample-sheet "$G/AKurg_06.10.22_b38.csv" \
  --qc-table "$PRS/validation/results/qc/samples/sample_qc_table.tsv" \
  --outdir "$PRS/validation/results/contamination" \
  --focal "${FOCAL_SAMPLE_IDS:-}"
