#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/22_prs_analysis.sbatch   (original step script: 52_prs_analysis.sbatch)
# Produces: Case-control statistics, stability, comparison with an external PRS, figures.
# Run: sbatch 05_PRS/22_prs_analysis.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_analysis
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# Step 52 - analysis and figures of the imputed PRS (see 05_PRS/22_prs_analysis.py) + one text file with the construction checks (bim steps, weight reports).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; R="$PRS/validation/results/prs_imputed"; mkdir -p "$R"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/05_PRS/22_prs_analysis.py" --scoredir "$PRS/prscs/scores_imputed" --weightdir "$PRS/prscs/weights_imputed" --clean-dir "$PRS/imputation/clean" \
  --labels "$PRS/metadata/dna_extraction_labels.tsv" --donor-key "$PRS/metadata/donor_key.tsv" --extprs "$PRS/metadata/extprs/extprs_prs_all.txt" \
  --bim-counts "$PRS/prscs/inputs_imputed/bim_counts.tsv" --grs-orientation "$PRS/prscs/weights_imputed/grs29_orientation.tsv" --outdir "$R"
{
  echo "### bim construction (step 49)"; cat "$PRS/prscs/inputs_imputed/bim_report.txt"
  echo; echo "### weights returned to chr:pos:ref:alt IDs (step 51)"; for s in V1 V2 V3 V1_r2ge0.8; do cat "$PRS/prscs/weights_imputed/${s}_weights.txt.report.txt"; done
  echo; echo "### GRS29 (step 51b)"; cat "$PRS/prscs/weights_imputed/grs29_report.txt" | head -n 4
} > "$R/checks_report.txt"
echo "done"; ls "$R" "$R/figures"
