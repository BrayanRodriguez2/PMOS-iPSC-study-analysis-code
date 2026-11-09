#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/08_external_pipeline_check.sbatch   (original step script: 66b_their_pipeline.sbatch)
# Produces: Runs the original analysis scripts on the original counts (requires original_analysis/).
# Run: sbatch 06_RNAseq_counts_and_DE_validation/08_external_pipeline_check.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=their_pipe
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 66b - the original analysis' pipeline EXACTLY as in Additional_file_3 (common.py, an1.py, an2.py, an_restore.py sections 1-4: PCA, sample correlation, variance partition,
# BMI-adjusted models) run on the original Normalized_counts.tabular, to confirm that the original code reproduces the original own numbers (separates "data" from "code" differences).
# The scripts are verbatim copies (DEG_validation/their_run/scripts/*.py); only the hard-coded absolute paths of the original session were replaced and an_restore.py was
# cut after section 4 (sections 5-8 need the September supplement workbook and the HipSci protein tables, not used here). Nothing of the original package is modified.
# Environment degval (python 3.11, numpy 2, pandas 3, scipy 1.17: newer than the original Python 3.13 stack of the time; any difference would show up in the logs).
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; W="$H/DEG_validation/their_run"; PK="$H/DEG_validation/package"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
python -I "${CODE_DIR}/06_RNAseq_counts_and_DE_validation/08_external_pipeline_check.py" "$PK" "$W"
cp "$PK/Normalized_counts.tabular" "$W/Normalized_counts.tabular"
mkdir -p "$W/pkg/5_Pipeline_CSC/scripts"; cp "$W/scripts/common.py" "$W/pkg/5_Pipeline_CSC/scripts/common.py"
cd "$W"
echo "######## an1.py"; python scripts/an1.py 2>&1 | tee "$W/log_an1.txt"
echo "######## an2.py"; python scripts/an2.py 2>&1 | tee "$W/log_an2.txt"
echo "######## an_restore.py (sections 1-4)"
python scripts/an_restore_s1_4.py "$W/pkg" "$W/Normalized_counts.tabular" "$W/pkg" "$W/restore_out" 2>&1 | tee "$W/log_an_restore.txt"
echo done
