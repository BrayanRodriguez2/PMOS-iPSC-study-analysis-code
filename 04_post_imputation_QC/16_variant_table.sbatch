#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/16_variant_table.sbatch   (original step script: 47c_variant_table.sbatch)
# Produces: Variant table of the clean dataset after the ER2 filter.
# Run: sbatch 04_post_imputation_QC/16_variant_table.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=er2_table
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%j.out
# Step 47c - variant table of the clean file AFTER removing the ER2 < 0.5 genotyped variants (python part of 47b: the second `conda activate` inside 47b
# failed, so it is done here in a fresh job). The pre-filter table stays as variant_table_pre_er2filter.tsv.gz.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
O="${PROJECT_DIR}/PRS/imputation/clean"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
cd "$O"
python -I - <<'PY'
import pandas as pd
t = pd.read_csv("variant_table_pre_er2filter.tsv.gz", sep="\t", dtype={"CHROM": str})
ids = set(l.strip() for l in open("removed_er2_lt_0.5_ids.txt"))
n = t[~t.ID.isin(ids)]
n.to_csv("variant_table.tsv.gz", sep="\t", index=False, compression="gzip", float_format="%.6g")
print(f"variant_table: {len(t):,} -> {len(n):,}; Typed {int((n.TYPE=='Typed').sum()):,}, Imputed {int((n.TYPE=='Imputed').sum()):,}; "
      f"R2>=0.8 {int((n.R2>=0.8).sum()):,} ({100*(n.R2>=0.8).mean():.2f} %); typed with ER2<0.5 left: {int(((n.TYPE=='Typed')&(n.ER2<0.5)).sum())}; "
      f"typed with ER2<0.8: {int(((n.TYPE=='Typed')&(n.ER2<0.8)).sum())}")
PY
cat counts_after_er2filter.tsv | tail -n 3
