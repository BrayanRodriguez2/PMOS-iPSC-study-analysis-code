#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/15_apply_er2_filter.sbatch   (original step script: 47b_apply_er2_filter.sbatch)
# Produces: Removes those variants from the clean dataset (previous files renamed, not deleted).
# Run: sbatch 04_post_imputation_QC/15_apply_er2_filter.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=er2_apply
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 47b - Remove from the clean imputed dataset the genotyped variants with ER2 < 0.5 (list of step 47a; decision of the analyst, 2026-11-04).
# NOTHING IS DELETED: the previous fileset and table are RENAMED to *_pre_er2filter and a new v2_imputed_clean.{pgen,pvar,psam} is written
# (plink2 --exclude). Checks: variants = before - removed, none of the removed IDs left, same 135 samples. Also writes the per-chromosome counts after
# the filter. Use as before: --keep keep_129.txt (main) / keep_135.txt (sensitivity).
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; O="$PRS/imputation/clean"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
T=$SLURM_CPUS_PER_TASK; MEM=$((SLURM_MEM_PER_NODE - 4000))
cd "$O"
[ -s removed_er2_lt_0.5_ids.txt ] || { echo "list missing (run 47a)"; exit 1; }
if [ -e v2_imputed_clean_pre_er2filter.pgen ]; then echo "pre_er2filter fileset already exists: refusing to overwrite"; exit 1; fi
for e in pgen pvar psam log; do mv "v2_imputed_clean.$e" "v2_imputed_clean_pre_er2filter.$e"; done
mv variant_table.tsv.gz variant_table_pre_er2filter.tsv.gz
before=$(grep -vc '^#' v2_imputed_clean_pre_er2filter.pvar); nrm=$(wc -l < removed_er2_lt_0.5_ids.txt)
plink2 --pfile v2_imputed_clean_pre_er2filter --exclude removed_er2_lt_0.5_ids.txt --make-pgen --threads "$T" --memory "$MEM" --out v2_imputed_clean
after=$(grep -vc '^#' v2_imputed_clean.pvar)
echo "variants before: $before | removed (list): $nrm | after: $after | expected: $((before - nrm))"
[ "$after" -eq $((before - nrm)) ] || { echo "COUNT MISMATCH"; exit 1; }
left=$(grep -v '^#' v2_imputed_clean.pvar | cut -f3 | sort | comm -12 - <(sort removed_er2_lt_0.5_ids.txt) | wc -l)
echo "removed IDs still present: $left (must be 0)"; [ "$left" -eq 0 ]
echo "samples: $(grep -vc '^#' v2_imputed_clean.psam) (same as before: $(grep -vc '^#' v2_imputed_clean_pre_er2filter.psam))"
for c in $(seq 1 22); do printf "chr%s\t%s\t%s\n" "$c" "$(grep -v '^#' v2_imputed_clean_pre_er2filter.pvar | awk -F'\t' -v c="$c" '$1==c' | wc -l)" "$(grep -v '^#' v2_imputed_clean.pvar | awk -F'\t' -v c="$c" '$1==c' | wc -l)"; done > counts_after_er2filter.tsv
sed -i '1i chrom\tbefore_er2_filter\tafter_er2_filter' counts_after_er2filter.tsv
# (the variant table of the filtered dataset is written by 04_post_imputation_QC/16_variant_table.sbatch; the original script ended with a second `conda activate` that failed in this environment)
cat counts_after_er2filter.tsv
