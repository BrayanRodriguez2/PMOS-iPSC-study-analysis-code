#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/12_clean_merge.sbatch   (original step script: 45b_clean_merge.sbatch)
# Produces: Merge of the per-chromosome clean files.
# Run: sbatch 04_post_imputation_QC/12_clean_merge.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=imp_merge
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 45b - Merge the 22 per-chromosome clean pgen files into ONE PLINK2 fileset with dosage (v2_imputed_clean.pgen/pvar/psam), compute MAF/HWE
# over the 129 QC samples for the variant table, and write the counts / variant tables (see 04_post_imputation_QC/13_clean_tables.py). Originals of TOPMed untouched.
# Output: PRS/imputation/clean/{v2_imputed_clean.{pgen,pvar,psam}, stats129.{afreq,hardy}, variant_table.tsv.gz, counts_clean.tsv, clean_report.txt}
# Use later:  plink2 --pfile v2_imputed_clean --keep keep_129.txt ... (main)   or   --keep keep_135.txt (sensitivity);  R2 >= 0.8:  --extract-if-info "R2 >= 0.8"
#             (R2 is in the pvar INFO if the merge kept it; otherwise filter with variant_table.tsv.gz).
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; O="$PRS/imputation/clean"; D="$PRS/imputation/topmed_imputation_v2"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
T=$SLURM_CPUS_PER_TASK; MEM=$((SLURM_MEM_PER_NODE - 4000))
for c in $(seq 1 22); do [ -s "$O/chr${c}_clean.pgen" ] && [ -s "$O/counts_chr$c.tsv" ] || { echo "chr$c clean file missing"; exit 1; }; done
: > "$O/merge_list.txt"; for c in $(seq 1 22); do echo "$O/chr${c}_clean" >> "$O/merge_list.txt"; done
plink2 --pmerge-list "$O/merge_list.txt" --make-pgen --threads "$T" --memory "$MEM" --out "$O/v2_imputed_clean"
echo "### final fileset"; ls -la "$O"/v2_imputed_clean.*; grep -vc '^#' "$O/v2_imputed_clean.pvar"; grep -vc '^#' "$O/v2_imputed_clean.psam" || true
plink2 --pfile "$O/v2_imputed_clean" --keep "$O/keep_129.txt" --freq --hardy --threads "$T" --memory "$MEM" --out "$O/stats129"
python -I "${CODE_DIR}/04_post_imputation_QC/13_clean_tables.py" --clean-dir "$O" --tables-dir "$D/info_tables"
