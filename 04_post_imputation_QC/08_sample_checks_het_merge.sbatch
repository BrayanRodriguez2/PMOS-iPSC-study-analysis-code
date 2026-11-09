#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/08_sample_checks_het_merge.sbatch   (original step script: 43b_het_merge.sbatch)
# Produces: Concordance, heterozygosity and outlier checks of the samples after imputation.
# Run: sbatch 04_post_imputation_QC/08_sample_checks_het_merge.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=imp_het
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --output=logs/%x_%j.out
# Step 43b - heterozygosity/F of the imputed best-guess genotypes (R2 >= 0.8, MAF >= 0.05) in the 66 windows of step 42 for all 140 samples,
# then merge with the concordance counts of 43a and the sample QC table; table, figure and report (see 04_post_imputation_QC/08_sample_checks.py).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; D="$PRS/imputation/topmed_imputation_v2"; O="$PRS/validation/results/imputation_qc"; W="$O/work_samples"
mkdir -p "$W"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
S="${CODE_DIR}"
python -I "$S/04_post_imputation_QC/07_freq_vs_1000g.py" windows --tables-dir "$D/info_tables" --out "$W/windows.tsv"
bcftools query -l "$D/unzipped/chr1.dose.vcf.gz" > "$W/samples140.txt"
awk 'NR>1 {print $2}' "$PRS/validation/results/qc/samples/keep_samples.txt" > "$W/samples129.txt"
: > "$W/het_gt.tsv"
for c in $(seq 1 22); do
  awk -F'\t' -v c="chr$c" '$1==c' "$W/windows.tsv" > "$W/win_chr$c.tsv"
  bcftools view -R "$W/win_chr$c.tsv" -i 'INFO/R2>=0.8 && INFO/MAF>=0.05' -Ou "$D/unzipped/chr$c.dose.vcf.gz" | bcftools query -f '[%GT\t]\n' | sed 's/\t$//' >> "$W/het_gt.tsv"
done
echo "variants (R2 >= 0.8, MAF >= 0.05) in the windows: $(wc -l < "$W/het_gt.tsv")"
python -I "$S/04_post_imputation_QC/08_sample_checks.py" het --gt "$W/het_gt.tsv" --samples "$W/samples140.txt" --keep129 "$W/samples129.txt" --out "$O/het_imputed_windows.tsv"
python -I "$S/04_post_imputation_QC/08_sample_checks.py" merge --counts-dir "$W" --het "$O/het_imputed_windows.tsv" --qc-table "$PRS/validation/results/qc/samples/sample_qc_table.tsv" --outdir "$O"
