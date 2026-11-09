#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/08_sample_checks_concordance.sbatch   (original step script: 43a_concordance.sbatch)
# Produces: Concordance, heterozygosity and outlier checks of the samples after imputation.
# Run: sbatch 04_post_imputation_QC/08_sample_checks_concordance.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=imp_conc
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=1-22%8
#SBATCH --output=logs/%x_%A_%a.out
# Step 43a - per-sample concordance between the original genotype (v2 upload VCF) and the imputed dosage at the GENOTYPED sites, one chromosome per task
# (see 04_post_imputation_QC/08_sample_checks.py). Input: PRS/topmed/upload_v2/chrN.vcf.gz, unzipped/chrN.dose.vcf.gz (+ .tbi). Output: imputation_qc/work_samples/concord_chrN.tsv
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
C=$SLURM_ARRAY_TASK_ID
PRS="${PROJECT_DIR}/PRS"; D="$PRS/imputation/topmed_imputation_v2"; W="$PRS/validation/results/imputation_qc/work_samples"
mkdir -p "$W"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
bcftools query -l "$PRS/topmed/upload_v2/chr$C.vcf.gz" > "$W/samples.txt"
bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' "$PRS/topmed/upload_v2/chr$C.vcf.gz" > "$W/orig_chr$C.tsv"
bcftools view -i 'INFO/TYPED=1' -Ou "$D/unzipped/chr$C.dose.vcf.gz" | bcftools query -f '%POS\t%REF\t%ALT[\t%DS]\n' > "$W/dose_typed_chr$C.tsv"
echo "chr$C: original sites $(wc -l < "$W/orig_chr$C.tsv") | typed sites in the dose VCF $(wc -l < "$W/dose_typed_chr$C.tsv")"
python -I "${CODE_DIR}/04_post_imputation_QC/08_sample_checks.py" concord --orig "$W/orig_chr$C.tsv" --dose "$W/dose_typed_chr$C.tsv" --samples "$W/samples.txt" --out "$W/concord_chr$C.tsv"
# the two intermediate text files (orig_chr*.tsv, dose_typed_chr*.tsv, ~20 MB each) are kept: nothing is deleted without the analyst's OK
