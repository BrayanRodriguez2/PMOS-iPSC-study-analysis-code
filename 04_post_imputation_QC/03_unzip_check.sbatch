#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/03_unzip_check.sbatch   (original step script: 39_unzip_check.sbatch)
# Produces: Decrypts/unzips each chromosome (password read from a file) and checks the VCFs (array job).
# Run: sbatch 04_post_imputation_QC/03_unzip_check.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=imp_unzip
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --array=1-22%8
#SBATCH --output=logs/%x_%A_%a.out
# ---------------------------------------------------------------------------
# Step 39 - Decrypt (AES256 zips) and check the TOPMed imputation results. Array task N = chromosome N.
# The zip password is read ONLY from ${TOPMED_PW_FILE} (chmod 600, put there by the analyst) inside the 7z command; it is never printed,
# logged or copied (no `set -x`; 7z output goes to /dev/null). Note: as for any `7z -p...` call, the argument is briefly visible to other
# users of the same compute node in `ps`.
# Each task: 7z x chr_N.zip -> unzipped/chrN.dose.vcf.gz + chrN.info.gz (the zips only hold these 2 files: no empiricalDose, no index),
#   gzip -t of both, tabix index of the dose VCF (my own derived file, needed for region queries), number of variants in the VCF
#   (from the index) and in the info table, the sample list of the VCF (md5 + comparison with the 140 samples of the v2 upload VCFs).
# The original zips are never modified or deleted.
# Output: PRS/imputation/topmed_imputation_v2/unzipped/{chrN.dose.vcf.gz(.tbi), chrN.info.gz, check_chrN.tsv}
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
C=$SLURM_ARRAY_TASK_ID
PRS="${PROJECT_DIR}/PRS"; D="$PRS/imputation/topmed_imputation_v2"; U="$D/unzipped"
mkdir -p "$U"; cd "$D"
eval "${P7ZIP_SETUP:-true}"
[ -s "${TOPMED_PW_FILE}" ] || { echo "password file missing"; exit 1; }
7z x -p"$(cat "${TOPMED_PW_FILE}")" -o"$U" -y "chr_$C.zip" > /dev/null
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
V="$U/chr$C.dose.vcf.gz"; I="$U/chr$C.info.gz"
[ -s "$V" ] && [ -s "$I" ] || { echo "chr$C: missing or empty output"; exit 1; }
gzip -t "$I"; gzip -t "$V"
bcftools index -f -t "$V"
nv=$(bcftools index -n "$V"); ni=$(( $(zcat "$I" | wc -l) - 1 ))
bcftools query -l "$V" > "$U/samples_chr$C.txt"
ns=$(wc -l < "$U/samples_chr$C.txt"); sm=$(md5sum < "$U/samples_chr$C.txt" | cut -d' ' -f1)
bcftools query -l "$PRS/topmed/upload_v2/chr1.vcf.gz" | md5sum | cut -d' ' -f1 > "$U/.v2_samples.md5.$C"
same=$([ "$sm" = "$(cat "$U/.v2_samples.md5.$C")" ] && echo yes || echo NO); rm -f "$U/.v2_samples.md5.$C"
printf "chr%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$C" "$(stat -c %s "$V")" "$(stat -c %s "$I")" "$nv" "$ni" "$ns" "$same" > "$U/check_chr$C.tsv"
cat "$U/check_chr$C.tsv"
if [ "$C" -eq 22 ]; then
  echo "### VCF header (FORMAT/INFO lines) and first records, chr22"
  bcftools view -h "$V" | grep -E '^##(FORMAT|INFO|source|fileformat)' | cut -c1-200
  bcftools view -H "$V" | sed -n '1,2p' | cut -f1-12 | cut -c1-250
  echo "### info header and first rows"; zcat "$I" | sed -n '1,4p'
fi
