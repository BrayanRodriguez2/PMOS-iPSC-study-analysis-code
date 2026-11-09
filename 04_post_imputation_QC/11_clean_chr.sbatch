#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/11_clean_chr.sbatch   (original step script: 45a_clean_chr.sbatch)
# Produces: Per-chromosome clean PLINK2 dosage files (R2/MAF/HWE filters).
# Run: sbatch 04_post_imputation_QC/11_clean_chr.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=imp_clean
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=1-22%6
#SBATCH --output=logs/%x_%A_%a.out
# ---------------------------------------------------------------------------
# Step 45a - Clean imputed dataset for the PRS (PLINK2 pgen/pvar/psam WITH DOSAGE), one chromosome per array task.
# Input : unzipped/chrN.dose.vcf.gz (TOPMed r3 imputation of the v2 VCFs; R2 >= 0.3 filter already applied by the server; the originals are untouched)
# Samples kept in the file: the 135 of the call-rate >= 0.97 sensitivity set (superset of the 129 of the main QC), so the 129 (main) or the 135
#   (sensitivity) are chosen later with --keep keep_129.txt / keep_135.txt; no second file is needed.
# Variant filters, applied in this order (counts per step written to clean/counts_chrN.tsv):
#   0 delivered by the server (R2 >= 0.3)      1 SNPs only (A/C/G/T), biallelic        2 MAF >= 0.01 in the 129 QC samples (plink2 --maf)
#   3 HWE exact p >= 1e-6 in the 129 (hard calls; plink2 --hwe 1e-6 0)                 4 unique IDs  chr:pos:ref:alt  (--rm-dup force-first as a guard)
#   autosomes only (the server output has only chr1-22). MAF and HWE are computed on the 129 although the file holds 135 samples.
# IDs: --set-all-var-ids '@:#:$r:$a' -> e.g. 1:12345:A:G (no "chr"); the original rsID and the R2/ER2/MAF of the server go to clean/variant_table (step 45b).
# Dosage: --vcf ... dosage=DS keeps the imputed dosage (hard calls come from the phased GT). Sample IDs: --double-id (FID = IID, as in the array data).
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
C=$SLURM_ARRAY_TASK_ID
PRS="${PROJECT_DIR}/PRS"; D="$PRS/imputation/topmed_imputation_v2"; O="$PRS/imputation/clean"; W="$O/work"
mkdir -p "$W"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
T=$SLURM_CPUS_PER_TASK; MEM=$((SLURM_MEM_PER_NODE - 2000))
[ -s "$O/keep_129.txt" ] && [ -s "$O/keep_135.txt" ] || { echo "keep lists missing (run 04_post_imputation_QC/10_make_keep_lists.sh first)"; exit 1; }
n0=$(bcftools index -n "$D/unzipped/chr$C.dose.vcf.gz" 2>/dev/null || true)
echo "[$(date +%T)] chr$C: import (SNPs only, dosage=DS, 135 samples)"
plink2 --vcf "$D/unzipped/chr$C.dose.vcf.gz" dosage=DS --double-id --keep "$O/keep_135.txt" --update-sex "$O/sex_all_female.txt" --snps-only just-acgt --max-alleles 2 \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 20 --make-pgen --threads "$T" --memory "$MEM" --out "$W/s1_chr$C"
n1=$(grep -vc '^#' "$W/s1_chr$C.pvar")
echo "[$(date +%T)] chr$C: MAF >= 0.01 in the 129"
plink2 --pfile "$W/s1_chr$C" --keep "$O/keep_129.txt" --maf 0.01 --write-snplist --threads "$T" --memory "$MEM" --out "$W/s2_chr$C"
n2=$(wc -l < "$W/s2_chr$C.snplist")
echo "[$(date +%T)] chr$C: HWE p >= 1e-6 in the 129"
plink2 --pfile "$W/s1_chr$C" --keep "$O/keep_129.txt" --extract "$W/s2_chr$C.snplist" --hwe 1e-6 0 --write-snplist --threads "$T" --memory "$MEM" --out "$W/s3_chr$C"
n3=$(wc -l < "$W/s3_chr$C.snplist")
echo "[$(date +%T)] chr$C: final pgen (135 samples, filtered variants, unique IDs)"
plink2 --pfile "$W/s1_chr$C" --extract "$W/s3_chr$C.snplist" --rm-dup force-first --make-pgen --threads "$T" --memory "$MEM" --out "$O/chr${C}_clean"
n4=$(grep -vc '^#' "$O/chr${C}_clean.pvar")
printf "chr%s\t%s\t%s\t%s\t%s\t%s\n" "$C" "${n0:-NA}" "$n1" "$n2" "$n3" "$n4" > "$O/counts_chr$C.tsv"
cat "$O/counts_chr$C.tsv"
if [ "$C" -eq 22 ]; then echo "### pvar header and first rows (chr22)"; grep '^#' "$O/chr22_clean.pvar" | tail -n 4 | cut -c1-200; grep -v '^#' "$O/chr22_clean.pvar" | sed -n '1,3p' | cut -c1-200; fi
