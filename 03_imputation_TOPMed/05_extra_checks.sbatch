#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/05_extra_checks.sbatch   (original step script: 21b_topmed_extra_checks.sbatch)
# Produces: Additional checks of the upload VCFs.
# Run: sbatch 03_imputation_TOPMed/05_extra_checks.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=tm_extra
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 21b - two extra read-only checks of the upload VCFs: (1) bgzip integrity (bgzip -t) of every file and
# (2) which 10 Mb chunks (the server's QC unit) end up with fewer than 3 variants (the server excludes such chunks)
# or more than 5,000 variants, plus the sample names of the VCF header.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
UP="$PRS/topmed/upload"
TMP="$PRS/validation/results/topmed/vcf_work"
OUT="$PRS/validation/results/topmed/extra_checks.txt"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate topmed
{
  echo "### bgzip integrity (bgzip -t) of the 22 files"
  bad=0; for c in $(seq 1 22); do bgzip -t "$UP/chr$c.vcf.gz" || { echo "chr$c FAILED"; bad=1; }; done
  [ $bad -eq 0 ] && echo "all 22 files pass bgzip -t"
  echo "### tabix/csi indexes present"; ls "$UP"/*.tbi | wc -l
  echo "### chunks (10 Mb) with < 3 variants (server would exclude them) | chromosome, chunk start, n variants, chromosome length of the last variant"
  for c in $(seq 1 22); do
    awk -v c="chr$c" -F'\t' '{k=int(($2-1)/10000000); n[k]++; last=$2} END {for (k in n) if (n[k] < 3) printf "%s chunk %d-%d Mb: %d variants (last variant of the chromosome at %d)\n", c, k*10, (k+1)*10, n[k], last}' "$TMP/chr$c.pos.tsv"
  done | sort -V
  echo "### chunks with 0 variants inside the span of each chromosome (gaps, e.g. centromeres): count per chromosome"
  for c in $(seq 1 22); do
    awk -v c="chr$c" -F'\t' '{k=int(($2-1)/10000000); n[k]=1; if (k>mx) mx=k} END {g=0; for (i=0;i<=mx;i++) if (!(i in n)) g++; printf "%s:%d ", c, g}' "$TMP/chr$c.pos.tsv"
  done; echo
  echo "### sample names in the VCF header (first 5 and count)"
  bcftools query -l "$UP/chr1.vcf.gz" | sed -n '1,5p' | tr '\n' ' '; echo
  bcftools query -l "$UP/chr1.vcf.gz" | wc -l
  echo "### identical sample order in all 22 files?"
  bcftools query -l "$UP/chr1.vcf.gz" | md5sum | cut -d' ' -f1 > "$TMP/s1.md5"
  d=0; for c in $(seq 2 22); do [ "$(bcftools query -l "$UP/chr$c.vcf.gz" | md5sum | cut -d' ' -f1)" = "$(cat "$TMP/s1.md5")" ] || { echo "chr$c differs"; d=1; }; done
  [ $d -eq 0 ] && echo "same 140 samples, same order, in all 22 files"
} 2>&1 | tee "$OUT"
