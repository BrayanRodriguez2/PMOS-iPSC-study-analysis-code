#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/04_unzip_summary.sh   (original step script: 39b_unzip_summary.sh)
# Produces: Summary of the per-chromosome unzip checks.
# Run: bash 04_post_imputation_QC/04_unzip_summary.sh <args>   (light job; source config/config.sh first)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
# Step 39b - one table with the per-chromosome checks of the decrypted TOPMed files (from the check_chrN.tsv written by step 39). Light.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
U="${PROJECT_DIR}/PRS/imputation/topmed_imputation_v2"
{ printf "chrom\tdose_vcf_bytes\tinfo_bytes\tvariants_in_vcf\tinfo_lines_minus_1\tsamples\tsame_samples_order_as_v2\n"
  for c in $(seq 1 22); do cat "$U/unzipped/check_chr$c.tsv"; done; } > "$U/unzip_check.tsv"
cat "$U/unzip_check.tsv"
awk -F'\t' 'NR>1 {v+=$4; if ($6!=140 || $7!="yes") bad++} END {print "total variants:", v, "| chromosomes with a sample problem:", bad+0}' "$U/unzip_check.tsv"
echo "note: info_lines_minus_1 = variants + 12 header lines (11 ## lines + #CHROM) in every file"
