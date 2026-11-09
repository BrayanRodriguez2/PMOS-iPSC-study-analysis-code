#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/21_prscs_log_counts.sh   (original step script: 51d_prscs_log_counts.sh)
# Produces: SNPs seen by PRS-CS in each version (from the job logs).
# Run: bash 05_PRS/21_prscs_log_counts.sh <args>   (light job; source config/config.sh first)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
# Step 51d - SNPs seen by PRS-CS in each version = sum over chromosomes of the log line "N common SNPs in the reference, sumstats, and validation set"
# (reads the SLURM logs only; imputed V1 job 3366877, V2 3366878; V3 = array of step 15, logs prscs_V3_*). Output -> results/prs_imputed/prscs_snps_seen.tsv
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
L="${PROJECT_DIR}/PRS/validation/logs"; R="${PROJECT_DIR}/PRS/validation/results/prs_imputed"; mkdir -p "$R"
{
printf "version\tlog_glob\tchromosomes_with_log\tcommon_snps_ref_sumstats_bim\n"
for spec in "V1_imputed:prscs_imp_V1_3366877_*.out" "V2_imputed:prscs_imp_V2_3366878_*.out" "V3_step15:prscs_V3_*_*.out" "V1_array_step15:prscs_V1_*_*.out"; do
  name=${spec%%:*}; glob=${spec#*:}
  n=0; tot=0
  for f in $L/$glob; do
    [ -e "$f" ] || continue
    c=$(tr '\r' '\n' < "$f" | grep -m1 'common SNPs in the reference' | awk '{print $2}')
    [ -n "$c" ] && { n=$((n+1)); tot=$((tot+c)); }
  done
  printf "%s\t%s\t%s\t%s\n" "$name" "$glob" "$n" "$tot"
done
} | tee "$R/prscs_snps_seen.tsv"
