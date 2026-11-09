#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/01_fetch_download_script.sh   (original step script: 38a_imputation_fetch_script.sh)
# Produces: Fetches (without running) the download script of the imputation results for inspection.
# Run: bash 04_post_imputation_QC/01_fetch_download_script.sh <args>   (light job; source config/config.sh first)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
# Step 38a - Fetch the TOPMed "download all" script of the full imputation job (v2 VCFs: 140 samples x 510,412 SNPs; panel r3, hg38,
# Eagle, R2 >= 0.3, AES256 results) into PRS/imputation/topmed_imputation_v2/ WITHOUT running it, and describe it, so it can be read
# before it runs. The user gave the command `curl -sL <share url> | bash`; saving it first is a safety step, the share URL is the analyst's.
# Light job (one small file) -> login node. Nothing is uploaded anywhere. URL query strings (signed tokens) are never printed.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
D="${PROJECT_DIR}/PRS/imputation/topmed_imputation_v2"
mkdir -p "$D"; cd "$D"
echo "### folder content before"; ls -la
curl -sL "${1:?private share URL of the imputation job (given by the TOPMed server)}" -o download_all.sh
echo "### size / lines / first bytes"; stat -c %s download_all.sh; wc -l < download_all.sh; head -c 300 download_all.sh | tr -c '[:print:]\n' '?' | sed -E 's/\?[^ ]*//g'; echo
echo "### commands used (first word of each non-comment line)"
{ grep -v '^\s*#' download_all.sh || true; } | awk '{print $1}' | sort | uniq -c | sort -rn | sed -n '1,10p'
echo "### lines (query strings cut, so no tokens are printed)"
sed -E 's/\?[^ "'"'"']*//g' download_all.sh | cut -c1-230 | sed -n '1,70p'
