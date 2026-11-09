#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/06_qc_only_download.sh   (original step script: 32_topmed_qc_download.sh)
# Produces: Download of the server QC-only report (link passed as an argument).
# Run: python 03_imputation_TOPMed/06_qc_only_download.sh <args>   (or from a SLURM script; see the usage line below)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Step 32 - Download the results of the TOPMed Imputation Server job "imputation_ipsc" (mode: Quality Control Only,
# panel TOPMed r3 2023, build hg38; upload = the 22 v1 VCFs of step 21: 140 samples x 510,607 SNPs).
# The links are the private share links given by the analyst (they expire after a few days; do not share them).
# Light job (5 small files, plain wget) -> run on the login node; nothing is uploaded anywhere.
# Output: PRS/imputation/topmed_qc_only/{qc_report.txt, quality-control.html, chunks-excluded.txt, snps-excluded.txt,
#         snps-typed-only.txt, md5sums_downloaded.txt}; the downloaded files are DATA (never executed; the html is only parsed as text).
# Checks: not empty, not an error page (no "<html" in the txt files), line counts and sizes are printed.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
D="${PROJECT_DIR}/PRS/imputation/topmed_qc_only"
mkdir -p "$D"; cd "$D"
wget -nv --tries=2 --timeout=60 "${TOPMED_URL_QC_REPORT:?private share URL of qc_report.txt (given by the TOPMed server)}"
wget -nv --tries=2 --timeout=60 "${TOPMED_URL_QUALITY_CONTROL:?private share URL of quality-control.html (given by the TOPMed server)}"
wget -nv --tries=2 --timeout=60 "${TOPMED_URL_CHUNKS_EXCLUDED:?private share URL of chunks-excluded.txt (given by the TOPMed server)}"
wget -nv --tries=2 --timeout=60 "${TOPMED_URL_SNPS_EXCLUDED:?private share URL of snps-excluded.txt (given by the TOPMed server)}"
wget -nv --tries=2 --timeout=60 "${TOPMED_URL_SNPS_TYPED_ONLY:?private share URL of snps-typed-only.txt (given by the TOPMed server)}"
echo "### sizes and line counts"
for f in qc_report.txt quality-control.html chunks-excluded.txt snps-excluded.txt snps-typed-only.txt; do
  [ -s "$f" ] || { echo "EMPTY or missing: $f"; exit 1; }
  printf "%-24s %10s bytes %8s lines\n" "$f" "$(stat -c %s "$f")" "$(wc -l < "$f")"
done
echo "### error-page check (txt files must not look like HTML)"
for f in qc_report.txt chunks-excluded.txt snps-excluded.txt snps-typed-only.txt; do
  if grep -q -i -m1 '<html' "$f"; then echo "WARNING: $f looks like an HTML page"; else echo "$f: plain text OK"; fi
done
grep -q -i -m1 '<html' quality-control.html && echo "quality-control.html: HTML OK" || echo "WARNING: quality-control.html has no <html tag"
md5sum qc_report.txt quality-control.html chunks-excluded.txt snps-excluded.txt snps-typed-only.txt | tee md5sums_downloaded.txt
echo "### heads"
for f in qc_report.txt chunks-excluded.txt snps-excluded.txt snps-typed-only.txt; do echo "--- $f"; sed -n '1,4p' "$f" | cut -c1-220; done
