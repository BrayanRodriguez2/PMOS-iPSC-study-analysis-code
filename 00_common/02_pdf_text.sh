#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 00_common/02_pdf_text.sh   (original step script: 60_pdf_text.sh)
# Produces: Text of PDF files with pdftotext (figure book, supplementary information).
# Run: bash 00_common/02_pdf_text.sh <args>   (light job; source config/config.sh first)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
# Step 60a - text of the verification PDFs (figure book, Additional_file_2) with pdftotext from the existing "pdftools" env (light, login node).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate pdftools
D="${PROJECT_DIR}/DEG_validation/package"; mkdir -p "$D"
for f in "$@"; do
  pdftotext -layout "$D/$f.pdf" "$D/$f.txt"; echo "$f: $(pdfinfo "$D/$f.pdf" | grep Pages) | $(wc -l < "$D/$f.txt") lines"
done
