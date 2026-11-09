#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/06_dump_package_tables.sh   (original step script: 62_xlsx_dump.sh)
# Produces: TSV dumps of the supplementary tables and of the clinical workbook.
# Run: bash 06_RNAseq_counts_and_DE_validation/06_dump_package_tables.sh <args>   (light job; source config/config.sh first)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
# Step 62a - dump the Excel files of the verification package (Additional_file_1 = supplementary tables S1-S15, the corrected clinical workbook) to TSV
# with 00_common/01_xlsx_to_tsv.py (standard library only). Read-only copies of the package live in DEG_validation/package.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
D="${PROJECT_DIR}/DEG_validation/package"; S="${CODE_DIR}"
bash "$S/probe_py.sh" 00_common/01_xlsx_to_tsv.py "$D/Additional_file_1.xlsx" "$D/af1_tsv"
bash "$S/probe_py.sh" 00_common/01_xlsx_to_tsv.py "$D/Clinical_corrected.xlsx" "$D/clinical_tsv"
