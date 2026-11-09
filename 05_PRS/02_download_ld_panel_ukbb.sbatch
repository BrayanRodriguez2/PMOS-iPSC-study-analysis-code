#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/02_download_ld_panel_ukbb.sbatch   (original step script: 13_download_ukbb_panel.sbatch)
# Produces: PRS-CS UK Biobank EUR LD panel.
# Run: sbatch 05_PRS/02_download_ld_panel_ukbb.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_ukbbref
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 13 - Download the PRS-CS UK Biobank European LD reference panel (main panel;
# 1000G EUR from step 03 is kept as sensitivity analysis).
# Source: PRS-CS README (https://github.com/getian107/PRScs), "UKBB EUR reference"
#         ldblk_ukbb_eur.tar.gz (~6.25 GB, Dropbox). Contains snpinfo_ukbb_hm3 and
#         ldblk_ukbb_chr*.hdf5.
# Output: PRS/reference/ldblk_ukbb_eur/ (+ tarball and its md5, kept).
# Network-bound job: 1 CPU / 4 GB; `main` because it may exceed 2 h on a slow link.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
REF_DIR="${PROJECT_DIR}/PRS/reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

URL="https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=1"
wget -c --progress=dot:giga -O ldblk_ukbb_eur.tar.gz "$URL"
ls -la ldblk_ukbb_eur.tar.gz
md5sum ldblk_ukbb_eur.tar.gz | tee ldblk_ukbb_eur.tar.gz.md5   # Dropbox publishes no checksum
gzip -t ldblk_ukbb_eur.tar.gz && echo "gzip integrity OK"
tar -zxf ldblk_ukbb_eur.tar.gz
ls ldblk_ukbb_eur | head -30
find ldblk_ukbb_eur -name 'snpinfo*' -exec ls -la {} \; -exec head -3 {} \;
