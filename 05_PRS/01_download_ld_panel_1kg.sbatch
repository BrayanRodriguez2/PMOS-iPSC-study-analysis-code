#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/01_download_ld_panel_1kg.sbatch   (original step script: 03_download_ld_panel.sbatch)
# Produces: PRS-CS 1000 Genomes EUR LD panel.
# Run: sbatch 05_PRS/01_download_ld_panel_1kg.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_ldref
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 03 - Download the PRS-CS European LD reference panel (1000 Genomes EUR).
#
# Source: PRS-CS README (https://github.com/getian107/PRScs), link "EUR reference"
#         ldblk_1kg_eur.tar.gz (~4.56 GB, Dropbox). The panel is HapMap3-restricted
#         and contains `snpinfo_1kg_hm3` (SNP list: CHR SNP BP A1 A2 MAF) plus one
#         ldblk_1kg_chr*.hdf5 per chromosome.
# Output: PRS/reference/ldblk_1kg_eur/  (+ the original .tar.gz and its checksums)
# Resources: network-bound, 1 CPU / 4 GB is enough; `main` because it can exceed
#            the 2 h limit of `testing` on a slow link. wget -c lets it resume.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail

REF_DIR="${PROJECT_DIR}/PRS/reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

URL="https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=1"
wget -c --progress=dot:giga -O ldblk_1kg_eur.tar.gz "$URL"

ls -la ldblk_1kg_eur.tar.gz
# Dropbox publishes no checksum: record ours so the copy can be re-verified later.
md5sum ldblk_1kg_eur.tar.gz | tee ldblk_1kg_eur.tar.gz.md5
gzip -t ldblk_1kg_eur.tar.gz && echo "gzip integrity OK"

tar -zxf ldblk_1kg_eur.tar.gz
echo "### extracted content:"
ls -la ldblk_1kg_eur | head -40
echo "### snpinfo file(s):"
find ldblk_1kg_eur -name 'snpinfo*' -exec ls -la {} \; -exec head -3 {} \;
