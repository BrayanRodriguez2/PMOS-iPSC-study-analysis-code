#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/03_snp_map.sbatch   (original step script: 23b_snp_map.sbatch)
# Produces: Locus table of the array manifest (BPM).
# Run: sbatch 02_contamination_identity/03_snp_map.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=snpmap
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 23b - Part A of the contamination diagnostic: locus table of the .bpm manifest (name, chromosome, position, in
# manifest order = the order of every per-locus array inside the GTC files). `bcftools +idat2gtc --snp-map` writes it.
# IlluminaBeadArrayFiles is not on PyPI, so step 24 reads the GTC files with a small documented parser of its own
# and uses this table instead of Illumina's BPM reader.
# Output: PRS/intensity/bpm_snp_map.csv ; prints header, first rows, loci per chromosome.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
G=${RAW_GENOTYPE_DIR}
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
mkdir -p "$PRS/intensity"
bcftools +idat2gtc --bpm "$G/ESTchip-2_GSAv3-MD_20032937X364826_A3.bpm" --egt "$G/NG_1_23_26012021_v1.egt" \
  --snp-map "$PRS/intensity/bpm_snp_map.csv" 2>&1 | sed -n '1,5p' || true
ls -la "$PRS/intensity/"
echo "### header + first rows"; sed -n '1,4p' "$PRS/intensity/bpm_snp_map.csv"
echo "### number of lines"; wc -l < "$PRS/intensity/bpm_snp_map.csv"
