#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/01_idat_env.sbatch   (original step script: 22_idat_env.sbatch)
# Produces: Conda environment 'idat' (bcftools + gtc2vcf plugin, Python).
# Run: sbatch 02_contamination_identity/01_idat_env.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=idat_env
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 22 - Diagnostic of the excess-heterozygosity samples (the focal samples): environment for the array
# intensity analysis (Part A: BAF / LRR from IDAT files).
#
# Creates the conda env "idat" (conda-forge + bioconda only, strict priority):
#   - bcftools + bcftools-gtc2vcf-plugin (G. Genovese, https://github.com/freeseek/gtc2vcf): provides the
#     bcftools plugins +idat2gtc (IDAT -> GTC using the .bpm manifest and .egt cluster file) and +gtc2vcf
#     (GTC -> VCF with GT, BAF, LRR, ...). Bioconda build 1.22.
#   - verifyidintensity (Jun et al. 2012, AJHG 91:839): intensity-based contamination estimator (tried as a
#     second opinion next to our own BAF-regression estimator).
#   - python 3.11 + numpy, pandas, scipy, matplotlib for the summaries and figures.
# Then prints the usage of the two plugins so the next step uses the exact option names of this version.
# Nothing is written outside PRS/env and the log.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
mkdir -p "$PRS/env"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx "idat"; then
  conda create -y -n idat --override-channels --strict-channel-priority -c conda-forge -c bioconda \
    python=3.11 bcftools bcftools-gtc2vcf-plugin verifyidintensity numpy pandas scipy matplotlib
fi
conda env export -n idat --no-builds > "$PRS/env/environment_idat.yml"   # before activating (see step 01)
conda activate idat
echo "### versions"
bcftools --version | sed -n '1,2p'
echo "### plugins"
bcftools plugin -l | grep -E 'gtc2vcf|idat2gtc' || echo "plugins NOT listed (check BCFTOOLS_PLUGINS=$BCFTOOLS_PLUGINS)"
echo "### +idat2gtc usage"; bcftools +idat2gtc 2>&1 | sed -n '1,60p' || true
echo "### +gtc2vcf usage"; bcftools +gtc2vcf 2>&1 | sed -n '1,90p' || true
echo "### verifyIDintensity usage"; (verifyIDintensity --help 2>&1 || true) | sed -n '1,60p'
