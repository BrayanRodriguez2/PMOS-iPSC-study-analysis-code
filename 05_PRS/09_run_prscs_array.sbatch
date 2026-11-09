#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/09_run_prscs_array.sbatch   (original step script: 15_run_prscs.sbatch)
# Produces: PRS-CS per chromosome (array job; versions V1, V2, V3 on the un-imputed array).
# Run: sbatch 05_PRS/09_run_prscs_array.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prscs
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=1-22
#SBATCH --output=logs/%x_%A_%a.out
# ---------------------------------------------------------------------------
# Step 15 - Part C.3: PRS-CS (Ge et al. 2019, Nat Commun 10:1776), one array task per
# chromosome, phi fixed at 1e-4 (the value used in the paper, i.e. no automatic phi).
#
# Submit one array per version (chosen with the VERSION environment variable):
#   sbatch --export=ALL,VERSION=V1 scripts/05_PRS/09_run_prscs_array.sbatch   # main: UKBB panel + QC array .bim
#   sbatch --export=ALL,VERSION=V2 scripts/05_PRS/09_run_prscs_array.sbatch   # sensitivity: 1000G panel + QC array .bim
#   sbatch --export=ALL,VERSION=V3 scripts/05_PRS/09_run_prscs_array.sbatch   # article-like: UKBB panel + panel SNP list as bim
# Test first on one chromosome: sbatch --partition=testing --time=02:00:00 --array=22 --export=ALL,VERSION=V1 ...
#
# PRS-CS parameters: --phi=1e-4 --n_gwas=<median of n> --a=1 --b=0.5 (defaults) --n_iter=1000
#   --n_burnin=500 --thin=5 (defaults) --seed=2026 (fixed, for reproducibility).
# Resources: one MCMC chain per chromosome, mostly single-threaded linear algebra; 2 CPUs
#   (BLAS threads) and 16 GB cover the largest chromosome blocks of the UKBB panel.
# Output: prscs/out/<VERSION>/<VERSION>_pst_eff_a1_b0.5_phi1e-04_chr<N>.txt  (CHR SNP BP A1 A2 BETA)
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
V="${VERSION:?set VERSION=V1|V2|V3 with --export=ALL,VERSION=...}"
PRS="${PROJECT_DIR}/PRS"
IN="$PRS/prscs/inputs"
case "$V" in
  V1) REF="$PRS/reference/ldblk_ukbb_eur"; BIM="$PRS/genotypes/plink_qc/b38_qc_rsid" ;;
  V2) REF="$PRS/reference/ldblk_1kg_eur";  BIM="$PRS/genotypes/plink_qc/b38_qc_rsid" ;;
  V3) REF="$PRS/reference/ldblk_ukbb_eur"; BIM="$IN/V3_panel_ukbb" ;;
  *) echo "unknown VERSION $V"; exit 1 ;;
esac
CHR="${SLURM_ARRAY_TASK_ID:?run as an array job}"
OUT="$PRS/prscs/out/$V"
mkdir -p "$OUT"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# PRS-CS runs in env "prscs" (numpy 1.26), NOT "prs": see 05_PRS/08_setup_prscs_env.sbatch for why.
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prscs

echo "### $V chr$CHR | ref $REF | bim $BIM | n_gwas $(cat "$IN/n_gwas.txt")"
cd "$PRS/software/PRScs"
python PRScs.py \
  --ref_dir="$REF" \
  --bim_prefix="$BIM" \
  --sst_file="$IN/sumstats_prscs.txt" \
  --n_gwas="$(cat "$IN/n_gwas.txt")" \
  --phi=1e-4 \
  --chrom="$CHR" \
  --seed=2026 \
  --out_dir="$OUT/$V"
ls -la "$OUT" | grep "chr${CHR}\.txt" || true
