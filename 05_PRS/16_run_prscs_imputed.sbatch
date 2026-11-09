#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/16_run_prscs_imputed.sbatch   (original step script: 50_run_prscs_imputed.sbatch)
# Produces: PRS-CS per chromosome on the imputed data (array job; V1, V2).
# Run: sbatch 05_PRS/16_run_prscs_imputed.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prscs_imp
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=1-22
#SBATCH --output=logs/%x_%A_%a.out
# ---------------------------------------------------------------------------
# Step 50 - PRS-CS on the IMPUTED dataset: same program, summary statistics, n_gwas, phi and seed as step 15 (copy of 05_PRS/09_run_prscs_array.sbatch with only the
# .bim changed). One array task per chromosome. Choose the version with --export=ALL,VERSION=V1|V2:
#   V1 (main)         UKBB EUR LD panel, bim = panel SNPs found in the clean imputed file (rsID of the panel; step 49: prscs/inputs_imputed/V1_imputed.bim)
#   V2 (stability)    1000G EUR LD panel, bim = V2_imputed.bim
#   V3 (article-like) is NOT re-run: its bim is the whole UKBB panel, independent of the genotypes, so the PRS-CS output of step 15 (prscs/out/V3, same
#                     sumstats / n_gwas / phi / seed) is re-used as is.
# Parameters: --phi=1e-4 --n_gwas=412081 (median of n; prscs/inputs/n_gwas.txt) --a=1 --b=0.5 (defaults) --n_iter=1000 --n_burnin=500 --thin=5 --seed=2026.
# Output: prscs/out_imputed/<VERSION>/<VERSION>_pst_eff_a1_b0.5_phi1e-04_chr<N>.txt  (CHR SNP BP A1 A2 BETA)
# Env "prscs" (numpy 1.26; see 14b). Resources as step 15 (2 CPUs / 16 GB per chromosome).
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
V="${VERSION:?set VERSION=V1|V2 with --export=ALL,VERSION=...}"
PRS="${PROJECT_DIR}/PRS"
IN="$PRS/prscs/inputs"
case "$V" in
  V1) REF="$PRS/reference/ldblk_ukbb_eur"; BIM="$PRS/prscs/inputs_imputed/V1_imputed" ;;
  V2) REF="$PRS/reference/ldblk_1kg_eur";  BIM="$PRS/prscs/inputs_imputed/V2_imputed" ;;
  *) echo "unknown VERSION $V (V3 is re-used from step 15)"; exit 1 ;;
esac
CHR="${SLURM_ARRAY_TASK_ID:?run as an array job}"
OUT="$PRS/prscs/out_imputed/$V"
mkdir -p "$OUT"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prscs
echo "### $V (imputed) chr$CHR | ref $REF | bim $BIM | n_gwas $(cat "$IN/n_gwas.txt")"
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
