#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/03_sample_qc.sbatch   (original step script: 09_qc_samples.sbatch)
# Produces: Sample QC (call rate, heterozygosity, sex, kinship) and the list of excluded samples.
# Run: sbatch 01_genotype_QC/03_sample_qc.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_qc_samples
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 09 - SAMPLE quality control of the 140 donor genotypes (b38, GRCh38).
#
# Metrics (plink2), each computed on the set of variants appropriate for it:
#   * call rate  : autosomal variants with variant call rate >= 95 % (so that a few
#                  bad probes do not drive the sample statistic).
#   * heterozygosity (F): autosomal, call rate >= 98 %, MAF >= 5 %, HWE p >= 1e-6,
#                  LD-pruned (indep-pairwise 50 5 0.2) -> standard for inbreeding / contamination.
#   * sex        : chrX F statistic (LD-pruned X variants, MAF >= 5 %) with plink2 --check-sex.
#                  Cut-offs max-female-xf=0.2 / min-male-xf=0.8 are plink's defaults for a
#                  F-statistic; chrY call rate is reported too (not used as a hard filter because
#                  array Y probes give noisy calls in women).
#   * relatedness: KING-robust kinship on the same pruned autosomal set.
# The decision logic (thresholds, reasons) is in 01_genotype_QC/03_sample_qc.py and is fully
# traceable in results/qc/samples/sample_qc_table.tsv.
# Inputs are only READ; outputs go to PRS/validation/results/qc/samples/.
# Resources: 140 samples x 777 k variants is small; 4 CPUs / 16 GB on `testing` is plenty.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail

PRS="${PROJECT_DIR}/PRS"
RAW="$PRS/genotypes/plink_bin/b38"
OUT="$PRS/validation/results/qc/samples"
mkdir -p "$OUT"

eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prs
T="${SLURM_CPUS_PER_TASK:-1}"

echo "### 1/6 LD-pruned autosomal set (call rate >= 98 %, MAF >= 5 %, HWE p >= 1e-6)"
plink2 --bfile "$RAW" --chr 1-22 --geno 0.02 --maf 0.05 --hwe 1e-6 0 \
       --indep-pairwise 50 5 0.2 --threads "$T" --out "$OUT/auto_prune"

echo "### 2/6 sample call rate (autosomal variants with call rate >= 95 %)"
plink2 --bfile "$RAW" --chr 1-22 --geno 0.05 --missing sample-only \
       --threads "$T" --out "$OUT/autosomal"

echo "### 3/6 autosomal heterozygosity (F)"
plink2 --bfile "$RAW" --extract "$OUT/auto_prune.prune.in" --het \
       --threads "$T" --out "$OUT/het"

echo "### 4/6 KING-robust kinship (pairs >= 0.0884 reported)"
plink2 --bfile "$RAW" --extract "$OUT/auto_prune.prune.in" \
       --make-king-table --king-table-filter 0.0884 \
       --threads "$T" --out "$OUT/king"

echo "### 5/6 sex: chrX F statistic + chrY call rate"
plink2 --bfile "$RAW" --chr X --geno 0.05 --maf 0.05 --indep-pairwise 50 5 0.2 \
       --threads "$T" --out "$OUT/x_prune"
plink2 --bfile "$RAW" --chr X --extract "$OUT/x_prune.prune.in" \
       --check-sex max-female-xf=0.2 min-male-xf=0.8 \
       --threads "$T" --out "$OUT/sexcheck"
plink2 --bfile "$RAW" --chr Y --missing sample-only --threads "$T" --out "$OUT/chrY"

echo "### 6/6 decision"
python -I "${CODE_DIR}/01_genotype_QC/03_sample_qc.py" \
  --fam "$RAW.fam" \
  --smiss "$OUT/autosomal.smiss" \
  --ysmiss "$OUT/chrY.smiss" \
  --het "$OUT/het.het" \
  --sexcheck "$OUT/sexcheck.sexcheck" \
  --kin0 "$OUT/king.kin0" \
  --donor-key "$PRS/metadata/donor_key.tsv" \
  --outdir "$OUT"
