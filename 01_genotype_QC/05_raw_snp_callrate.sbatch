#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/05_raw_snp_callrate.sbatch   (original step script: 10c_raw_snp_callrate.sbatch)
# Produces: Call rate per SNP of the raw array in the QC-pass samples.
# Run: sbatch 01_genotype_QC/05_raw_snp_callrate.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_rawvmiss
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 10c - per-variant call rate of the RAW array (all 777,016 variants), computed on the
# 129 samples that pass sample QC (step 09), i.e. exactly the population on which filter S2
# (call rate >= 98 %) of step 10 was applied. Needed to draw the removed tail in figure 05 of
# step 10b, which the table written by step 10 (b38_qc.vmiss, kept variants only) cannot show.
# Only a missingness table is produced (plink2 --missing variant-only): no filtering, no new dataset.
# Consistency check: the number of variants with missing rate <= 0.02 must equal the 746,167
# variants that step 10 kept after S2 (qc_counts.tsv).
# Input (read-only): genotypes/plink_bin/b38 + results/qc/samples/keep_samples.txt
# Output: results/qc/variants/raw_129.vmiss
# Resources: small; 2 CPUs / 8 GB on `testing`, ~1 min.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs

OUT="$PRS/validation/results/qc/variants/raw_129"
plink2 --bfile "$PRS/genotypes/plink_bin/b38" \
       --keep "$PRS/validation/results/qc/samples/keep_samples.txt" \
       --missing variant-only \
       --threads "${SLURM_CPUS_PER_TASK:-1}" --out "$OUT"

echo "### consistency check against step 10 (S2 kept 746,167 variants)"
awk 'NR>1 {n++; if ($5 <= 0.02) ok++} END {printf "variants: %d | missing rate <= 0.02: %d | > 0.02: %d\n", n, ok, n-ok}' "$OUT.vmiss"
grep -P '^S2' "$PRS/validation/results/qc/variants/qc_counts.tsv"
echo "### number of missing calls per variant (among 129 samples), first levels"
awk 'NR>1 {c[$3]++} END {for (m in c) print m, c[m]}' "$OUT.vmiss" | sort -n | head -n 12
