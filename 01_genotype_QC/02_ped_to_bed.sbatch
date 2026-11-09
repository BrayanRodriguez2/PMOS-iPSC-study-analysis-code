#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/02_ped_to_bed.sbatch   (original step script: 02_ped_to_bed.sbatch)
# Produces: GenomeStudio PED/MAP export -> sorted PLINK binary dataset.
# Run: sbatch 01_genotype_QC/02_ped_to_bed.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_ped2bed
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 02 - Convert the donor genotypes (PLINK text .ped/.map) to binary
#           PLINK (.bed/.bim/.fam) with plink2 (env "prs").
#
# Input  (READ-ONLY, original data - never written to):
#   ${RAW_GENOTYPE_DIR}/PLINK_100626_0356/b38.{ped,map}
#   Illumina GSAv3 array, GRCh38, 140 donors, 777,016 variants.
# Output (own area): PRS/genotypes/plink_bin/b38.{bed,bim,fam,log}
# No filtering is applied here (no QC, no sample/variant removal): this is a pure
# format conversion so that the .bim reflects the raw array content.
# Skipped if the output already exists.
# Resources: 435 MB text file; 4 CPUs / 16 GB is generous, `testing` is enough.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail

SRC_PREFIX="${RAW_GENOTYPE_DIR}/PLINK_100626_0356/b38"
OUT_DIR="${PROJECT_DIR}/PRS/genotypes/plink_bin"
mkdir -p "$OUT_DIR"

if [ -s "$OUT_DIR/b38.bed" ] && [ -s "$OUT_DIR/b38.bim" ] && [ -s "$OUT_DIR/b38.fam" ]; then
  echo "### binary PLINK files already exist in $OUT_DIR - nothing to do"
  exit 0
fi

eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prs

plink2 --version

# The original .map is NOT sorted by chromosome ("split chromosome" error when
# writing .bed directly), so it is done in two steps:
#   (1) import + sort variants (chr, position) into a temporary .pgen,
#   (2) convert the sorted .pgen to .bed/.bim/.fam.
# Sorting only changes the variant ORDER; no variant or sample is added/dropped.
TMP="$OUT_DIR/tmp_sorted"
plink2 --pedmap "$SRC_PREFIX" \
       --make-pgen --sort-vars \
       --threads "${SLURM_CPUS_PER_TASK:-1}" \
       --memory 12000 \
       --out "$TMP"
plink2 --pfile "$TMP" \
       --make-bed \
       --threads "${SLURM_CPUS_PER_TASK:-1}" \
       --memory 12000 \
       --out "$OUT_DIR/b38"
rm -f "$TMP".pgen "$TMP".pvar "$TMP".psam   # own temporary intermediates only

echo "### output"
ls -la "$OUT_DIR"
echo "variants: $(wc -l < "$OUT_DIR/b38.bim") | samples: $(wc -l < "$OUT_DIR/b38.fam")"
