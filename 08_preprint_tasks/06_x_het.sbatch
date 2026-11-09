#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/06_x_het.sbatch   (original step script: 83b_x_het.sbatch)
# Produces: Heterozygosity of chromosome X per sample (raw array).
# Run: sbatch 08_preprint_tasks/06_x_het.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=x_het
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# Step 83b - per-sample heterozygosity of chromosome X (fraction of heterozygous calls among called X genotypes) from the RAW array dataset (genotypes/plink_bin/b38, 140 samples, read-only):
# the column 'x_het' of the original analysis' Table S2 (all samples are genetically female, so ~0.33-0.36 is expected). Output: results/array_pca/x_het_nonPAR.scount (X without the pseudo-autosomal regions) and x_het_allX.scount (whole X incl. PAR).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
P="${PROJECT_DIR}/PRS"; O="$P/validation/results/array_pca"; mkdir -p "$O"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
plink2 --bfile "$P/genotypes/plink_bin/b38" --chr X --not-chr PAR1,PAR2 --sample-counts cols=hom,het,missing --out "$O/x_het" 2>&1 | tail -n 5 || true
if [ ! -s "$O/x_het.scount" ]; then plink2 --bfile "$P/genotypes/plink_bin/b38" --chr 23 --sample-counts cols=hom,het,missing --out "$O/x_het"; fi
mv "$O/x_het.scount" "$O/x_het_nonPAR.scount"; mv "$O/x_het.log" "$O/x_het_nonPAR.log"
plink2 --bfile "$P/genotypes/plink_bin/b38" --chr X --sample-counts cols=hom,het,missing --out "$O/x_het_allX" 2>&1 | tail -n 2
head -n 3 "$O/x_het_nonPAR.scount" "$O/x_het_allX.scount"; wc -l "$O/x_het_allX.scount"
