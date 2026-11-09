#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/15_prscs_bim_imputed.sbatch   (original step script: 49_prscs_bim_imputed.sbatch)
# Produces: PRS-CS .bim with panel rsIDs matched to the imputed data by position and alleles.
# Run: sbatch 05_PRS/15_prscs_bim_imputed.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=bim_imp
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
# Step 49 - PRS-CS bim (rsID of the panel, matched by GRCh38 position + alleles) and rsID <-> chr:pos:ref:alt map for the imputed clean file (see 05_PRS/15_prscs_bim_imputed.py).
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
python -I "${CODE_DIR}/05_PRS/15_prscs_bim_imputed.py" --variant-table "$PRS/imputation/clean/variant_table.tsv.gz" \
  --snpinfo-ukbb "$PRS/reference/ldblk_ukbb_eur/snpinfo_ukbb_hm3" --snpinfo-1kg "$PRS/reference/ldblk_1kg_eur/snpinfo_1kg_hm3" \
  --sumstats-gz "${PROJECT_DIR}/gwas/GCST90570579_harmonised/GCST90570579.h.tsv.gz" --sumstats-prscs "$PRS/prscs/inputs/sumstats_prscs.txt" \
  --outdir "$PRS/prscs/inputs_imputed"
wc -l "$PRS"/prscs/inputs_imputed/*.bim
