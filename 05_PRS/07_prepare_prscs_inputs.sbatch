#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/07_prepare_prscs_inputs.sbatch   (original step script: 14_prepare_prscs_inputs.sbatch)
# Produces: Summary statistics in PRS-CS format, n_gwas and the whole-panel .bim.
# Run: sbatch 05_PRS/07_prepare_prscs_inputs.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_prep
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 14 - Part C.1/C.2: inputs for PRS-CS.
#   (1) sumstats in PRS-CS format (SNP A1 A2 BETA SE) + median n -> prscs/inputs/
#   (2) "bim" for version V3: the full SNP list of the UKBB LD panel (snpinfo) written as a
#       PLINK .bim (CHR SNP 0 BP A1 A2). With it PRS-CS returns weights for ALL HapMap3 SNPs that
#       are in the sumstats, independent of what the array contains (article-like PRS).
#   (3) record how many SNPs each version can use BEFORE running (rsID overlap of sumstats /
#       panel / bim), so that the PRS-CS logs can be compared with the expectation.
# Versions:  V1 main  = UKBB panel + QC array .bim   (genotypes/plink_qc/b38_qc_rsid.bim)
#            V2 sens. = 1000G panel + QC array .bim
#            V3 article-like = UKBB panel + panel SNP list as bim; weights then applied to the array
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
IN="$PRS/prscs/inputs"
mkdir -p "$IN"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs

python -I "${CODE_DIR}/05_PRS/07_prepare_prscs_inputs.py" \
  --sumstats "${PROJECT_DIR}/gwas/GCST90570579_harmonised/GCST90570579.h.tsv.gz" --outdir "$IN"

# V3 pseudo-bim from the UKBB panel: snpinfo columns CHR SNP BP A1 A2 MAF -> bim CHR SNP 0 BP A1 A2
awk 'NR>1 {print $1"\t"$2"\t0\t"$3"\t"$4"\t"$5}' "$PRS/reference/ldblk_ukbb_eur/snpinfo_ukbb_hm3" > "$IN/V3_panel_ukbb.bim"
wc -l "$IN/V3_panel_ukbb.bim"

echo "### expected rsID overlap per version"
S="$IN/sumstats_prscs.txt"
for V in "UKBB:$PRS/reference/ldblk_ukbb_eur/snpinfo_ukbb_hm3" "1KG:$PRS/reference/ldblk_1kg_eur/snpinfo_1kg_hm3"; do
  name=${V%%:*}; file=${V#*:}
  # rsIDs in (sumstats & panel) and in (sumstats & panel & QC array bim)
  awk -v name="$name" 'FILENAME==ARGV[1] && FNR>1 {s[$1]=1; next}
       FILENAME==ARGV[2] && FNR>1 {if ($2 in s) {p[$2]=1; np++}; next}
       {if ($2 in p) nb++}
       END {print name " panel: sumstats&panel=" np " | sumstats&panel&QC-array=" nb}' \
      "$S" "$file" "$PRS/genotypes/plink_qc/b38_qc_rsid.bim"
done
