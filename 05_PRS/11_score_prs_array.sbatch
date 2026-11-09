#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/11_score_prs_array.sbatch   (original step script: 16_score_prs.sbatch)
# Produces: plink2 --score of V1/V2/V3 on the array data and summary.
# Run: sbatch 05_PRS/11_score_prs_array.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_score
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 16 - Part C.4/C.5: scores with plink2 --score, standardisation, comparison and figures.
#
# For each version V1/V2/V3:
#   1. 05_PRS/10_build_weights_array.py: merge chr1-22 PRS-CS outputs, align A1 to the ARRAY alleles
#      (strand flip when needed) and keep only variants present in the array -> weights file.
#   2. plink2 --score weights 1 2 3 header cols=+scoresums list-variants
#        * columns: 1 = variant ID (rsID), 2 = scored allele (the PRS-CS effect allele A1, aligned
#          to the array), 3 = weight (PRS-CS posterior-mean BETA).
#        * the score is SCORE1_SUM = sum over variants of BETA x (number of copies of A1).
#          (SCORE1_AVG = sum / number of scored alleles is also written, but SUM is used.)
#        * MISSING genotypes: plink2 default = mean imputation, i.e. a missing call contributes
#          BETA x (2 x ALT allele frequency). Frequencies are fixed with --read-freq from the
#          129 samples that PASS sample QC (--freq on b38_qc_rsid), so that the imputation does
#          not depend on whether QC-failing samples are scored too. ('no-mean-imputation' would
#          instead drop missing calls and divide by fewer alleles.)
#        * no 'center' / 'variance-standardize' option: genotypes enter as raw dosages.
#        * dataset: b38_qc_rsid_all140 (same final variants, all 140 samples) so that donors failing
#          sample QC can be reported, clearly flagged. Standardisation uses QC-pass samples only.
#   3. 05_PRS/12_summarize_prs_array.py: z-scores, V1-V2-V3 correlations, histograms with the 14 donors.
# Env: "prs" (+ matplotlib, installed once below before activating, as in step 01).
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
Q="$PRS/genotypes/plink_qc"
W="$PRS/prscs/weights"
S="$PRS/prscs/scores"
R="$PRS/validation/results/prs"
mkdir -p "$W" "$S" "$R"

eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda run -n prs python -c "import matplotlib" 2>/dev/null; then
  conda install -y -n prs --override-channels --strict-channel-priority -c conda-forge -c bioconda matplotlib
  conda env export -n prs --no-builds > "$PRS/env/environment.yml"
fi
conda activate prs
T="${SLURM_CPUS_PER_TASK:-1}"

echo "### allele frequencies of the QC cohort (129 samples) for mean imputation"
plink2 --bfile "$Q/b38_qc_rsid" --freq --threads "$T" --out "$S/qc_cohort"

for V in V1 V2 V3; do
  echo "### $V"
  python -I "${CODE_DIR}/05_PRS/10_build_weights_array.py" \
    --prefix "$PRS/prscs/out/$V/$V" --bim "$Q/b38_qc_rsid.bim" --out "$W/${V}_weights_array.txt"
  plink2 --bfile "$Q/b38_qc_rsid_all140" \
         --read-freq "$S/qc_cohort.afreq" \
         --score "$W/${V}_weights_array.txt" 1 2 3 header cols=+scoresums list-variants \
         --threads "$T" --out "$S/$V"
done

python -I "${CODE_DIR}/05_PRS/12_summarize_prs_array.py" \
  --scoredir "$S" \
  --qc-table "$PRS/validation/results/qc/samples/sample_qc_table.tsv" \
  --donor-key "$PRS/metadata/donor_key.tsv" \
  --logdir "$PRS/validation/logs" \
  --weightdir "$W" \
  --outdir "$R"
