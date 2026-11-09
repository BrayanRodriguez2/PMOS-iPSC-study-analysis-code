#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/18_score_imputed.sbatch   (original step script: 51_score_imputed.sbatch)
# Produces: plink2 --score on imputed dosages (V1, V2, V3, V1 R2>=0.8, GRS29).
# Run: sbatch 05_PRS/18_score_imputed.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_score_imp
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 51 - weights back to the clean-file IDs and scoring with plink2 --score on the imputed DOSAGES (adapted from step 05_PRS/11_score_prs_array.sbatch).
#   Sets (env SETS, default all):
#     V1          PRS-CS, UKBB EUR panel, bim = panel SNPs in the clean file  (prscs/out_imputed/V1; main)
#     V2          PRS-CS, 1000G EUR panel                                      (prscs/out_imputed/V2)
#     V3          PRS-CS, UKBB panel, whole panel as bim (output of step 15, re-used: prscs/out/V3), applied to the panel SNPs present in the clean file
#     V1_r2ge0.8  the V1 weights restricted to clean variants with server R2 >= 0.8 (same PRS-CS run, fewer scored SNPs)
#     V1_chr22    V1 weights of chromosome 22 only (to verify the scoring independently against the VCF dosages, step 51c)
#     GRS29       29 ST2 lead SNPs: risk allele and effect (51b)
#   plink2 --pfile v2_imputed_clean --keep keep_135.txt (all 135 clean samples; the analysis then uses 129 / 133) --score <weights> 1 2 3 header
#          cols=+scoresums list-variants --read-freq stats129.afreq
#     columns: 1 = variant ID (chr:pos:ref:alt), 2 = scored allele (the effect allele in the clean file's alleles), 3 = weight; plink2 uses the DOSAGES of the pgen
#     (DS field of the imputation); SCORE1_SUM = sum(weight x dosage of the scored allele) is the score used; imputed dosages have no missing values, but
#     the allele frequencies of the 129 QC samples (--read-freq) would replace any missing genotype by plink2's mean imputation; ALLELE_CT in the .sscore shows
#     whether every sample contributed all alleles. No centring/standardisation inside plink2 (z-scores are made in step 52).
# Output: prscs/weights_imputed/<SET>_weights.txt(+.report.txt), prscs/scores_imputed/<SET>.sscore(.vars,.log)
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; C="$PRS/imputation/clean"; SC="${CODE_DIR}"
W="$PRS/prscs/weights_imputed"; S="$PRS/prscs/scores_imputed"; mkdir -p "$W" "$S"
SETS="${SETS:-V1 V2 V3 V1_r2ge0.8 V1_chr22 GRS29}"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
T="${SLURM_CPUS_PER_TASK:-1}"; MEM=$((SLURM_MEM_PER_NODE - 4000))
for s in $SETS; do
  echo "######## $s"
  case "$s" in
    V1)         python -I "$SC/05_PRS/17_build_weights_imputed.py" --prefix "$PRS/prscs/out_imputed/V1/V1" --map "$PRS/prscs/inputs_imputed/map_ukbb.tsv.gz" --out "$W/V1_weights.txt" --label V1 ;;
    V2)         python -I "$SC/05_PRS/17_build_weights_imputed.py" --prefix "$PRS/prscs/out_imputed/V2/V2" --map "$PRS/prscs/inputs_imputed/map_1kg.tsv.gz"  --out "$W/V2_weights.txt" --label V2 ;;
    V3)         python -I "$SC/05_PRS/17_build_weights_imputed.py" --prefix "$PRS/prscs/out/V3/V3" --map "$PRS/prscs/inputs_imputed/map_ukbb.tsv.gz" --out "$W/V3_weights.txt" --label V3 ;;
    V1_r2ge0.8) python -I "$SC/05_PRS/17_build_weights_imputed.py" --prefix "$PRS/prscs/out_imputed/V1/V1" --map "$PRS/prscs/inputs_imputed/map_ukbb.tsv.gz" --out "$W/V1_r2ge0.8_weights.txt" --min-r2 0.8 --label V1_r2ge0.8 ;;
    V1_chr22)   awk -F'\t' 'NR==1 || $1 ~ /^22:/' "$W/V1_weights.txt" > "$W/V1_chr22_weights.txt"; echo "chr22 weights: $(($(wc -l < "$W/V1_chr22_weights.txt") - 1))" ;;
    GRS29)      cp "$W/grs29_weights.txt" "$W/GRS29_weights.txt" ;;
    *) echo "unknown set $s"; exit 1 ;;
  esac
  plink2 --pfile "$C/v2_imputed_clean" --keep "$C/keep_135.txt" --read-freq "$C/stats129.afreq" \
         --score "$W/${s}_weights.txt" 1 2 3 header cols=+scoresums list-variants \
         --threads "$T" --memory "$MEM" --out "$S/$s"
  n_w=$(( $(wc -l < "$W/${s}_weights.txt") - 1 )); n_s=$(wc -l < "$S/$s.sscore.vars"); n_ind=$(( $(wc -l < "$S/$s.sscore") - 1 ))
  echo "$s: weights in file $n_w | variants scored (sscore.vars) $n_s | samples scored $n_ind"
  [ "$n_w" -eq "$n_s" ] || { echo "WARNING: number of scored variants differs from the weight file"; }
  [ "$n_ind" -eq 135 ] || { echo "ERROR: samples scored != 135"; exit 1; }
  awk 'NR>1 {ct[$0 ~ /x/ ? 0 : $3]++} END {for (k in ct) print "  ALLELE_CT value " k ": " ct[k] " samples"}' "$S/$s.sscore"
done
echo "done"
