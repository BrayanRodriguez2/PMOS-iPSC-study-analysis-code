#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/04_snp_qc_dedup_rename.sbatch   (original step script: 10_qc_variants.sbatch)
# Produces: SNP QC, duplicate handling and renaming to rsIDs; QC PLINK datasets.
# Run: sbatch 01_genotype_QC/04_snp_qc_dedup_rename.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_qc_variants
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 10 - VARIANT QC, exclusions, de-duplication and rsID renaming.
# Input : genotypes/plink_bin/b38 (raw binary copy) + the sample keep-list of step 09.
# Filters applied SEQUENTIALLY (so the effect of each one is counted):
#   S1 keep samples that pass step 09
#   S2 variant call rate >= 98 %           (--geno 0.02)
#   S3 remove monomorphic variants         (--mac 1)
#   S4 Hardy-Weinberg exact test p < 1e-6  (--hwe 1e-6 0; k=0 is required by this plink2)
#   S5 MAF >= 1 %                          (--maf 0.01; extra filter for the PRS)
#   S6 autosomes only (1-22)               -> genotypes/plink_qc/b38_qc  (original IDs)
#   then 01_genotype_QC/04_snp_qc_dedup_rename.py: A palindromic recovered rsIDs out, B 508 incompatible out,
#   C duplicates by rsID, D duplicates by chr:pos+alleles (best call rate kept)
#   and rename recovered variants to their rsID -> genotypes/plink_qc/b38_qc_rsid
# Output: results/qc/variants/qc_counts.tsv (samples and variants after each filter).
# HWE is applied to the whole cohort (case/control status of the 140 donors is not known here);
# with n ~ 138 the 1e-6 threshold is conservative.
# Intermediate .pgen files (genotypes/qc_work) are own temporary files and are removed at the end.
# Resources: small data; 4 CPUs / 16 GB on `testing`.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail

PRS="${PROJECT_DIR}/PRS"
RAW="$PRS/genotypes/plink_bin/b38"
KEEP="$PRS/validation/results/qc/samples/keep_samples.txt"
Q="$PRS/genotypes/plink_qc"
W="$PRS/genotypes/qc_work"
RES="$PRS/validation/results/qc/variants"
mkdir -p "$Q" "$W" "$RES"

eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prs
T="${SLURM_CPUS_PER_TASK:-1}"
COUNTS="$RES/qc_counts.tsv"
printf "step\tsamples\tvariants\tremoved_in_step\n" > "$COUNTS"
PREV=""

add_count() {  # label n_samples n_variants
  local removed="NA"
  if [ -n "$PREV" ]; then removed=$((PREV - $3)); fi
  printf "%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$removed" >> "$COUNTS"
  PREV="$3"
}
count_pgen() { add_count "$1" "$(grep -vc '^#' "$2.psam")" "$(grep -vc '^#' "$2.pvar")"; }

add_count "S0 raw (b38, as converted)" "$(wc -l < "$RAW.fam")" "$(wc -l < "$RAW.bim")"

plink2 --bfile "$RAW" --keep "$KEEP" --make-pgen --threads "$T" --out "$W/s1"
count_pgen "S1 samples passing QC (step 09)" "$W/s1"
plink2 --pfile "$W/s1" --geno 0.02 --make-pgen --threads "$T" --out "$W/s2"
count_pgen "S2 variant call rate >= 98%" "$W/s2"
plink2 --pfile "$W/s2" --mac 1 --make-pgen --threads "$T" --out "$W/s3"
count_pgen "S3 remove monomorphic" "$W/s3"
plink2 --pfile "$W/s3" --hwe 1e-6 0 --make-pgen --threads "$T" --out "$W/s4"
count_pgen "S4 HWE p >= 1e-6" "$W/s4"
plink2 --pfile "$W/s4" --maf 0.01 --make-pgen --threads "$T" --out "$W/s5"
count_pgen "S5 MAF >= 0.01" "$W/s5"
plink2 --pfile "$W/s5" --chr 1-22 --make-bed --threads "$T" --out "$Q/b38_qc"
add_count "S6 autosomes only = b38_qc" "$(wc -l < "$Q/b38_qc.fam")" "$(wc -l < "$Q/b38_qc.bim")"

echo "### per-variant call rate of b38_qc (used to choose among duplicates)"
plink2 --bfile "$Q/b38_qc" --missing variant-only --threads "$T" --out "$RES/b38_qc"

echo "### exclusions, duplicates and rsID renaming"
python -I "${CODE_DIR}/01_genotype_QC/04_snp_qc_dedup_rename.py" \
  --bim "$Q/b38_qc.bim" \
  --vmiss "$RES/b38_qc.vmiss" \
  --recovered "$PRS/validation/results/overlap/recovered_rsids.tsv" \
  --incompatible "$PRS/validation/results/overlap/incompatible_alleles.tsv" \
  --outdir "$RES"

N0=$(wc -l < "$Q/b38_qc.bim")
for tag in A_after_palindromic B_after_incompatible C_after_dedup_rsid D_final; do
  n=$(wc -l < "$RES/keep_$tag.txt")
  add_count "10.$tag" "$(wc -l < "$Q/b38_qc.fam")" "$n"
done

plink2 --bfile "$Q/b38_qc" --extract "$RES/keep_D_final.txt" --make-pgen --threads "$T" --out "$W/s7"
# --update-name file: column 1 = old ID, column 2 = new ID
plink2 --pfile "$W/s7" --update-name "$RES/rename_map.tsv" --make-bed --threads "$T" --out "$Q/b38_qc_rsid"
add_count "FINAL b38_qc_rsid (renamed copy)" "$(wc -l < "$Q/b38_qc_rsid.fam")" "$(wc -l < "$Q/b38_qc_rsid.bim")"

# Same final variant set but ALL 140 samples (also those failing sample QC): used only to
# report, clearly flagged, the PRS of donors that fail QC (e.g. 2 of the 14 iPSC donors).
plink2 --bfile "$RAW" --extract "$RES/keep_D_final.txt" --make-pgen --threads "$T" --out "$W/s8"
plink2 --pfile "$W/s8" --update-name "$RES/rename_map.tsv" --make-bed --threads "$T" --out "$Q/b38_qc_rsid_all140"

echo "### sanity: number of rs<digits> IDs in the final .bim and duplicated IDs"
grep -c -E '	rs[0-9]+	' "$Q/b38_qc_rsid.bim" || true
cut -f2 "$Q/b38_qc_rsid.bim" | sort | uniq -d | wc -l

echo "### counts table"; cat "$COUNTS"
rm -f "$W"/s*.pgen "$W"/s*.pvar "$W"/s*.psam
