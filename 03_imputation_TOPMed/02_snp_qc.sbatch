#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/02_snp_qc.sbatch   (original step script: 19_topmed_snp_qc.sbatch)
# Produces: SNP QC for imputation (call rate, monomorphic, duplicated positions).
# Run: sbatch 03_imputation_TOPMed/02_snp_qc.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=tm_snpqc
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 19 - SNP-only QC of the RAW array for the TOPMed Imputation Server. ALL 140 samples are kept:
# no sample QC (exclusions will be decided after imputation).
# Input (read-only): genotypes/plink_bin/b38 (raw binary copy of the donor genotypes, GRCh38)
# Filters, applied SEQUENTIALLY so that each loss is counted (results/topmed/snp_qc_counts.tsv):
#   T1 autosomes 1-22 only (drops X, XY/PAR, Y, MT and unplaced "0")
#   T2 SNPs only: biallelic A/C/G/T (--snps-only just-acgt; drops indels, missing alleles)
#   T3 variant call rate >= 0.95 in the 140 samples (--geno 0.05)
#   T4 no monomorphic variants (--mac 1)
#   T5 no duplicated positions (03_imputation_TOPMed/02_snp_qc.py: redundant copies -> keep best call rate;
#      conflicting allele sets -> drop all)
# Why these: the server's own QC removes invalid alleles, duplicates, indels, monomorphic sites, allele
# mismatches and SNP call rate < 90 % (docs: statgen.github.io/tis-v2-docs/imputation/getting-started/);
# doing it here makes the losses visible, and 0.95 is stricter than the server's 0.90. The documentation
# recommends no other pre-filters (no HWE, no MAF): none is applied.
# Output: genotypes/topmed/b38_snpqc.{bed,bim,fam}; intermediate .pgen files are removed at the end.
# Resources: 777 k x 140; 4 CPUs / 16 GB on `testing`.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
RAW="$PRS/genotypes/plink_bin/b38"
OUTD="$PRS/genotypes/topmed"
W="$OUTD/work"
RES="$PRS/validation/results/topmed"
mkdir -p "$OUTD" "$W" "$RES"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate topmed
T="${SLURM_CPUS_PER_TASK:-1}"
COUNTS="$RES/snp_qc_counts.tsv"
printf "step\tsamples\tvariants\tremoved_in_step\n" > "$COUNTS"
PREV=""
add_count() {  # label n_samples n_variants
  local removed="NA"
  if [ -n "$PREV" ]; then removed=$((PREV - $3)); fi
  printf "%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$removed" >> "$COUNTS"
  PREV="$3"
}
count_pgen() { add_count "$1" "$(grep -vc '^#' "$2.psam")" "$(grep -vc '^#' "$2.pvar")"; }

add_count "T0 raw array (b38)" "$(wc -l < "$RAW.fam")" "$(wc -l < "$RAW.bim")"
plink2 --bfile "$RAW" --chr 1-22 --make-pgen --threads "$T" --out "$W/t1"
count_pgen "T1 autosomes 1-22" "$W/t1"
plink2 --pfile "$W/t1" --snps-only just-acgt --make-pgen --threads "$T" --out "$W/t2"
count_pgen "T2 biallelic A/C/G/T SNPs only" "$W/t2"
plink2 --pfile "$W/t2" --geno 0.05 --make-pgen --threads "$T" --out "$W/t3"
count_pgen "T3 variant call rate >= 0.95 (140 samples)" "$W/t3"
plink2 --pfile "$W/t3" --mac 1 --make-pgen --threads "$T" --out "$W/t4"
count_pgen "T4 not monomorphic" "$W/t4"

plink2 --pfile "$W/t4" --missing variant-only --threads "$T" --out "$W/t4"
python -I "${CODE_DIR}/03_imputation_TOPMed/02_snp_qc.py" --pvar "$W/t4.pvar" --vmiss "$W/t4.vmiss" --out "$RES/dup_positions_remove.txt"
plink2 --pfile "$W/t4" --exclude "$RES/dup_positions_remove.txt" --make-bed --threads "$T" --out "$OUTD/b38_snpqc"
add_count "T5 no duplicated positions = b38_snpqc" "$(wc -l < "$OUTD/b38_snpqc.fam")" "$(wc -l < "$OUTD/b38_snpqc.bim")"

echo "### counts"; cat "$COUNTS"
echo "### chromosome content of b38_snpqc"; cut -f1 "$OUTD/b38_snpqc.bim" | sort -n | uniq -c | tr '\n' ' '; echo
rm -f "$W"/t[1-4].pgen "$W"/t[1-4].pvar "$W"/t[1-4].psam
