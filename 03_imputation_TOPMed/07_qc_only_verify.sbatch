#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/07_qc_only_verify.sbatch   (original step script: 33_qc_only_verify.sbatch)
# Produces: Independent verification of the server QC-only statistics.
# Run: sbatch 03_imputation_TOPMed/07_qc_only_verify.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=qc_verify
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 33 - Independent verification of the TOPMed QC-only job (see 03_imputation_TOPMed/07_qc_only_verify.py).
# Input : PRS/imputation/topmed_qc_only/{quality-control.html, snps-excluded.txt} (step 32, data only),
#         PRS/topmed/upload/chr*.vcf.gz (the v1 files that were uploaded; 140 samples), sample QC table (step 09)
# Output: PRS/validation/results/topmed_qc_only/{freq_mismatch_446_verified.tsv, allele_mismatch_182_verified.tsv,
#         freq_mismatch_446_carriers_per_sample.tsv, verify_report.txt} (+ work files)
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
Q="$PRS/imputation/topmed_qc_only"; UP="$PRS/topmed/upload"; OUT="$PRS/validation/results/topmed_qc_only"
mkdir -p "$OUT"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat      # bcftools 1.22 + python (numpy, pandas, scipy)
echo "### v1 site table (ID CHROM POS REF ALT) of the 22 uploaded files"
: > "$OUT/v1_sites.tsv"
for c in $(seq 1 22); do bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' "$UP/chr$c.vcf.gz" >> "$OUT/v1_sites.tsv"; done
wc -l "$OUT/v1_sites.tsv"
python -I "${CODE_DIR}/03_imputation_TOPMed/07_qc_only_verify.py" parse --html "$Q/quality-control.html" --snps-excluded "$Q/snps-excluded.txt" \
  --v1-sites "$OUT/v1_sites.tsv" --outdir "$OUT"
echo "### genotypes of the target positions (140 samples)"
bcftools query -l "$UP/chr1.vcf.gz" > "$OUT/work_samples.txt"
: > "$OUT/work_gt.tsv"
for c in $(seq 1 22); do
  bcftools view -T "$OUT/work_targets.tsv" -Ou "$UP/chr$c.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' >> "$OUT/work_gt.tsv"
done
wc -l "$OUT/work_gt.tsv"
python -I "${CODE_DIR}/03_imputation_TOPMed/07_qc_only_verify.py" analyze --outdir "$OUT" --gt "$OUT/work_gt.tsv" --samples "$OUT/work_samples.txt" \
  --qc-table "$PRS/validation/results/qc/samples/sample_qc_table.tsv"
