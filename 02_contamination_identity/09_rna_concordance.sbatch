#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/09_rna_concordance.sbatch   (original step script: 28_rna_concordance.sbatch)
# Produces: Line x donor genotype concordance (identity check).
# Run: sbatch 02_contamination_identity/09_rna_concordance.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=rna_conc
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 28 - Part B of the contamination diagnostic: 14 lines x 140 donors genotype concordance (see 02_contamination_identity/09_rna_concordance.py).
# This wrapper extracts the donors' array genotypes at the exonic array SNPs (bcftools, from the TOPMed-prepared VCF of
# step 20: 140 donors, no sample QC) and runs the python analysis.
# Input : PRS/rnaseq/{rna_samples.tsv, genotypes/*.calls.tsv.gz, sites/array_exon_sites.regions.tsv.gz}, array VCF
# Output: PRS/validation/results/rna_identity/
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
VCF="$PRS/genotypes/topmed/aligned/all_chr.vcf.gz"
R="$PRS/rnaseq"; OUT="$PRS/validation/results/rna_identity"; mkdir -p "$OUT/work"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate rnaseq
bcftools query -l "$VCF" > "$OUT/work/donor_samples.txt"
echo "donors: $(wc -l < "$OUT/work/donor_samples.txt")"
bcftools view -T "$R/sites/array_exon_sites.regions.tsv.gz" -Ou "$VCF" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > "$OUT/work/donor_gt.tsv"
echo "exonic array SNPs: $(wc -l < "$OUT/work/donor_gt.tsv")"
python -I "${CODE_DIR}/02_contamination_identity/09_rna_concordance.py" \
  --rna-samples "$R/rna_samples.tsv" --calls-dir "$R/genotypes" \
  --donor-gt "$OUT/work/donor_gt.tsv" --donor-samples "$OUT/work/donor_samples.txt" \
  --qc-table "$PRS/validation/results/qc/samples/sample_qc_table.tsv" \
  --outdir "$OUT" --min-dp 10 --focal-donors "${FOCAL_SAMPLE_IDS:-}"
