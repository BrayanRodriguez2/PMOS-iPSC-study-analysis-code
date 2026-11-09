#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/08_rna_genotypes.sbatch   (original step script: 27b_rna_genotypes.sbatch)
# Produces: Genotypes of each iPSC line called from its RNA-seq BAM at those SNPs.
# Run: sbatch 02_contamination_identity/08_rna_genotypes.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=rna_gt
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-13%7
#SBATCH --output=logs/%x_%A_%a.out
# ---------------------------------------------------------------------------
# Step 27b - Part B of the contamination diagnostic: genotype each iPSC line from its RNA-seq BAM at the array SNP positions
# (step 27a), with `bcftools mpileup | bcftools call -m -C alleles` restricted to those positions and to the array's
# REF/ALT alleles (so the call is "which of the two known alleles are in the reads").
#   mpileup: -q 20 (mapping quality; STAR gives 255 to unique reads, <= 3 to multi-mappers) -Q 20 (base quality),
#            -I (no indel calling), -a FORMAT/AD,FORMAT/DP, --max-depth 2000; reference = the STAR genome FASTA.
# Output: PRS/rnaseq/genotypes/<fastq_name>.calls.tsv.gz  with CHROM POS REF ALT GT DP AD  (one row per covered SNP;
#         the >= 10 reads filter is applied in step 28, where each line is compared with the 140 donors).
# Array task i = i-th line of PRS/rnaseq/rna_samples.tsv (same order as step 26).
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
FA="$PRS/reference/star_gencode_v50/GRCh38.primary_assembly.genome.fa"
S="$PRS/rnaseq/sites"
OUT="$PRS/rnaseq/genotypes"; mkdir -p "$OUT"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate rnaseq
NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$PRS/rnaseq/rna_samples.tsv" | cut -f1)   # +2: header line
BAM="$PRS/rnaseq/star/$NAME/$NAME.bam"
echo "task $SLURM_ARRAY_TASK_ID: $NAME  $BAM"; bcftools --version | sed -n '1p'
bcftools mpileup -f "$FA" -T "$S/array_exon_sites.regions.tsv.gz" -q 20 -Q 20 -I --max-depth 2000 \
    -a FORMAT/AD,FORMAT/DP -Ou "$BAM" \
  | bcftools call -m -C alleles -T "$S/array_exon_sites.alleles.tsv.gz" -Ou \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%AD]\n' \
  | bgzip -c > "$OUT/$NAME.calls.tsv.gz"
echo "rows: $(zcat "$OUT/$NAME.calls.tsv.gz" | wc -l)"
echo "rows with DP>=10: $(zcat "$OUT/$NAME.calls.tsv.gz" | awk -F'\t' '$6>=10' | wc -l)"
