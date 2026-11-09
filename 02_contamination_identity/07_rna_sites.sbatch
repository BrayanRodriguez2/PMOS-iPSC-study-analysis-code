#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/07_rna_sites.sbatch   (original step script: 27a_rna_sites.sbatch)
# Produces: Array SNPs that can be genotyped from RNA-seq (inside exons).
# Run: sbatch 02_contamination_identity/07_rna_sites.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=rna_sites
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 27a - Part B of the contamination diagnostic: list of array SNPs that can be genotyped from RNA-seq.
# The SNPs are those of the TOPMed-prepared array VCF (step 20: 510,607 autosomal biallelic SNPs, GRCh38, forward strand,
# REF = GRCh38 base, 140 donors) that lie inside a GENCODE v50 exon (RNA reads can only cover transcribed sequence).
# The depth filter (>= 10 reads) is applied later, per iPSC line.
# Output: PRS/rnaseq/sites/{exons.merged.tsv, array_exon_sites.regions.tsv.gz (CHROM POS),
#         array_exon_sites.alleles.tsv.gz (CHROM POS REF,ALT; bcftools call -C alleles format)} (+ .tbi)
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
VCF="$PRS/genotypes/topmed/aligned/all_chr.vcf.gz"
GTF="$PRS/reference/star_gencode_v50/gencode.v50.primary_assembly.annotation.gtf"
S="$PRS/rnaseq/sites"; mkdir -p "$S"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate rnaseq

echo "### merged exon intervals (all transcripts, chromosomes 1-22)"
awk -F'\t' '$3=="exon" && $1 ~ /^chr([0-9]+)$/ {print $1"\t"$4"\t"$5}' "$GTF" \
  | sort -k1,1V -k2,2n \
  | awk -F'\t' 'BEGIN{OFS="\t"} {if ($1==c && $2<=e+1) {if ($3>e) e=$3} else {if (c!="") print c,s,e; c=$1; s=$2; e=$3}} END{print c,s,e}' > "$S/exons.merged.tsv"
wc -l "$S/exons.merged.tsv"

[ -f "$VCF.tbi" ] || [ -f "$VCF.csi" ] || bcftools index -t "$VCF"
echo "### array SNPs inside exons"
bcftools view -G -R "$S/exons.merged.tsv" "$VCF" -Ou \
  | bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' | sort -u -k1,1V -k2,2n > "$S/sites.alleles.tsv"
wc -l "$S/sites.alleles.tsv"
cut -f1,2 "$S/sites.alleles.tsv" | bgzip -c > "$S/array_exon_sites.regions.tsv.gz"
bgzip -c "$S/sites.alleles.tsv" > "$S/array_exon_sites.alleles.tsv.gz"
tabix -f -s1 -b2 -e2 "$S/array_exon_sites.regions.tsv.gz"
tabix -f -s1 -b2 -e2 "$S/array_exon_sites.alleles.tsv.gz"
rm -f "$S/sites.alleles.tsv"
ls -la "$S"
