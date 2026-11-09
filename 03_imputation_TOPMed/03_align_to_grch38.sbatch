#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/03_align_to_grch38.sbatch   (original step script: 20_topmed_align.sbatch)
# Produces: Strand/allele alignment to GRCh38 (bcftools +fixref, frequency-based strand resolution, frequency filter).
# Run: sbatch 03_imputation_TOPMed/03_align_to_grch38.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=tm_align
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 20 - Strand and allele alignment of the SNP-QC dataset to GRCh38 for the TOPMed Imputation Server.
#
# Recommended tool (server docs): W. Rayner's HRC/1000G/TOPMed pre-imputation check. NOT VIABLE here: its
# download site answers HTTP 403 to automated clients (and the page returned 404 in a browser on 2026-11-02),
# and I do not circumvent that. Fallback requested: bcftools +fixref with the GRCh38 FASTA.
#
# Order of operations (each loss is counted in results/topmed/alignment_counts.tsv):
#   A  palindromic A/T and C/G SNPs: strand resolved by allele frequency vs 1000G EUR (03_imputation_TOPMed/03_align_to_grch38_resolve_strand.py);
#      removed if study MAF > 0.4, not in the reference, or discordant frequency (> 0.2); the rest are
#      complemented with plink2 --flip when the array is on the reverse strand.
#   B  export VCF 4.2 with "chr" prefix (required by the server for GRCh38) and run bcftools +fixref
#      -m flip-all -d with hg38.fa: swaps/flips REF/ALT (and recodes GT) so that REF = GRCh38 base; sites whose
#      alleles match the FASTA neither directly nor after flip/swap are discarded (-d).
#      flip-all is safe for the ambiguous sites because their strand was already resolved in A: fixref then only
#      swaps REF/ALT when needed.
#   C  allele-frequency discordance |study - 1000G EUR| > 0.2 (03_imputation_TOPMed/03_align_to_grch38_freq_filter.py): removed.
# Output: genotypes/topmed/aligned/all_chr.vcf.gz (+ .csi), all 140 samples; work files in aligned/work.
# 1000G EUR is a SUBSTITUTE for TOPMed frequencies (not downloadable by script).
# Resources: 450 k variants x 140 samples; 4 CPUs / 16 GB on `testing`.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
IN="$PRS/genotypes/topmed/b38_snpqc"
AL="$PRS/genotypes/topmed/aligned"
W="$AL/work"
RES="$PRS/validation/results/topmed"
FA="$PRS/reference/grch38/hg38.fa"
SITES="$PRS/reference/grch38/1000G_GRCh38.sites.vcf.gz"
mkdir -p "$W" "$RES"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate topmed
T="${SLURM_CPUS_PER_TASK:-1}"
CNT="$RES/alignment_counts.tsv"
printf "step\tvariants\tremoved_in_step\n" > "$CNT"
PREV=""
add() {  # label n
  local rem="NA"; if [ -n "$PREV" ]; then rem=$((PREV - $2)); fi
  printf "%s\t%s\t%s\n" "$1" "$2" "$rem" >> "$CNT"; PREV="$2"
}
nvar() { bcftools view -H "$1" | wc -l; }

add "A0 SNP-QC dataset (b38_snpqc)" "$(wc -l < "$IN.bim")"

echo "### study frequencies + reference frequencies at our positions"
plink2 --bfile "$IN" --freq --threads "$T" --out "$W/study"
awk '{print $1"\t"$4}' "$IN.bim" | sort -k1,1n -k2,2n -u > "$W/positions.tsv"
bcftools query -R "$W/positions.tsv" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/EUR_AF\n' "$SITES" > "$W/ref_sites.tsv"
wc -l "$W/positions.tsv" "$W/ref_sites.tsv"

echo "### A: palindromic SNPs"
python -I "${CODE_DIR}/03_imputation_TOPMed/03_align_to_grch38_resolve_strand.py" --bim "$IN.bim" --afreq "$W/study.afreq" \
       --ref-sites "$W/ref_sites.tsv" --outdir "$RES"
FLAGS=""
[ -s "$RES/palindromic_remove.txt" ] && FLAGS="$FLAGS --exclude $RES/palindromic_remove.txt"
[ -s "$RES/palindromic_flip.txt" ] && FLAGS="$FLAGS --flip $RES/palindromic_flip.txt"
plink2 --bfile "$IN" $FLAGS --make-bed --threads "$T" --out "$W/pal"
add "A. palindromic SNPs removed (MAF>0.4 / not in reference / discordant)" "$(wc -l < "$W/pal.bim")"

echo "### B: VCF 4.2 with chr prefix + bcftools +fixref"
plink2 --bfile "$W/pal" --export vcf-4.2 bgz id-paste=iid --output-chr chrM --threads "$T" --out "$W/pal"
bcftools index -f -t "$W/pal.vcf.gz"
echo "--- fixref stats BEFORE (strand/REF check against hg38.fa)"
bcftools +fixref "$W/pal.vcf.gz" -- -f "$FA" 2>&1 | tail -n 14
bcftools +fixref "$W/pal.vcf.gz" -Oz -o "$W/fixed.vcf.gz" -- -f "$FA" -m flip-all -d 2>&1 | tail -n 5
bcftools index -f -t "$W/fixed.vcf.gz"
echo "--- INFO/FIXREF tags (what was changed)"
bcftools query -f '%INFO/FIXREF\n' "$W/fixed.vcf.gz" | sort | uniq -c
echo "--- fixref stats AFTER"
bcftools +fixref "$W/fixed.vcf.gz" -- -f "$FA" 2>&1 | tail -n 14
add "B. REF/ALT not matching hg38 after flip/swap discarded (fixref -d)" "$(nvar "$W/fixed.vcf.gz")"

echo "### C: allele-frequency discordance vs 1000G EUR"
bcftools +fill-tags "$W/fixed.vcf.gz" -Ou -- -t AF | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' > "$W/fixed_af.tsv"
python -I "${CODE_DIR}/03_imputation_TOPMed/03_align_to_grch38_freq_filter.py" --study "$W/fixed_af.tsv" --ref-sites "$W/ref_sites.tsv" --outdir "$RES"
if [ -s "$RES/freq_remove_sites.tsv" ]; then
  bcftools view "$W/fixed.vcf.gz" -T "^$RES/freq_remove_sites.tsv" -Ou | bcftools sort -Ou | bcftools annotate -x INFO -Oz -o "$AL/all_chr.vcf.gz"
else
  bcftools sort "$W/fixed.vcf.gz" -Ou | bcftools annotate -x INFO -Oz -o "$AL/all_chr.vcf.gz"
fi
bcftools index -f "$AL/all_chr.vcf.gz"
add "C. allele-frequency discordance vs 1000G EUR (|diff|>0.2) removed" "$(nvar "$AL/all_chr.vcf.gz")"

echo "### counts table"; cat "$CNT"
echo "samples in the VCF: $(bcftools query -l "$AL/all_chr.vcf.gz" | wc -l)"
bcftools view -h "$AL/all_chr.vcf.gz" | sed -n '1,2p'   # sed, not head: no SIGPIPE under pipefail
