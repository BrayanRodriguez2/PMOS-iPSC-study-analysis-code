#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/08_frequency_vs_1000g.sbatch   (original step script: 34_1000g_compare.sbatch)
# Produces: Allele frequencies of discordant SNPs against 1000 Genomes EUR/FIN.
# Run: sbatch 03_imputation_TOPMed/08_frequency_vs_1000g.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=g1000_cmp
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 34 - Compare the 446 mismatched-frequency SNPs, the 182 allele-mismatch SNPs and the 29 lead SNPs with 1000 Genomes
# phase 3 (GRCh38): EUR (503 samples, includes FIN) and FIN (99). See 03_imputation_TOPMed/08_frequency_vs_1000g.py.
# Remote reading: only the target positions (~650) are read from the per-chromosome 1000G VCFs with bcftools + their .tbi index
# (https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/); no
# whole-chromosome file is downloaded. 1000G contigs are named "1", "2", ... (no "chr" prefix).
# Also uses the local 1000G sites file (INFO/EUR_AF, 2 decimals) = exactly what the old step-20 filter used, to explain why it missed the 88.
# Output: PRS/validation/results/topmed_qc_only/{compare_1000g_446.tsv, compare_1000g_182.tsv, compare_1000g_report.txt, g1000_targets.tsv, g1000_gt.tsv}
# Network job: `main`, 2 CPUs / 8 GB.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; OUT="$PRS/validation/results/topmed_qc_only"; W="$OUT/g1000_work"
LD="$PRS/validation/results/loci"
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"
mkdir -p "$W"; cd "$W"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
cp "$LD/ld_work/eur_samples.txt" "$W/eur_samples.txt"        # 503 EUR samples (includes the 99 FIN), from step 12
cp "$LD/ld_work/1kg.panel" "$W/1kg.panel"                     # sample / pop / super_pop / gender
python -I "${CODE_DIR}/03_imputation_TOPMed/08_frequency_vs_1000g.py" --outdir "$OUT" --make-targets "$OUT/g1000_targets.tsv" --leads "$LD/lead_snps_lookup.tsv"

echo "### local sites file (old filter's view): EUR_AF at the target positions"
bcftools query -R "$OUT/g1000_targets.tsv" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/EUR_AF\n' "$PRS/reference/grch38/1000G_GRCh38.sites.vcf.gz" > "$OUT/sites_eur_af.tsv"
wc -l "$OUT/sites_eur_af.tsv"

echo "### remote genotypes (EUR 503 samples) at the target positions"
URL1="$BASE/ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
bcftools view -h -S eur_samples.txt --force-samples "$URL1" | tail -n 1 | cut -f10- | tr '\t' '\n' > "$OUT/g1000_samples.txt"
echo "samples in the extraction: $(wc -l < "$OUT/g1000_samples.txt")"
: > "$OUT/g1000_gt.tsv"
for c in $(seq 1 22); do
  awk -F'\t' -v c="$c" '$1==c' "$OUT/g1000_targets.tsv" > "$W/targets_chr$c.tsv"
  [ -s "$W/targets_chr$c.tsv" ] || continue
  url="$BASE/ALL.chr${c}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  bcftools view -S eur_samples.txt --force-samples -R "$W/targets_chr$c.tsv" -m2 -M2 -v snps "$url" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' >> "$OUT/g1000_gt.tsv"
  echo "chr$c: $(wc -l < "$W/targets_chr$c.tsv") target positions, records so far: $(wc -l < "$OUT/g1000_gt.tsv")"
done
python -I "${CODE_DIR}/03_imputation_TOPMed/08_frequency_vs_1000g.py" --outdir "$OUT" --gt "$OUT/g1000_gt.tsv" --samples "$OUT/g1000_samples.txt" \
  --panel "$W/1kg.panel" --sites-eur-af "$OUT/sites_eur_af.tsv"
