#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/07_freq_vs_1000g.sbatch   (original step script: 42_freq_vs_1000g.sbatch)
# Produces: Imputed allele frequencies against 1000 Genomes EUR/FIN in sampled windows.
# Run: sbatch 04_post_imputation_QC/07_freq_vs_1000g.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=freq_1kg
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 42 - Post-imputation QC (d): imputed frequency (129 QC samples) vs 1000G EUR / FIN in common variants, in 66 windows of 300 kb
# (see 04_post_imputation_QC/07_freq_vs_1000g.py for the design and its limitation). Network job: the 1000G genotypes of the windows are read remotely
# (bcftools + .tbi, only the windows) from
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ (contigs "1".."22").
# Input : unzipped/chrN.dose.vcf.gz (+ .tbi, step 39), info_tables (step 40), results/qc/samples/keep_samples.txt (129), step-34 files
#         (eur_samples.txt, g1000_samples.txt, 1kg.panel)
# Output: PRS/validation/results/imputation_qc/{freq_vs_1000g_*.tsv, freq_vs_1000g_report.txt, figures/F4_*} ; work files in imputation_qc/work_freq
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; D="$PRS/imputation/topmed_imputation_v2"; O="$PRS/validation/results/imputation_qc"; W="$O/work_freq"
Q="$PRS/validation/results/topmed_qc_only"; G="$Q/g1000_work"
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"
mkdir -p "$W"; cd "$W"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
S="${CODE_DIR}/04_post_imputation_QC/07_freq_vs_1000g.py"
python -I "$S" windows --tables-dir "$D/info_tables" --out "$W/windows.tsv"
awk 'NR>1 {print $2}' "$PRS/validation/results/qc/samples/keep_samples.txt" > "$W/samples129.txt"
echo "samples in the dose extraction: $(wc -l < "$W/samples129.txt")"
: > "$W/dose_windows.tsv"; : > "$W/g1000_windows_gt.tsv"
for c in $(seq 1 22); do
  awk -F'\t' -v c="chr$c" '$1==c' "$W/windows.tsv" > "$W/win_chr$c.tsv"
  bcftools view -S "$W/samples129.txt" -R "$W/win_chr$c.tsv" -Ou "$D/unzipped/chr$c.dose.vcf.gz" \
    | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/R2\t%INFO/MAF\t%INFO/IMPUTED\t%INFO/TYPED[\t%DS]\n' >> "$W/dose_windows.tsv"
  sed 's/^chr//' "$W/win_chr$c.tsv" > "$W/win_g1000_chr$c.tsv"
  url="$BASE/ALL.chr${c}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  bcftools view -S "$G/eur_samples.txt" --force-samples -R "$W/win_g1000_chr$c.tsv" -m2 -M2 -v snps "$url" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' >> "$W/g1000_windows_gt.tsv"
  echo "[$(date +%T)] chr$c: dose rows $(wc -l < "$W/dose_windows.tsv") | 1000G rows $(wc -l < "$W/g1000_windows_gt.tsv")"
done
python -I "$S" compare --dose "$W/dose_windows.tsv" --g1000 "$W/g1000_windows_gt.tsv" --g1000-samples "$Q/g1000_samples.txt" --panel "$G/1kg.panel" --outdir "$O"
