#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/01_setup_references.sbatch   (original step script: 18_topmed_setup_refs.sbatch)
# Produces: GRCh38 reference FASTA and reference-panel resources for the pre-imputation checks.
# Run: sbatch 03_imputation_TOPMed/01_setup_references.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=tm_refs
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 18 - TOPMed Imputation Server preparation: environment and reference files.
#
# (1) conda env "topmed" (conda-forge + bioconda only, strict priority): python 3.11, bcftools (with the
#     +fixref plugin; bgzip/tabix come with it) and plink2.
# (2) GRCh38 reference FASTA with "chr" contig names (UCSC hg38.fa.gz, ~0.98 GB, md5 checked against the
#     file published by UCSC) -> uncompressed hg38.fa, needed by bcftools +fixref (the .fai index is
#     built automatically on first use).
# (3) 1000 Genomes phase 3 GRCh38 sites-only VCF (EBI, 932 MB) with INFO/EUR_AF: used as the ONLY
#     programmatically available stand-in for reference allele frequencies (TOPMed/Bravo frequencies and
#     Rayner's TOPMed site list are not downloadable by script: HTTP 403). It is used to resolve the strand
#     of palindromic SNPs and to flag allele-frequency discordance. It is NOT the TOPMed panel.
# Output: PRS/reference/grch38/{hg38.fa, hg38.fa.gz, md5sum.txt, 1000G_GRCh38.sites.vcf.gz(.tbi)}
# Resources: network + gunzip, 2 CPUs / 8 GB; `main` (may exceed 2 h on a slow link). wget -c resumes.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
REF="$PRS/reference/grch38"
mkdir -p "$REF" "$PRS/env"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx "topmed"; then
  conda create -y -n topmed --override-channels --strict-channel-priority -c conda-forge -c bioconda \
    python=3.11 bcftools plink2
fi
conda env export -n topmed --no-builds > "$PRS/env/environment_topmed.yml"   # before activating (see step 01)
conda activate topmed
bcftools --version | head -n 1; plink2 --version; bgzip --version | head -n 1

cd "$REF"
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips
wget -c -nv -O hg38.fa.gz "$UCSC/hg38.fa.gz"
wget -nv -O md5sum.txt "$UCSC/md5sum.txt"
echo "### md5 of hg38.fa.gz vs UCSC"
grep ' hg38.fa.gz$' md5sum.txt
md5sum hg38.fa.gz
[ "$(grep ' hg38.fa.gz$' md5sum.txt | cut -d' ' -f1)" = "$(md5sum hg38.fa.gz | cut -d' ' -f1)" ] && echo "md5 OK" || { echo "md5 MISMATCH"; exit 1; }
[ -s hg38.fa ] || gzip -dc hg38.fa.gz > hg38.fa
ls -la hg38.fa
grep -c '^>' hg38.fa | sed 's/^/contigs in FASTA: /'
grep '^>chr[0-9XYM]*$' hg38.fa | head -n 30 | tr '\n' ' '; echo

B1KG=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL
S=ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz
wget -c -nv -O 1000G_GRCh38.sites.vcf.gz "$B1KG/$S"
wget -c -nv -O 1000G_GRCh38.sites.vcf.gz.tbi "$B1KG/$S.tbi"
md5sum 1000G_GRCh38.sites.vcf.gz | tee 1000G_GRCh38.sites.vcf.gz.md5
# (the two previews below end in `head`: SIGPIPE under `pipefail` must not fail the job)
bcftools view -h 1000G_GRCh38.sites.vcf.gz | grep -E '^##INFO=<ID=(AF|EUR_AF)' || true
bcftools view -H 1000G_GRCh38.sites.vcf.gz | head -n 2 | cut -c1-200 || true
