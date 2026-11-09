#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/05_ld_proxies.sbatch   (original step script: 12_ld_proxies.sbatch)
# Produces: LD proxies of the lead SNPs and locus coverage on the array.
# Run: sbatch 05_PRS/05_ld_proxies.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_ldproxy
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 12 - Part B.2: coverage of the 29 lead SNPs by the QC array, with LD proxies.
#
# LD reference: 1000 Genomes phase 3 GRCh38 (EBI, 20190312 biallelic SNV+INDEL, shapeit2 phased),
# EUR super-population (503 samples) -- an EXTERNAL panel; the donors are not used for LD.
#   https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
#   https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
# For each lead SNP only the +/- 500 kb window is read remotely (bcftools + tabix index), so no
# whole-chromosome VCF is downloaded. A new small env "ldtools" (conda-forge+bioconda: bcftools) is
# created; the analysis env "prs" is not modified. r2 is computed in 05_PRS/05_ld_proxies.py.
# Resources: ~29 remote slices of 1 Mb, then numpy: 2 CPUs / 16 GB, `testing` partition.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
RES="$PRS/validation/results/loci"
GT="$RES/ld_gt"
mkdir -p "$GT" "$RES/ld_work"
WIN=500000
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"

eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda env list | awk '{print $1}' | grep -qx "ldtools"; then
  conda create -y -n ldtools --override-channels --strict-channel-priority -c conda-forge -c bioconda bcftools
fi

# EUR sample list (503 samples)
cd "$RES/ld_work"
curl -sS -L --max-time 60 -o 1kg.panel "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
awk '$3=="EUR" {print $1}' 1kg.panel > eur_samples.txt
echo "EUR samples: $(wc -l < eur_samples.txt)"

conda activate ldtools
bcftools --version | head -n 1
# bcftools writes the downloaded .tbi next to the cwd: stay in ld_work
tail -n +2 "$RES/lead_snps_lookup.tsv" | while IFS=$'\t' read -r rsid chr pos rest; do
  [ "$chr" = "NA" ] && { echo "skip $rsid (not in sumstats)"; continue; }
  start=$((pos - WIN)); [ $start -lt 1 ] && start=1
  end=$((pos + WIN))
  url="$BASE/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  out="$GT/${rsid}.gt.tsv.gz"
  [ -s "$out" ] && { echo "have $rsid"; continue; }
  echo "== $rsid chr$chr:$start-$end"
  # NOTE: contigs in these VCFs are named "22", NOT "chr22" (a "chr" prefix returns an empty slice silently)
  ( bcftools view -S eur_samples.txt --force-samples -r "${chr}:${start}-${end}" -m2 -M2 -v snps "$url" \
      | bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' | gzip > "$out" ) \
    || { echo "FAILED slice for $rsid (chr$chr)"; rm -f "$out"; continue; }
  # an empty slice still yields a 20-byte gzip stream: treat it as a failure
  if [ "$(zcat "$out" | head -c 1 | wc -c)" -eq 0 ]; then echo "EMPTY slice for $rsid"; rm -f "$out"; fi
done
conda deactivate

echo "### r2 and coverage"
conda activate prs
python -I "${CODE_DIR}/05_PRS/05_ld_proxies.py" \
  --leads "$RES/lead_snps_lookup.tsv" \
  --gtdir "$GT" \
  --qc-bim "$PRS/genotypes/plink_qc/b38_qc_rsid.bim" \
  --raw-bim "$PRS/genotypes/plink_bin/b38.bim" \
  --panel-snpinfo "$PRS/reference/ldblk_ukbb_eur/snpinfo_ukbb_hm3" \
  --outdir "$RES"
