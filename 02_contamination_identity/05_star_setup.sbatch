#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/05_star_setup.sbatch   (original step script: 25_star_setup.sbatch)
# Produces: Conda environment 'rnaseq', GENCODE v50 references and STAR genome index.
# Run: sbatch 02_contamination_identity/05_star_setup.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=star_idx
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 25 - Part B of the contamination diagnostic (identity check against the RNA-seq): environment, GENCODE references
# and STAR genome index. The index is reusable for the RNA-seq pipeline of the project.
#
# (1) conda env "rnaseq" (conda-forge + bioconda, strict priority): STAR (2.7.11b at the time of writing), samtools,
#     bcftools, python 3.11 + numpy/pandas/scipy/matplotlib (the last ones for step 28).
# (2) References (GENCODE human release 50, GRCh38; https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/):
#     GRCh38.primary_assembly.genome.fa.gz (chromosomes 1-22, X, Y, M + unlocalized/unplaced scaffolds, "chr" names)
#     gencode.v50.primary_assembly.annotation.gtf.gz (comprehensive annotation, primary assembly). md5 checked against
#     the MD5SUMS file published in the same folder.
# (3) STAR index: --sjdbOverhang 80 (read length 81 bp -> overhang = read length - 1), GTF = GENCODE v50.
# Output: PRS/reference/star_gencode_v50/{GRCh38.primary_assembly.genome.fa, gencode.v50.primary_assembly.annotation.gtf,
#         MD5SUMS, star_index/}; PRS/env/environment_rnaseq.yml; versions in the log.
# Resources: index generation ~30-60 min with 12 threads, ~35 GB RAM for the human genome; 12 CPUs / 48 GB on `main`.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
REF="$PRS/reference/star_gencode_v50"
BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50
mkdir -p "$REF" "$PRS/env"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx "rnaseq"; then
  conda create -y -n rnaseq --override-channels --strict-channel-priority -c conda-forge -c bioconda \
    python=3.11 star samtools bcftools numpy pandas scipy matplotlib
fi
conda env export -n rnaseq --no-builds > "$PRS/env/environment_rnaseq.yml"   # before activating (see step 01)
conda activate rnaseq
echo "### versions"; STAR --version; samtools --version | sed -n '1p'; bcftools --version | sed -n '1p'

cd "$REF"
for f in GRCh38.primary_assembly.genome.fa.gz gencode.v50.primary_assembly.annotation.gtf.gz; do
  wget -c -nv -O "$f" "$BASE/$f"
done
wget -nv -O MD5SUMS "$BASE/MD5SUMS"
echo "### md5 vs GENCODE MD5SUMS"
for f in GRCh38.primary_assembly.genome.fa.gz gencode.v50.primary_assembly.annotation.gtf.gz; do
  exp=$(grep " $f\$" MD5SUMS | cut -d' ' -f1); got=$(md5sum "$f" | cut -d' ' -f1)
  echo "$f expected=$exp got=$got"
  [ "$exp" = "$got" ] || { echo "MD5 MISMATCH for $f"; exit 1; }
done
[ -f GRCh38.primary_assembly.genome.fa ] || gunzip -k GRCh38.primary_assembly.genome.fa.gz
[ -f gencode.v50.primary_assembly.annotation.gtf ] || gunzip -k gencode.v50.primary_assembly.annotation.gtf.gz
echo "### GTF header"; grep '^##' gencode.v50.primary_assembly.annotation.gtf | sed -n '1,6p'

mkdir -p star_index
STAR --runMode genomeGenerate --runThreadN "$SLURM_CPUS_PER_TASK" \
     --genomeDir star_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v50.primary_assembly.annotation.gtf \
     --sjdbOverhang 80 \
     --outFileNamePrefix star_index/gen_
echo "### index done"; du -sh star_index; ls star_index | sed -n '1,20p'
