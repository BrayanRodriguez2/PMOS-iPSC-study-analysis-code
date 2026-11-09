#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/02_idat_to_gtc.sbatch   (original step script: 23_idat_to_gtc.sbatch)
# Produces: IDAT -> GTC (genotype calls, BAF, Log R Ratio) with the project manifest and cluster file.
# Run: sbatch 02_contamination_identity/02_idat_to_gtc.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=idat2gtc
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 23 - Part A of the contamination diagnostic: IDAT -> GTC (genotype calls + B allele frequency + Log R Ratio).
#
# Why this route: the project folder only has the raw IDAT files (2 per sample: _Grn and _Red), the .bpm manifest
# and the .egt cluster file (the GenomeStudio "Data" folder holds only a project database, not per-sample GTCs).
# `bcftools +idat2gtc` (Genovese, https://github.com/freeseek/gtc2vcf, v. 2024-09-27, bioconda bcftools-gtc2vcf-plugin
# 1.22) re-implements Illumina's normalisation + GenCall using the SAME .bpm and .egt, so BAF and LRR are computed
# with the project's own cluster definitions. The GTC files are then read directly in python (step 24) instead of
# going through a VCF, because no CSV manifest is available for this custom chip (ESTchip-2_GSAv3-MD).
#
# Inputs (read-only): ${RAW_GENOTYPE_DIR}/
#     <chip>/<chip>_R??C??_{Grn,Red}.idat  (7 chips, 144 sample positions), *.bpm, *.egt
# Output: PRS/intensity/gtc/<chip>/<chip>_R??C??.gtc   (one GTC per sample position)
# Resources: ~1-3 min per chip expected; 4 chips in parallel, 4 CPUs / 16 GB on `testing`.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
G=${RAW_GENOTYPE_DIR}
BPM="$G/ESTchip-2_GSAv3-MD_20032937X364826_A3.bpm"
EGT="$G/NG_1_23_26012021_v1.egt"
OUT="$PRS/intensity/gtc"
mkdir -p "$OUT"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
bcftools --version | sed -n '1p'

run_chip() {
  chip="$1"
  mkdir -p "$OUT/$chip"
  echo "[$(date +%T)] start $chip"
  bcftools +idat2gtc --bpm "$BPM" --egt "$EGT" --idats "$G/$chip" --output "$OUT/$chip" 2> "$OUT/$chip.log" \
    && echo "[$(date +%T)] done $chip: $(ls "$OUT/$chip"/*.gtc | wc -l) GTC files" \
    || { echo "[$(date +%T)] FAILED $chip (see $OUT/$chip.log)"; tail -n 5 "$OUT/$chip.log"; return 1; }
}
export -f run_chip; export OUT G BPM EGT

ls "$G" | grep -E '^[0-9]{12}$' | xargs -P 4 -I{} bash -c 'run_chip {}'

echo "### GTC files produced"; find "$OUT" -name '*.gtc' | wc -l
du -sh "$OUT"
echo "### one log"; sed -n '1,8p' "$OUT/206149980019.log"
