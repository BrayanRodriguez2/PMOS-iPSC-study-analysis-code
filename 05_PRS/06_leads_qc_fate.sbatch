#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/06_leads_qc_fate.sbatch   (original step script: 12b_leads_qc_fate.sbatch)
# Produces: What happened to the lead SNPs in the SNP QC.
# Run: sbatch 05_PRS/06_leads_qc_fate.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_leadfate
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 12b - Which QC filter removed the lead SNPs that ARE on the raw array but not in the QC
# dataset (from locus_coverage.tsv: lead_on_raw_array != no and lead_on_qc_array == no)?
# For each of them, in the samples that pass sample QC: call rate, MAF, HWE p, monomorphic?
# Filters of step 10: call rate >= 98 %, monomorphic, HWE p < 1e-6, MAF < 0.01.
# Read-only on the genotypes; output to results/loci/lead_qc_fate.txt
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
RES="$PRS/validation/results/loci"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs

awk -F'\t' 'NR>1 && $7!="no" && $7!="NA" && $8=="no" {print $7}' "$RES/locus_coverage.tsv" > "$RES/ld_work/leads_removed_by_qc_ids.txt"
echo "raw-array IDs of leads removed by QC:"; cat "$RES/ld_work/leads_removed_by_qc_ids.txt"
plink2 --bfile "$PRS/genotypes/plink_bin/b38" \
       --keep "$PRS/validation/results/qc/samples/keep_samples.txt" \
       --extract "$RES/ld_work/leads_removed_by_qc_ids.txt" \
       --freq --missing variant-only --hardy \
       --out "$RES/ld_work/leads_removed_by_qc"
{
  echo "## ALT_FREQS / call rate / HWE of lead SNPs removed by QC (129 QC-pass samples)"
  echo "## columns afreq:"; head -n 1 "$RES/ld_work/leads_removed_by_qc.afreq"; cat "$RES/ld_work/leads_removed_by_qc.afreq" | tail -n +2
  echo "## vmiss:"; head -n 1 "$RES/ld_work/leads_removed_by_qc.vmiss"; tail -n +2 "$RES/ld_work/leads_removed_by_qc.vmiss"
  echo "## hardy:"; head -n 1 "$RES/ld_work/leads_removed_by_qc.hardy"; tail -n +2 "$RES/ld_work/leads_removed_by_qc.hardy"
} | tee "$RES/lead_qc_fate.txt"
