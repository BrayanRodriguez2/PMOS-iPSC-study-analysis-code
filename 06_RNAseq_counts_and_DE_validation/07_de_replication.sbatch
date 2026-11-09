#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/07_de_replication.sbatch   (original step script: 65_deg_run.sbatch)
# Produces: Independent replication of the expression analyses (OLS, PCA, permutations, variance partition, chr20).
# Run: sbatch 06_RNAseq_counts_and_DE_validation/07_de_replication.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=deg_rep
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 65 - clinical table (65a), gene positions (GENCODE v50) and the independent replication (06_RNAseq_counts_and_DE_validation/07_de_replication.py) on three inputs:
#   ours_s2  our counts, stranded (-s 2, the correct orientation)   [main]
#   ours_s0  our counts counted unstranded (-s 0)                    [to see whether the counting mode explains differences with the original file]
#   theirs   the original Normalized_counts.tabular                          [our code on the original data: separates data from code differences]
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; D="$H/DEG_validation"; S="${CODE_DIR}"; P="$D/package"; R="$D/results"
GTF="$H/PRS/reference/star_gencode_v50/gencode.v50.primary_assembly.annotation.gtf"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
mkdir -p "$R/clinical" "$R/analysis"
python -I "$S/06_RNAseq_counts_and_DE_validation/05_clinical.py" --clinical "$P/clinical_tsv/sheet03_Matching_Final_sent_from_Oulu.tsv" --s1 "$P/af1_tsv/sheet02_S1_donors.tsv" --rna-samples "$H/PRS/rnaseq/rna_samples.tsv" --out "$R/clinical"
awk -F'\t' 'BEGIN{OFS="\t"; print "gene_id","gene_name","chr","start","end","gene_type"} $3=="gene" {id="";nm="";ty=""; if (match($9,/gene_id "[^"]+"/)) id=substr($9,RSTART+9,RLENGTH-10); if (match($9,/gene_name "[^"]+"/)) nm=substr($9,RSTART+11,RLENGTH-12); if (match($9,/gene_type "[^"]+"/)) ty=substr($9,RSTART+11,RLENGTH-12); print id,nm,$1,$4,$5,ty}' "$GTF" > "$R/gene_pos.tsv"
echo "gene positions: $(($(wc -l < "$R/gene_pos.tsv") - 1))"
for spec in "ours_s2:$R/step1/norm_counts_symbol.tsv.gz:--figs" "ours_s0:$R/step1_s0/norm_counts_symbol.tsv.gz:" "theirs:$P/Normalized_counts.tabular:"; do
  tag=${spec%%:*}; rest=${spec#*:}; f=${rest%%:*}; fig=${rest#*:}
  echo "################ $tag"
  python -I "$S/06_RNAseq_counts_and_DE_validation/07_de_replication.py" --norm "$f" --tag "$tag" --meta "$R/clinical/donor_meta.tsv" --s4 "$P/af1_tsv/sheet05_S4_DE.tsv" --s15 "$P/af1_tsv/sheet19_S15_top40_genes.tsv" \
    --gene-pos "$R/gene_pos.tsv" --out "$R/analysis" $fig
done
echo done
