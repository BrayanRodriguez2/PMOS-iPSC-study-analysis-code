#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/02_geneset_score.sbatch   (original step script: 72_geneset_score.sbatch)
# Produces: Replication of the GWAS-gene-set score test with expression-matched random sets.
# Run: sbatch 07_gene_set_analyses/02_geneset_score.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=geneset_score
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# Step 72 - replication of Fig. 5B / Table S7b (see 07_gene_set_analyses/02_geneset_score.py) on our stranded counts (06_RNAseq_counts_and_DE_validation) and, as a check, on the original Normalized_counts.tabular.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; D="$H/DEG_validation"; G="$H/GSEA_validation"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
for spec in "ours_s2:$D/results/step1/norm_counts_symbol.tsv.gz" "theirs:$D/package/Normalized_counts.tabular"; do
  tag=${spec%%:*}; f=${spec#*:}
  python -I "${CODE_DIR}/07_gene_set_analyses/02_geneset_score.py" --norm "$f" --tag "$tag" --meta "$D/results/clinical/donor_meta.tsv" \
    --s7a "$G/package_scrt/sheet10_S7a_GWAS_genes.tsv" --s7b "$G/package_scrt/sheet11_S7b_geneset_tests.tsv" --out "$G/results/fig5B"
done

# sensitivity: our counts without AMH (the only gene whose expression status differs from the original besides CYP3A5, which is below the filter here)
python -I "${CODE_DIR}/07_gene_set_analyses/02_geneset_score.py" --norm "$D/results/step1/norm_counts_symbol.tsv.gz" --tag ours_s2_noAMH --exclude-genes AMH --meta "$D/results/clinical/donor_meta.tsv" \
  --s7a "$G/package_scrt/sheet10_S7a_GWAS_genes.tsv" --s7b "$G/package_scrt/sheet11_S7b_geneset_tests.tsv" --out "$G/results/fig5B"
