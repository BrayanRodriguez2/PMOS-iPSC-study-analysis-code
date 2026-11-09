#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/01_setup_gsea_env.sbatch   (original step script: 72a_setup_gsea_env.sbatch)
# Produces: Conda environment 'gsea' (R, fgsea, msigdbr).
# Run: sbatch 07_gene_set_analyses/01_setup_gsea_env.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=gsea_env
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 72a - conda environment "gsea" for the gene-set analyses: R 4.3 + fgsea + msigdbr (MSigDB Hallmark; the conda-forge build embeds the data) + ggplot2 + data.table.
# Separate from "degval" so that the DESeq2/limma environment of 06_RNAseq_counts_and_DE_validation is not touched. Strict priority conda-forge > bioconda.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
if conda env list | awk '{print $1}' | grep -qx gsea; then echo "env gsea already exists"; else
  conda create -y -n gsea --override-channels --strict-channel-priority -c conda-forge -c bioconda r-base=4.3 bioconductor-fgsea r-msigdbr r-ggplot2 r-data.table
fi
conda activate gsea
Rscript -e 'for (p in c("fgsea","msigdbr","ggplot2","data.table")) cat(p, as.character(packageVersion(p)), "\n"); m <- msigdbr::msigdbr(species = "Homo sapiens", category = "H"); cat("Hallmark rows:", nrow(m), "sets:", length(unique(m$gs_name)), "\n"); if ("db_version" %in% names(m)) cat("db_version", unique(m$db_version), "\n")'
echo done
