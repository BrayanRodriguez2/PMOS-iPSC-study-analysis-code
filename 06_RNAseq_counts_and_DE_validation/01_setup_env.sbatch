#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/01_setup_env.sbatch   (original step script: 61_setup_deg_env.sbatch)
# Produces: Conda environment 'degval' (featureCounts, R, DESeq2, limma, edgeR, Python).
# Run: sbatch 06_RNAseq_counts_and_DE_validation/01_setup_env.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=deg_env
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 61 - conda environment "degval" for the expression-validation task: subread (featureCounts), RSeQC (infer_experiment), R + DESeq2 + limma + edgeR,
# python (numpy, scipy, pandas, matplotlib). Strict channel priority conda-forge > bioconda, as for the other environments. The exported YAML is kept in PRS/env/.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
if conda env list | awk '{print $1}' | grep -qx degval; then echo "env degval already exists"; else
  conda create -y -n degval --override-channels --strict-channel-priority -c conda-forge -c bioconda \
    python=3.11 numpy scipy pandas matplotlib r-base=4.3 bioconductor-deseq2 bioconductor-limma bioconductor-edger subread rseqc samtools
fi
conda activate degval
echo "### versions"; featureCounts -v 2>&1 | head -n 2; R --version | head -n 1
Rscript -e 'for (p in c("DESeq2","limma","edgeR")) cat(p, as.character(packageVersion(p)), "\n")'
python -c "import numpy, scipy, pandas, matplotlib; print('numpy', numpy.__version__, 'scipy', scipy.__version__, 'pandas', pandas.__version__)"
infer_experiment.py --version 2>&1 | head -n 1
conda env export -n degval --no-builds > "${PROJECT_DIR}/PRS/env/environment_degval.yml"
echo done
