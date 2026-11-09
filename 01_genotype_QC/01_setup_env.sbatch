#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/01_setup_env.sbatch   (original step script: 01_setup_env.sbatch)
# Produces: Conda environment 'prs' (PLINK 1.9/2, Python, R).
# Run: sbatch 01_genotype_QC/01_setup_env.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_env
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 01 - Build the conda environment "prs" and install PRS-CS.
#
# What it does:
#   1. Loads the RHEL9 miniconda module (miniconda3/25.5.1; NOT the centos7 ones).
#   2. Creates env "prs" using ONLY conda-forge + bioconda (--override-channels,
#      strict channel priority) with python 3.11, scipy, h5py and plink2.
#   3. Exports the env to PRS/env/environment.yml (--no-builds).
#   4. Verifies plink2, scipy and h5py.
#   5. Clones PRS-CS (https://github.com/getian107/PRScs), records the commit
#      and checks that PRScs.py starts.
# Resources: small (conda solving is single-threaded); `testing` partition is enough.
# Idempotent: re-running skips steps that are already done.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail

PRS_DIR="${PROJECT_DIR}/PRS"
mkdir -p "$PRS_DIR/env" "$PRS_DIR/software"

eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "### conda: $(conda --version) | base: $(conda info --base)"

if conda env list | awk '{print $1}' | grep -qx "prs"; then
  echo "### env 'prs' already exists - skipping creation"
else
  # --override-channels + --strict-channel-priority: ignore the local .condarc / defaults channels
  conda create -y -n prs \
    --override-channels --strict-channel-priority \
    -c conda-forge -c bioconda \
    python=3.11 scipy h5py plink2
fi

# Export BEFORE activating: once "prs" is active its python shadows the module's
# python and the module's `conda` executable fails with "No module named 'conda'".
conda env export -n prs --no-builds > "$PRS_DIR/env/environment.yml"
echo "### exported $PRS_DIR/env/environment.yml"
conda activate prs

echo "### verification"
plink2 --version
python -c "import sys, scipy, h5py; print('python', sys.version.split()[0], '| scipy', scipy.__version__, '| h5py', h5py.__version__)"

# --- PRS-CS ---------------------------------------------------------------
if [ ! -d "$PRS_DIR/software/PRScs/.git" ]; then
  git clone https://github.com/getian107/PRScs.git "$PRS_DIR/software/PRScs"
fi
cd "$PRS_DIR/software/PRScs"
git rev-parse HEAD > "$PRS_DIR/software/PRScs_COMMIT.txt"
git log -1 --format='commit %H%ndate   %ad%nsubject %s' >> "$PRS_DIR/software/PRScs_COMMIT.txt"
cat "$PRS_DIR/software/PRScs_COMMIT.txt"
echo "### PRScs.py --help"
python PRScs.py --help || true   # PRS-CS prints its usage and may exit non-zero
