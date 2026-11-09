#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/08_setup_prscs_env.sbatch   (original step script: 14b_setup_prscs_env.sbatch)
# Produces: Conda environment 'prscs' (NumPy 1.26) and PRS-CS checkout.
# Run: sbatch 05_PRS/08_setup_prscs_env.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_env2
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 14b - Environment "prscs": PRS-CS needs an OLD numpy.
# Why: PRS-CS (commit fbc9ec0, 2018-era code) calls float() on 1-element numpy arrays in gigrnd.py;
# with the numpy installed in env "prs" (2.x) this raises
#   TypeError: only 0-dimensional arrays can be converted to Python scalars
# (seen in the chr22 test, job 3276963). Instead of patching the PRS-CS source (which would no
# longer match the commit recorded in software/PRScs_COMMIT.txt) PRS-CS is run in this separate env:
# python 3.11, numpy 1.26, scipy 1.13, h5py (conda-forge only, strict priority).
# Env "prs" (python 3.11, scipy, h5py, plink2) is left untouched.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda env list | awk '{print $1}' | grep -qx "prscs"; then
  conda create -y -n prscs --override-channels --strict-channel-priority -c conda-forge \
    python=3.11 "numpy=1.26" "scipy=1.13" h5py
fi
conda env export -n prscs --no-builds > "$PRS/env/environment_prscs.yml"   # before activating (see step 01)
conda activate prscs
python -c "import sys, numpy, scipy, h5py; print('python', sys.version.split()[0], '| numpy', numpy.__version__, '| scipy', scipy.__version__, '| h5py', h5py.__version__)"
cd "$PRS/software/PRScs" && python PRScs.py --help | head -n 3
