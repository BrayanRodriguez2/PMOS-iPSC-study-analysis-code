#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Example configuration. Copy to config/config.sh, edit, and `source` it in the shell from which you submit jobs
# (SLURM exports the environment to the jobs). All paths below are placeholders; nothing here is a real path or secret.
# Do NOT commit your edited copy (config.sh is listed in .gitignore).
# ---------------------------------------------------------------------------

# ---- Directories ------------------------------------------------------------
export PROJECT_DIR="/path/to/project"                    # working directory (layout described in README.md)
export CODE_DIR="/path/to/PMOS_iPSC_code_release"        # root of this repository (scripts are called as ${CODE_DIR}/<folder>/<script>)
export RAW_GENOTYPE_DIR="/path/to/raw_genotyping"        # GenomeStudio export (b38 .ped/.map), sample-sheet CSV, IDAT folders, .bpm and .egt files
export FASTQ_DIR="/path/to/rnaseq_fastq"                 # one sub-folder per iPSC line: <line>/<line>_S*_L001_R1_001.fastq.gz

# ---- Software set-up (site specific) ------------------------------------------
# Command(s) that make `conda` available in a job. Use 'true' if conda is already on the PATH.
export CONDA_SETUP='module purge; module load miniconda3'
# Command(s) that make 7z available (only 04_post_imputation_QC/03_unzip_check.sbatch). Use 'true' if 7z is on the PATH.
export P7ZIP_SETUP='module purge; module load p7zip'

# ---- TOPMed Imputation Server (secrets: never commit) -----------------------------
export TOPMED_PW_FILE="${HOME}/.topmed_pw"               # file with the password of the encrypted results (chmod 600)
# Private share links given by the server for each job (used by 03_imputation_TOPMed/06_qc_only_download.sh):
export TOPMED_URL_QC_REPORT="PASTE_PRIVATE_URL"
export TOPMED_URL_QUALITY_CONTROL="PASTE_PRIVATE_URL"
export TOPMED_URL_CHUNKS_EXCLUDED="PASTE_PRIVATE_URL"
export TOPMED_URL_SNPS_EXCLUDED="PASTE_PRIVATE_URL"
export TOPMED_URL_SNPS_TYPED_ONLY="PASTE_PRIVATE_URL"
# (the link of the full imputation job is passed as the first argument of 04_post_imputation_QC/01_fetch_download_script.sh)

# ---- Sample-ID lists (comma separated; NOT included in the repository) -----------------
# Analyses that annotate or exclude specific samples read their IDs from these variables. Leave empty if not needed;
# set them with the IDs held by the data owners. IDs are the array sample IDs (IID of the PLINK .fam files).
export FOCAL_SAMPLE_IDS=""              # samples examined in the contamination/identity diagnostics (02_contamination_identity)
export COMPARISON_SAMPLE_IDS=""         # extra samples drawn next to the focal ones in 02_contamination_identity/11_het_baf_vs_rna.py
export ANNOTATE_SAMPLE_IDS=""           # samples labelled in the IBS0 figure (02_contamination_identity/12_ibs0_donor_profile.py)
export CONTAMINATED_SAMPLE_IDS=""       # samples with a contamination signal, highlighted in 03_.../07_qc_only_verify.py and 04_.../08_sample_checks.py
export SENSITIVITY_EXCLUDE_IDS=""       # samples removed from the 135-sample sensitivity set (05_PRS/22_prs_analysis.py)
export LOWLEVEL_SAMPLE_IDS=""           # sample(s) with a low-level contamination estimate, excluded in one sensitivity analysis (05_PRS/22_prs_analysis.py)
export KEY_CHECK_SAMPLE_IDS=""          # donors printed in the array-ID <-> line key check (08_preprint_tasks/02_array_pca_report.py)
export EXTRA_PLOT_SAMPLE_IDS=""         # samples marked in the all-140 PCA figure (08_preprint_tasks/02_array_pca_report.py)
export EXCLUDED_FROM_135_IDS=""         # samples missing from the 135-sample sensitivity set (08_preprint_tasks/07_update_S2_S3.py)
export RELATED_PAIR_SAMPLE_IDS=""       # samples of a related pair omitted from the genotype-PC panel of Fig 2G (original_analysis/.../figs.py; IDs held by the authors)

# ---- SLURM -----------------------------------------------------------------------
# The .sbatch files contain the placeholders YOUR_SLURM_ACCOUNT and YOUR_PARTITION in their #SBATCH lines (SLURM does not expand variables there).
# Replace them once, e.g.:  sed -i 's/YOUR_SLURM_ACCOUNT/myaccount/; s/YOUR_PARTITION/myqueue/' ${CODE_DIR}/*/*.sbatch
# (the analysis used a short "testing" queue for light jobs and the main queue for long ones).
