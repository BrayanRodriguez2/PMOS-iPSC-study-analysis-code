#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 08_preprint_tasks/01_array_pca.sbatch   (original step script: 80_array_pca.sbatch)
# Produces: Array genotype PCA (LD-pruned, long-range LD regions removed).
# Run: sbatch 08_preprint_tasks/01_array_pca.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=array_pca
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 80 - genotype PCA of the donor cohort from the ARRAY data (not imputed): the project QC dataset b38_qc_rsid (129 samples passing the standard sample QC; 447,260 SNPs after the
# SNP QC of step 10, GRCh38) and, for comparison, the same SNPs in all 140 samples (b38_qc_rsid_all140).
#   SNPs: autosomes only, MAF >= 0.05 (computed in the 129), long-range LD regions removed (list below), LD pruning --indep-pairwise 200 50 0.2 (window in variants, step 50, r2 0.2);
#   the SNP set is chosen on the 129 and used for both PCAs. PCA: plink2 --pca 128 (all components, so that the % of variance is eigenvalue / sum of eigenvalues), standardised genotypes.
# Long-range LD regions (GRCh38, APPROXIMATE: the regions of Price et al. 2008 / Anderson et al. 2010 as lifted over for GRCh38, plus the chr17q21 inversion): MHC chr6:25.4-33.4 Mb, chr8 inversion
#   chr8:7.1-13.1 Mb, chr17:43-46.5 Mb, and the other classic regions listed in the file below. Margins are generous; they are not an exact liftover.
# Output: PRS/validation/results/array_pca/{pca129.eigenvec,.eigenval,pca140.eigenvec,.eigenval,pruned_snps.txt,ld_regions_hg38.txt,...}
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
P="${PROJECT_DIR}/PRS"; Q="$P/genotypes/plink_qc"; O="$P/validation/results/array_pca"; mkdir -p "$O"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate prs
cat > "$O/ld_regions_hg38.txt" <<'EOF'
1	47534328	51534328	LD1
2	133742429	137242430	LD2
2	182135273	189135274	LD3
3	47458510	49962567	LD4
3	83450849	86950850	LD5
5	98664296	101164296	LD6
5	129664307	132664308	LD7
5	136164311	139164312	LD8
6	25392021	33392022	MHC
6	56892041	63142042	LD10
6	139958307	142458308	LD11
7	55225450	66225451	LD12
8	7142478	13142478	chr8_inversion
8	110987771	113487771	LD14
11	87789108	90289108	LD15
12	109062195	111562196	LD16
17	43000000	46500000	chr17q21_inversion
20	33412194	35912195	LD18
EOF
T="${SLURM_CPUS_PER_TASK:-1}"
plink2 --bfile "$Q/b38_qc_rsid" --autosome --maf 0.05 --exclude range "$O/ld_regions_hg38.txt" --indep-pairwise 200 50 0.2 --threads "$T" --out "$O/prune"
cp "$O/prune.prune.in" "$O/pruned_snps.txt"
echo "pruned SNPs: $(wc -l < "$O/pruned_snps.txt")"
plink2 --bfile "$Q/b38_qc_rsid" --extract "$O/pruned_snps.txt" --pca 128 --threads "$T" --out "$O/pca129"
plink2 --bfile "$Q/b38_qc_rsid_all140" --extract "$O/pruned_snps.txt" --pca 139 --threads "$T" --out "$O/pca140"
ls -l "$O" | head -n 20
