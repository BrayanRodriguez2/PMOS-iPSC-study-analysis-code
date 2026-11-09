#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/02_featurecounts.sbatch   (original step script: 63_featurecounts.sbatch)
# Produces: Gene counts with featureCounts in three strand settings (-s 0/1/2).
# Run: sbatch 06_RNAseq_counts_and_DE_validation/02_featurecounts.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=fcounts
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 63 - gene-level counts of the 14 STAR BAMs (step 26: GRCh38 primary assembly, GENCODE v50, whole genome, coordinate-sorted, 14 single-end FASTQ files)
# with featureCounts (Subread) in THREE strand settings, to decide the library orientation from the fraction of reads assigned:
#   -s 0 unstranded | -s 1 stranded (read = same strand as the gene) | -s 2 reverse-stranded (read = opposite strand; dUTP / TruSeq stranded mRNA)
# Common settings: single-end (no -p), -t exon, -g gene_id (gene level), GENCODE v50 comprehensive annotation, default multi-mapping handling (multi-mapping
# reads, flagged NH > 1 by STAR, are NOT counted), no -O (reads overlapping several genes are not counted), --extraAttributes gene_name,gene_type to keep the symbols.
# Output: DEG_validation/counts/fc_s{0,1,2}.txt (+ .summary); STAR's own ReadsPerGene.out.tab is summarised too (columns 2-4 = unstranded / strand 1 / strand 2).
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
H="${PROJECT_DIR}"; O="$H/DEG_validation/counts"; mkdir -p "$O"
GTF="$H/PRS/reference/star_gencode_v50/gencode.v50.primary_assembly.annotation.gtf"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate degval
featureCounts -v 2>&1 | head -n 2
BAMS=$(for d in "$H"/PRS/rnaseq/star/*/; do l=$(basename "$d"); echo "$H/PRS/rnaseq/star/$l/$l.bam"; done)   # one STAR output folder per line (step 06 of 02_contamination_identity)
for S in 0 1 2; do
  echo "### featureCounts -s $S"
  featureCounts -T "$SLURM_CPUS_PER_TASK" -s "$S" -t exon -g gene_id --extraAttributes gene_name,gene_type -a "$GTF" -o "$O/fc_s$S.txt" $BAMS 2>&1 | tail -n 4
done
echo "### assigned % per strand setting (Assigned / all alignments counted)"
for S in 0 1 2; do
  awk -v s="$S" 'NR>1 {for (i=2;i<=NF;i++) tot[i]+=$i; if ($1=="Assigned") for (i=2;i<=NF;i++) a[i]=$i} END {m=0; for (i in a) {m+=a[i]/tot[i]; n++}; printf "s%s mean assigned fraction over %d samples: %.4f\n", s, n, m/n}' "$O/fc_s$S.txt.summary"
done
echo done
