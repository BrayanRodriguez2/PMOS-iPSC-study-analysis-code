#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/06_star_align.sbatch   (original step script: 26_star_align.sbatch)
# Produces: STAR alignment of the 14 RNA-seq FASTQ files (array job).
# Run: sbatch 02_contamination_identity/06_star_align.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=star_aln
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --array=0-13%4
#SBATCH --output=logs/%x_%A_%a.out
# ---------------------------------------------------------------------------
# Step 26 - Part B of the contamination diagnostic: align the 14 single-end RNA-seq FASTQ files (one per iPSC line) to GRCh38
# with STAR + GENCODE v50 (index from step 25). Array task i = i-th FASTQ (sorted by sample number S1..S14);
# at most 4 tasks run at the same time (each loads the ~30 GB index).
#
# Input (read-only): ${FASTQ_DIR}/HEL*/HEL*_S*_L001_R1_001.fastq.gz
# Output: PRS/rnaseq/star/<fastq_name>/{Aligned.sortedByCoord.out.bam(+.bai), ReadsPerGene.out.tab, Log.final.out, ...}
#         PRS/rnaseq/rna_samples.tsv  (fastq_name, ipsc_line, esf_id, group, fastq path)
# NAMING NOTE: the FASTQ "HEL362_3" is really the clone HEL362.1 (naming error confirmed by the lab). The file is NOT
# renamed; rna_samples.tsv carries the corrected line name (ipsc_line) next to the original fastq_name.
# Options: ENCODE-style STAR parameters (BySJout filtering, max 4 % mismatches over read length, intron 20 b - 1 Mb);
# `--quantMode GeneCounts` is added because the same alignments feed the RNA-seq pipeline (counts are not used here).
# Read group: ID = fastq_name, SM = corrected line name.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
IDX="$PRS/reference/star_gencode_v50/star_index"
FQDIR=${FASTQ_DIR}
OUT="$PRS/rnaseq"
mkdir -p "$OUT/star"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate rnaseq

# sample table (identical for every task; written once per task, idempotent)
KEY="$PRS/metadata/donor_key.tsv"        # ipsc_line  esf_id  group (header + 14 rows)
ls "$FQDIR"/HEL*/HEL*_S*_L001_R1_001.fastq.gz | sort -V > "$OUT/fastq_list.tmp.$$" && mv -f "$OUT/fastq_list.tmp.$$" "$OUT/fastq_list.txt"
{
  printf "fastq_name\tipsc_line\tesf_id\tgroup\tfastq\n"
  while read -r fq; do
    name=$(basename "$(dirname "$fq")")                 # e.g. HEL351_3
    line=$(echo "$name" | tr '_' '.')                    # HEL351.3
    [ "$name" = "HEL362_3" ] && line="HEL362.1"          # naming error of the FASTQ (see note above)
    donor=$(awk -F'\t' -v l="$line" '$1==l {print $2"\t"$3}' "$KEY")
    printf "%s\t%s\t%s\t%s\n" "$name" "$line" "$donor" "$fq"
  done < "$OUT/fastq_list.txt"
} > "$OUT/rna_samples.tmp.$$" && mv -f "$OUT/rna_samples.tmp.$$" "$OUT/rna_samples.tsv"   # atomic: 14 tasks write it

i=$SLURM_ARRAY_TASK_ID
FQ=$(sed -n "$((i + 1))p" "$OUT/fastq_list.txt")
NAME=$(basename "$(dirname "$FQ")")
LINE=$(awk -F'\t' -v n="$NAME" '$1==n {print $2}' "$OUT/rna_samples.tsv")
D="$OUT/star/$NAME"; mkdir -p "$D"
echo "task $i: $NAME ($LINE) <- $FQ"
STAR --version
TMPD="${TMPDIR:-/tmp}/star_${SLURM_JOB_ID}_$i"; rm -rf "$TMPD"
STAR --runThreadN "$SLURM_CPUS_PER_TASK" --genomeDir "$IDX" \
     --readFilesIn "$FQ" --readFilesCommand zcat \
     --outFileNamePrefix "$D/" --outTmpDir "$TMPD" \
     --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000 \
     --outSAMattrRGline ID:"$NAME" SM:"$LINE" PL:ILLUMINA \
     --outSAMattributes NH HI AS nM NM MD \
     --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 --alignIntronMax 1000000 --quantMode GeneCounts
rm -rf "$TMPD"
mv "$D/Aligned.sortedByCoord.out.bam" "$D/$NAME.bam"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$D/$NAME.bam"
echo "### mapping summary"; grep -E 'Number of input reads|Uniquely mapped reads (number|%)|% of reads mapped to multiple loci' "$D/Log.final.out"
