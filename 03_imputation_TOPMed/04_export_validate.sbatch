#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/04_export_validate.sbatch   (original step script: 21_topmed_export_validate.sbatch)
# Produces: Per-chromosome VCF export and validation for the TOPMed server.
# Run: sbatch 03_imputation_TOPMed/04_export_validate.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=tm_export
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 21 - Export one bgzip VCF per chromosome (chr1 ... chr22) for upload and validate them.
# Server rules (docs V2, "Data preparation"): bgzip, one file per chromosome, sorted, GRCh38 with "chr" prefix,
# VCF version <= 4.2, GT present, 20 <= samples <= 25,000, <= 20,000 variants per 10 Mb chunk.
# Input : genotypes/topmed/aligned/all_chr.vcf.gz (step 20; 140 samples)
# Output: PRS/topmed/upload/chr<N>.vcf.gz (+ .tbi) and upload_manifest.tsv (size, md5, variants, samples);
#         results/topmed/{vcf_checks.tsv, vcf_checks_report.txt, bcftools_validation.txt}
# Validation: bcftools (index/sort/REF check against hg38.fa with `norm --check-ref e`, per-sample missingness)
#   + 03_imputation_TOPMed/04_export_validate.py (contig names, order, alleles, variants per 10 Mb chunk). checkVCF.py is deliberately not run (see note below).
# Resources: small; 4 CPUs / 16 GB on `testing`.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
AL="$PRS/genotypes/topmed/aligned/all_chr.vcf.gz"
UP="$PRS/topmed/upload"
RES="$PRS/validation/results/topmed"
FA="$PRS/reference/grch38/hg38.fa"
TMP="$RES/vcf_work"
mkdir -p "$UP" "$RES" "$TMP"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate topmed
VAL="$RES/bcftools_validation.txt"
: > "$VAL"

echo "### export one VCF per chromosome"
for c in $(seq 1 22); do
  bcftools view -r "chr$c" "$AL" -Oz -o "$UP/chr$c.vcf.gz"
  # force the header to VCFv4.2 if bcftools wrote another version
  # (sed -n '1p' instead of head: head closes the pipe early and `pipefail` would turn the SIGPIPE into a failure)
  v=$(bcftools view -h "$UP/chr$c.vcf.gz" | sed -n '1p')
  if [ "$v" != "##fileformat=VCFv4.2" ]; then
    echo "chr$c header was '$v' -> rewritten to VCFv4.2" | tee -a "$VAL"
    bcftools view -h "$UP/chr$c.vcf.gz" | sed '1s/.*/##fileformat=VCFv4.2/' > "$TMP/hdr.txt"
    bcftools reheader -h "$TMP/hdr.txt" "$UP/chr$c.vcf.gz" -o "$TMP/chr$c.re.vcf.gz"
    mv "$TMP/chr$c.re.vcf.gz" "$UP/chr$c.vcf.gz"
  fi
  bcftools index -f -t "$UP/chr$c.vcf.gz"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$UP/chr$c.vcf.gz" > "$TMP/chr$c.pos.tsv"
done

echo "### bcftools validation (per file)" | tee -a "$VAL"
printf "file\tbytes\tmd5\tvariants\tsamples\theader\tref_check_vs_hg38\tmissing_genotypes\n" > "$UP/upload_manifest.tsv"
for c in $(seq 1 22); do
  f="$UP/chr$c.vcf.gz"
  nv=$(bcftools view -H "$f" | wc -l)
  ns=$(bcftools query -l "$f" | wc -l)
  hdr=$(bcftools view -h "$f" | sed -n '1p')
  # norm --check-ref e fails (non-zero exit) at the first record whose REF differs from the FASTA
  if bcftools norm -f "$FA" --check-ref e "$f" -Ou -o /dev/null 2>"$TMP/norm$c.err"; then rc="REF=hg38 for all records"; else rc="MISMATCH: $(head -n 1 "$TMP/norm$c.err")"; fi
  miss=$(bcftools stats "$f" | awk -F'\t' '$1=="SN" && /missing/ && !done {print $4; done=1}')
  sz=$(stat -c %s "$f"); md=$(md5sum "$f" | cut -d' ' -f1)
  printf "chr%s.vcf.gz\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$c" "$sz" "$md" "$nv" "$ns" "$hdr" "$rc" "${miss:-NA}" | tee -a "$UP/upload_manifest.tsv" >> "$VAL"
done
echo "--- GT is the only FORMAT field? (first record of chr1)" | tee -a "$VAL"
bcftools view -H "$UP/chr1.vcf.gz" | sed -n '1p' | cut -f1-12 | tee -a "$VAL"
bcftools view -h "$UP/chr1.vcf.gz" | grep -E '^##(fileformat|FORMAT|contig)' | sed -n '1,8p' | tee -a "$VAL"

echo "### per-sample call rate over all chromosomes (minimum)" | tee -a "$VAL"
# bcftools stats PSC columns: $4 nRefHom, $5 nNonRefHom, $6 nHets, ..., $14 nMissing
bcftools stats -s - "$AL" | awk -F'\t' '$1=="PSC" {n=$4+$5+$6+$14; cr=1-$14/n; if (min==""||cr<min) {min=cr; who=$3}; sum+=cr; k++}
     END {printf "samples: %d | min sample call rate: %.4f (sample %s) | mean: %.4f\n", k, min, who, sum/k}' | tee -a "$VAL"

echo "### structural checks (python)"
python -I "${CODE_DIR}/03_imputation_TOPMed/04_export_validate.py" --dir "$TMP" --outdir "$RES"

# NOTE: checkVCF.py (the validator named in the server docs) is NOT run here: it would have to be downloaded from a
# third-party GitHub repository without a pinned version and executed on the cluster. Its checks (VCF validity,
# REF vs FASTA, sorted records, contig names) are covered above by bcftools (index/sort, norm --check-ref e) and
# 03_imputation_TOPMed/04_export_validate.py. It can be run by hand from a pinned release if wanted:
#   checkVCF.py -r hg38.fa -o out chr22.vcf.gz   (the "Action Items" it prints can be ignored, per the docs)
echo "### manifest"; cat "$UP/upload_manifest.tsv" | cut -f1-5,7
du -sh "$UP"
