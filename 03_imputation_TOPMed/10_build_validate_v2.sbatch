#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/10_build_validate_v2.sbatch   (original step script: 36_build_validate_v2.sbatch)
# Produces: Cleaned v2 upload VCFs and their validation.
# Run: sbatch 03_imputation_TOPMed/10_build_validate_v2.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=build_v2
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 36 - Build and validate the v2 upload VCFs for the TOPMed Imputation Server.
# v2 = v1 (PRS/topmed/upload, 140 samples x 510,607 SNPs; kept untouched as a copy) minus the SNPs in exclude_v2.tsv (step 35):
#      182 TOPMed "Allele mismatch" + the frequency-discordant SNPs that 1000G FIN does not support. NO samples are removed.
# Output (NEW folders only; nothing is overwritten): PRS/topmed/upload_v2/chr1..22.vcf.gz (+ .tbi, upload_manifest_v2.tsv with md5)
#         PRS/validation/results/topmed_v2/{counts_v2.tsv, counts_by_chrom_v2.tsv, bcftools_validation_v2.txt, extra_checks_v2.txt, vcf_checks*.}
# Same checks as step 21/21b: bgzip -t, tabix index, "chr" contigs, VCF 4.2, GT only, REF = hg38 (bcftools norm --check-ref e),
# identical 140 samples in identical order in the 22 files, chunks of 10 Mb with < 3 variants (the server drops those) and > 20,000;
# plus: v2 sites == v1 sites minus the exclusion list, none of the excluded IDs left, genotypes of the kept sites identical to v1.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"
V1="$PRS/topmed/upload"; V2="$PRS/topmed/upload_v2"; RES="$PRS/validation/results/topmed_v2"
EX="$PRS/validation/results/topmed_qc_only/exclude_v2.tsv"; FA="$PRS/reference/grch38/hg38.fa"; W="$RES/work"
mkdir -p "$V2" "$RES" "$W"
if ls "$V2"/chr*.vcf.gz >/dev/null 2>&1; then echo "v2 files already exist in $V2: refusing to overwrite"; exit 1; fi
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
VAL="$RES/bcftools_validation_v2.txt"; : > "$VAL"

echo "### exclusion list"; tail -n +2 "$EX" | cut -f6 | sort | uniq -c | tee -a "$VAL"
tail -n +2 "$EX" | cut -f2,3 | sort -k1,1V -k2,2n > "$W/exclude_pos.tsv"
tail -n +2 "$EX" | cut -f1 | sort > "$W/exclude_ids.txt"
echo "positions to exclude: $(wc -l < "$W/exclude_pos.tsv") (unique: $(sort -u "$W/exclude_pos.tsv" | wc -l))" | tee -a "$VAL"

echo "### build v2"
for c in $(seq 1 22); do
  awk -F'\t' -v c="chr$c" '$1==c' "$W/exclude_pos.tsv" > "$W/exclude_chr$c.tsv"
  if [ -s "$W/exclude_chr$c.tsv" ]; then bcftools view -T "^$W/exclude_chr$c.tsv" "$V1/chr$c.vcf.gz" -Oz -o "$V2/chr$c.vcf.gz"
  else echo "chr$c: nothing to exclude" | tee -a "$VAL"; bcftools view "$V1/chr$c.vcf.gz" -Oz -o "$V2/chr$c.vcf.gz"; fi
  v=$(bcftools view -h "$V2/chr$c.vcf.gz" | sed -n '1p')
  if [ "$v" != "##fileformat=VCFv4.2" ]; then
    echo "chr$c header was '$v' -> rewritten to VCFv4.2" | tee -a "$VAL"
    bcftools view -h "$V2/chr$c.vcf.gz" | sed '1s/.*/##fileformat=VCFv4.2/' > "$W/hdr.txt"
    bcftools reheader -h "$W/hdr.txt" "$V2/chr$c.vcf.gz" -o "$W/chr$c.re.vcf.gz"; mv "$W/chr$c.re.vcf.gz" "$V2/chr$c.vcf.gz"
  fi
  bcftools index -f -t "$V2/chr$c.vcf.gz"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$V2/chr$c.vcf.gz" > "$W/chr$c.pos.tsv"
done

echo "### counts" | tee -a "$VAL"
printf "chrom\tv1_variants\texcluded_allele_mismatch\texcluded_freq_rule\texpected_v2\tv2_variants\tmatch\n" > "$RES/counts_by_chrom_v2.tsv"
for c in $(seq 1 22); do
  n1=$(bcftools index -n "$V1/chr$c.vcf.gz"); n2=$(bcftools index -n "$V2/chr$c.vcf.gz")
  ea=$(awk -F'\t' -v c="chr$c" 'NR>1 && $2==c && $6=="TOPMed_allele_mismatch"' "$EX" | wc -l)
  ef=$(awk -F'\t' -v c="chr$c" 'NR>1 && $2==c && $6!="TOPMed_allele_mismatch"' "$EX" | wc -l)
  printf "chr%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$c" "$n1" "$ea" "$ef" "$((n1 - ea - ef))" "$n2" "$([ $((n1 - ea - ef)) -eq "$n2" ] && echo OK || echo MISMATCH)" >> "$RES/counts_by_chrom_v2.tsv"
done
cat "$RES/counts_by_chrom_v2.tsv" | tee -a "$VAL"
T1=$(awk 'NR>1 {s+=$2} END {print s}' "$RES/counts_by_chrom_v2.tsv"); TA=$(awk 'NR>1 {s+=$3} END {print s}' "$RES/counts_by_chrom_v2.tsv")
TF=$(awk 'NR>1 {s+=$4} END {print s}' "$RES/counts_by_chrom_v2.tsv"); T2=$(awk 'NR>1 {s+=$6} END {print s}' "$RES/counts_by_chrom_v2.tsv")
printf "step\tvariants_after\tremoved\n" > "$RES/counts_v2.tsv"
printf "v1 (uploaded to the QC-only job)\t%s\t-\n" "$T1" >> "$RES/counts_v2.tsv"
printf -- "- TOPMed 'Allele mismatch'\t%s\t%s\n" "$((T1 - TA))" "$TA" >> "$RES/counts_v2.tsv"          # `--`: the format starts with "-"
printf -- "- frequency-discordant (> 0.2 vs TOPMed) not supported by 1000G FIN\t%s\t%s\n" "$((T1 - TA - TF))" "$TF" >> "$RES/counts_v2.tsv"
printf "v2 (final)\t%s\t%s\n" "$T2" "$((TA + TF))" >> "$RES/counts_v2.tsv"
cat "$RES/counts_v2.tsv" | tee -a "$VAL"

echo "### site-set check: v2 == v1 minus exclusion list; no excluded ID left; genotypes identical to v1 on kept sites" | tee -a "$VAL"
for c in $(seq 1 22); do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$V1/chr$c.vcf.gz" | sort > "$W/v1_$c.sorted"
  awk -F'\t' -v c="chr$c" '$1==c {print $1"\t"$2}' "$W/exclude_pos.tsv" | sort > "$W/ex_$c.sorted"
  awk -F'\t' 'NR==FNR {e[$1"\t"$2]=1; next} !(($1"\t"$2) in e)' "$W/ex_$c.sorted" "$W/v1_$c.sorted" | sort > "$W/expected_$c.sorted"
  sort "$W/chr$c.pos.tsv" > "$W/v2_$c.sorted"
  cmp -s "$W/expected_$c.sorted" "$W/v2_$c.sorted" || { echo "chr$c: v2 sites differ from v1 minus exclusions" | tee -a "$VAL"; exit 1; }
  left=$(bcftools query -f '%ID\n' "$V2/chr$c.vcf.gz" | sort | comm -12 - "$W/exclude_ids.txt" | wc -l)
  [ "$left" -eq 0 ] || { echo "chr$c: $left excluded IDs still present" | tee -a "$VAL"; exit 1; }
  if [ -s "$W/exclude_chr$c.tsv" ]; then g1=$(bcftools view -T "^$W/exclude_chr$c.tsv" "$V1/chr$c.vcf.gz" | bcftools query -f '[%GT ]\n' | md5sum | cut -d' ' -f1)
  else g1=$(bcftools query -f '[%GT ]\n' "$V1/chr$c.vcf.gz" | md5sum | cut -d' ' -f1); fi
  g2=$(bcftools query -f '[%GT ]\n' "$V2/chr$c.vcf.gz" | md5sum | cut -d' ' -f1)
  [ "$g1" = "$g2" ] || { echo "chr$c: genotypes differ from v1" | tee -a "$VAL"; exit 1; }
done
echo "OK: v2 sites = v1 sites minus exclusion list in all 22 chromosomes; no excluded ID left; genotypes identical to v1 on kept sites" | tee -a "$VAL"

echo "### bcftools validation (per file) + manifest" | tee -a "$VAL"
printf "file\tbytes\tmd5\tvariants\tsamples\theader\tref_check_vs_hg38\tmissing_genotypes\n" > "$V2/upload_manifest_v2.tsv"
for c in $(seq 1 22); do
  f="$V2/chr$c.vcf.gz"
  bgzip -t "$f" || { echo "chr$c bgzip -t FAILED" | tee -a "$VAL"; exit 1; }
  nv=$(bcftools view -H "$f" | wc -l); ns=$(bcftools query -l "$f" | wc -l); hdr=$(bcftools view -h "$f" | sed -n '1p')
  if bcftools norm -f "$FA" --check-ref e "$f" -Ou -o /dev/null 2>"$W/norm$c.err"; then rc="REF=hg38 for all records"; else rc="MISMATCH: $(sed -n '1p' "$W/norm$c.err")"; fi
  miss=$(bcftools stats "$f" | awk -F'\t' '$1=="SN" && /missing/ && !done {print $4; done=1}')
  printf "chr%s.vcf.gz\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$c" "$(stat -c %s "$f")" "$(md5sum "$f" | cut -d' ' -f1)" "$nv" "$ns" "$hdr" "$rc" "${miss:-NA}" | tee -a "$V2/upload_manifest_v2.tsv" >> "$VAL"
done
echo "--- FORMAT fields (first record of chr1) and contigs" | tee -a "$VAL"
bcftools view -H "$V2/chr1.vcf.gz" | sed -n '1p' | cut -f1-10 | tee -a "$VAL"
for c in $(seq 1 22); do bcftools view -H "$V2/chr$c.vcf.gz" | cut -f1 | sort -u; done | tr '\n' ' ' | tee -a "$VAL"; echo | tee -a "$VAL"
echo "--- only GT in FORMAT? distinct FORMAT columns across all 22 files:" | tee -a "$VAL"
for c in $(seq 1 22); do bcftools view -H "$V2/chr$c.vcf.gz" | cut -f9 | sort -u; done | sort -u | tee -a "$VAL"
echo "--- samples identical and in the same order in the 22 files, and the same as v1" | tee -a "$VAL"
bcftools query -l "$V1/chr1.vcf.gz" | md5sum | cut -d' ' -f1 > "$W/s1.md5"
d=0; for c in $(seq 1 22); do [ "$(bcftools query -l "$V2/chr$c.vcf.gz" | md5sum | cut -d' ' -f1)" = "$(cat "$W/s1.md5")" ] || { echo "chr$c samples differ" | tee -a "$VAL"; d=1; }; done
[ $d -eq 0 ] && echo "same 140 samples, same order, in all 22 files (identical to v1)" | tee -a "$VAL"
echo "--- per-sample call rate over all chromosomes" | tee -a "$VAL"
for c in $(seq 1 22); do bcftools stats -s - "$V2/chr$c.vcf.gz" | awk -F'\t' '$1=="PSC" {print $3"\t"$4+$5+$6+$14"\t"$14}'; done \
  | awk -F'\t' '{n[$1]+=$2; m[$1]+=$3} END {mn=2; for (s in n) {cr=1-m[s]/n[s]; if (cr<mn) {mn=cr; who=s}; sum+=cr; k++} printf "samples: %d | min call rate: %.4f (sample %s) | mean: %.4f\n", k, mn, who, sum/k}' | tee -a "$VAL"
echo "--- chunks of 10 Mb with < 3 variants / maximum per chunk" | tee "$RES/extra_checks_v2.txt"
for c in $(seq 1 22); do
  awk -v c="chr$c" -F'\t' '{k=int(($2-1)/10000000); n[k]++} END {for (k in n) if (n[k] < 3) printf "%s chunk %d-%d Mb: %d variants\n", c, k*10, (k+1)*10, n[k]}' "$W/chr$c.pos.tsv"
done | sort -V | tee -a "$RES/extra_checks_v2.txt"
python -I "${CODE_DIR}/03_imputation_TOPMed/04_export_validate.py" --dir "$W" --outdir "$RES" | tee -a "$VAL"
cp "$EX" "$RES/exclude_v2.tsv"
echo "### manifest"; cut -f1-5,7 "$V2/upload_manifest_v2.tsv"; du -sh "$V2"
