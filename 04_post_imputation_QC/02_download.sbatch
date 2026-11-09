#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/02_download.sbatch   (original step script: 38b_imputation_download.sbatch)
# Produces: Downloads the imputation results and checks their md5 sums.
# Run: sbatch 04_post_imputation_QC/02_download.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=imp_dl
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%j.out
# ---------------------------------------------------------------------------
# Step 38b - Download the full TOPMed imputation results (job with the v2 VCFs: 140 samples x 510,412 SNPs; panel TOPMed r3, hg38,
# Eagle phasing, rsq >= 0.3; AES256-encrypted zips) with the "download all" script the server generated (step 38a read it first: it
# only contains `echo` and `curl -L <url> -o <file>` lines, 48 files = 22 chr_N.zip + 22 chr_N.zip.md5 + qc_report.txt,
# quality-control.html, statistics/chunks-excluded.txt, statistics/snps-typed-only.txt), then checks every zip against its .md5.
# Output: PRS/imputation/topmed_imputation_v2/{chr_*.zip, chr_*.zip.md5, qc_report.txt, quality-control.html, statistics/, md5_check.tsv}
# The original downloads are never modified or deleted. Network job: `main`, 1 CPU / 2 GB. curl progress meters are filtered out of the log.
# ---------------------------------------------------------------------------
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
D="${PROJECT_DIR}/PRS/imputation/topmed_imputation_v2"
cd "$D"
echo "[$(date +%T)] start download (48 files)"
bash download_all.sh 2>&1 | tr '\r' '\n' | grep -E '^Downloading|curl: |rror' || true
echo "[$(date +%T)] download finished"
echo "### files"; ls -la | sed -n '1,60p'; du -sh .
echo "### md5 check (each chr_N.zip against chr_N.zip.md5)"
printf "file\tbytes\tmd5_expected\tmd5_computed\tresult\n" > md5_check.tsv
bad=0
for c in $(seq 1 22); do
  f="chr_$c.zip"
  [ -s "$f" ] && [ -s "$f.md5" ] || { printf "%s\tmissing\t-\t-\tMISSING\n" "$f" >> md5_check.tsv; bad=1; continue; }
  exp=$(awk '{print $1; exit}' "$f.md5" | tr 'A-F' 'a-f'); got=$(md5sum "$f" | cut -d' ' -f1)
  if [ "$exp" = "$got" ]; then r=OK; else r=MISMATCH; bad=1; fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$f" "$(stat -c %s "$f")" "$exp" "$got" "$r" >> md5_check.tsv
done
cat md5_check.tsv
for f in qc_report.txt quality-control.html statistics/chunks-excluded.txt statistics/snps-typed-only.txt; do [ -s "$f" ] || { echo "EMPTY/missing: $f"; bad=1; }; done
[ $bad -eq 0 ] && echo "ALL 22 ZIPS: MD5 OK" || { echo "SOME FILE FAILED - stop and tell the analyst"; exit 1; }
