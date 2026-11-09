#!/usr/bin/env bash
#SBATCH --job-name=08a_dos --time=00:20:00 --mem=8G
# ---------------------------------------------------------------------
# STEP 8a  Genotype dosages of the short-listed SNPs for the 14 iPSC donors
#   SNP list = 'snp' column of the candidate-pair file (default: 6b short list).
#   Output: $OUT/08_targeted/dosages.raw  (PLINK --export A, counted allele in header)
# ---------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.sh"
D=$OUT/08_targeted; mkdir -p "$D"; cd "$D"
PAIRS=${1:-$OUT/06_i2qtl/i2qtl_best_per_locus.tsv}
 
# donors with an iPSC line (array_iid column of the sample key, hel_line not empty)
$PYTHON - "$SAMPLE_KEY" <<'EOF'
import sys, pandas as pd
k = pd.read_csv(sys.argv[1]); k = k[k.hel_line.fillna("").str.strip() != ""]
pd.DataFrame({"FID": k.array_iid, "IID": k.array_iid}).to_csv("ipsc_donors.txt", sep="\t", index=False, header=False)
print(len(k), "iPSC donors")
EOF
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="snp") c=i; next} {print $c}' "$PAIRS" | sort -u > snps.txt
 
$PLINK2 --pfile "$OUT/03_imputed/imputed" --keep ipsc_donors.txt --extract snps.txt \
        --export A --out dosages
echo "$(($(head -1 dosages.raw | wc -w) - 6)) SNPs exported for $(($(wc -l < dosages.raw) - 1)) donors"
