#!/usr/bin/env bash
#SBATCH --job-name=01_qc --time=01:00:00 --mem=16G --cpus-per-task=4
# ---------------------------------------------------------------------
# STEP 1  Genotype QC of the GSA array data (all 140 women)
#   ped/map -> binary PLINK, sample & SNP QC, sex check,
#   heterozygosity, relatedness/duplicates, ancestry PCA.
#   Samples are only FLAGGED here; you decide what to remove (qc_flags.tsv).
# ---------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.sh"
D=$OUT/01_qc; mkdir -p "$D"; cd "$D"
 
# 1. text -> binary (much smaller/faster)
$PLINK --file "$RAW_PLINK" --make-bed --allow-extra-chr --out raw
 
# 2. keep chr 1-22 + X, drop SNPs with unknown position
$PLINK --bfile raw --chr 1-22,X --make-bed --out step1
 
# 3. sex check needs pseudo-autosomal region split (GRCh38)
$PLINK --bfile step1 --split-x b38 no-fail --make-bed --out step2
$PLINK --bfile step2 --check-sex --out sexcheck
 
# 4. missingness
$PLINK --bfile step2 --missing --out missing
 
# 5. SNP filters, then sample call rate
$PLINK --bfile step2 --geno 0.05 --maf "$MAF" --hwe 1e-6 --make-bed --out step3
$PLINK --bfile step3 --mind 0.05 --make-bed --out qc
 
# 6. LD-pruned autosomal set for het / kinship / PCA
$PLINK --bfile qc --chr 1-22 --indep-pairwise 200 50 0.2 --out prune
$PLINK --bfile qc --extract prune.prune.in --het --out het
$PLINK2 --bfile qc --extract prune.prune.in --make-king-table --king-table-filter 0.0884 --out king
$PLINK2 --bfile qc --extract prune.prune.in --pca 10 --out pca
 
# 7. collect flags (does not remove anything)
$PYTHON - <<'EOF'
import pandas as pd
f=[]
sx=pd.read_csv("sexcheck.sexcheck",sep=r"\s+")
for _,r in sx[sx.SNPSEX!=2].iterrows(): f.append((r.FID,r.IID,"sex_not_female",f"F={r.F:.2f}"))
mi=pd.read_csv("missing.imiss",sep=r"\s+")
for _,r in mi[mi.F_MISS>0.05].iterrows(): f.append((r.FID,r.IID,"call_rate<95%",f"{r.F_MISS:.3f}"))
h=pd.read_csv("het.het",sep=r"\s+"); m,s=h.F.mean(),h.F.std()
for _,r in h[(h.F-m).abs()>3*s].iterrows(): f.append((r.FID,r.IID,"heterozygosity_outlier",f"F={r.F:.3f}"))
try:
    k=pd.read_csv("king.kin0",sep=r"\s+")
    for _,r in k.iterrows():
        lab="duplicate" if r.KINSHIP>0.354 else ("1st_degree" if r.KINSHIP>0.177 else "2nd_degree")
        f.append((r["#FID1"],r.IID1,lab,f"with {r.IID2} kin={r.KINSHIP:.3f}"))
except FileNotFoundError: pass
pd.DataFrame(f,columns=["FID","IID","flag","detail"]).to_csv("qc_flags.tsv",sep="\t",index=False)
print(f"{len(f)} flags written to qc_flags.tsv")
EOF
 
echo "QC done. Review $D/qc_flags.tsv and $D/pca.eigenvec"
echo "To remove samples: put FID IID in remove.txt and run:"
echo "  $PLINK --bfile qc --remove remove.txt --make-bed --out qc_final"
echo "Otherwise:  cp qc.bed qc_final.bed; cp qc.bim qc_final.bim; cp qc.fam qc_final.fam"
