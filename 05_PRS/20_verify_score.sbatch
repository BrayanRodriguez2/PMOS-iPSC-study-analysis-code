#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 05_PRS/20_verify_score.sbatch   (original step script: 51c_verify_score.sbatch)
# Produces: Independent re-computation of the scores from the VCF dosages.
# Run: sbatch 05_PRS/20_verify_score.sbatch   (after: source config/config.sh; mkdir -p logs)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
#SBATCH --job-name=prs_verify
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# Step 51c - independent verification of the plink2 scoring against the DS dosages of the TOPMed VCF (not the PLINK2 file).
#   MODE=grs29 : the 29 GRS variants, queried one by one in the VCF of their chromosome (python re-computation of the total GRS)
#   MODE=chr22 : all V1 weights of chromosome 22 (see 05_PRS/20_verify_score.py); needs the set V1_chr22 scored in step 51
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
MODE="${MODE:-grs29}"
PRS="${PROJECT_DIR}/PRS"; C="$PRS/imputation/clean"; D="$PRS/imputation/topmed_imputation_v2/unzipped"; SC="${CODE_DIR}"
W="$PRS/prscs/weights_imputed"; S="$PRS/prscs/scores_imputed"
eval "${CONDA_SETUP:-true}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
awk 'NR>1 {print $2}' "$C/keep_135.txt" > "$S/samples_135.txt"
if [ "$MODE" = "chr22" ]; then
  echo "### V1 chromosome 22: every weight of the chromosome recomputed from the VCF dosages"
  bcftools query -S "$S/samples_135.txt" -f '%POS\t%REF\t%ALT[\t%SAMPLE,%DS]\n' "$D/chr22.dose.vcf.gz" | python -I "$SC/05_PRS/20_verify_score.py" 22 "$W/V1_chr22_weights.txt" "$S/V1_chr22.sscore"
else
  echo "### GRS29: python re-computation from the VCF dosages of the 29 variants"
  python -I - <<'PY'
import os, subprocess
import pandas as pd
H = os.environ["PROJECT_DIR"] + "/PRS"
w = pd.read_csv(f"{H}/prscs/weights_imputed/GRS29_weights.txt", sep="\t")
samples = [l.strip() for l in open(f"{H}/prscs/scores_imputed/samples_135.txt")]
tot = {s: 0.0 for s in samples}; n = 0
for r in w.itertuples():
    c, pos, ref, alt = r.ID.split(":")
    out = subprocess.run(["bcftools", "query", "-S", f"{H}/prscs/scores_imputed/samples_135.txt", "-r", f"chr{c}:{pos}-{pos}", "-f", "%POS\t%REF\t%ALT[\t%SAMPLE,%DS]\n",
                          f"{H}/imputation/topmed_imputation_v2/unzipped/chr{c}.dose.vcf.gz"], capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        f = line.split("\t")
        if f[0] == pos and f[1] == ref and f[2] == alt:
            n += 1
            for cell in f[3:]:
                s, ds = cell.rsplit(",", 1); d = float(ds)
                tot[s] += r.BETA * (d if r.A1 == alt else 2.0 - d)
ps = pd.read_csv(f"{H}/prscs/scores_imputed/GRS29.sscore", sep="\t"); ps.columns = [x.lstrip("#") for x in ps.columns]
ps = ps.set_index("IID"); ps.index = ps.index.astype(str); ps["python_sum"] = pd.Series(tot)
d = (ps.SCORE1_SUM - ps.python_sum).abs()
print(f"GRS29: variants found in the VCF {n} of {len(w)} | samples {int(ps.python_sum.notna().sum())} | max |plink2 - python| = {d.max():.2e} | r = {ps.SCORE1_SUM.corr(ps.python_sum):.8f}")
print("GRS29 VERIFICATION", "PASSED" if (n == len(w) and d.max() < 1e-3) else "FAILED")
PY
fi
