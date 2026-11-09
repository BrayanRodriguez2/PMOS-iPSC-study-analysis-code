#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/20_verify_score.py   (original step script: 51c_verify_score.py)
# Produces: Independent re-computation of the scores from the VCF dosages.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 51c - independent check of the plink2 --score results: the score of a weight set is recomputed in python from the DS (dosage) field of the TOPMed
VCF of one chromosome (reading the VCF with bcftools query, NOT the PLINK2 file) and compared with SCORE1_SUM of plink2 for the same samples.
Reads from stdin the output of:  bcftools query -S <samples> -f '%POS\\t%REF\\t%ALT[\\t%SAMPLE,%DS]\\n' chr<N>.dose.vcf.gz
Usage: python 05_PRS/20_verify_score.py <chrom> <weights.txt (ID A1 BETA, header)> <plink2 .sscore> < query_output
The weight file must hold only variants of <chrom>. Score = sum(beta x dosage of the scored allele); dosage of the scored allele = DS if scored = ALT else 2 - DS.
PASS when every sample agrees within 1e-3 (plink2 stores dosages with ~1e-4 resolution) and all weight variants were found in the VCF.
"""
import sys
import numpy as np
import pandas as pd

chrom, wfile, sfile = sys.argv[1], sys.argv[2], sys.argv[3]
w = pd.read_csv(wfile, sep="\t")
wd = {r.ID: (r.A1, float(r.BETA)) for r in w.itertuples()}
ps = pd.read_csv(sfile, sep="\t")
ps.columns = [c.lstrip("#") for c in ps.columns]
score = {}
found = 0
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    vid = f"{chrom}:{f[0]}:{f[1]}:{f[2]}"
    if vid not in wd:
        continue
    a1, beta = wd[vid]
    found += 1
    for cell in f[3:]:
        s, ds = cell.rsplit(",", 1)
        d = float(ds)
        score[s] = score.get(s, 0.0) + beta * (d if a1 == f[2] else 2.0 - d)
mine = pd.Series(score, name="python_sum")
cmp_ = ps.assign(IID=ps.IID.astype(str)).set_index("IID")[["SCORE1_SUM"]].join(mine)
diff = (cmp_.SCORE1_SUM - cmp_.python_sum).abs()
print(f"weight variants: {len(wd):,} | found in the VCF: {found:,} | samples compared: {int(cmp_.python_sum.notna().sum())} of {len(cmp_)}")
print(f"max |plink2 - python| = {diff.max():.2e} | correlation = {cmp_.SCORE1_SUM.corr(cmp_.python_sum):.8f} | score range plink2: {cmp_.SCORE1_SUM.min():.4f} .. {cmp_.SCORE1_SUM.max():.4f}")
ok = found == len(wd) and cmp_.python_sum.notna().all() and diff.max() < 1e-3
print("VERIFICATION", "PASSED" if ok else "FAILED")
sys.exit(0 if ok else 1)
