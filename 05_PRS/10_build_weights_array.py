#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/10_build_weights_array.py   (original step script: 16_build_weights.py)
# Produces: Merges PRS-CS output and aligns the effect allele to the array alleles.
# Run: python 05_PRS/10_build_weights_array.py <args>   (or from a SLURM script; see the usage line below)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 16a - merge the per-chromosome PRS-CS outputs of one version into a single weight file
for `plink2 --score`, with the effect allele aligned to the alleles of the donor array.

PRS-CS output (one file per chromosome, no header): CHR SNP BP A1 A2 BETA
  BETA is the posterior-mean effect of allele A1 (A1 is the allele PRS-CS reports; see below).
plink2 --score needs the scored allele to be one of the two alleles of the variant in the genotype
file, written on the SAME strand as the genotypes. The array and the LD panel can use opposite
strands for the same SNP (e.g. array A/G vs panel T/C), so for every weight:
  * A1 is one of the array alleles                      -> kept as is
  * complement(A1) is one of the array alleles          -> A1 is replaced by its complement
    (strand flip; BETA unchanged because it is the effect of the same biological allele)
  * neither (allele mismatch / different SNP)           -> weight dropped and counted
Strand-palindromic SNPs (A/T, C/G) cannot be distinguished by this rule: they are checked
explicitly and counted (the HapMap3 PRS-CS panels contain none, so none are expected).
Only variants present in the array .bim (the final QC dataset) are written, i.e. V3 weights (whole
HapMap3 panel) are restricted to the SNPs the array actually has.

Inputs : --prefix <dir>/<VERSION> (files <prefix>_pst_eff_a1_b0.5_phi1e-04_chr<N>.txt), --bim, --out
Outputs: <out> (header "SNP A1 BETA") and <out>.report.txt
Standard library only.
"""
import argparse
import glob
import os

COMP = str.maketrans("ACGT", "TGCA")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--bim", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    alleles = {}
    with open(args.bim) as fh:
        for line in fh:
            f = line.split()
            alleles[f[1]] = (f[4], f[5])

    files = sorted(glob.glob(f"{args.prefix}_pst_eff_a1_b0.5_phi1e-04_chr*.txt"))
    chroms = sorted(int(f.rsplit("chr", 1)[1].split(".")[0]) for f in files)
    n_read = n_in_array = n_kept = n_flipped = n_mismatch = n_pal = n_zero = 0
    with open(args.out, "w") as out:
        out.write("SNP\tA1\tBETA\n")
        for fn in files:
            with open(fn) as fh:
                for line in fh:
                    f = line.split()
                    if len(f) < 6:
                        continue
                    n_read += 1
                    snp, a1, beta = f[1], f[3].upper(), f[5]
                    if snp not in alleles:
                        continue
                    n_in_array += 1
                    arr = {x.upper() for x in alleles[snp]}
                    if a1 in arr:
                        use = a1
                    elif a1.translate(COMP) in arr:
                        use = a1.translate(COMP)
                        n_flipped += 1
                    else:
                        n_mismatch += 1
                        continue
                    if len(arr) == 2 and arr in ({"A", "T"}, {"C", "G"}):
                        n_pal += 1
                    if float(beta) == 0.0:
                        n_zero += 1
                    out.write(f"{snp}\t{use}\t{beta}\n")
                    n_kept += 1
    rep = [f"chromosomes with output: {chroms} ({len(chroms)} files)",
           f"weights read: {n_read:,} | in the array .bim: {n_in_array:,} | written for scoring: {n_kept:,}",
           f"strand-flipped to array alleles: {n_flipped:,} | dropped for allele mismatch: {n_mismatch:,} | "
           f"palindromic among written: {n_pal:,} | weights exactly 0: {n_zero:,}"]
    print("\n".join(rep))
    with open(args.out + ".report.txt", "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
