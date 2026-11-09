#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/17_build_weights_imputed.py   (original step script: 51_build_weights_imputed.py)
# Produces: PRS-CS weights returned to chr:pos:ref:alt IDs.
# Run: python 05_PRS/17_build_weights_imputed.py <args>   (or from a SLURM script; see the usage line below)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 51a - PRS-CS weights -> weights for the imputed clean file (adapted from 05_PRS/10_build_weights_array.py, which aligned the weights to the array .bim).

PRS-CS output (one file per chromosome, no header): CHR SNP BP A1 A2 BETA; SNP = rsID of the LD panel, A1 = effect allele in the PANEL's strand and BETA its
posterior-mean effect. The map of step 49 (rsID -> clean ID chr:pos:ref:alt, strand = direct | complement) returns every weight to the clean file:
  scored allele = A1 if A1 is REF or ALT of the clean variant, else complement(A1) if that is (otherwise the weight is dropped and counted). (The strand of the
  PRS-CS output alleles is that of the .bim given to PRS-CS: clean alleles for V1/V2, panel alleles for V3; a first version of this script assumed the panel strand
  for V1/V2 and wrongly dropped the 538 strand-complement SNPs; corrected before any result was used.)
Only PRS-CS weights whose rsID is in the map (panel SNP present in the clean file) are written: for V3 (bim = the whole panel) the weights of the panel SNPs
absent from the clean file are therefore left out (counted). Optional --min-r2 keeps only clean variants with server R2 >= threshold (sensitivity).
Output: <out> (header "ID A1 BETA"; plink2 --score columns 1 2 3) and <out>.report.txt. Standard library only.
"""
import argparse, glob, gzip, os

COMP = str.maketrans("ACGT", "TGCA")
ap = argparse.ArgumentParser()
ap.add_argument("--prefix", required=True); ap.add_argument("--map", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--min-r2", type=float, default=None); ap.add_argument("--label", default="")
a = ap.parse_args()

mp = {}
with gzip.open(a.map, "rt") as fh:
    h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        mp[f[ix["rsID"]]] = (f[ix["ID"]], f[ix["REF"]], f[ix["ALT"]], f[ix["strand"]], f[ix["panel_A1"]], f[ix["panel_A2"]], float(f[ix["R2"]]), f[ix["TYPE"]])
files = sorted(glob.glob(f"{a.prefix}_pst_eff_a1_b0.5_phi1e-04_chr*.txt"))
chroms = sorted(int(x.rsplit("chr", 1)[1].split(".")[0]) for x in files)
n_read = n_in_map = n_kept = n_mismatch = n_r2 = n_zero = n_comp = 0
ids = set()
with open(a.out, "w") as out:
    out.write("ID\tA1\tBETA\n")
    for fn in files:
        for line in open(fn):
            f = line.split()
            if len(f) < 6:
                continue
            n_read += 1
            snp, a1, beta = f[1], f[3].upper(), f[5]
            if snp not in mp:
                continue
            n_in_map += 1
            cid, ref, alt, strand, pa1, pa2, r2, typ = mp[snp]
            # The strand of the PRS-CS output alleles is the one of the .bim given to PRS-CS: the clean file's alleles for V1/V2 (A1 is REF or ALT as is) and the
            # panel's alleles for V3 (whole-panel bim; A1 may be the complement of REF/ALT). Strand-ambiguous (A/T, C/G) SNPs are not in the map, so the rule below is
            # unambiguous: A1 itself if it is REF or ALT, otherwise its complement if that one is.
            if a1 in (ref, alt):
                use = a1
            elif a1.translate(COMP) in (ref, alt):
                use = a1.translate(COMP); n_comp += 1
            else:
                n_mismatch += 1; continue
            if a.min_r2 is not None and r2 < a.min_r2:
                n_r2 += 1; continue
            if cid in ids:
                raise SystemExit(f"duplicated clean ID in the weights: {cid}")
            ids.add(cid)
            n_zero += float(beta) == 0.0
            out.write(f"{cid}\t{use}\t{beta}\n"); n_kept += 1
rep = [f"[{a.label}] chromosomes with output: {len(chroms)} ({chroms[0]}-{chroms[-1]})" if chroms else f"[{a.label}] NO PRS-CS OUTPUT FILES",
       f"[{a.label}] weights read from PRS-CS: {n_read:,} | with a clean-file ID (panel SNP present in the clean file): {n_in_map:,} | not in the clean file: {n_read - n_in_map:,}",
       f"[{a.label}] dropped for allele mismatch: {n_mismatch:,} | dropped for R2 < {a.min_r2}: {n_r2:,} | WRITTEN for scoring: {n_kept:,} (A1 replaced by its complement: {n_comp:,}; weights exactly 0: {n_zero:,})"]
print("\n".join(rep))
open(a.out + ".report.txt", "w").write("\n".join(rep) + "\n")
