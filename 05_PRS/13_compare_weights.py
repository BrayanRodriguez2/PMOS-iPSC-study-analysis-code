#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/13_compare_weights.py   (original step script: 17_compare_weights.py)
# Produces: Comparison of the weights between PRS-CS versions.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 17 - Why do V1, V2 and V3 correlate only 0.71-0.87? Compare the WEIGHTS directly.

The three weight files (prscs/weights/V*_weights_array.txt: SNP A1 BETA, already aligned to the array
alleles) refer to (almost) the same 169 k array SNPs. For the SNPs present in all three it reports
  * Pearson / Spearman correlation of the weights (sign-aware, all aligned to the array allele),
  * share of SNPs with the same sign,
  * the variance of the PRS contributed by the 29 lead loci regions is NOT computed here;
    instead the concentration of the weights is shown: fraction of sum(beta^2) carried by the top 1 %
    and top 10 % SNPs in each version (a PRS dominated by a few large weights is more fragile).
Outputs (--out): weight_comparison.txt. Standard library + numpy.
"""
import argparse

import numpy as np


def load(path):
    d = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            s, a, b = line.split()
            d[s] = (a, float(b))
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--weightdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    w = {v: load(f"{args.weightdir}/{v}_weights_array.txt") for v in ("V1", "V2", "V3")}
    rep = []

    def log(m=""):
        print(m)
        rep.append(m)

    for v in w:
        log(f"[{v}] SNPs with a weight: {len(w[v]):,}")
    common = sorted(set(w["V1"]) & set(w["V2"]) & set(w["V3"]))
    log(f"[common] SNPs with a weight in V1, V2 and V3: {len(common):,}")
    a1_diff = sum(1 for s in common if not (w["V1"][s][0] == w["V2"][s][0] == w["V3"][s][0]))
    log(f"[alleles] SNPs where the scored allele differs between versions: {a1_diff:,} (weights are comparable only if 0; "
        f"otherwise the sign of one of them is flipped below)")
    X = {}
    for v in w:
        X[v] = np.array([w[v][s][1] for s in common])
    # sign-align to V1's scored allele when the allele differs (beta of the other allele = -beta)
    for v in ("V2", "V3"):
        flip = np.array([1.0 if w[v][s][0] == w["V1"][s][0] else -1.0 for s in common])
        X[v] = X[v] * flip
    from scipy import stats
    for a, b in (("V1", "V2"), ("V1", "V3"), ("V2", "V3")):
        r = stats.pearsonr(X[a], X[b])[0]
        rho = stats.spearmanr(X[a], X[b])[0]
        same = float(np.mean(np.sign(X[a]) == np.sign(X[b])))
        slope = float(np.polyfit(X[a], X[b], 1)[0])
        log(f"[weights {a}-{b}] Pearson r = {r:.3f} | Spearman rho = {rho:.3f} | same sign = {100 * same:.1f}% | slope {b}~{a} = {slope:.3f}")
    for v in w:
        b2 = np.sort(X[v] ** 2)[::-1]
        tot = b2.sum()
        n = len(b2)
        log(f"[concentration {v}] sum(beta^2) = {tot:.4g} | top 1% SNPs carry {100 * b2[:n // 100].sum() / tot:.1f}% | "
            f"top 10% carry {100 * b2[:n // 10].sum() / tot:.1f}% | max |beta| = {np.sqrt(b2[0]):.3g} | "
            f"median |beta| = {np.median(np.abs(X[v])):.3g}")
    with open(args.out, "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
