#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/02_snp_qc.py   (original step script: 19_dup_positions.py)
# Produces: SNP QC for imputation (call rate, monomorphic, duplicated positions).
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 19 (python part) - variants that share a genomic position ("sin duplicados de posicion").

The TOPMed Imputation Server removes duplicated positions itself, but it is cleaner (and countable) to do it
before. For every chr:position carried by more than one variant of the SNP-QC dataset:
  * all copies have the SAME allele set (redundant probes for the same SNP): keep ONE, the one with the best
    call rate (ties: the first in the file), drop the others;
  * the copies have DIFFERENT allele sets (e.g. different alternative alleles at the same site): which probe is
    right cannot be decided here, so ALL copies are dropped.
Inputs : --pvar (plink2 .pvar of the dataset), --vmiss (plink2 .vmiss: per-variant missingness), --out (list of IDs
         to remove, one per line) ; also writes <out>.report.txt.
Standard library only.
"""
import argparse
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pvar", required=True)
    ap.add_argument("--vmiss", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    fmiss = {}
    with open(args.vmiss) as fh:
        h = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for line in fh:
            r = dict(zip(h, line.rstrip("\n").split("\t")))
            fmiss[r["ID"]] = float(r["F_MISS"])
    groups = defaultdict(list)
    n = 0
    with open(args.pvar) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            groups[(f[0], int(f[1]))].append((f[2], frozenset((f[3], f[4])), n))
            n += 1
    remove, pos_dup, same_set, diff_set = [], 0, 0, 0
    drop_same = drop_diff = 0
    for key, g in groups.items():
        if len(g) < 2:
            continue
        pos_dup += 1
        if len({s for _, s, _ in g}) == 1:
            same_set += 1
            g.sort(key=lambda x: (fmiss[x[0]], x[2]))
            for vid, _, _ in g[1:]:
                remove.append(vid)
                drop_same += 1
        else:
            diff_set += 1
            for vid, _, _ in g:
                remove.append(vid)
                drop_diff += 1
    with open(args.out, "w") as out:
        out.write("\n".join(remove) + ("\n" if remove else ""))
    rep = [f"variants in dataset: {n:,} | positions carried by more than one variant: {pos_dup:,} "
           f"(same allele set: {same_set:,}; different allele sets: {diff_set:,})",
           f"variants removed: {len(remove):,} = {drop_same:,} redundant copies (best call rate kept) + {drop_diff:,} in conflicting positions (all copies dropped)"]
    print("\n".join(rep))
    with open(args.out + ".report.txt", "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
