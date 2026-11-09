#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/03_align_to_grch38_resolve_strand.py   (original step script: 20_resolve_strand.py)
# Produces: Strand/allele alignment to GRCh38 (bcftools +fixref, frequency-based strand resolution, frequency filter).
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 20 (python part A) - strand of the strand-AMBIGUOUS (palindromic A/T and C/G) SNPs.

Why: for A/T and C/G SNPs the alleles look identical on both strands, so neither the FASTA check of
bcftools +fixref nor the TOPMed server (which ignores A/T and C/G in its strand test) can tell whether the array
reports them on the forward or on the reverse strand. The only information left is the allele FREQUENCY.
Reference of frequencies: 1000 Genomes phase 3 GRCh38, EUR super-population (INFO/EUR_AF of the sites-only VCF).
SUBSTITUTE for TOPMed: TOPMed/Bravo frequencies are not downloadable by script (HTTP 403); EUR_AF has only 2 decimals.

Rules, per ambiguous variant (study allele frequency from plink2 --freq on all 140 samples):
  1. study MAF > 0.4        -> REMOVED. Orientation cannot be told apart by frequency when p ~ 0.5.
  2. no record of the 1000G sites file at that chr:pos with the same allele pair -> REMOVED (strand unresolvable).
  3. otherwise compare the study frequency p of its ALT allele with the reference frequency q of the same allele:
        d_same = |p - q|     (array already on the forward strand)
        d_flip = |p - (1-q)| (array on the reverse strand: its ALT is the other base on the forward strand)
     orientation = the smaller one; if that distance is > 0.2 -> REMOVED (discordant frequency).
     If orientation is "flip" the variant goes to the strand-flip list (plink2 --flip complements its alleles).
Non-ambiguous SNPs are not touched here: their strand follows unambiguously from the reference FASTA and is
fixed by bcftools +fixref.

Inputs : --bim (SNP-QC dataset), --afreq (plink2 .afreq of the same dataset), --ref-sites (TSV: CHROM POS REF ALT EUR_AF,
         chromosome names without "chr"), --outdir.
Outputs: palindromic_remove.txt, palindromic_flip.txt, palindromic_decisions.tsv, palindromic_report.txt
Standard library only.
"""
import argparse
import os
from collections import Counter, defaultdict

PAL = ({"A", "T"}, {"C", "G"})
MAF_MAX = 0.4
D_MAX = 0.2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("bim", "afreq", "ref_sites", "outdir"):
        ap.add_argument("--" + a.replace("_", "-"), required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    ref = defaultdict(list)
    with open(args.ref_sites) as fh:
        for line in fh:
            c, p, r, a, f = line.rstrip("\n").split("\t")
            if f in (".", ""):
                continue
            ref[(c.replace("chr", ""), int(p))].append((r, a, float(f)))
    freq = {}
    with open(args.afreq) as fh:
        h = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for line in fh:
            r = dict(zip(h, line.rstrip("\n").split("\t")))
            freq[r["ID"]] = (r["REF"], r["ALT"], float(r["ALT_FREQS"]))
    stats = Counter()
    remove, flip, rows = [], [], []
    with open(args.bim) as fh:
        for line in fh:
            f = line.split()
            chrom, vid, pos = f[0], f[1], int(f[3])
            if vid not in freq:
                stats["no_frequency"] += 1
                continue
            r0, a0, p = freq[vid]
            if {r0, a0} not in PAL:
                stats["non_ambiguous"] += 1
                continue
            stats["ambiguous"] += 1
            maf = min(p, 1 - p)
            if maf > MAF_MAX:
                stats["removed_maf_gt_0.4"] += 1
                remove.append(vid)
                rows.append((vid, chrom, pos, r0, a0, p, "", "", "", "", "remove_maf_gt_0.4"))
                continue
            hit = [x for x in ref.get((chrom, pos), []) if {x[0], x[1]} == {r0, a0}]
            if not hit:
                stats["removed_not_in_reference"] += 1
                remove.append(vid)
                rows.append((vid, chrom, pos, r0, a0, p, "", "", "", "", "remove_not_in_reference"))
                continue
            rr, ra, q_ra = hit[0]
            q = q_ra if ra == a0 else 1 - q_ra  # reference frequency of the study ALT allele (forward strand)
            d_same, d_flip = abs(p - q), abs(p - (1 - q))
            best, dec = (d_same, "same") if d_same <= d_flip else (d_flip, "flip")
            if best > D_MAX:
                stats["removed_freq_discordant"] += 1
                remove.append(vid)
                rows.append((vid, chrom, pos, r0, a0, p, q, d_same, d_flip, best, "remove_freq_discordant"))
                continue
            stats["kept_" + dec] += 1
            if dec == "flip":
                flip.append(vid)
            rows.append((vid, chrom, pos, r0, a0, p, q, d_same, d_flip, best, "keep_" + dec))
    with open(os.path.join(args.outdir, "palindromic_remove.txt"), "w") as out:
        out.write("\n".join(remove) + ("\n" if remove else ""))
    with open(os.path.join(args.outdir, "palindromic_flip.txt"), "w") as out:
        out.write("\n".join(flip) + ("\n" if flip else ""))
    with open(os.path.join(args.outdir, "palindromic_decisions.tsv"), "w") as out:
        out.write("id\tchr\tpos\tstudy_ref\tstudy_alt\tstudy_alt_freq\tref_freq_same_allele\td_same\td_flip\tbest_d\tdecision\n")
        for r in rows:
            out.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
    rep = [f"variants with a frequency: {stats['non_ambiguous'] + stats['ambiguous']:,} | strand-ambiguous (A/T, C/G): {stats['ambiguous']:,} "
           f"({100 * stats['ambiguous'] / max(stats['non_ambiguous'] + stats['ambiguous'], 1):.1f} %) | non-ambiguous: {stats['non_ambiguous']:,}",
           f"ambiguous removed: MAF > {MAF_MAX}: {stats['removed_maf_gt_0.4']:,} | not in the 1000G reference: {stats['removed_not_in_reference']:,} | "
           f"frequency discordant (best distance > {D_MAX}): {stats['removed_freq_discordant']:,}",
           f"ambiguous kept: array on forward strand ('same'): {stats['kept_same']:,} | on reverse strand ('flip', alleles complemented): {stats['kept_flip']:,}"]
    print("\n".join(rep))
    with open(os.path.join(args.outdir, "palindromic_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
