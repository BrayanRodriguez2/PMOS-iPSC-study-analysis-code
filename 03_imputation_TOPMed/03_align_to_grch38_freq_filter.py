#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/03_align_to_grch38_freq_filter.py   (original step script: 20_freq_filter.py)
# Produces: Strand/allele alignment to GRCh38 (bcftools +fixref, frequency-based strand resolution, frequency filter).
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 20 (python part B) - allele-frequency discordance of the ALIGNED variants against the reference.

Run after bcftools +fixref, i.e. when REF/ALT already match the GRCh38 FASTA (forward strand). For every variant
that has a record with the same allele pair in the 1000 Genomes GRCh38 sites file, the study frequency of the ALT
allele (INFO/AF computed over the 140 samples by bcftools +fill-tags) is compared with the EUR frequency of that
same allele (INFO/EUR_AF, 2 decimals); variants with |difference| > 0.2 are REMOVED (the usual pre-imputation
threshold). Reference = 1000G EUR, a SUBSTITUTE for the TOPMed frequencies, which are not downloadable by script.
Not removed (only counted): variants without a reference record ("not assessed") and variants whose allele pair at
that position differs from the reference one ("allele pair differs": e.g. another ALT allele at the same site).
Inputs : --study (TSV from bcftools query: CHROM POS ID REF ALT AF; chromosome names with "chr"),
         --ref-sites (TSV: CHROM POS REF ALT EUR_AF; names without "chr"), --outdir.
Outputs: freq_remove_sites.tsv (CHROM POS of the variants to remove, with the "chr" prefix, for bcftools -T ^),
         freq_removed.tsv (details), freq_report.txt
Standard library only.
"""
import argparse
import os
from collections import Counter, defaultdict

D_MAX = 0.2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("study", "ref_sites", "outdir"):
        ap.add_argument("--" + a.replace("_", "-"), required=True)
    args = ap.parse_args()
    ref = defaultdict(list)
    with open(args.ref_sites) as fh:
        for line in fh:
            c, p, r, a, f = line.rstrip("\n").split("\t")
            if f in (".", ""):
                continue
            ref[(c.replace("chr", ""), int(p))].append((r, a, float(f)))
    st = Counter()
    removed = []
    diffs = []
    with open(args.study) as fh:
        for line in fh:
            c, p, vid, r, a, af = line.rstrip("\n").split("\t")
            st["variants"] += 1
            recs = ref.get((c.replace("chr", ""), int(p)))
            if not recs:
                st["not_assessed_no_reference_record"] += 1
                continue
            q = None
            for rr, ra, f in recs:
                if (rr, ra) == (r, a):
                    q = f
                    break
                if (rr, ra) == (a, r):
                    q = 1 - f
                    break
            if q is None:
                st["not_assessed_allele_pair_differs"] += 1
                continue
            st["assessed"] += 1
            d = abs(float(af) - q)
            diffs.append(d)
            if d > D_MAX:
                removed.append((c, int(p), vid, r, a, float(af), q, d))
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "freq_remove_sites.tsv"), "w") as out:
        for c, p, *_ in removed:
            out.write(f"{c}\t{p}\n")
    with open(os.path.join(args.outdir, "freq_removed.tsv"), "w") as out:
        out.write("chr\tpos\tid\tref\talt\tstudy_alt_freq\tref_freq_same_allele\tabs_diff\n")
        for r in removed:
            out.write("\t".join(f"{x:.3f}" if isinstance(x, float) else str(x) for x in r) + "\n")
    diffs.sort()
    med = diffs[len(diffs) // 2] if diffs else float("nan")
    rep = [f"aligned variants: {st['variants']:,} | assessed against the 1000G EUR frequency: {st['assessed']:,} | "
           f"no reference record: {st['not_assessed_no_reference_record']:,} | allele pair differs from the reference record: "
           f"{st['not_assessed_allele_pair_differs']:,}",
           f"|study - reference| frequency difference among assessed: median {med:.3f} | 95th percentile "
           f"{diffs[int(0.95 * len(diffs))] if diffs else float('nan'):.3f} | REMOVED (> {D_MAX}): {len(removed):,}"]
    print("\n".join(rep))
    with open(os.path.join(args.outdir, "freq_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
