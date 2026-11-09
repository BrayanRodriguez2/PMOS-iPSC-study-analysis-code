#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/04_export_validate.py   (original step script: 21_vcf_checks.py)
# Produces: Per-chromosome VCF export and validation for the TOPMed server.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 21 (python part) - structural checks of the per-chromosome VCFs against the TOPMed Imputation Server rules
(docs: statgen.github.io/tis-v2-docs/imputation/data-preparation/ and .../pipeline-overview/).

Input : one TSV per chromosome made with  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'  (file name chr<N>.pos.tsv).
Checks per chromosome:
  * every record is on the expected contig name "chr<N>" (GRCh38 input needs the "chr" prefix);
  * positions strictly increasing (sorted, and therefore no duplicated position);
  * REF and ALT are single A/C/G/T bases and REF != ALT (no indels, no invalid alleles);
  * the server cuts each chromosome into 10 Mb chunks and refuses chunks with more than 20,000 variants (FAQ: it is
    meant for arrays) or fewer than 3 variants: counts per chunk are reported (the last, partial chunk of a chromosome
    may legitimately have few variants).
Outputs (--outdir): vcf_checks.tsv (one row per chromosome) and vcf_checks_report.txt.
Standard library only.
"""
import argparse
import os
from collections import Counter

CHUNK = 10_000_000


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="directory with chr<N>.pos.tsv")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    rows, rep = [], []
    for c in range(1, 23):
        path = os.path.join(args.dir, f"chr{c}.pos.tsv")
        n = bad_contig = unsorted = bad_allele = same = 0
        last = 0
        chunks = Counter()
        first = None
        with open(path) as fh:
            for line in fh:
                ch, pos, ref, alt = line.rstrip("\n").split("\t")
                pos = int(pos)
                n += 1
                first = pos if first is None else first
                if ch != f"chr{c}":
                    bad_contig += 1
                if pos <= last:
                    unsorted += 1
                last = pos
                if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                    bad_allele += 1
                if ref == alt:
                    same += 1
                chunks[(pos - 1) // CHUNK] += 1
        mx = max(chunks.values()) if chunks else 0
        full = [v for k, v in chunks.items()]
        rows.append((f"chr{c}", n, first, last, len(chunks), mx, min(full) if full else 0, bad_contig, unsorted, bad_allele, same))
    head = ["contig", "variants", "first_pos", "last_pos", "n_10Mb_chunks_with_variants", "max_variants_in_a_10Mb_chunk",
            "min_variants_in_a_10Mb_chunk", "wrong_contig_name", "not_strictly_increasing", "invalid_alleles", "ref_equals_alt"]
    with open(os.path.join(args.outdir, "vcf_checks.tsv"), "w") as out:
        out.write("\t".join(head) + "\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    tot = sum(r[1] for r in rows)
    ok = all(r[7] == r[8] == r[9] == r[10] == 0 for r in rows)
    over = [r[0] for r in rows if r[5] > 20000]
    rep.append(f"chromosomes: {len(rows)} | total variants: {tot:,} | structural problems (wrong contig / unsorted-or-duplicated / invalid alleles / REF=ALT): "
               f"{'NONE' if ok else 'FOUND - see vcf_checks.tsv'}")
    rep.append(f"max variants in any 10 Mb chunk: {max(r[5] for r in rows):,} (server limit 20,000) | chromosomes over the limit: {over or 'none'}")
    rep.append("smallest chunk (chunks with at least one variant): " + ", ".join(f"{r[0]}={r[6]}" for r in rows if r[6] < 100) or "none < 100")
    print("\n".join(rep))
    with open(os.path.join(args.outdir, "vcf_checks_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
