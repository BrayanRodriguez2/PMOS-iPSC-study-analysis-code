#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/03_count_overlap.py   (original step script: 04_count_overlap.py)
# Produces: Overlap of array SNPs with the PRS-CS panel (by rsID, with allele checks).
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 04 - rsID overlap between GWAS summary statistics, the PRS-CS LD panel and
the donor genotyping array, plus rsID recovery for unlabeled array SNPs.

Inputs
  --sumstats  GWAS Catalog harmonised file (GCST90570579.h.tsv.gz, GRCh38, GWAS-SSF)
              columns used: chromosome, base_pair_location, effect_allele,
              other_allele, rsid
  --snpinfo   PRS-CS LD reference SNP list (snpinfo_1kg_hm3): CHR SNP BP A1 A2 MAF
  --bim       PLINK .bim of the donor array (GRCh38; col 4 = GRCh38 position).
              The .bim is only READ; it is never modified.

What it computes
  (a) = rsIDs in the summary statistics
  (b) = rsIDs in the LD panel (snpinfo)
  (c) = rsIDs in the .bim (variant IDs of the form rs<digits>)
  Sizes of a, b, c and the intersections a&b, a&c, b&c, a&b&c.

  Recovery: .bim variants WITHOUT rsID (e.g. "1:103380393") are looked up in the
  summary statistics by chromosome + GRCh38 position + alleles. Alleles match in
  any order (A1/A2 vs effect/other) and on either strand (complement). Strand-
  palindromic SNPs (A/T, C/G) match trivially, so they are counted separately
  because their strand cannot be verified from alleles alone. If several different
  rsIDs match the same variant it is flagged "ambiguous" and NOT recovered.
  Then the a&b&c intersection is recomputed with the recovered rsIDs added to (c).

  Extra checks (not requested but needed to interpret the numbers):
   - build sanity check: for .bim rsIDs also present in the sumstats, how often
     the .bim position equals the sumstats GRCh38 position;
   - allele compatibility inside a&b&c (sumstats vs .bim vs panel alleles, any
     order/strand), because PRS-CS drops SNPs whose alleles cannot be reconciled.

Memory: keeps the set of sumstats rsIDs (~10 M strings) -> ask for >= 16 GB.
Only the Python standard library is used.

Outputs (in --outdir)
  overlap_summary.tsv       table of list sizes and intersections
  overlap_report.txt        human-readable report with all counts
  recovered_rsids.tsv       one row per recovered .bim variant (bim ID -> rsID)
  incompatible_alleles.tsv  rsIDs in a&b&c whose alleles cannot be reconciled (added 2026-10-30)
"""
import argparse
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

RS_RE = re.compile(r"^rs\d+$")
COMP = str.maketrans("ACGT", "TGCA")
PALINDROMIC = {frozenset("AT"), frozenset("CG")}


def norm_chr(c):
    """Normalise chromosome labels: strip 'chr', map 23/24/26 to X/Y/MT."""
    c = str(c)
    if c.lower().startswith("chr"):
        c = c[3:]
    return {"23": "X", "24": "Y", "26": "MT", "M": "MT"}.get(c, c)


def complement(allele):
    return allele.upper().translate(COMP)


def is_acgt(allele):
    return bool(allele) and set(allele.upper()) <= set("ACGT")


def allele_match(pair_a, pair_b):
    """Compare two allele pairs ignoring order and strand.
    Returns 'direct', 'flip' or None. Palindromic pairs are reported as 'direct'
    (their strand is undeterminable; callers check palindromic status separately)."""
    sa = frozenset(x.upper() for x in pair_a)
    sb = frozenset(x.upper() for x in pair_b)
    if sa == sb:
        return "direct"
    if all(is_acgt(x) for x in pair_b) and sa == frozenset(complement(x) for x in pair_b):
        return "flip"
    return None


def is_palindromic(pair):
    return frozenset(x.upper() for x in pair) in PALINDROMIC


def read_snpinfo(path):
    """PRS-CS panel SNP list -> dict rsID -> (A1, A2). Header: CHR SNP BP A1 A2 MAF."""
    panel = {}
    n_lines = 0
    with open(path) as fh:
        header = fh.readline().split()
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.split()
            if len(f) < 5:
                continue
            n_lines += 1
            panel[f[idx["SNP"]]] = (f[idx["A1"]], f[idx["A2"]])
    return panel, n_lines, header


def read_bim(path):
    """PLINK .bim -> list of (chr, id, pos, a1, a2)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            rows.append((norm_chr(f[0]), f[1], int(f[3]), f[4], f[5]))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sumstats", required=True)
    ap.add_argument("--snpinfo", required=True)
    ap.add_argument("--bim", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rep = []  # report lines

    def log(msg=""):
        print(msg, flush=True)
        rep.append(msg)

    # ---- (b) LD panel ------------------------------------------------------
    panel, n_panel_lines, panel_header = read_snpinfo(args.snpinfo)
    B = {k for k in panel if RS_RE.match(k)}
    log(f"[panel] snpinfo lines: {n_panel_lines:,} | unique SNP IDs: {len(panel):,} | rs-IDs: {len(B):,} | header: {panel_header}")

    # ---- (c) array .bim ----------------------------------------------------
    bim = read_bim(args.bim)
    bim_rs = {}      # rsID -> (chr, pos, a1, a2)
    unlabeled = []   # variants whose ID is not rs<digits>
    n_rs_like_other = 0
    for chrom, vid, pos, a1, a2 in bim:
        if RS_RE.match(vid):
            bim_rs[vid] = (chrom, pos, a1, a2)
        else:
            unlabeled.append((chrom, vid, pos, a1, a2))
            if vid.lower().startswith("rs"):
                n_rs_like_other += 1
    C = set(bim_rs)
    log(f"[bim] variants: {len(bim):,} | with rs<digits> ID: {len(bim_rs):,} (unique {len(C):,}) | without: {len(unlabeled):,} | start with 'rs' but not rs<digits>: {n_rs_like_other:,}")
    want = defaultdict(list)  # (chr, pos) -> indices in `unlabeled`
    for i, (chrom, vid, pos, a1, a2) in enumerate(unlabeled):
        want[(chrom, pos)].append(i)

    # ---- (a) summary statistics (single streaming pass) ---------------------
    A = set()
    sst_alleles = defaultdict(list)  # rsID -> [(ea, oa)] only for rsIDs in the panel (for allele check)
    cand = defaultdict(list)         # (chr, pos) -> [(ea, oa, rsid)] only at positions of unlabeled variants
    pos_check = Counter()            # build sanity check on .bim rsIDs
    n_rows = n_no_rsid = 0
    with gzip.open(args.sumstats, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        need = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "rsid"]
        missing = [c for c in need if c not in ix]
        if missing:
            sys.exit(f"ERROR: sumstats missing columns {missing}; header = {header}")
        ic, ip, ie, io, ir = (ix[c] for c in need)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            n_rows += 1
            rsid = f[ir]
            chrom = norm_chr(f[ic])
            pos = int(f[ip])
            ea, oa = f[ie], f[io]
            if RS_RE.match(rsid):
                A.add(rsid)
                if rsid in panel:
                    sst_alleles[rsid].append((ea, oa))
                if rsid in bim_rs:
                    pos_check["same_position" if bim_rs[rsid][:2] == (chrom, pos) else "different_position"] += 1
            else:
                n_no_rsid += 1
            if (chrom, pos) in want:
                cand[(chrom, pos)].append((ea, oa, rsid))
    log(f"[sumstats] rows: {n_rows:,} | rows without valid rsID: {n_no_rsid:,} | unique rs-IDs: {len(A):,}")
    tot = sum(pos_check.values())
    log(f"[build check] .bim rsIDs found in sumstats: {tot:,} | same chr:pos: {pos_check['same_position']:,} "
        f"({100 * pos_check['same_position'] / max(tot, 1):.2f}%) | different: {pos_check['different_position']:,}")

    # ---- intersections by rsID ---------------------------------------------
    AB, AC, BC = A & B, A & C, B & C
    ABC = A & B & C
    log("")
    log("[intersections by rsID]")
    log(f"  |a| sumstats {len(A):,} | |b| panel {len(B):,} | |c| bim {len(C):,}")
    log(f"  a&b {len(AB):,} | a&c {len(AC):,} | b&c {len(BC):,} | a&b&c {len(ABC):,}")

    # ---- recovery of rsIDs for unlabeled .bim variants ---------------------
    status = Counter()
    recovered = []  # (bim_id, chr, pos, a1, a2, rsid, ea, oa, match, palindromic)
    for chrom, vid, pos, a1, a2 in unlabeled:
        if a1 in ("0", ".", "") or a2 in ("0", ".", ""):
            status["skipped_missing_allele_in_bim"] += 1
            continue
        hits = cand.get((chrom, pos))
        if not hits:
            status["no_variant_at_position_in_sumstats"] += 1
            continue
        matches = []
        for ea, oa, rsid in hits:
            m = allele_match((a1, a2), (ea, oa))
            if m and RS_RE.match(rsid):
                matches.append((rsid, ea, oa, m))
        if not matches:
            status["position_found_but_alleles_incompatible"] += 1
            continue
        rs_set = {m[0] for m in matches}
        if len(rs_set) > 1:
            status["ambiguous_several_rsids"] += 1
            continue
        rsid, ea, oa, m = matches[0]
        pal = is_palindromic((a1, a2))
        status["recovered_palindromic_strand_unverifiable" if pal else "recovered_nonpalindromic"] += 1
        recovered.append((vid, chrom, pos, a1, a2, rsid, ea, oa, m, "yes" if pal else "no"))
    log("")
    log("[recovery of rsIDs for .bim variants without rsID]")
    for k in ("recovered_nonpalindromic", "recovered_palindromic_strand_unverifiable",
              "no_variant_at_position_in_sumstats", "position_found_but_alleles_incompatible",
              "ambiguous_several_rsids", "skipped_missing_allele_in_bim"):
        log(f"  {k}: {status[k]:,}")
    log(f"  total unlabeled variants: {len(unlabeled):,}")

    R_all = {r[5] for r in recovered}
    R_np = {r[5] for r in recovered if r[9] == "no"}
    log(f"  recovered unique rsIDs (all): {len(R_all):,} | of which already in (c): {len(R_all & C):,} | NEW to (c): {len(R_all - C):,}")
    C_all = C | R_all
    C_np = C | R_np
    ABC_all = A & B & C_all
    ABC_np = A & B & C_np
    log(f"  new rsIDs that are also in the panel (b): {len((R_all - C) & B):,} (non-palindromic only: {len((R_np - C) & B):,})")
    log(f"  a&b&c after recovery (all recovered): {len(ABC_all):,} | (non-palindromic recovered only): {len(ABC_np):,} | before: {len(ABC):,}")

    # ---- allele compatibility inside a&b&c ---------------------------------
    # Map every rsID of the (c + recovered) set to the .bim alleles
    bim_alleles = {rs: (v[2], v[3]) for rs, v in bim_rs.items()}
    for r in recovered:
        bim_alleles.setdefault(r[5], (r[3], r[4]))

    def compat(rs_set, bad=None):
        """Count allele compatibility; if `bad` is a list, also collect the incompatible rsIDs
        (added 2026-10-30 so that step 10 can exclude them; counts are unchanged)."""
        ok = Counter()
        for rs in rs_set:
            if rs not in sst_alleles or rs not in bim_alleles:
                ok["missing_alleles"] += 1
                continue
            ba, pa = bim_alleles[rs], panel[rs]
            good = any(allele_match(ba, s) and allele_match(pa, s) for s in sst_alleles[rs])
            if good:
                ok["compatible_palindromic" if is_palindromic(ba) else "compatible_nonpalindromic"] += 1
            else:
                ok["alleles_incompatible"] += 1
                if bad is not None:
                    bad.append((rs, "/".join(ba), "/".join(pa), ";".join("/".join(s) for s in sst_alleles[rs])))
        return ok

    incompatible = []
    ok0, ok1 = compat(ABC), compat(ABC_all, incompatible)
    log("")
    log("[allele compatibility (sumstats vs .bim vs panel; any order/strand)]")
    for label, ok, n in (("before recovery", ok0, len(ABC)), ("after recovery ", ok1, len(ABC_all))):
        log(f"  {label}: a&b&c = {n:,} | compatible non-palindromic: {ok['compatible_nonpalindromic']:,} | "
            f"compatible palindromic: {ok['compatible_palindromic']:,} | incompatible: {ok['alleles_incompatible']:,} | missing: {ok['missing_alleles']:,}")

    # ---- write outputs ------------------------------------------------------
    with open(os.path.join(args.outdir, "overlap_summary.tsv"), "w") as out:
        out.write("list\tn_rsids\n")
        rows = [("a_sumstats", len(A)), ("b_ld_panel", len(B)), ("c_bim", len(C)),
                ("a_and_b", len(AB)), ("a_and_c", len(AC)), ("b_and_c", len(BC)), ("a_and_b_and_c", len(ABC)),
                ("c_bim_plus_recovered", len(C_all)), ("a_and_b_and_c_after_recovery", len(ABC_all)),
                ("a_and_b_and_c_after_recovery_nonpalindromic_only", len(ABC_np))]
        for k, v in rows:
            out.write(f"{k}\t{v}\n")
    with open(os.path.join(args.outdir, "incompatible_alleles.tsv"), "w") as out:
        out.write("rsid\tbim_alleles\tpanel_alleles\tsumstats_alleles_effect/other\n")
        for row in sorted(incompatible):
            out.write("\t".join(row) + "\n")
    with open(os.path.join(args.outdir, "recovered_rsids.tsv"), "w") as out:
        out.write("bim_id\tchr\tpos_grch38\tbim_a1\tbim_a2\trsid\tsumstats_effect_allele\tsumstats_other_allele\tmatch\tpalindromic\n")
        for r in recovered:
            out.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(args.outdir, "overlap_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")
    log(f"\nwritten to {args.outdir}")


if __name__ == "__main__":
    main()
