#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/17_er2_missing_variants.py   (original step script: 71_er2_missing7.py)
# Produces: Why a few kept variants have no ER2 (absent from the server output).
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 71 - Why do 7 variants of the 75 / 358 subsets have no ER2?  (follow-up of the post-imputation QC; READ-ONLY, nothing is modified)

Context: step 41 reported "of the 75: 72 found, of the 358: 354 found" because it only keeps TYPE == "Typed" rows with a
non-missing ER2 and matches the subsets by the exact key CHROM|POS|REF|ALT. This script re-matches all 75 + 358 variants and
classifies the ones that were not found, to tell apart four possible causes:
  - typed_but_ER2_missing         : exact allele match, TYPE == Typed, ER2 empty       (ER2 cannot be evaluated)
  - exact_allele_match_not_typed  : exact allele match, but TYPE != Typed              (imputed only, never genotyped by the server)
  - same_position_other_alleles   : the position exists in info.gz with other REF/ALT   (allele matching problem, e.g. swapped/split)
  - position_absent_from_info     : the position is not in the TOPMed output           (dropped by the server QC / not in the panel)
For every problem variant it also checks: was it in the v2 upload VCF, is it in the clean variant table (after / before the
ER2 < 0.5 filter), and was it removed by the ER2 filter.

Inputs : info_tables/chr*.tsv.gz (step 40), keep_fin_variants_88.tsv (the 75 kept) and compare_1000g_446.tsv (358 = |diff vs TOPMed| <= 0.2),
         imputation/clean/variant_table{,_pre_er2filter}.tsv.gz, removed_er2_lt_0.5_ids.txt, topmed/upload_v2/chr*.vcf.gz
Outputs: <outdir>/er2_missing_variants.tsv, er2_status_counts.tsv, er2_missing_report.txt
"""
import argparse
import gzip
import os

import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for name in ["keep75", "compare446", "info-dir", "clean-table", "pre-table", "removed-ids", "vcf-dir", "outdir"]:
    ap.add_argument("--" + name, required=True)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
LOG = []


def say(msg):
    print(msg, flush=True)
    LOG.append(msg)


def norm_chr(s):
    """'chr12' -> '12' so that info.gz (chrN), the variant tables (N) and the VCFs can be compared."""
    return s.astype(str).str.replace(r"^chr", "", regex=True)


# ---- 1. Rebuild the two subsets exactly as step 41 does ------------------------------------------------------------
k75 = pd.read_csv(a.keep75, sep="\t", dtype={"CHROM": str})
c446 = pd.read_csv(a.compare446, sep="\t", dtype={"CHROM": str})
over = c446["abs_diff_topmed_gt_0_2"].astype(str).str.lower().eq("true")  # True = |diff vs TOPMed| > 0.2 (the 88)
k358 = c446[~over]
cols = ["CHROM", "POS", "REF", "ALT"]
t = pd.concat([k75[cols].assign(subset="75 kept (FIN-supported)"),
               k358[cols].assign(subset="358 discordant kept")], ignore_index=True)
t["chrom"] = norm_chr(t.CHROM)
t["POS"] = t.POS.astype(int)
say(f"subset sizes: 75 -> {int((t.subset.str.startswith('75')).sum())} | 358 -> {int((t.subset.str.startswith('358')).sum())}")

# ---- 2. Look the variants up in the TOPMed info tables (by position, then by alleles) -------------------------------
hits = []
for c, g in t.groupby("chrom"):
    f = os.path.join(a.info_dir, f"chr{c}.tsv.gz")
    if not os.path.exists(f):
        say(f"WARNING: missing info table {f}")
        continue
    d = pd.read_csv(f, sep="\t", dtype={"CHROM": str}, usecols=["CHROM", "POS", "ID", "REF", "ALT", "R2", "ER2", "TYPE"])
    d = d[d.POS.isin(set(g.POS))].copy()
    d["chrom"] = c
    hits.append(d)
hits = pd.concat(hits, ignore_index=True)

recs = []
for r in t.itertuples(index=False):
    h = hits[(hits.chrom == r.chrom) & (hits.POS == r.POS)]
    ex = h[(h.REF == r.REF) & (h.ALT == r.ALT)]
    info_txt = "; ".join(f"{x.ID}|{x.REF}>{x.ALT}|{x.TYPE}|R2={x.R2}|ER2={x.ER2}" for x in h.itertuples())
    if len(ex):
        e = ex.iloc[0]
        if e.TYPE == "Typed" and pd.notna(e.ER2):
            status = "OK_typed_with_ER2"
        elif e.TYPE == "Typed":
            status = "typed_but_ER2_missing"
        else:
            status = "exact_allele_match_not_typed"
    elif len(h):
        status = "same_position_other_alleles"
    else:
        status = "position_absent_from_info"
    recs.append({"subset": r.subset, "chrom": r.chrom, "POS": r.POS, "REF": r.REF, "ALT": r.ALT,
                 "status": status, "info_rows_at_position": info_txt})
res = pd.DataFrame(recs)
counts = res.groupby(["subset", "status"]).size().rename("n").reset_index()
counts.to_csv(os.path.join(a.outdir, "er2_status_counts.tsv"), sep="\t", index=False)
say("status counts (OK_typed_with_ER2 must reproduce 72 and 354):\n" + counts.to_string(index=False))

bad = res[res.status != "OK_typed_with_ER2"].copy()
say(f"variants without ER2 in the step-41 sense: {len(bad)}")

# ---- 3. Problem variants: were they uploaded to TOPMed (v2 VCFs)? ---------------------------------------------------
def scan_vcf(c, positions):
    out = {}
    f = os.path.join(a.vcf_dir, f"chr{c}.vcf.gz")
    if not os.path.exists(f):
        return out
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t", 5)
            if int(p[1]) in positions:
                out.setdefault(int(p[1]), []).append(f"{p[2]}|{p[3]}>{p[4]}")
    return out


bad["in_v2_upload_vcf"] = ""
for c, g in bad.groupby("chrom"):
    found = scan_vcf(c, set(g.POS))
    for i, r in g.iterrows():
        bad.loc[i, "in_v2_upload_vcf"] = "; ".join(found.get(r.POS, [])) or "NO"

# ---- 4. Problem variants: are they in the clean tables, and were they removed by the ER2 < 0.5 filter? ---------------
removed = set(open(a.removed_ids).read().split())


def scan_table(path, label):
    found = {}
    pos_by_chr = {c: set(g.POS) for c, g in bad.groupby("chrom")}
    for ch in pd.read_csv(path, sep="\t", dtype={"CHROM": str}, usecols=["CHROM", "POS", "ID", "REF", "ALT", "R2", "ER2", "TYPE"],
                          chunksize=1_000_000):
        ch["chrom"] = norm_chr(ch.CHROM)
        ch = ch[[p in pos_by_chr.get(c, ()) for c, p in zip(ch.chrom, ch.POS)]]
        for x in ch.itertuples():
            found.setdefault((x.chrom, int(x.POS)), []).append(f"{x.ID}|{x.REF}>{x.ALT}|{x.TYPE}|ER2={x.ER2}")
    bad[label] = [("; ".join(found.get((r.chrom, r.POS), [])) or "NO") for r in bad.itertuples()]
    return found


scan_table(a.pre_table, "in_clean_table_before_er2_filter")
scan_table(a.clean_table, "in_clean_table_now")
bad["in_removed_er2_list"] = ["YES" if any(i.split("|")[0] in removed for i in v.split("; ")) else "no"
                              for v in bad.in_clean_table_before_er2_filter]

bad.to_csv(os.path.join(a.outdir, "er2_missing_variants.tsv"), sep="\t", index=False)
say("\nProblem variants:")
for r in bad.itertuples():
    say(f" {r.subset} | chr{r.chrom}:{r.POS} {r.REF}>{r.ALT} | status={r.status}\n"
        f"    info.gz rows at position : {r.info_rows_at_position or 'NONE'}\n"
        f"    v2 upload VCF            : {r.in_v2_upload_vcf}\n"
        f"    clean table (pre filter) : {r.in_clean_table_before_er2_filter}\n"
        f"    clean table (now)        : {r.in_clean_table_now}\n"
        f"    removed by ER2<0.5 list  : {r.in_removed_er2_list}")
with open(os.path.join(a.outdir, "er2_missing_report.txt"), "w") as fh:
    fh.write("\n".join(LOG) + "\n")
