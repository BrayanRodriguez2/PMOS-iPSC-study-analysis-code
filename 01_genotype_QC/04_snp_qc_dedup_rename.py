#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 01_genotype_QC/04_snp_qc_dedup_rename.py   (original step script: 10_dedup_rename.py)
# Produces: SNP QC, duplicate handling and renaming to rsIDs; QC PLINK datasets.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 10 (python part) - exclusions, de-duplication and rsID renaming of the QC'd variants.

Works on the .bim of the QC dataset (b38_qc: after sample QC, call rate, monomorphic, HWE, MAF,
autosomes). Nothing is written to the original files: it only produces ID lists / a rename map
that step 10 (bash) feeds to plink2.

Order of operations (each step produces its own keep-list so that the counts can be reported):
  A. Exclude variants whose rsID was RECOVERED by position+alleles (step 04) AND that are
     strand-palindromic (A/T, C/G): their strand cannot be verified from the alleles.
  B. Exclude variants whose rsID is in incompatible_alleles.tsv (the 508 rsIDs whose alleles
     in the .bim / LD panel / sumstats cannot be reconciled, from step 04).
  C. De-duplicate by FINAL rsID (plain rs ID or recovered rsID): when several array probes carry
     the same rsID keep the one with the best call rate (ties: plain rs ID before recovered, then the
     first in the .bim).
  D. De-duplicate by chr:position(GRCh38)+allele set (any ID, including variants without rsID):
     same rule.
  Finally every kept variant whose rsID was recovered is renamed to that rsID (rename_map.tsv).
Why dedup comes after recovery: duplicates by rsID only become visible once the recovered rsIDs
are assigned (e.g. "GSA-rs123" and "rs123" are the same SNP).

Inputs : --bim, --vmiss (plink2 .vmiss of the QC dataset), --recovered (recovered_rsids.tsv),
         --incompatible (incompatible_alleles.tsv), --outdir.
Outputs: keep_A_*.txt ... keep_D_*.txt (variant IDs of the QC dataset), rename_map.tsv (old<TAB>new),
         variant_decisions.tsv (every variant and what happened), dedup_report.txt.
Standard library only.
"""
import argparse
import os
import re
from collections import defaultdict

RS_RE = re.compile(r"^rs\d+$")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("bim", "vmiss", "recovered", "incompatible", "outdir"):
        ap.add_argument("--" + a, required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rep = []

    def log(m=""):
        print(m, flush=True)
        rep.append(m)

    # ---- load ----------------------------------------------------------------
    callrate = {}
    with open(args.vmiss) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            callrate[r["ID"]] = 1.0 - float(r["F_MISS"])
    variants = []  # dicts in .bim order
    with open(args.bim) as fh:
        for i, line in enumerate(fh):
            f = line.split()
            variants.append({"i": i, "chr": f[0], "id": f[1], "pos": int(f[3]), "a1": f[4], "a2": f[5]})
    recovered = {}
    with open(args.recovered) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            recovered[f[0]] = (f[5], f[9])  # bim_id -> (rsid, palindromic yes/no)
    incompatible = set()
    with open(args.incompatible) as fh:
        next(fh)
        for line in fh:
            incompatible.add(line.split("\t")[0])
    log(f"[input] QC variants: {len(variants):,} | recovered-rsID table: {len(recovered):,} | incompatible rsIDs: {len(incompatible):,}")

    # ---- final rsID per variant ------------------------------------------------
    for v in variants:
        v["call"] = callrate[v["id"]]
        if RS_RE.match(v["id"]):
            v["final"], v["kind"] = v["id"], "plain"
        elif v["id"] in recovered:
            rs, pal = recovered[v["id"]]
            v["final"], v["kind"] = rs, ("recovered_palindromic" if pal == "yes" else "recovered")
        else:
            v["final"], v["kind"] = v["id"], "no_rsid"
        v["decision"] = "kept"
    kinds = defaultdict(int)
    for v in variants:
        kinds[v["kind"]] += 1
    log(f"[kinds] {dict(kinds)}")

    def write_keep(tag, label):
        kept = [v for v in variants if v["decision"] == "kept"]
        with open(os.path.join(args.outdir, f"keep_{tag}.txt"), "w") as out:
            for v in kept:
                out.write(v["id"] + "\n")
        log(f"[{label}] variants remaining: {len(kept):,}")
        return len(kept)

    counts = [("QC dataset (input of this step)", len(variants))]

    # ---- A. recovered palindromic -----------------------------------------------
    for v in variants:
        if v["kind"] == "recovered_palindromic":
            v["decision"] = "excluded_recovered_palindromic"
    counts.append(("A. exclude recovered palindromic (A/T, C/G)", write_keep("A_after_palindromic", "A")))

    # ---- B. incompatible alleles ---------------------------------------------------
    for v in variants:
        if v["decision"] == "kept" and v["kind"] in ("plain", "recovered") and v["final"] in incompatible:
            v["decision"] = "excluded_incompatible_alleles"
    counts.append(("B. exclude incompatible alleles (sumstats/panel/.bim)", write_keep("B_after_incompatible", "B")))

    def rank(v):
        # best call rate first; ties: plain rs ID before recovered / no-rsID; then .bim order
        return (-v["call"], 0 if v["kind"] == "plain" else 1, v["i"])

    # ---- C. duplicates by final rsID ----------------------------------------------------
    groups = defaultdict(list)
    for v in variants:
        if v["decision"] == "kept" and v["kind"] in ("plain", "recovered"):
            groups[v["final"]].append(v)
    n_groups = 0
    for rs, g in groups.items():
        if len(g) > 1:
            n_groups += 1
            g.sort(key=rank)
            for loser in g[1:]:
                loser["decision"] = f"excluded_duplicate_rsid_kept_{g[0]['id']}"
    log(f"[C] rsIDs carried by >1 probe: {n_groups:,}")
    counts.append(("C. de-duplicate by rsID (best call rate)", write_keep("C_after_dedup_rsid", "C")))

    # ---- D. duplicates by chr:pos + allele set -------------------------------------------
    groups = defaultdict(list)
    for v in variants:
        if v["decision"] == "kept":
            groups[(v["chr"], v["pos"], frozenset((v["a1"], v["a2"])))].append(v)
    n_groups = 0
    for key, g in groups.items():
        if len(g) > 1:
            n_groups += 1
            g.sort(key=rank)
            for loser in g[1:]:
                loser["decision"] = f"excluded_duplicate_pos_alleles_kept_{g[0]['id']}"
    log(f"[D] chr:pos+alleles shared by >1 probe: {n_groups:,}")
    counts.append(("D. de-duplicate by chr:pos+alleles (best call rate)", write_keep("D_final", "D")))

    # ---- rename map -------------------------------------------------------------------------
    n_ren = 0
    with open(os.path.join(args.outdir, "rename_map.tsv"), "w") as out:
        for v in variants:
            if v["decision"] == "kept" and v["kind"] == "recovered" and v["final"] != v["id"]:
                out.write(f"{v['id']}\t{v['final']}\n")
                n_ren += 1
    with open(os.path.join(args.outdir, "variant_decisions.tsv"), "w") as out:
        out.write("bim_id\tchr\tpos_grch38\ta1\ta2\tcall_rate\tkind\tfinal_id\tdecision\n")
        for v in variants:
            out.write(f"{v['id']}\t{v['chr']}\t{v['pos']}\t{v['a1']}\t{v['a2']}\t{v['call']:.4f}\t{v['kind']}\t{v['final']}\t{v['decision']}\n")
    kept = [v for v in variants if v["decision"] == "kept"]
    n_rs = sum(1 for v in kept if v["kind"] in ("plain", "recovered"))
    log(f"[final] kept {len(kept):,} | renamed to recovered rsID: {n_ren:,} | with an rsID: {n_rs:,} | without rsID: {len(kept) - n_rs:,}")
    log("[summary of exclusions] " + str({k: sum(1 for v in variants if v['decision'].startswith(k))
                                          for k in ("excluded_recovered_palindromic", "excluded_incompatible_alleles",
                                                    "excluded_duplicate_rsid", "excluded_duplicate_pos_alleles")}))
    with open(os.path.join(args.outdir, "dedup_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
