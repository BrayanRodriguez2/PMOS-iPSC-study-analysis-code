#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/05_ld_proxies.py   (original step script: 12_ld_proxies.py)
# Produces: LD proxies of the lead SNPs and locus coverage on the array.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 12 - Part B.2: is each lead SNP covered by the donor array, and if not, is there a proxy?

LD reference: 1000 Genomes phase 3, GRCh38 (EBI "20190312_biallelic_SNV_and_INDEL"), the 503
EUR samples (integrated_call_samples_v3.20130502.ALL.panel). The donors themselves (n~129) are NOT
used for LD (too few samples) -- only the external European panel.
For each lead SNP, genotypes of every biallelic SNP within +/- 500 kb were extracted with bcftools
(step 12 .sbatch, one text file per locus, columns POS REF ALT then one GT per sample) and the squared
Pearson correlation of allele dosages (r2, unphased, equals plink2 --r2-unphased) between the lead and
each SNP with MAF >= 0.01 is computed here (numpy). Missing genotypes are mean-imputed.

"On the array" = the lead (or a candidate proxy) matches a variant of the QC dataset b38_qc_rsid at
the same chr:GRCh38 position with the same alleles (any order, either strand), or the same rsID.
"Usable by PRS-CS" = on the QC array AND its rsID is in the PRS-CS LD panel snpinfo (HapMap3 SNPs):
PRS-CS can only give a weight to such variants when run on the array without imputation.

Coverage definitions reported per locus:
  genotyped            lead on the QC array (also: on the raw array before QC)
  proxy_r2_ge_0.8      best proxy (r2 >= 0.8) among QC-array SNPs
  usable_*             same, but restricted to the SNPs usable by PRS-CS
Outputs (--outdir): locus_coverage.tsv and locus_coverage_report.txt
"""
import argparse
import gzip
import os
import re

import numpy as np

COMP = str.maketrans("ACGT", "TGCA")


def same_snp(p, q):
    """Same biallelic SNP ignoring allele order and strand."""
    sp, sq = frozenset(x.upper() for x in p), frozenset(x.upper() for x in q)
    if sp == sq:
        return True
    return all(len(x) == 1 and x.upper() in "ACGT" for x in q) and sp == frozenset(x.upper().translate(COMP) for x in q)


def load_bim(path):
    by_pos, ids = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            by_pos.setdefault((f[0], int(f[3])), []).append((f[1], f[4], f[5]))
            ids[f[1]] = (f[0], int(f[3]))
    return by_pos, ids


def on_array(by_pos, ids, rsid, chrom, pos, a, b):
    """Return the array variant ID matching (rsid) or (chrom,pos,alleles), else None."""
    if rsid in ids:
        return rsid
    for vid, a1, a2 in by_pos.get((chrom, pos), []):
        if "." in (a1, a2) or "0" in (a1, a2):
            continue
        if same_snp((a1, a2), (a, b)):
            return vid
    return None


def read_gt(path):
    """Text from `bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'` -> (pos, ref, alt, dosage matrix)."""
    pos, ref, alt, rows = [], [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f[2]) != 1 or len(f[1]) != 1:
                continue  # SNPs only
            g = []
            for t in f[3:]:
                if "." in t:
                    g.append(np.nan)
                else:
                    g.append(sum(1 for c in re.split(r"[|/]", t) if c == "1"))
            pos.append(int(f[0]))
            ref.append(f[1])
            alt.append(f[2])
            rows.append(g)
    return np.array(pos), ref, alt, np.array(rows, dtype=float)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--leads", required=True, help="lead_snps_lookup.tsv from step 11")
    ap.add_argument("--gtdir", required=True, help="directory with <rsid>.gt.tsv.gz from step 12 bash")
    ap.add_argument("--qc-bim", required=True)
    ap.add_argument("--raw-bim", required=True)
    ap.add_argument("--panel-snpinfo", required=True, help="PRS-CS snpinfo of the MAIN (UKBB) panel")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--window", type=int, default=500_000)
    ap.add_argument("--r2", type=float, default=0.8)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    qc_pos, qc_ids = load_bim(args.qc_bim)
    raw_pos, raw_ids = load_bim(args.raw_bim)
    panel = set()
    with open(args.panel_snpinfo) as fh:
        next(fh)
        for line in fh:
            panel.add(line.split()[1])
    rep = []

    def log(m=""):
        print(m, flush=True)
        rep.append(m)

    leads = []
    with open(args.leads) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            leads.append(dict(zip(h, line.rstrip("\n").split("\t"))))
    cols = ["rsid", "chr", "pos_grch38", "p_value_deposited", "lead_in_1kg_eur", "lead_maf_1kg_eur", "lead_on_raw_array",
            "lead_on_qc_array", "lead_in_prscs_panel", "lead_usable_by_prscs", "n_snps_window_1kg", "n_proxies_r2_ge_thr_on_qc_array",
            "best_proxy_qc_array", "best_proxy_pos", "best_proxy_r2", "best_proxy_usable_by_prscs", "best_usable_proxy_r2",
            "best_usable_proxy", "covered_genotyped", "covered_genotyped_or_proxy", "covered_usable_by_prscs"]
    table = []
    for L in leads:
        rs, chrom = L["rsid"], L["chr"]
        row = dict.fromkeys(cols, "NA")
        row.update(rsid=rs, chr=chrom, pos_grch38=L["pos_grch38"], p_value_deposited=L["p_value"])
        if chrom == "NA":
            table.append(row)
            continue
        pos0, a, b = int(L["pos_grch38"]), L["effect_allele"], L["other_allele"]
        raw_hit = on_array(raw_pos, raw_ids, rs, chrom, pos0, a, b)
        qc_hit = on_array(qc_pos, qc_ids, rs, chrom, pos0, a, b)
        row["lead_on_raw_array"] = raw_hit or "no"
        row["lead_on_qc_array"] = qc_hit or "no"
        row["lead_in_prscs_panel"] = "yes" if rs in panel else "no"
        row["lead_usable_by_prscs"] = "yes" if (qc_hit and (rs in panel or qc_hit in panel)) else "no"
        gt = os.path.join(args.gtdir, f"{rs}.gt.tsv.gz")
        best, best_u, n_prox = None, None, 0
        loaded = None
        if os.path.exists(gt) and os.path.getsize(gt) > 0:
            loaded = read_gt(gt)
            if loaded[3].ndim < 2 or loaded[3].shape[0] == 0:  # empty slice
                loaded = None
        if loaded is not None:
            pos, ref, alt, G = loaded
            n = G.shape[1]
            varmean = np.nanmean(G, axis=1)  # per-variant mean over samples, used to impute missing genotypes
            G = np.where(np.isnan(G), varmean[:, None], G)
            maf = np.minimum(G.mean(axis=1) / 2, 1 - G.mean(axis=1) / 2)
            li = [i for i in range(len(pos)) if pos[i] == pos0 and same_snp((ref[i], alt[i]), (a, b))]
            row["n_snps_window_1kg"] = len(pos)
            if li:
                li = li[0]
                row["lead_in_1kg_eur"] = "yes"
                row["lead_maf_1kg_eur"] = f"{maf[li]:.3f}"
                z = (G - G.mean(axis=1, keepdims=True))
                sd = z.std(axis=1)
                ok = sd > 0
                zl = z[li] / (sd[li] if sd[li] > 0 else 1)
                r = np.zeros(len(pos))
                r[ok] = (z[ok] @ zl) / (n * sd[ok])
                r2 = r ** 2
                for i in range(len(pos)):
                    if i == li or maf[i] < 0.01 or abs(int(pos[i]) - pos0) > args.window:
                        continue
                    hit = None
                    for vid, a1, a2 in qc_pos.get((chrom, int(pos[i])), []):
                        if "." in (a1, a2) or "0" in (a1, a2):
                            continue
                        if same_snp((a1, a2), (ref[i], alt[i])):
                            hit = vid
                            break
                    if hit is None:
                        continue
                    usable = hit in panel
                    if r2[i] >= args.r2:
                        n_prox += 1
                    if best is None or r2[i] > best[0]:
                        best = (r2[i], hit, int(pos[i]), usable)
                    if usable and (best_u is None or r2[i] > best_u[0]):
                        best_u = (r2[i], hit, int(pos[i]))
            else:
                row["lead_in_1kg_eur"] = "no"
        else:
            row["lead_in_1kg_eur"] = "not_extracted"
        row["n_proxies_r2_ge_thr_on_qc_array"] = n_prox
        if best:
            row.update(best_proxy_qc_array=best[1], best_proxy_pos=best[2], best_proxy_r2=f"{best[0]:.3f}",
                       best_proxy_usable_by_prscs="yes" if best[3] else "no")
        if best_u:
            row.update(best_usable_proxy=best_u[1], best_usable_proxy_r2=f"{best_u[0]:.3f}")
        geno = row["lead_on_qc_array"] != "no"
        prox = best is not None and best[0] >= args.r2
        uprox = best_u is not None and best_u[0] >= args.r2
        row["covered_genotyped"] = "yes" if geno else "no"
        row["covered_genotyped_or_proxy"] = "yes" if (geno or prox) else "no"
        row["covered_usable_by_prscs"] = "yes" if (row["lead_usable_by_prscs"] == "yes" or uprox) else "no"
        table.append(row)

    with open(os.path.join(args.outdir, "locus_coverage.tsv"), "w") as out:
        out.write("\t".join(cols) + "\n")
        for row in table:
            out.write("\t".join(str(row[c]) for c in cols) + "\n")
    ev = [r for r in table if r["lead_in_1kg_eur"] in ("yes", "no")]
    cnt = lambda k: sum(1 for r in table if r[k] == "yes")
    log(f"[loci] leads: {len(table)} | evaluated against 1000G EUR: {len(ev)} | lead absent from the 1000G SNP set: "
        f"{sum(1 for r in table if r['lead_in_1kg_eur'] == 'no')} | genotype slice not extracted: "
        f"{sum(1 for r in table if r['lead_in_1kg_eur'] == 'not_extracted')}")
    log(f"[coverage, r2 >= {args.r2}, +/-{args.window // 1000} kb, 1000G EUR GRCh38]")
    log(f"  lead on the RAW array: {sum(1 for r in table if r['lead_on_raw_array'] not in ('no', 'NA'))} of {len(table)}")
    log(f"  lead on the QC array (genotyped): {cnt('covered_genotyped')} of {len(table)}")
    log(f"  genotyped OR proxy: {cnt('covered_genotyped_or_proxy')} of {len(table)}")
    log(f"  usable by PRS-CS (lead or proxy in array & HapMap3 panel): {cnt('covered_usable_by_prscs')} of {len(table)}")
    with open(os.path.join(args.outdir, "locus_coverage_report.txt"), "w") as out:
        out.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
