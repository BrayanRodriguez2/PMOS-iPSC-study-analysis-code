#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/08_frequency_vs_1000g.py   (original step script: 34_1000g_compare.py)
# Produces: Allele frequencies of discordant SNPs against 1000 Genomes EUR/FIN.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 34 - Compare the 446 "mismatched frequency" SNPs (and the 182 allele-mismatch SNPs, and the 29 lead SNPs) with 1000 Genomes
phase 3 (GRCh38 release 20190312, biallelic SNV/INDEL): EUR (503 samples, INCLUDES the 99 FIN) and, separately, FIN (99 Finns).
Diagnostic only; nothing is modified.

Frequencies of OUR ALT allele in 1000G are obtained from the genotypes of the 1000G samples at the same position, matching the allele
pair (same REF/ALT -> AF of ALT; swapped -> 1 - AF; any other pair, or no record -> not available). Genotype extraction from the remote
VCFs is done by the sbatch wrapper (bcftools + tabix index, only the target positions are read).
Why the earlier filter (03_imputation_TOPMed/03_align_to_grch38_freq_filter.py, |study - 1000G EUR_AF| > 0.2) did not remove the 88: the old filter used INFO/EUR_AF of the
1000G *sites* file rounded to 2 decimals and only compared variants with a record with the same allele pair; variants without such a
record were only counted ("not assessed"). Here every one of the 446 is classified into: no 1000G record / allele pair differs /
record present with |study - EUR| <= 0.2 / record present with |study - EUR| > 0.2.

Inputs : results/topmed_qc_only/freq_mismatch_446_verified.tsv, allele_mismatch_182_verified.tsv (step 33), the old-filter table
         (EUR_AF of the sites file at the 446+182 positions: sites_eur_af.tsv), 1000G genotypes (g1000_gt.tsv + g1000_samples.txt + the
         panel file with the population of each sample), lead_snps_lookup.tsv
Outputs: results/topmed_qc_only/compare_1000g_446.tsv, compare_1000g_182.tsv, compare_1000g_report.txt, leads_1000g_positions.tsv
"""
import argparse, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--outdir", required=True)
ap.add_argument("--make-targets", help="write the target positions (CHROM without 'chr', POS) of the 446 + 182 SNPs and the 29 leads, then exit")
ap.add_argument("--leads", help="lead_snps_lookup.tsv (columns rsid, chr, pos_grch38)")
ap.add_argument("--gt"); ap.add_argument("--samples"); ap.add_argument("--panel"); ap.add_argument("--sites-eur-af")
a = ap.parse_args()
if a.make_targets:
    pos = []
    for f in ("freq_mismatch_446_verified.tsv", "allele_mismatch_182_verified.tsv"):
        d = pd.read_csv(os.path.join(a.outdir, f), sep="\t", dtype={"CHROM": str})
        pos.append(pd.DataFrame({"CHROM": d.CHROM.str.replace("chr", "", regex=False), "POS": d.POS.astype(int)}))
    ld = pd.read_csv(a.leads, sep="\t", dtype={"chr": str}).dropna(subset=["pos_grch38"])
    pos.append(pd.DataFrame({"CHROM": ld.chr.astype(str), "POS": ld.pos_grch38.astype(int)}))
    t = pd.concat(pos).drop_duplicates()
    t["k"] = t.CHROM.astype(int); t = t.sort_values(["k", "POS"]).drop(columns="k")
    t.to_csv(a.make_targets, sep="\t", index=False, header=False)
    print(f"{len(t)} target positions written ({len(ld)} leads)")
    raise SystemExit(0)
GT2D = {"0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2, "0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}
samples = [l.strip() for l in open(a.samples) if l.strip()]
panel = pd.read_csv(a.panel, sep="\t").set_index("sample")
fin = [i for i, s in enumerate(samples) if panel.loc[s, "pop"] == "FIN"]
eur = [i for i, s in enumerate(samples) if panel.loc[s, "super_pop"] == "EUR"]
print(f"1000G samples read: {len(samples)} | EUR: {len(eur)} | FIN: {len(fin)}", flush=True)
g = pd.read_csv(a.gt, sep="\t", header=None, dtype=str)
g.columns = ["CHROM", "POS", "REF", "ALT"] + samples
g["CHROM"] = g.CHROM.str.replace("chr", "", regex=False); g["POS"] = g.POS.astype(int)
G = np.vstack([g[s].map(GT2D).fillna(-1).astype(np.int8).to_numpy() for s in samples]).T
print(f"1000G records read: {len(g):,}", flush=True)
rec = {}
for i, r in enumerate(g.itertuples()):
    rec.setdefault((r.CHROM, r.POS), []).append(i)


def af_of(idx, cols):
    d = G[idx][cols]; ok = d >= 0
    return (d.clip(min=0) * ok).sum() / (2 * ok.sum()) if ok.sum() else np.nan, int(ok.sum())


def classify(chrom, pos, ref, alt):
    """-> (status, AF_EUR_of_our_ALT, AF_FIN_of_our_ALT, n_FIN_called, 1000G ref/alt)"""
    hits = rec.get((str(chrom).replace("chr", ""), int(pos)), [])
    if not hits: return "no_1000G_record", np.nan, np.nan, 0, ""
    for i in hits:
        r0, a0 = g.REF.iat[i], g.ALT.iat[i]
        if (r0, a0) == (ref, alt) or (r0, a0) == (alt, ref):
            sw = (r0, a0) == (alt, ref)
            e, _ = af_of(i, eur); f, nf = af_of(i, fin)
            return ("record_same_alleles" if not sw else "record_same_alleles_swapped"), (1 - e if sw else e), (1 - f if sw else f), nf, f"{r0}/{a0}"
    return "allele_pair_differs", np.nan, np.nan, 0, ";".join(f"{g.REF.iat[i]}/{g.ALT.iat[i]}" for i in hits)


old = pd.read_csv(a.sites_eur_af, sep="\t", header=None, names=["CHROM", "POS", "REF", "ALT", "EUR_AF"], dtype={"CHROM": str, "POS": int, "EUR_AF": str})
old["CHROM"] = old.CHROM.str.replace("chr", "", regex=False)
oldrec = {}
for r in old.itertuples():
    oldrec.setdefault((r.CHROM, r.POS), []).append(r)


def old_filter_view(chrom, pos, ref, alt, study_alt):
    """What the step-20 filter saw: (status, EUR_AF of our ALT (2 decimals), |diff|)"""
    rs = oldrec.get((str(chrom).replace("chr", ""), int(pos)), [])
    if not rs: return "no_record", np.nan, np.nan
    for r in rs:
        try: e = float(r.EUR_AF.split(",")[0])
        except ValueError: continue
        if (r.REF, r.ALT) == (ref, alt): return "same_pair", e, abs(study_alt - e)
        if (r.REF, r.ALT) == (alt, ref): return "same_pair", 1 - e, abs(study_alt - (1 - e))
    return "allele_pair_differs", np.nan, np.nan


m = pd.read_csv(os.path.join(a.outdir, "freq_mismatch_446_verified.tsv"), sep="\t", dtype={"CHROM": str})
rows = []
for r in m.itertuples():
    st, e, f, nf, p = classify(r.CHROM, r.POS, r.REF, r.ALT)
    ost, oe, od = old_filter_view(r.CHROM, r.POS, r.REF, r.ALT, r.alt_freq_study_vcf)
    rows.append(dict(ID=r.ID, CHROM=r.CHROM, POS=r.POS, REF=r.REF, ALT=r.ALT, study_alt_freq=r.alt_freq_study_vcf,
                     topmed_alt_freq=r.alt_freq_panel_report, diff_vs_topmed=r.alt_freq_study_vcf - r.alt_freq_panel_report,
                     abs_diff_topmed_gt_0_2=abs(r.alt_freq_study_vcf - r.alt_freq_panel_report) > 0.2,
                     g1000_status=st, g1000_refalt=p, eur_af_genotypes=e, abs_diff_eur=abs(r.alt_freq_study_vcf - e),
                     fin_af=f, n_fin_called=nf, abs_diff_fin=abs(r.alt_freq_study_vcf - f),
                     old_filter_status=ost, old_filter_eur_af_2dec=oe, old_filter_abs_diff=od,
                     n_carriers_here=r.n_carriers, hwe_p=r.hwe_p, chisq=r.chisq))
c = pd.DataFrame(rows)
c.to_csv(os.path.join(a.outdir, "compare_1000g_446.tsv"), sep="\t", index=False, float_format="%.5g")
e182 = pd.read_csv(os.path.join(a.outdir, "allele_mismatch_182_verified.tsv"), sep="\t", dtype={"CHROM": str})
rows = []
for r in e182.itertuples():
    st, e, f, nf, p = classify(r.CHROM, r.POS, r.REF, r.ALT)
    rows.append(dict(ID=r.ID, CHROM=r.CHROM, POS=r.POS, REF=r.REF, ALT=r.ALT, panel_ref=r.panel_ref, panel_alt=r.panel_alt,
                     g1000_status=st, g1000_refalt=p, fin_af=f, eur_af_genotypes=e))
pd.DataFrame(rows).to_csv(os.path.join(a.outdir, "compare_1000g_182.tsv"), sep="\t", index=False, float_format="%.5g")

b = c[c.abs_diff_topmed_gt_0_2]
L = [f"446 SNPs classified against 1000G | all 446: status {c.g1000_status.value_counts().to_dict()}",
     f"88 SNPs with |study - TOPMed| > 0.2 (recomputed: {len(b)})",
     f"  1000G record status: {b.g1000_status.value_counts().to_dict()}",
     f"  what the old step-20 filter saw: {b.old_filter_status.value_counts().to_dict()}",
     f"  with a same-pair record: |study - EUR| (genotypes, 503 EUR): median {b.abs_diff_eur.median():.3f}; > 0.2: {int((b.abs_diff_eur > 0.2).sum())}; <= 0.2: {int((b.abs_diff_eur <= 0.2).sum())}",
     f"  old filter, EUR_AF 2 decimals, |diff|: median {b.old_filter_abs_diff.median():.3f}; > 0.2: {int((b.old_filter_abs_diff > 0.2).sum())}",
     f"  FIN (99): |study - FIN| median {b.abs_diff_fin.median():.3f}; <= 0.2: {int((b.abs_diff_fin <= 0.2).sum())}; > 0.2: {int((b.abs_diff_fin > 0.2).sum())}; "
     f"no FIN frequency (no record / allele pair differs): {int(b.abs_diff_fin.isna().sum())}",
     f"  median ALT freq: study {b.study_alt_freq.median():.3f} | TOPMed {b.topmed_alt_freq.median():.3f} | EUR {b.eur_af_genotypes.median():.3f} | FIN {b.fin_af.median():.3f}",
     "", f"Among all 446: FIN available {int(c.fin_af.notna().sum())}; |study - FIN| <= 0.2: {int((c.abs_diff_fin <= 0.2).sum())}; "
         f"median study {c.study_alt_freq.median():.3f} | TOPMed {c.topmed_alt_freq.median():.4f} | EUR {c.eur_af_genotypes.median():.3f} | FIN {c.fin_af.median():.3f}",
     f"Among the other 358 (|diff| <= 0.2 vs TOPMed): median study {c[~c.abs_diff_topmed_gt_0_2].study_alt_freq.median():.3f}, FIN {c[~c.abs_diff_topmed_gt_0_2].fin_af.median():.3f}, "
     f"EUR {c[~c.abs_diff_topmed_gt_0_2].eur_af_genotypes.median():.3f}; FIN/EUR frequency ratio median {(c.fin_af / c.eur_af_genotypes.replace(0, np.nan)).median():.2f}"]
open(os.path.join(a.outdir, "compare_1000g_report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
