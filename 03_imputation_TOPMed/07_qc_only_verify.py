#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 03_imputation_TOPMed/07_qc_only_verify.py   (original step script: 33_qc_only_verify.py)
# Produces: Independent verification of the server QC-only statistics.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 33 - Independent check of the TOPMed Imputation Server "Quality Control Only" job (job "imputation_ipsc": the 22 v1 VCFs of
step 21, 140 samples x 510,607 SNPs, panel TOPMed r3 2023, hg38). Nothing is modified; the downloaded files are only parsed as text.

Two sub-commands (run in this order by 03_imputation_TOPMed/07_qc_only_verify.sbatch):
  parse   : reads snps-excluded.txt and the 'Mismatched frequencies' lines of quality-control.html, joins them to the v1 site table
            (ID, CHROM, POS, REF, ALT of the uploaded VCFs) and writes the target positions to extract genotypes from the VCFs.
  analyze : with the genotypes of those positions (140 samples) recomputes, for the 446 frequency-mismatch SNPs, the allele frequency,
            AA/AB/BB counts, exact HWE p-value (Wigginton et al. 2005), number of carriers and carriers per sample, and checks the
            182 "Allele mismatch" SNPs. Also tests whether the 11 sample-QC-excluded donors (and the 5 contaminated ones) carry more
            of these variants than the other donors.
Server report format (quality-control.html, pandoc rendering of the job log), one line per SNP:
  ## Mismatched frequencies for 'ID' f[A1,A2] = [studyA1,studyA2] vs [panelA1,panelA2], chisq X
where [A1,A2] follows the REF,ALT order of the uploaded VCF (verified here against the VCF); the ALT frequency is the second number.
"""
import argparse, html, os, re
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)
p1 = sub.add_parser("parse")
p1.add_argument("--html", required=True); p1.add_argument("--snps-excluded", required=True)
p1.add_argument("--v1-sites", required=True, help="ID CHROM POS REF ALT of the v1 VCFs"); p1.add_argument("--outdir", required=True)
p2 = sub.add_parser("analyze")
p2.add_argument("--outdir", required=True); p2.add_argument("--gt", required=True, help="CHROM POS ID REF ALT + 140 GT columns")
p2.add_argument("--samples", required=True); p2.add_argument("--qc-table", required=True)
a = ap.parse_args()
GT2D = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/0": 1, "1/1": 2, "1|1": 2}
RX = re.compile(r"Mismatched frequencies for '([^']+)' f\[([^\]]+)\] = \[([^\]]+)\] vs \[([^\]]+)\], chisq ([0-9.eE+-]+)")


def hwe_exact(obs_hets, obs_hom1, obs_hom2):
    """Exact two-sided HWE test, Wigginton, Cutler & Abecasis 2005 (AJHG 76:887) - the algorithm used by PLINK/VCFtools."""
    obs_homc, obs_homr = max(obs_hom1, obs_hom2), min(obs_hom1, obs_hom2)
    rare = 2 * obs_homr + obs_hets; n = obs_hets + obs_homc + obs_homr
    if n == 0: return 1.0
    probs = np.zeros(rare + 1)
    mid = int(rare * (2 * n - rare) / (2 * n))
    if (mid % 2) != (rare % 2): mid += 1
    probs[mid] = 1.0; total = 1.0
    h, hr = mid, (rare - mid) // 2; hc = n - h - hr
    while h > 1:
        probs[h - 2] = probs[h] * h * (h - 1.0) / (4.0 * (hr + 1.0) * (hc + 1.0)); total += probs[h - 2]; hr += 1; hc += 1; h -= 2
    h, hr = mid, (rare - mid) // 2; hc = n - h - hr
    while h <= rare - 2:
        probs[h + 2] = probs[h] * 4.0 * hr * hc / ((h + 2.0) * (h + 1.0)); total += probs[h + 2]; hr -= 1; hc -= 1; h += 2
    probs /= total
    return float(min(1.0, probs[probs <= probs[obs_hets] + 1e-12].sum()))


if a.cmd == "parse":
    os.makedirs(a.outdir, exist_ok=True)
    sites = pd.read_csv(a.v1_sites, sep="\t", header=None, names=["ID", "CHROM", "POS", "REF", "ALT"], dtype={"CHROM": str, "POS": int})
    print(f"v1 sites: {len(sites):,} | unique IDs: {sites.ID.is_unique}", flush=True)
    rows = []
    for line in open(a.html, encoding="utf-8", errors="replace"):
        m = RX.search(html.unescape(line))
        if m:
            i, al, fs, fp, chi = m.groups()
            al = al.split(","); fs = [float(x) for x in fs.split(",")]; fp = [float(x) for x in fp.split(",")]
            rows.append(dict(ID=i, a1=al[0], a2=al[1], study_f1=fs[0], study_f2=fs[1], panel_f1=fp[0], panel_f2=fp[1], chisq=float(chi)))
    mm = pd.DataFrame(rows)
    print(f"'Mismatched frequencies' lines in the html: {len(mm)} | unique IDs: {mm.ID.is_unique}", flush=True)
    mm = mm.merge(sites, on="ID", how="left")
    print("IDs of the 446 not found in the v1 VCFs:", int(mm.CHROM.isna().sum()), flush=True)
    ex = pd.read_csv(a.snps_excluded, sep="\t", dtype={"CHROM": str, "POS": int})
    print(f"snps-excluded.txt: {len(ex)} rows | reasons: {ex.INFO.str.replace(r'Reference Panel: .*', '', regex=True).value_counts().to_dict()}", flush=True)
    ex.columns = ["ID", "CHROM_x", "POS_x", "REF_x", "ALT_x", "INFO"]
    ex = ex.merge(sites, on="ID", how="left")
    same = ((ex.CHROM == ex.CHROM_x) & (ex.POS == ex.POS_x) & (ex.REF == ex.REF_x) & (ex.ALT == ex.ALT_x))
    print("excluded IDs not found in v1:", int(ex.CHROM.isna().sum()), "| position/alleles identical to v1:", int(same.sum()), flush=True)
    mm.to_csv(os.path.join(a.outdir, "freq_mismatch_446_parsed.tsv"), sep="\t", index=False)
    ex.to_csv(os.path.join(a.outdir, "allele_mismatch_182_parsed.tsv"), sep="\t", index=False)
    t = pd.concat([mm[["CHROM", "POS"]], ex[["CHROM", "POS"]]]).dropna().astype({"POS": int}).drop_duplicates()
    t.to_csv(os.path.join(a.outdir, "work_targets.tsv"), sep="\t", index=False, header=False)
    print(f"target positions for genotype extraction: {len(t)}", flush=True)

else:
    mm = pd.read_csv(os.path.join(a.outdir, "freq_mismatch_446_parsed.tsv"), sep="\t", dtype={"CHROM": str})
    ex = pd.read_csv(os.path.join(a.outdir, "allele_mismatch_182_parsed.tsv"), sep="\t", dtype={"CHROM": str, "CHROM_x": str})
    samples = [l.strip() for l in open(a.samples) if l.strip()]
    g = pd.read_csv(a.gt, sep="\t", header=None, dtype=str)
    g.columns = ["CHROM", "POS", "ID", "REF", "ALT"] + samples
    G = np.vstack([g[s].map(GT2D).fillna(-1).astype(np.int8).to_numpy() for s in samples]).T      # sites x samples
    gi = {k: i for i, k in enumerate(g.ID)}
    qc = pd.read_csv(a.qc_table, sep="\t", dtype={"IID": str}).set_index("IID")
    rows = []
    for r in mm.itertuples():
        d = G[gi[r.ID]]
        n0, n1, n2 = int((d == 0).sum()), int((d == 1).sum()), int((d == 2).sum()); n = n0 + n1 + n2
        rows.append(dict(ID=r.ID, CHROM=r.CHROM, POS=r.POS, REF=r.REF, ALT=r.ALT, report_a1=r.a1, report_a2=r.a2,
                         alleles_order_matches_vcf=(r.a1 == r.REF and r.a2 == r.ALT),
                         alt_freq_study_report=r.study_f2, alt_freq_panel_report=r.panel_f2, chisq=r.chisq,
                         alt_freq_study_vcf=(n1 + 2 * n2) / (2 * n) if n else np.nan, n_called=n, n_AA=n0, n_AB=n1, n_BB=n2,
                         hwe_p=hwe_exact(n1, n0, n2), n_carriers=n1 + n2, diff_study_minus_panel=r.study_f2 - r.panel_f2))
    t = pd.DataFrame(rows)
    t["abs_diff"] = t.diff_study_minus_panel.abs(); t["gt_0.2"] = t.abs_diff > 0.2
    t.to_csv(os.path.join(a.outdir, "freq_mismatch_446_verified.tsv"), sep="\t", index=False, float_format="%.6g")
    Gm = G[np.array([gi[i] for i in t.ID])]                         # 446 x 140
    cps = (Gm >= 1).sum(0)
    per_sample = pd.DataFrame({"IID": samples, "n_carried_of_446": cps,
                               "excluded_sample_qc": [qc.loc[s, "excluded"] if s in qc.index else "NA" for s in samples]})
    per_sample.to_csv(os.path.join(a.outdir, "freq_mismatch_446_carriers_per_sample.tsv"), sep="\t", index=False)
    cont = [x for x in __import__("os").environ.get("CONTAMINATED_SAMPLE_IDS", "").split(",") if x]
    L = [f"446 SNPs parsed from the html: {len(t)} | all found in the v1 VCFs: {bool(t.POS.notna().all())}",
         f"allele order in the report == (REF, ALT) of the VCF: {int(t.alleles_order_matches_vcf.sum())}/{len(t)}",
         f"ALT frequency in the report vs recomputed from the VCF: max |difference| = {(t.alt_freq_study_report - t.alt_freq_study_vcf).abs().max():.4f}",
         f"study ALT freq > panel ALT freq in: {int((t.diff_study_minus_panel > 0).sum())}/{len(t)} SNPs",
         f"median ALT freq: study {t.alt_freq_study_report.median() * 100:.2f} % | panel {t.alt_freq_panel_report.median() * 100:.3f} %",
         f"panel ALT freq < 1 %: {int((t.alt_freq_panel_report < 0.01).sum())}/{len(t)}",
         f"|study - panel| > 0.2: {int(t['gt_0.2'].sum())} SNPs",
         f"HWE exact p < 1e-6: {int((t.hwe_p < 1e-6).sum())}/{len(t)} " + (f"(ID {', '.join(t[t.hwe_p < 1e-6].ID)}; p = {', '.join(f'{x:.2g}' for x in t[t.hwe_p < 1e-6].hwe_p)})" if (t.hwe_p < 1e-6).any() else ""),
         f"carriers per SNP: median {t.n_carriers.median():.0f}, min {t.n_carriers.min()}, max {t.n_carriers.max()} (of 140 samples)",
         f"SNPs per sample carried (of 446): median {np.median(cps):.0f}, min {cps.min()}, max {cps.max()}; samples carrying >= 1: {int((cps > 0).sum())}/140"]
    ex_s = per_sample[per_sample.excluded_sample_qc == "yes"].n_carried_of_446; ok_s = per_sample[per_sample.excluded_sample_qc == "no"].n_carried_of_446
    L.append(f"11 sample-QC-excluded donors: median {ex_s.median():.0f} vs 129 others: median {ok_s.median():.0f} (Mann-Whitney p = {mannwhitneyu(ex_s, ok_s, alternative='two-sided')[1]:.3g})")
    cs = per_sample[per_sample.IID.isin(cont)].n_carried_of_446; os_ = per_sample[~per_sample.IID.isin(cont)].n_carried_of_446
    L.append(f"5 contaminated donors {cont}: {cs.tolist()} (median {cs.median():.0f}) vs the other 135: median {os_.median():.0f} (p = {mannwhitneyu(cs, os_, alternative='two-sided')[1]:.3g})")
    ex["panel_ref"] = ex.INFO.str.extract(r"Reference Panel: ([ACGT]+)/")[0]; ex["panel_alt"] = ex.INFO.str.extract(r"Reference Panel: [ACGT]+/([ACGT]+)")[0]
    same = (ex.CHROM == ex.CHROM_x) & (ex.POS == ex.POS_x) & (ex.REF == ex.REF_x) & (ex.ALT == ex.ALT_x)
    L += ["", f"snps-excluded.txt: {len(ex)} SNPs, all 'Allele mismatch': {bool(ex.INFO.str.startswith('Allele mismatch').all())}; "
              f"present in the v1 VCFs with identical position and alleles: {int(same.sum())}/{len(ex)}",
          f"  panel REF equals our REF: {int((ex.panel_ref == ex.REF_x).sum())}/{len(ex)} | panel lists the same base twice (X/X: no ALT in the panel at that site): "
          f"{int((ex.panel_ref == ex.panel_alt).sum())}/{len(ex)} | our ALT == panel ALT: {int((ex.ALT_x == ex.panel_alt).sum())}",
          "  ID prefixes: " + str(ex.ID.str.extract(r"^([A-Za-z]+)")[0].value_counts().head(6).to_dict())]
    ge = [gi[i] for i in ex.ID if i in gi]
    if ge:
        De = G[ge]; nalt = (De == 1).sum(1) + 2 * (De == 2).sum(1); tot = 2 * (De >= 0).sum(1); maf = np.minimum(nalt, tot - nalt) / tot
        L.append(f"  182 allele-mismatch SNPs in the VCF: median MAF in our samples {np.median(maf) * 100:.1f} %, MAF >= 5 %: {int((maf >= 0.05).sum())}")
    open(os.path.join(a.outdir, "verify_report.txt"), "w").write("\n".join(L) + "\n")
    ex.to_csv(os.path.join(a.outdir, "allele_mismatch_182_verified.tsv"), sep="\t", index=False)
    print("\n".join(L))
