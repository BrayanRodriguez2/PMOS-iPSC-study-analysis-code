#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/19_grs29.py   (original step script: 51b_grs29.py)
# Produces: Weights and allele-orientation table of the 29-SNP genetic risk score.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 51b - weights and allele-orientation table of the 29-SNP genetic risk score (GRS) from the lead SNPs of Supplementary Table 2 (ST2) of the article
(Moolhuijsen et al., Nat Genet 2026; "PCOS signals identified in the full 3-strata model").

For each of the 29 leads:
  * ST2 gives the PCOS RISK-INCREASING allele, the other allele, its frequency and the effect (beta of the risk allele, all positive) on the forward strand
    (GRCh37 coordinates in ST2; GRCh38 position from the harmonised GCST90570579 summary statistics via the rsID, table lead_snps_lookup.tsv of step 11).
  * The variant is looked up in the clean imputed file by chromosome + GRCh38 position + allele SET (clean variant_table.tsv.gz). TOPMed/clean alleles are
    on the forward strand (REF = GRCh38 base), the same as the GWAS, so a lead matches DIRECTLY ({risk, other} = {REF, ALT}); a strand-complement match is
    accepted only for non-palindromic SNPs. Palindromic leads (A/T, C/G) can only be matched directly; their orientation is checked with the frequency.
  * Scored allele = the risk allele written in the clean file's alleles; weight = ST2 effect (positive). GRS = sum(effect x dosage of the risk allele)
    with imputed dosages, computed with plink2 --score in step 51.
Orientation checks written to the table (one row per lead):
  risk allele freq. in ST2 vs. in the 129 samples (from the clean file, AF129 of the ALT allele converted to the risk allele) vs. the effect-allele frequency
  in the deposited sumstats (2-strata); the sign of the deposited beta relative to the ST2 risk allele (+ means the deposited effect allele is the risk
  allele with beta > 0, or the other allele with beta < 0); R2 / ER2 / typed-or-imputed; flags for palindromic SNPs and for |freq difference| > 0.25.
Inputs : ST2 tsv (step 48a: sheet03_ST2*.tsv), lead_snps_lookup.tsv, clean/variant_table.tsv.gz
Outputs: --outdir/grs29_weights.txt (ID A1 BETA, header), grs29_orientation.tsv, grs29_report.txt
"""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--st2", required=True); ap.add_argument("--lookup", required=True); ap.add_argument("--variant-table", required=True); ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
comp = lambda s: "".join(COMP.get(x, "N") for x in s)
aset = lambda x, y: min(x, y) + max(x, y)

# ---- ST2 (header row = the line that starts with 'SNP')
raw = [l.rstrip("\n").split("\t") for l in open(a.st2, encoding="utf-8")]
hi = next(i for i, r in enumerate(raw) if r and r[0].strip() == "SNP")
hdr = [h.strip() for h in raw[hi]]
rows = [r for r in raw[hi + 1:] if r and r[0].strip().startswith("rs")]
st2 = pd.DataFrame([r + [""] * (len(hdr) - len(r)) for r in rows], columns=hdr)
col = lambda key: next(c for c in st2.columns if c.lower().startswith(key.lower()))
st2 = pd.DataFrame({"rsid": st2["SNP"].str.strip(), "st2_variant_grch37": st2["Variant"], "gene": st2[col("Nearest")],
                    "risk_allele": st2[col("PCOS risk increasing")].str.strip().str.upper(), "other_allele": st2[col("Other allele")].str.strip().str.upper(),
                    "st2_risk_freq": st2[col("Risk increasing allele freq")].astype(float), "effect": st2["Effect"].astype(float),
                    "se": st2["Standard Error"].astype(float), "p_st2": st2["P-value"].astype(float)})
assert len(st2) == 29 and st2.rsid.is_unique, f"ST2 rows: {len(st2)}"
lk = pd.read_csv(a.lookup, sep="\t", dtype={"chr": str}).rename(columns={"effect_allele": "ss_effect_allele", "other_allele": "ss_other_allele", "beta": "ss_beta", "eaf": "ss_eaf", "p_value": "ss_p"})
lk = lk[["rsid", "chr", "pos_grch38", "ss_effect_allele", "ss_other_allele", "ss_beta", "ss_eaf", "ss_p"]]
t = st2.merge(lk, on="rsid", how="left")
assert t.pos_grch38.notna().all(), "lead without GRCh38 position"
t["chr"] = t.chr.astype(int); t["pos_grch38"] = t.pos_grch38.astype(int)

vt = pd.read_csv(a.variant_table, sep="\t", dtype={"CHROM": int}, usecols=["CHROM", "POS", "ID", "REF", "ALT", "rsID", "R2", "ER2", "TYPE", "AF129"])
vt["aset"] = [aset(r, al) for r, al in zip(vt.REF, vt.ALT)]
out = []
for r in t.itertuples():
    pal = comp(r.risk_allele) == r.other_allele
    s_direct = aset(r.risk_allele, r.other_allele)
    s_comp = aset(comp(r.risk_allele), comp(r.other_allele))
    cand = vt[(vt.CHROM == r.chr) & (vt.POS == r.pos_grch38) & ((vt.aset == s_direct) | ((vt.aset == s_comp) & (not pal)))]
    row = r._asdict(); row.pop("Index")
    row["palindromic"] = pal
    if len(cand) == 0:
        row.update(status="NOT in the clean file", clean_id="", ref="", alt="", match="", scored_allele="", rsid_in_clean="", R2=float("nan"), ER2=float("nan"), type="", cohort_risk_freq=float("nan"))
    else:
        c = cand.iloc[0]
        if len(cand) > 1:
            row["status"] = "AMBIGUOUS (%d records)" % len(cand)
        else:
            row["status"] = "ok"
        direct = c.aset == s_direct
        scored = r.risk_allele if direct else comp(r.risk_allele)
        assert scored in (c.REF, c.ALT)
        f_alt = float(c.AF129)
        row.update(clean_id=c.ID, ref=c.REF, alt=c.ALT, match="direct" if direct else "complement", scored_allele=scored, rsid_in_clean=c.rsID, R2=float(c.R2),
                   ER2=float(c.ER2) if pd.notna(c.ER2) else float("nan"), type=c.TYPE, cohort_risk_freq=f_alt if scored == c.ALT else 1 - f_alt)
    out.append(row)
o = pd.DataFrame(out)
# deposited sumstats orientation: beta of the deposited effect allele; sign relative to the ST2 risk allele
def ss_sign(r):
    if r.ss_effect_allele == r.risk_allele or comp(r.ss_effect_allele) == r.risk_allele and not r.palindromic:
        return "+" if r.ss_beta > 0 else "-"
    return "+" if r.ss_beta < 0 else "-"
o["sumstats_sign_vs_ST2_risk"] = o.apply(ss_sign, axis=1)
def ss_freq_risk(r):
    return r.ss_eaf if (r.ss_effect_allele == r.risk_allele or (comp(r.ss_effect_allele) == r.risk_allele and not r.palindromic)) else 1 - r.ss_eaf
o["sumstats_risk_freq"] = o.apply(ss_freq_risk, axis=1)
o["freq_diff_cohort_vs_ST2"] = o.cohort_risk_freq - o.st2_risk_freq
o["flag"] = ""
o.loc[o.freq_diff_cohort_vs_ST2.abs() > 0.25, "flag"] += "freq_diff>0.25;"
o.loc[o.palindromic, "flag"] += "palindromic(orientation by frequency);"
o.loc[o.sumstats_sign_vs_ST2_risk != "+", "flag"] += "sumstats_sign_differs;"
o.loc[o.status != "ok", "flag"] += "NOT_USED;"
cols = ["rsid", "chr", "pos_grch38", "gene", "st2_variant_grch37", "risk_allele", "other_allele", "effect", "se", "st2_risk_freq", "clean_id", "ref", "alt", "match", "scored_allele",
        "type", "R2", "ER2", "cohort_risk_freq", "freq_diff_cohort_vs_ST2", "sumstats_risk_freq", "ss_effect_allele", "ss_other_allele", "ss_beta", "sumstats_sign_vs_ST2_risk",
        "palindromic", "rsid_in_clean", "status", "flag"]
o[cols].to_csv(os.path.join(a.outdir, "grs29_orientation.tsv"), sep="\t", index=False, float_format="%.4g")
use = o[o.status == "ok"]
with open(os.path.join(a.outdir, "grs29_weights.txt"), "w") as f:
    f.write("ID\tA1\tBETA\n")
    for r in use.itertuples():
        f.write(f"{r.clean_id}\t{r.scored_allele}\t{r.effect}\n")
L = [f"leads in ST2: {len(o)} | found in the clean file: {len(use)} | direct matches: {int((use.match == 'direct').sum())}, complement: {int((use.match == 'complement').sum())} | palindromic leads: {int(o.palindromic.sum())}",
     f"typed {int((use.type == 'Typed').sum())}, imputed {int((use.type == 'Imputed').sum())} | R2 min {use.R2.min():.3f} | sum of weights {use.effect.sum():.3f}",
     f"sign of the deposited (2-strata) beta agrees with the ST2 risk allele in {int((o.sumstats_sign_vs_ST2_risk == '+').sum())}/29",
     f"|cohort - ST2 risk-allele freq| > 0.25: {int((o.freq_diff_cohort_vs_ST2.abs() > 0.25).sum())} | correlation cohort vs ST2 risk freq: {use.cohort_risk_freq.corr(use.st2_risk_freq):.3f}",
     o[["rsid", "gene", "risk_allele", "other_allele", "effect", "st2_risk_freq", "cohort_risk_freq", "scored_allele", "match", "type", "R2", "sumstats_sign_vs_ST2_risk", "flag"]].round(3).to_string(index=False)]
open(os.path.join(a.outdir, "grs29_report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
