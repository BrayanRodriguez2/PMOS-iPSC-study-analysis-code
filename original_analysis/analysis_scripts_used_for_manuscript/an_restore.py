"""v3 (2 Nov 2026): restores the Sept 2026 analyses (correlation, volcano, variance partition,
BMI-adjusted DE, BMI x PMOS interaction, BMI-stratified donors, HipSci protein cross-reference,
literature candidate screen) on the SAME 16,489-gene set and models as the rest of the manuscript.
Usage: python3 an_restore.py <package_dir> <counts.tabular> <sept_package_dir> <out_dir>"""
import sys, os, itertools
import numpy as np, pandas as pd
from scipy import stats
PKG, COUNTS, SEP, OUT = sys.argv[1:5]
sys.path.insert(0, PKG + "/5_Pipeline_CSC/scripts")
from common import read_counts, ols, bh          # hel_name maps HEL362_3 -> HEL362.1
os.makedirs(OUT, exist_ok=True)
R = PKG + "/4_Analysis_results"
 
cnt = read_counts(COUNTS)
key = pd.read_csv(f"{R}/metadata/sample_key_14_ipsc_donors.csv").set_index("hel_line")
assert set(key.index) <= set(cnt.columns), set(key.index) - set(cnt.columns)
ph = pd.read_csv(f"{R}/metadata/donor_phenotypes.tsv", sep="\t")
lines = list(key.sort_values(["group", "bmi"]).index); cnt = cnt[lines]
expd = cnt[(cnt >= 10).sum(1) >= 7]; lg = np.log2(expd + 1)
n = len(lines); G = len(lg)
grp = (key.loc[lines, "group"] == "PMOS").astype(float).values
xist = np.log2(cnt.loc["XIST", lines] + 1).values
bmi = key.loc[lines, "bmi"].values.astype(float); age = key.loc[lines, "age"].values.astype(float)
obese = (bmi > 30).astype(float)
assert ((bmi < 25) | (bmi > 30)).all()
Y = lg.values
log = open(f"{OUT}/log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n")
 
# ---------- vectorised OLS for all genes ----------
def fit(X, Y=Y):
    X1 = np.column_stack([np.ones(n), X]); XtXi = np.linalg.pinv(X1.T @ X1)
    B = Y @ X1 @ XtXi.T; E = Y - B @ X1.T; df = n - X1.shape[1]
    s2 = (E ** 2).sum(1) / df; se = np.sqrt(s2[:, None] * np.diag(XtXi)[None, :])
    t = B / se; p = 2 * stats.t.sf(np.abs(t), df)
    return B[:, 1:], p[:, 1:], (E ** 2).sum(1)
 
# 1. reproduce current numbers
b0, p0, _ = fit(grp[:, None]); b1, p1, _ = fit(np.column_stack([grp, xist]))
say("genes", G, "| PMOS p<0.05", (p0[:, 0] < .05).sum(), "| +XIST", (p1[:, 0] < .05).sum(),
    "| q<0.1", (bh(p0[:, 0]) < .1).sum())
s4 = pd.read_csv(f"{R}/transcriptome/S4_DE_PMOS_vs_control.tsv", sep="\t").set_index("gene").loc[lg.index]
say("max |p diff| vs S4:", np.abs(s4.p.values - p0[:, 0]).max())
v = lg.var(1).sort_values(ascending=False).index[:1000]
Z = lg.loc[v].sub(lg.loc[v].mean(1), axis=0).T.values
U, Sv, Vt = np.linalg.svd(Z, full_matrices=False); ve = Sv ** 2 / np.sum(Sv ** 2)
say("PC var % (1000 genes):", np.round(ve[:5] * 100, 1))
Za = lg.sub(lg.mean(1), axis=0).T.values; Ua, Sa, _ = np.linalg.svd(Za, full_matrices=False)
vea = Sa ** 2 / np.sum(Sa ** 2)
pd.DataFrame((Ua * Sa)[:, :5], index=lines, columns=[f"PC{i}" for i in range(1, 6)]).assign(
    var_pct=None).drop(columns="var_pct").to_csv(f"{OUT}/pcs_allgenes.tsv", sep="\t")
say("PC var % (all 16,489 genes):", np.round(vea[:5] * 100, 1))
 
# 2. sample correlation
C = np.corrcoef(Y.T); iu = np.triu_indices(n, 1)
def blockstat(lab):
    same = lab[:, None] == lab[None, :]
    return C[iu][same[iu]].mean() - C[iu][~same[iu]].mean()
obs = blockstat(grp)
null = []
for c in itertools.combinations(range(n), 7):
    z = np.zeros(n); z[list(c)] = 1; null.append(blockstat(z))
null = np.array(null); pc = (np.sum(null >= obs)) / len(null)
obs_b = blockstat(obese); null_b = []
for c in itertools.combinations(range(n), 7):
    z = np.zeros(n); z[list(c)] = 1; null_b.append(blockstat(z))
pb = np.sum(np.array(null_b) >= obs_b) / len(null_b)
xe = (xist < 8).astype(float); obs_x = blockstat(xe)
null_x = []
for c in itertools.combinations(range(n), 3):
    z = np.zeros(n); z[list(c)] = 1; null_x.append(blockstat(z))
px = np.sum(np.array(null_x) >= obs_x) / len(null_x)
corr_sum = pd.DataFrame([
    dict(grouping="PMOS vs control", within_minus_between_r=obs, exact_p=pc, n_splits=len(null)),
    dict(grouping="BMI >30 vs <25", within_minus_between_r=obs_b, exact_p=pb, n_splits=len(null_b)),
    dict(grouping="XIST-eroded vs XIST-expressing", within_minus_between_r=obs_x, exact_p=px, n_splits=len(null_x))])
same = grp[:, None] == grp[None, :]
say(f"Pearson r range {C[iu].min():.3f}-{C[iu].max():.3f}; within-group mean {C[iu][same[iu]].mean():.4f}, between {C[iu][~same[iu]].mean():.4f}")
say(corr_sum.to_string(index=False))
pd.DataFrame(C, index=lines, columns=lines).to_csv(f"{OUT}/sample_correlation.tsv", sep="\t", float_format="%.4f")
corr_sum.to_csv(f"{OUT}/sample_correlation_tests.tsv", sep="\t", index=False)
 
# 3. variance partition: Type II SS as % of total, PMOS + BMI + age + XIST, with permutation null
SST = ((Y - Y.mean(1, keepdims=True)) ** 2).sum(1)
def vp(Xd):
    _, _, sse_full = fit(Xd); out = {}
    for j in range(Xd.shape[1]):
        _, _, sse_red = fit(np.delete(Xd, j, 1)); out[j] = (sse_red - sse_full) / SST
    out["res"] = sse_full / SST
    return out
facs = ["PMOS status", "BMI", "Age", "XIST (log2)"]
Xd = np.column_stack([grp, bmi, age, xist])
VP = vp(Xd)
vpt = pd.DataFrame({facs[j]: VP[j] for j in range(4)}, index=lg.index); vpt["Residual"] = VP["res"]
rng = np.random.default_rng(2026); NP = 1000
nullmed = np.zeros((NP, 4))
for k in range(NP):
    pm = rng.permutation(n); V = vp(Xd[pm])
    nullmed[k] = [np.median(V[j]) for j in range(4)]
rows = []
for j, f in enumerate(facs):
    m = np.median(VP[j])
    rows.append(dict(factor=f, median_pct=100 * m, null_median_pct=100 * np.median(nullmed[:, j]),
                     null_95pct=100 * np.percentile(nullmed[:, j], 95),
                     perm_p=(np.sum(nullmed[:, j] >= m) + 1) / (NP + 1),
                     genes_above_20pct=int((VP[j] > .2).sum())))
rows.append(dict(factor="Residual", median_pct=100 * np.median(VP["res"])))
vps = pd.DataFrame(rows)
# the September model (PMOS + BMI + age, no XIST) for comparison
VP3 = vp(np.column_stack([grp, bmi, age]))
for j, f in enumerate(["PMOS status", "BMI", "Age"]):
    vps.loc[vps.factor == f, "median_pct_model_without_XIST"] = 100 * np.median(VP3[j])
vps.loc[vps.factor == "Residual", "median_pct_model_without_XIST"] = 100 * np.median(VP3["res"])
say(vps.round(3).to_string(index=False))
vpt.to_csv(f"{OUT}/variance_partition_per_gene.tsv.gz", sep="\t", float_format="%.5f")
vps.to_csv(f"{OUT}/variance_partition_summary.tsv", sep="\t", index=False)
np.save(f"{OUT}/vp_nullmed.npy", nullmed)
 
# 4. BMI-adjusted DE and interaction
bA, pA, _ = fit(np.column_stack([grp, obese]))
bC, pC, _ = fit(np.column_stack([grp, bmi]))
bI, pI, _ = fit(np.column_stack([grp, obese, grp * obese]))
u = p0[:, 0] < .05; a = pA[:, 0] < .05; c = pC[:, 0] < .05
qI = bh(pI[:, 2])
bmi_sum = pd.DataFrame([
    dict(model="PMOS status only (main analysis)", genes_p05=int(u.sum()), genes_FDR10=int((bh(p0[:, 0]) < .1).sum())),
    dict(model="PMOS + BMI category (lean <25 / obese >30)", genes_p05=int(a.sum()), genes_FDR10=int((bh(pA[:, 0]) < .1).sum()),
         overlap_with_main=int((u & a).sum()), only_in_this_model=int((a & ~u).sum()), lost_from_main=int((u & ~a).sum())),
    dict(model="PMOS + BMI continuous (sensitivity)", genes_p05=int(c.sum()), genes_FDR10=int((bh(pC[:, 0]) < .1).sum()),
         overlap_with_main=int((u & c).sum()), only_in_this_model=int((c & ~u).sum()), lost_from_main=int((u & ~c).sum())),
    dict(model="PMOS x BMI category interaction term", genes_p05=int((pI[:, 2] < .05).sum()),
         genes_FDR10=int((qI < .1).sum()), genes_FDR05=int((qI < .05).sum()))])
bmi_sum["expected_by_chance"] = round(0.05 * G)
say(bmi_sum.to_string(index=False))
say("interaction FDR<0.1 genes:", list(lg.index[qI < .1]))
bmide = pd.DataFrame({"gene": lg.index,
    "log2FC_unadjusted": b0[:, 0], "p_unadjusted": p0[:, 0],
    "log2FC_adjBMIcat": bA[:, 0], "p_adjBMIcat": pA[:, 0], "q_adjBMIcat": bh(pA[:, 0]),
    "log2FC_adjBMIcont": bC[:, 0], "p_adjBMIcont": pC[:, 0],
    "interaction_beta": bI[:, 2], "interaction_p": pI[:, 2], "interaction_q": qI}).sort_values("p_adjBMIcat")
bmide.to_csv(f"{OUT}/bmi_adjusted_DE.tsv.gz", sep="\t", index=False, float_format="%.6g")
bmi_sum.to_csv(f"{OUT}/bmi_adjusted_summary.tsv", sep="\t", index=False)
 
# 5. BMI-stratified donor characteristics
ph["g"] = np.where(ph.Group == "CTRL", "Control", "PMOS")
vars_ = [("AGE", "Age (years)"), ("BMI", "BMI (kg/m²)"), ("FGscore", "Ferriman–Gallwey score"), ("AMH", "AMH (ng/mL)"),
         ("Testo", "Testosterone (pmol/L)"), ("DHT", "DHT (pmol/L)"), ("SHBG", "SHBG (nmol/L)"), ("FAI", "Free androgen index"),
         ("LH", "LH (IU/L)"), ("HbA1c", "HbA1c (mmol/mol)"), ("Gluk_mM", "Fasting glucose (mmol/L)"),
         ("Insuliini_mIU_L", "Fasting insulin (mIU/L)"), ("HOMA_IR", "HOMA-IR")]
def fmt(v):
    f = lambda z: f"{z:.3g}" if abs(z) < 100 else f"{z:.0f}"
    return f"{f(np.median(v))} ({f(v.min())}–{f(v.max())})"
strat = []
for sname, mask in [("Lean (BMI <25)", ph.BMI < 25), ("Obese (BMI >30)", ph.BMI > 30)]:
    d = ph[mask]
    for c, lab in vars_:
        a_ = d[d.g == "Control"][c].astype(float); b_ = d[d.g == "PMOS"][c].astype(float)
        strat.append(dict(stratum=sname, characteristic=lab, control=f"{fmt(a_.values)} [n={len(a_)}]",
                          PMOS=f"{fmt(b_.values)} [n={len(b_)}]",
                          p_MannWhitney=stats.mannwhitneyu(a_, b_, alternative="two-sided").pvalue))
strat = pd.DataFrame(strat); strat.to_csv(f"{OUT}/bmi_stratified_donors.tsv", sep="\t", index=False)
say(strat[strat.characteristic.isin(["AMH (ng/mL)", "Testosterone (pmol/L)", "Free androgen index", "SHBG (nmol/L)"])].to_string(index=False))
 
# 6. HipSci protein cross-reference (values from the Sept supplementary table; Mirauta et al.)
hp = pd.read_excel(f"{SEP}/tables_and_supplementary/PMOS_iPSC_supplementary_FINAL.xlsx", sheet_name="TableS2_Confirmed_Protein_Genes")
hp = hp[["gene_symbol", "uniprot_id", "hipsci_mean_transcript", "hipsci_mean_protein", "cv_protein",
         "hipsci_rna_protein_corr", "pqtl", "eqtl"]].drop_duplicates("gene_symbol").set_index("gene_symbol")
mc = expd.mean(1)
prot = pd.DataFrame({"gene": lg.index, "mean_norm_count": mc.values, "log2FC_PMOS_vs_control": b0[:, 0],
                     "p": p0[:, 0], "q_BH": bh(p0[:, 0]), "p_adjXIST": p1[:, 0]}).set_index("gene")
prot = prot.join(hp, how="left"); prot["protein_detected_HipSci"] = prot.hipsci_mean_protein.notna()
prot["expression_decile"] = pd.qcut(prot.mean_norm_count.rank(method="first"), 10, labels=range(1, 11)).astype(int)
dec = prot.groupby("expression_decile").agg(genes=("p", "size"), protein_detected=("protein_detected_HipSci", "sum"),
      median_count=("mean_norm_count", "median")).reset_index()
dec["pct_detected"] = 100 * dec.protein_detected / dec.genes
cand = prot[(prot.p < .05) & prot.protein_detected_HipSci & (prot.hipsci_rna_protein_corr >= .3)].sort_values("p")
old30 = pd.read_excel(f"{SEP}/tables_and_supplementary/PMOS_iPSC_supplementary_FINAL.xlsx", sheet_name="TableS2_Best_Supported_Cand").gene_symbol
say(f"protein detected: {prot.protein_detected_HipSci.sum()} of {G} ({100*prot.protein_detected_HipSci.mean():.1f}%); "
    f"candidates p<0.05 & r>=0.3: {len(cand)}; overlap with Sept 30: {len(set(cand.index)&set(old30))}")
say(dec.round(1).to_string(index=False))
prot.reset_index().sort_values("p").to_csv(f"{OUT}/protein_crossref.tsv.gz", sep="\t", index=False, float_format="%.6g")
cand.reset_index().to_csv(f"{OUT}/protein_candidates.tsv", sep="\t", index=False, float_format="%.6g")
dec.to_csv(f"{OUT}/protein_by_decile.tsv", sep="\t", index=False)
gw = pd.read_csv(f"{R}/transcriptome/S7_PCOS_GWAS_genes_expression.tsv", sep="\t")
extra = ["RPS26", "SUOX", "IKZF4", "CDK2", "PA2G4", "IRX3", "RPGRIP1L"]
loc = pd.DataFrame({"gene": list(gw.gene) + extra,
                    "source": ["consensus gene (Table S7)"] * len(gw) + ["eQTL-tested neighbour (Table S6)"] * len(extra)})
loc = loc.drop_duplicates("gene").set_index("gene").join(hp[["hipsci_mean_protein", "hipsci_rna_protein_corr", "pqtl", "eqtl"]], how="left")
loc["expressed_iPSC"] = loc.index.isin(lg.index); loc["protein_detected_HipSci"] = loc.hipsci_mean_protein.notna()
loc.reset_index().to_csv(f"{OUT}/protein_pcos_locus_genes.tsv", sep="\t", index=False)
e = loc[loc.source.str.startswith("consensus") & loc.expressed_iPSC]
say(f"PCOS consensus genes expressed: {len(e)}, protein-detected: {e.protein_detected_HipSci.sum()}: {list(e.index[e.protein_detected_HipSci])}")
say("eQTL-tested neighbours:", loc.loc[extra, ["protein_detected_HipSci", "pqtl", "eqtl"]].to_string())
say("ERBB3:", loc.loc["ERBB3", ["protein_detected_HipSci"]].to_dict())
 
# 7. literature candidate screen (60 genes from Sept Table S1)
lit = pd.read_excel(f"{SEP}/tables_and_supplementary/PMOS_iPSC_supplementary_FINAL.xlsx", sheet_name="TableS1_GWAS_PathobiologyScreen", header=3)
lit = lit[lit.Category.isin(["GWAS locus", "Pathobiology"])][["Gene", "Category"]]
notexp = ["ADIPOQ", "C9orf3", "CYP17A1", "CYP21A2", "FSHB", "FSHR", "INS", "LINC02905", "MC4R", "SUMO1P1", "TNF"]
lit = pd.concat([lit, pd.DataFrame({"Gene": notexp, "Category": "not expressed in Sept filter"})]).drop_duplicates("Gene")
ix = pd.Series(range(G), index=lg.index)
qs0 = bh(p0[:, 0]); qs1 = bh(p1[:, 0])
rows = []
for gname, cat in lit.values:
    if gname in ix:
        i = ix[gname]
        rows.append(dict(gene=gname, category=cat, expressed=True, mean_norm_count=mc.iloc[i], log2FC_PMOS_vs_control=b0[i, 0],
                         p=p0[i, 0], q_BH_transcriptome=qs0[i], log2FC_adjXIST=b1[i, 0], p_adjXIST=p1[i, 0],
                         protein_detected_HipSci=gname in hp.index))
    else:
        rows.append(dict(gene=gname, category=cat, expressed=False,
                         mean_norm_count=cnt.loc[gname].mean() if gname in cnt.index else np.nan,
                         protein_detected_HipSci=gname in hp.index))
lit = pd.DataFrame(rows).sort_values(["expressed", "p"], ascending=[False, True])
lit.to_csv(f"{OUT}/literature_candidates.tsv", sep="\t", index=False, float_format="%.6g")
say(f"literature genes: {len(lit)}, expressed (16,489 filter): {lit.expressed.sum()}, nominal p<0.05: {(lit.p<0.05).sum()}, min q: {lit.q_BH_transcriptome.min():.3f}")
say(lit[lit.p < .05][["gene", "category", "log2FC_PMOS_vs_control", "p"]].to_string(index=False))
 
# 8. top 40 genes for heatmap; old annotation carried where gene still in top list
top = pd.DataFrame({"gene": lg.index, "log2FC": b0[:, 0], "p": p0[:, 0], "mean_norm_count": mc.values}).sort_values("p").head(40)
ann = pd.read_excel(f"{SEP}/tables_and_supplementary/PMOS_iPSC_manuscript_tables_FINAL.xlsx", sheet_name="Table5_GO_Annotation").set_index("Gene")
top["functional_category_Sept"] = top.gene.map(ann["Functional category"])
top.to_csv(f"{OUT}/top40_genes.tsv", sep="\t", index=False, float_format="%.6g")
say("top-40 genes with Sept annotation:", top.functional_category_Sept.notna().sum())
lg.loc[top.gene].to_csv(f"{OUT}/top40_log2.tsv", sep="\t")
# volcano data
pd.DataFrame({"gene": lg.index, "log2FC": b0[:, 0], "p": p0[:, 0]}).to_csv(f"{OUT}/volcano.tsv.gz", sep="\t", index=False)
pd.DataFrame({"line": lines, "group": key.loc[lines, "group"].values, "bmi": bmi, "obese": obese, "xist_log2": xist}).to_csv(f"{OUT}/line_annot.tsv", sep="\t", index=False)
log.close()
