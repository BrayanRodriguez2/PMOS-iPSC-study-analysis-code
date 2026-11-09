#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 02_contamination_identity/12_ibs0_donor_profile.py   (original step script: 31_ibs0_donor_profile.py)
# Produces: IBS0 profile of focal donors against the rest of the cohort.
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 31 - Independent signature of excess heterozygosity in the Part B matrix (diagnostic only; no QC data modified).

A heterozygous genotype can never be "IBS0" (0/0 against 1/1) with anyone. A donor whose array calls have too many false hets
(DNA mixture) therefore looks CLOSER to every RNA line, i.e. its COLUMN of the IBS0 matrix (step 28) is lower than the others, even
though none of the lines is its owner. For every donor we average IBS0 over the 13 or 14 lines that are NOT its own, and compare with
the other donors (robust z, median/MAD). The same average is computed for the concordance. Expect: low IBS0 in the donors with F < 0
(the focal samples, ...), high IBS0 for donors with a deficit of hets (F > 0).
Input : results/rna_identity/ibs0_matrix.tsv, concordance_matrix.tsv (labels "HELxxx.x (donor N, GROUP)"), sample_qc_table.tsv
Output: results/rna_identity/ibs0_donor_profile.tsv, figures/B4_ibs0_donor_profile.{png,pdf}
"""
import argparse, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--identity-dir", required=True); ap.add_argument("--qc-table", required=True)
a = ap.parse_args()
ib = pd.read_csv(os.path.join(a.identity_dir, "ibs0_matrix.tsv"), sep="\t", index_col=0)
cc = pd.read_csv(os.path.join(a.identity_dir, "concordance_matrix.tsv"), sep="\t", index_col=0)
own = [re.search(r"donor (\d+)", i).group(1) for i in ib.index]
ib.columns = ib.columns.astype(str); cc.columns = cc.columns.astype(str)
qc = pd.read_csv(a.qc_table, sep="\t", dtype={"IID": str}).set_index("IID")
rows = []
for d in ib.columns:
    keep = [i for i, o in enumerate(own) if o != d]
    rows.append(dict(donor=d, mean_ibs0_other_lines=ib[d].to_numpy()[keep].mean(), mean_conc_other_lines=cc[d].to_numpy()[keep].mean(),
                     n_lines=len(keep), het_F=qc.loc[d, "het_F"], call_rate=qc.loc[d, "call_rate"], excluded=qc.loc[d, "excluded"]))
t = pd.DataFrame(rows)
m = t.mean_ibs0_other_lines.median(); s = 1.4826 * (t.mean_ibs0_other_lines - m).abs().median()
t["ibs0_robust_z"] = (t.mean_ibs0_other_lines - m) / s
t.to_csv(os.path.join(a.identity_dir, "ibs0_donor_profile.tsv"), sep="\t", index=False, float_format="%.5g")
print("correlation (Spearman) mean IBS0 vs het_F over 140 donors:", round(t[["mean_ibs0_other_lines", "het_F"]].corr(method="spearman").iloc[0, 1], 3))
print("lowest IBS0 (closest to every line):\n", t.sort_values("mean_ibs0_other_lines").head(10).round(4).to_string(index=False))
print("highest IBS0:\n", t.sort_values("mean_ibs0_other_lines").tail(6).round(4).to_string(index=False))
fig, ax = plt.subplots(figsize=(6.4, 4.8))
col = np.where(t.excluded == "yes", "#c0392b", "#2c7fb8")
ax.scatter(t.het_F, t.mean_ibs0_other_lines * 100, c=col, s=22)
for r in t.itertuples():
    if r.donor in {x for x in __import__("os").environ.get("ANNOTATE_SAMPLE_IDS", "").split(",") if x}:
        ax.annotate(r.donor, (r.het_F, r.mean_ibs0_other_lines * 100), fontsize=8, xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("heterozygosity F of the donor (array, step 09)"); ax.set_ylabel("mean IBS0 vs the lines of OTHER donors (%)")
ax.set_title("Excess hets make a donor look 'closer' to every RNA line\n(red = fails sample QC)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(a.identity_dir, "figures", "B4_ibs0_donor_profile.png"), dpi=170)
fig.savefig(os.path.join(a.identity_dir, "figures", "B4_ibs0_donor_profile.pdf")); plt.close(fig)
