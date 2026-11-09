#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/08_external_pipeline_check.py   (original step script: 66a_their_inputs.py)
# Produces: Runs the original analysis scripts on the original counts (requires original_analysis/).
# Run: called by the SLURM script of the same number; usage in the docstring below
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 66a - lay out the inputs that the original analysis' own scripts (Additional_file_3: an1.py, an2.py, an_restore.py, common.py) expect, using ONLY files of the original package:
  * their_run/Normalized_counts.tabular   (the original Galaxy size-factor-normalised counts; copy of the file she sent, unchanged)
  * their_run/pheno/sample_key_14_ipsc_donors.csv + donor_phenotypes.tsv   (rebuilt from the original Table S1 = sheet S1_donors of Additional_file_1; columns as the scripts read them:
    array_iid, hel_line, group (CTRL/PMOS), bmi, age, amh, testosterone, fai, fg_score, homa_ir)
  * their_run/pkg/4_Analysis_results/transcriptome/S4_DE_PMOS_vs_control.tsv   (the original Table S4, so that an_restore.py can compare)
The original scripts are used verbatim (see 66b); only their hard-coded paths were replaced. This reproduces the original numbers from the original inputs: a check of the original code, not an independent analysis.
Usage: python 06_RNAseq_counts_and_DE_validation/08_external_pipeline_check.py <package_dir> <their_run_dir>
"""
import os, sys
import pandas as pd

pkg, run = sys.argv[1:3]
os.makedirs(f"{run}/pheno", exist_ok=True); os.makedirs(f"{run}/ms/tables", exist_ok=True); os.makedirs(f"{run}/scripts", exist_ok=True)
os.makedirs(f"{run}/pkg/4_Analysis_results/metadata", exist_ok=True); os.makedirs(f"{run}/pkg/4_Analysis_results/transcriptome", exist_ok=True)
os.makedirs(f"{run}/pkg/5_Pipeline_CSC", exist_ok=True)
s1 = pd.read_csv(f"{pkg}/af1_tsv/sheet02_S1_donors.tsv", sep="\t")
s1 = s1[s1.iPSC_line_RNAseq.astype(str).str.startswith("HEL")]
key = pd.DataFrame({"array_iid": s1.Genotype_array_ID.astype(int).astype(str), "hel_line": s1.iPSC_line_RNAseq,
                    "group": s1.Group.map({"Control": "CTRL", "PMOS": "PMOS"}), "bmi": s1.BMI, "age": s1.AGE, "amh": s1.AMH, "testosterone": s1.Testo,
                    "fai": s1.FAI, "fg_score": s1.FGscore, "homa_ir": s1.HOMA_IR})
assert key.group.notna().all() and len(key) == 14
key.to_csv(f"{run}/pheno/sample_key_14_ipsc_donors.csv", index=False)
key.to_csv(f"{run}/pkg/4_Analysis_results/metadata/sample_key_14_ipsc_donors.csv", index=False)
ph = s1.rename(columns={"Group": "G"}).assign(Group=key.group.values)
ph.to_csv(f"{run}/pheno/donor_phenotypes.tsv", sep="\t", index=False)
ph.to_csv(f"{run}/pkg/4_Analysis_results/metadata/donor_phenotypes.tsv", sep="\t", index=False)
s4 = pd.read_csv(f"{pkg}/af1_tsv/sheet05_S4_DE.tsv", sep="\t")
s4.to_csv(f"{run}/pkg/4_Analysis_results/transcriptome/S4_DE_PMOS_vs_control.tsv", sep="\t", index=False)
print(f"key: {len(key)} donors ({(key.group == 'PMOS').sum()} PMOS) | S4 rows: {len(s4)}")
