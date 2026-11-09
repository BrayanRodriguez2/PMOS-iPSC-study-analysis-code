#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Release file: 04_post_imputation_QC/10_make_keep_lists.sh   (original step script: 45_make_keep_lists.sh)
# Produces: Sample lists (main QC set and sensitivity set) for PLINK --keep.
# Run: bash 04_post_imputation_QC/10_make_keep_lists.sh <args>   (light job; source config/config.sh first)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
# Step 45 (prep) - keep lists for the clean imputed file: keep_129.txt (main QC: copy of results/qc/samples/keep_samples.txt, FID IID) and
# keep_135.txt (sensitivity: call rate >= 0.97, het 3 SD -> the samples excluded in that row of callrate_sensitivity.tsv are removed from the
# 140 donors of the v2 VCFs). Both in PLINK format (FID IID, FID = IID). Checks that the 129 are a subset of the 135.
: "${PROJECT_DIR:?Set PROJECT_DIR and CODE_DIR first: source config/config.sh}"
set -eo pipefail
PRS="${PROJECT_DIR}/PRS"; O="$PRS/imputation/clean"; mkdir -p "$O"
S="$PRS/validation/results/qc"
eval "${CONDA_SETUP:-true}"; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate idat
# NOTE: in the array PLINK data FID is the row number and IID the eSF ID (they differ after the missing eSF 37), but the imputed VCFs only have
# the eSF ID, so the clean file uses FID = IID (--double-id): the 129 list is rebuilt from the IID column only.
{ echo "#FID IID"; awk -F'\t' 'NR>1 {print $2"\t"$2}' "$S/samples/keep_samples.txt"; } > "$O/keep_129.txt"
bcftools query -l "$PRS/topmed/upload_v2/chr1.vcf.gz" > "$O/work_all140.txt"
ex=$(awk -F'\t' 'NR>1 && $1==3 && $2==0.97 {print $5}' "$S/sensitivity/callrate_sensitivity.tsv" | grep -o '[0-9]\+(' | tr -d '(' | tr '\n' ' ')
echo "samples excluded in the call-rate 0.97 / het 3 SD sensitivity row: $ex"
{ echo "#FID IID"; for s in $(cat "$O/work_all140.txt"); do case " $ex " in *" $s "*) ;; *) echo "$s $s";; esac; done; } | sed 's/ /\t/' > "$O/keep_135.txt"
echo "keep_129: $(( $(wc -l < "$O/keep_129.txt") - 1 )) | keep_135: $(( $(wc -l < "$O/keep_135.txt") - 1 ))"
echo "129 not in 135: $(comm -23 <(tail -n +2 "$O/keep_129.txt" | cut -f2 | sort) <(tail -n +2 "$O/keep_135.txt" | cut -f2 | sort) | wc -l) (must be 0)"
echo "135 not in 129 (the 6 samples that only the sensitivity set includes): $(comm -13 <(tail -n +2 "$O/keep_129.txt" | cut -f2 | sort) <(tail -n +2 "$O/keep_135.txt" | cut -f2 | sort) | tr '\n' ' ')"
{ echo "#FID IID SEX"; for s in $(cat "$O/work_all140.txt"); do echo "$s	$s	2"; done; } > "$O/sex_all_female.txt"     # the 140 donors are genetically female (contexto 2026-11-03)
echo "sex file lines: $(( $(wc -l < "$O/sex_all_female.txt") - 1 ))"
head -n 3 "$O/keep_129.txt" | cat -A | head -3
