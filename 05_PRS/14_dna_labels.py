#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 05_PRS/14_dna_labels.py   (original step script: 48b_dna_labels.py)
# Produces: Case/control labels from the DNA-extraction sheet (TSV dump).
# Run: python 05_PRS/14_dna_labels.py <args>   (or from a SLURM script; see the usage line below)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 48b - case/control labels of the 140 donors from the DNA-extraction sheet of the lab (DNA_extraction_blood_141_samples_iPSC_090626.xlsx,
sheet 'Oulu Bank for genotyping', converted to TSV by 00_common/01_xlsx_to_tsv.py). The sheet holds three side-by-side blocks of columns
(Tube, Group, ng/ul, 260/280, 260/230); Tube = donor / eSF ID (1-141 without 37), Group = PMOS or CTRL.
Tubes 22-30 carry an asterisk ("22*" ... "30*") and the lab note "Samples 22-30 maybe problems with Nanodrop". The asterisk and the note refer to the
DNA MEASUREMENT (anomalous 260/280), NOT to the diagnosis: these 9 donors HAVE a group and passed the genotype QC.
CORRECTION 2026-11-05: the first version of this script treated tubes 22-30 as unlabelled. Now the asterisk is stripped from the tube number and the
group of 22-30 is used in the case-control comparison (22, 23, 24, 26, 27, 28 = PMOS; 25, 29, 30 = CTRL; confirmed by the analyst).
The column nanodrop_flag_22_30 is kept for transparency only (yes = the tube had an asterisk); it does NOT change group_cc.
Tubes repeated in the sheet (a second table lists the donors sent for reprogramming) must carry the same group, otherwise the script stops.
Output: TSV  esf_id, group_sheet, group_cc (= group_sheet), nanodrop_flag_22_30.  Standard library only (runs locally or on the cluster).
Usage: python 05_PRS/14_dna_labels.py <sheet01.tsv> <out.tsv> [keep_129.txt]   (optional keep list: also prints the PMOS/CTRL counts within those samples)
"""
import re, sys

rows, dups = {}, []
for line in open(sys.argv[1], encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    for i, cell in enumerate(f):
        m = re.fullmatch(r"(\d+)\*?", cell.strip())   # tube number; the optional asterisk is dropped by the group (1)
        if m and i + 1 < len(f) and f[i + 1].strip() in ("PMOS", "CTRL"):
            tube, grp = int(m.group(1)), f[i + 1].strip()
            if tube in rows:   # the lower table of the sheet repeats the donors sent for reprogramming: accept the repeat ONLY if the group agrees
                assert rows[tube][0] == grp, f"tube {tube}: conflicting groups {rows[tube][0]} / {grp}"
                dups.append(tube)
                continue
            rows[tube] = (grp, "*" in cell)
flagged = sorted(t for t, (g, star) in rows.items() if star)
assert flagged == list(range(22, 31)), f"unexpected asterisk tubes: {flagged}"   # the sheet marks exactly 22-30

with open(sys.argv[2], "w") as out:
    out.write("esf_id\tgroup_sheet\tgroup_cc\tnanodrop_flag_22_30\n")
    for t in sorted(rows):
        g, star = rows[t]
        out.write(f"{t}\t{g}\t{g}\t{'yes' if star else 'no'}\n")      # group_cc = group_sheet for ALL tubes (asterisk is about the DNA reading)

cnt = lambda ts: f"PMOS {sum(rows[t][0] == 'PMOS' for t in ts)}, CTRL {sum(rows[t][0] == 'CTRL' for t in ts)}"
print(f"tubes: {len(rows)} | all labelled: {cnt(rows)} | tubes with asterisk (22-30): {flagged} -> {cnt(flagged)} | "
      f"missing tubes 1-141: {[t for t in range(1, 142) if t not in rows]} | repeated tubes with the same group (ignored): {sorted(dups)}")
if len(sys.argv) > 3:
    keep = {int(l.split()[1]) for l in open(sys.argv[3]) if l.strip() and not l.startswith("#")}
    print(f"samples in {sys.argv[3]}: {len(keep)} | labelled: {cnt([t for t in keep if t in rows])} | without label: {sorted(t for t in keep if t not in rows)}")
