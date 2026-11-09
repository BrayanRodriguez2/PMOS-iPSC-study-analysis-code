#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Release file: 00_common/01_xlsx_to_tsv.py   (original step script: 48_xlsx_to_tsv.py)
# Produces: Dumps every sheet of an .xlsx to TSV (standard library only); used to read supplementary tables.
# Run: python 00_common/01_xlsx_to_tsv.py <args>   (or from a SLURM script; see the usage line below)
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
"""
Step 48a - dump every sheet of the paper's Supplementary Tables workbook (xlsx) to a TSV, with the standard library only
(no openpyxl in the cluster environments). Used to read Supplementary Table 2 (ST2: the 29 lead SNPs with risk allele and effect).
Reads xl/workbook.xml, xl/_rels/workbook.xml.rels, xl/sharedStrings.xml and the sheet XMLs; inline strings and shared strings are handled; empty cells
are written as empty fields; cell references (A1, B1, ...) give the column so gaps are kept.
Usage: python 00_common/01_xlsx_to_tsv.py <in.xlsx> <outdir>
"""
import os, re, sys, zipfile
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "pr": "http://schemas.openxmlformats.org/package/2006/relationships"}


def col_index(ref):
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n - 1


def main():
    xlsx, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    z = zipfile.ZipFile(xlsx)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels")).findall("pr:Relationship", NS)}
    for i, sh in enumerate(wb.find("m:sheets", NS).findall("m:sheet", NS), 1):
        name = sh.get("name")
        target = rels[sh.get("{%s}id" % NS["r"])]
        path = target.lstrip("/") if target.startswith("/") else "xl/" + target
        rows = []
        for ev, el in ET.iterparse(z.open(path)):
            if el.tag != "{%s}row" % NS["m"]:
                continue
            cells = {}
            for c in el.findall("m:c", NS):
                t = c.get("t")
                v = c.find("m:v", NS)
                if t == "s" and v is not None:
                    val = shared[int(v.text)]
                elif t == "inlineStr":
                    val = "".join(x.text or "" for x in c.iter("{%s}t" % NS["m"]))
                else:
                    val = v.text if v is not None else ""
                cells[col_index(c.get("r"))] = val.replace("\t", " ").replace("\n", " ")
            if cells:
                rows.append([cells.get(k, "") for k in range(max(cells) + 1)])
            el.clear()
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
        fn = os.path.join(outdir, f"sheet{i:02d}_{safe}.tsv")
        with open(fn, "w") as out:
            for r in rows:
                out.write("\t".join(r) + "\n")
        print(f"{i:02d} '{name}': {len(rows)} rows -> {fn}")


if __name__ == "__main__":
    main()
