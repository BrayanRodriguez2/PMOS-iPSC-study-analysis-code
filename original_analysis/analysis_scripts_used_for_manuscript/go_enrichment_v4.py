"""v4 (4 Nov 2026): GO / pathway enrichment for the HipSci protein cross-reference (Table S13e).
 
Gene sets (HGNC symbols, from S13d of an_protein_v4.py; protein-coding expressed genes only):
  detected     expressed genes with protein detected in HipSci     vs background: all expressed protein-coding genes
  notdetected  expressed genes without detected protein             vs background: all expressed protein-coding genes
  candidates   49 protein-supported PMOS candidates (S13c rule)     vs background: protein detected and r >= 0.3
  nominal      nominal PMOS p < 0.05 with protein detected          vs background: all protein-detected genes
g:Profiler REST API (https://biit.cs.ut.ee/gprofiler/api/gost/profile/), organism hsapiens,
sources GO:BP, GO:MF, GO:CC, KEGG, REAC; g:SCS multiple-testing correction, threshold 0.05;
domain_scope = custom (background as listed). Kolberg et al., Nucleic Acids Res 2023;51:W207-W212.
 
NOTE: in the analysis environment the identical requests were sent from a browser session (the
analysis server had no direct internet access); this script reproduces them and needs internet.
g:Profiler database versions change over time, so term counts may differ slightly on re-run.
 
Usage: python3 go_enrichment_v4.py <v4_results_dir>
"""
import sys, json
import pandas as pd, requests
V4 = sys.argv[1]
API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
d = pd.read_csv(f"{V4}/S13d_protein_crossref_all_genes.tsv.gz", sep="\t")
pc = d[d.protein_coding]
det = pc[pc.protein_detected_HipSci]
elig = det[det.hipsci_rna_protein_corr >= .3]
U = lambda s: sorted(set(s.hgnc_symbol))
sets = {
    "detected": (U(det), U(pc)),
    "notdetected": (U(pc[~pc.protein_detected_HipSci]), U(pc)),
    "candidates": (U(elig[elig.p < .05]), U(elig)),
    "nominal": (U(det[det.p < .05]), U(det)),
}
 
def gost(query, background):
    body = dict(organism="hsapiens", query=query, sources=["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
                user_threshold=0.05, significance_threshold_method="g_SCS", domain_scope="custom",
                background=background, no_evidences=True)
    r = requests.post(API, json=body, timeout=300); r.raise_for_status()
    return r.json()
 
rows, meta = [], {}
for name, (q, bg) in sets.items():
    j = gost(q, bg)
    meta[name] = dict(query=len(q), background=len(bg), significant_terms=len(j["result"]),
                      gprofiler_version=j["meta"]["version"])
    for t in j["result"]:
        rows.append(dict(set=name, source=t["source"], term_id=t["native"], term_name=t["name"],
                         p_adj_gSCS=t["p_value"], intersection=t["intersection_size"], term_size=t["term_size"],
                         query_size=t["query_size"], background=t["effective_domain_size"]))
    print(name, meta[name])
df = pd.DataFrame(rows)
df["fold_enrichment"] = (df.intersection / df.term_size) / (df.query_size / df.background)
df.to_csv(f"{V4}/GO_enrichment_rerun.tsv", sep="\t", index=False)
json.dump(meta, open(f"{V4}/GO_enrichment_rerun_meta.json", "w"), indent=1)
