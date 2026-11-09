// v4 (4 Nov 2026): the g:Profiler request exactly as executed in the browser (page origin biit.cs.ut.ee).
// Gene lists were decoded in the page from a 2-bit code over the sorted HGNC protein-coding list
// (0 = not expressed, 1 = expressed/no protein, 2 = protein detected, 3 = protein detected and PMOS p < 0.05),
// giving 13,240 background, 8,726 detected, 4,514 not detected and 233 nominal genes.
// The 49-candidate set was queried against its 1,962-gene eligible background in the same way.
// go_enrichment_v4.py reproduces all four requests from Python.
async function gost(q, bg) {
  const r = await fetch('https://biit.cs.ut.ee/gprofiler/api/gost/profile/', {
    method: 'POST', headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({organism: 'hsapiens', query: q, sources: ['GO:BP','GO:MF','GO:CC','KEGG','REAC'],
      user_threshold: 0.05, significance_threshold_method: 'g_SCS', domain_scope: 'custom',
      background: bg, no_evidences: true})});
  return (await r.json()).result;
}
// detected vs background, not detected vs background, nominal vs detected:
// await gost(window.__det, window.__bgPC); await gost(window.__nd, window.__bgPC); await gost(window.__nom, window.__det);
