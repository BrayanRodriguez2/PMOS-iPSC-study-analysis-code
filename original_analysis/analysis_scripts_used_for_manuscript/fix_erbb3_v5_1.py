"""v5.1 (4 Nov 2026): correct the ERBB3 diagnosis-adjustment statements found by independent verification.
The reported beta 0.27, p 0.03 came from ERBB3 ~ dosage + diagnosis (no XIST). With the main-model XIST
covariate plus diagnosis: beta 0.17, p 0.19; controls only rho 0.62, p 0.19 (n = 6).
Usage: python3 fix_erbb3_v5_1.py <unpacked_dir> <highlight: magenta|none>"""
import sys, re
from xml.sax.saxutils import escape
D, HLV = sys.argv[1], sys.argv[2]
p = f"{D}/word/document.xml"; x = open(p, encoding="utf-8").read()
HL = '<w:highlight w:val="magenta"/>' if HLV == "magenta" else ""
def sub(old, new):
    global x
    eo = escape(old)
    ms = [m for m in re.finditer(r"<w:r>(?:<w:rPr>((?:(?!</w:rPr>|<w:r>).)*)</w:rPr>)?<w:t( xml:space=\"preserve\")?>([^<]*)</w:t></w:r>", x, flags=re.S) if eo in m.group(3)]
    assert len(ms) == 1, (old[:50], len(ms)); m = ms[0]
    rpr = m.group(1) or ""; a, b = m.group(3).split(eo, 1)
    keep = lambda s: f'<w:r><w:rPr>{rpr}</w:rPr><w:t xml:space="preserve">{s}</w:t></w:r>' if s else ""
    nr = re.sub(r'<w:highlight w:val="\w+"/>', "", rpr) + HL
    x = x[:m.start()] + keep(a) + f'<w:r><w:rPr>{nr}</w:rPr><w:t xml:space="preserve">{escape(new)}</w:t></w:r>' + keep(b) + x[m.end():]
sub("All seven PMOS donors were heterozygous at rs10876864, so we refitted the model with diagnosis as a covariate. The ERBB3 association persisted (β = 0.27, p = 0.03) while diagnosis did not (p = 0.22), and the direction held among controls alone (ρ = 0.62, n = 6).",
    "All seven PMOS donors were heterozygous at rs10876864, so genotype and diagnosis are partly confounded. With diagnosis instead of XIST as covariate the association remained (β = 0.27, p = 0.03; diagnosis p = 0.22), but with both covariates it weakened (β = 0.17, p = 0.19), and among controls alone the direction was the same but not significant (ρ = 0.62, p = 0.19, n = 6). In these donors the genotype effect therefore cannot be separated from diagnosis.")
sub("The signal held after adjustment for diagnosis, which matters because all PMOS donors happened to be heterozygous.",
    "Because all PMOS donors happened to be heterozygous, the signal could not be separated from diagnosis: it remained after adjustment for diagnosis alone but weakened when XIST and diagnosis were both included.")
open(p, "w", encoding="utf-8").write(x); print("fixed", D)
