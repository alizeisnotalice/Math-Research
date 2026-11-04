# Eighth-order paper: assembly and actual PDF verification

The final paper is `paper/main.tex`, `sections/global.tex`, `sections/sharp.tex`, `sections/eighth.tex`, and `sources.bib`. The final PDF has 53 A4 pages. The result identifies a positive two-sided eighth-order scale above the full P7 baseline; neither an exact eighth-order coefficient nor existence of the normalized eighth-order limit is asserted.

The unchanged foundations are global/sharp/bib. The main statements and all-m transfer were reviewed separately from the complete new module. An additional adversarial review checks the new seventh/eighth-prefix route. The layout review binds the mathematical module to the final source by exact reversible replacement blocks. All four reports are in `reviews/048-*.md`.

`layout-replacements.json` retains the exact twelve old/new blocks, base SHA, and final SHA; private local target paths from the working patch specification were removed. The blocks themselves are unchanged. Inverting them on the final eighth.tex bytes recovers the fully reviewed mathematical source with SHA `f19499528e7fe5d322cf2a2d0333d7743c2fd6ffa3ac059eaa5798adc34b98cd`. This archive does not repeat the complete frozen module.

From this directory, verify the exact reconstruction without writing any source:

```sh
python3 - <<'PYVERIFY'
from pathlib import Path
import hashlib,json
spec=json.loads(Path('layout-replacements.json').read_text())
s=Path('../../paper/sections/eighth.tex').read_text()
assert hashlib.sha256(s.encode()).hexdigest()==spec['final_sha256']
for block in reversed(spec['replacements']):
    assert s.count(block['new'])==1
    s=s.replace(block['new'],block['old'],1)
assert hashlib.sha256(s.encode()).hexdigest()==spec['base_sha256']
print('exact inverse SHA verified')
PYVERIFY
```

Compile from the project `paper` directory:

```sh
latexmk -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

The actual local TeX Live run exited 0. The final XeLaTeX log has no formal warnings, missing characters, over/underfull boxes, or undefined references. Early clean-build passes requested reruns for references; latexmk resolved them. The first 52-page pass, its warning log, and its first render were preserved outside the repository before fixing the twelve layout blocks. `pdftotext` was unavailable (exit 127); no text extraction was used as evidence.

The final PDF was actually rasterized page by page with Poppler `pdftoppm` and a Python 15-worker pool: 53 successful page renderers, about 1.120 seconds wall time and 6.601 seconds summed renderer CPU. The main agent inspected all three page-contact images and 150-dpi pages 2, 3, 50, 51, 53. The worker count is not a claim of sustained full CPU utilization. `render-summary.json` binds the final PDF hash, source hashes, actual aggregate timing, and final-log checks. Reproducible build output may differ in PDF timestamp metadata.

Large page PNGs, intermediate TeX output, repeated old symbolic batches, and replaced paper copies remain outside Git. The mathematical computations needed by this paper are the already archived, actually executed source-cost, sharp settlement, feedback, and eight-layer packages; finite algebra does not replace the geometric proofs or all-input supremum argument.
