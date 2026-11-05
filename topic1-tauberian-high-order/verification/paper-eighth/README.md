# Exact eighth-order paper: R058 assembly and actual PDF verification

The final editable paper is `paper/main.tex`, `sections/global.tex`, `sections/sharp.tex`, `sections/eighth.tex`, and `sources.bib`. The PDF has 67 A4 pages. It proves the all-input, all-positive-radius exact eighth-order asymptotic for every fixed n>=2, with Jhat_n = (n-1)n^3(2609+1971n+81n^2)/26244 + (16/27)1_(n=2). In particular Jhat_2=17638/6561. The ordinary coefficients through j=8 and the all-m uniform transfer are included. There is no dimension-uniform claim.

Mathematical reviews are `reviews/058-main.md`, `058-module.md`, and `058-module-adversarial.md`. The separate `058-layout.md` binds the final source to those mathematical freezes by exact reversible blocks. Global/sharp/bib remain unchanged. The old coarse eighth-order proof remains as the noncircular bootstrap for the competitive-input stability theorem; the weaker old Theta main conclusion has been removed.

`layout-replacements.json` records eight exact old/new blocks for two files: formula wrapping, one CJK text-font wrapper, two terminology edits, and explicitly stating measurability in the local lemma. The latter was already used in all applications and was approved as prose clarification by both module reviewers; no new geometric or radius hypothesis is introduced. Independent forward and inverse byte checks passed. The reviewed mathematical hashes are main `413a02da01cccd08f33e25709b5ac579f19d092dbfd4aa965804e675c22f68ed` and module `56e4c2c04839f508aa6b546d41f69c168571181a243540b1cfe62a22e5eb5ae9`.

From this directory, verify the exact reconstruction without writing any source:

```sh
python3 - <<'PYVERIFY'
from pathlib import Path
import hashlib,json
spec=json.loads(Path('layout-replacements.json').read_text())
for name,info in spec['files'].items():
    s=(Path('../../paper')/name).read_bytes().decode()
    assert hashlib.sha256(s.encode()).hexdigest()==info['final_sha256']
    for block in reversed(spec['replacements']):
        if block['file']==name:
            assert s.count(block['new'])==1
            s=s.replace(block['new'],block['old'],1)
    assert hashlib.sha256(s.encode()).hexdigest()==info['base_sha256']
print('both exact inverse SHAs verified')
PYVERIFY
```

Compile from `paper`:

```sh
latexmk -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

The actual local XeLaTeX build exited 0. The final log has zero formal warnings, over/underfull boxes, missing characters, undefined references, or fatal errors. There are 211 unique labels and 240 resolved references. First and second 67-page passes, including their original source/PDF/logs, remain outside Git; their four/two and one/zero overfull/missing issues were fixed before final acceptance. The replacement old 53-page project and PDF were preserved separately with exact hashes.

Poppler `pdftoppm` rendered all 67 final pages in a Python pool with 15 actual workers: 1.401 seconds wall time, 8.471 seconds summed renderer CPU. A second pool with 15 actual workers rendered 15 distinct meaningful pages at 150 dpi, selected from the final aux labels plus the first and final pages: 0.451 seconds wall, 3.590 seconds summed CPU. These worker counts do not assert sustained full CPU use. The main agent viewed all three contact sheets and 150-dpi pages 1,2,50,55,59,63,66,67, plus full 85-dpi pages 3,54,64,65. No clipping or missing glyphs was observed. `render-summary.json` binds all source/PDF hashes, actual aggregate timing and inspection scope. PDF timestamps may differ on rebuild.

The exact scalar kernels needed for the new theorem were already actually run and archived in the normal-eighth-packet, whole-height-complement, sharp-inside-source, two-scale-lower and global-eighth-assembly packages. R058 reviewers additionally retained their original narrow Fraction scripts/stdout outside Git; the reports state their exact scope and hashes. These computations support algebra and transcription; the geometric and all-input proof is in the paper. No old large symbolic batch was rerun for utilization.

Large PNGs, temporary build output, original reviewer scripts with private workspace bindings, and replaced paper copies stay outside Git. The user's instruction is to pause after this delivery and confirmed GitHub synchronization; no next research round was started.
