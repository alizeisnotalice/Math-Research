# Actual PDF visual check

Final PDF: `tauberian.pdf`, SHA256 `5f84aaa0a3d18fad3e9a86eea5a275bf9c5222a13d8171f2d10a7e47d6ff98a1`.
73 pages, 703251 bytes.

The first complete compiled version (6126f3...) was actually rendered at 100 dpi on
all 73 pages. All five numbered 4x4 contact sheets were viewed. Pages
1, 3, 5, 6, 60, 64, 65, 66, 68, 69, 70, 71, 72 were independently rendered at
150 dpi and each image was viewed. No clipping, text collisions, missing glyphs,
broken formula placement, unexpected blank page or abnormal page spacing was observed.

The parent-requested final removal and two pointer corrections were then actually
compiled again. The final PDF was rendered on all 73 pages at 100 dpi into
`render-delivery/`. Exact image/pixel comparison against the previously viewed
version found only pages 36, 47, 72, 73 changed; the other 69 pages are byte-identical
PNG files. Those four changed pages were rendered again at 150 dpi into
`render-delivery-keys/` and each was viewed. All four are clear and fit the page.

The final log contains zero warnings, overfull/underfull boxes, undefined references,
undefined citations or missing-character notices. Raster ink bounds are clear of
all page edges. `pdf-delivery-check.json` stores page counts, text counts and ink bounds.

This is a layout check, separate from the nonauthor mathematical proof review.
No mathematical source was changed after the three documented layout edits and
three final trim/pointer edits. The full inverse chain is executable as
`python3 evidence/verify_inverse.py`, without writing or reverting any source file.
