# R056 sharp inside-source scalar checks

Frozen source: `scratch/r056-sharp-small-source/sharp_inside_cost.md`, SHA-256 `18fa60b96a6b2a2be7092d06dc59b4741a049e4214941323a33f846dad12ffc0`. The formal reviews are `scratch/r056-review/audit.md` (SHA-256 `a87f991f05c0e79ac661cf89fd388be0f3ffaf51a9c98d738b92ae90e0cc1ab4`) and `scratch/r056-review-adversarial/audit.md` (SHA-256 `2e5f13ea480534ea1d3ef440fc8fa18f267059de0da455a858b3f2e20e8785bd`). Public proof path: `proofs/056-sharp-inside-source.md`.

## Recorded Wolfram checks

The preserved Wolfram 15.0.0 ARM log records these scalar results:

- For two oriented one-dimensional intervals with positive radii and center separation greater than `Max[r1,r2]`, the positive-overlap violation set is `False`.
- The two-end nominal coefficient is `8/(27(1+v)^3)` and its residual against that expression is zero.
- The fresh-derivative identity and final exponent identity have residual zero.
- The inverse-square-root and fresh-monotonicity violation sets are `False` on their stated domains.

Reproduce with `wolframscript -file sharp_local.wl`. The script and output were copied byte for byte; no computation was rerun.

## Scope boundary

These are scalar and one-dimensional interval checks. They do not certify a full-point measurable cover, nonintersection of `Pi` boxes inside or across angular sectors, or the same-family quantifiers. The computations do not establish the geometric conclusion on their own.

## Files and SHA-256

- `sharp_local.wl` — 919 bytes — SHA-256 `e7209ffe22d3b9f70b51c189f60cc8c38ceb0c88a47fba0dc0db7578c8e388e7`
- `sharp_local.log` — 353 bytes — SHA-256 `9440677862d9eba7594760f74ec79b350729ce4f7c7fa55ea2d9a3d5d12c0aed`
