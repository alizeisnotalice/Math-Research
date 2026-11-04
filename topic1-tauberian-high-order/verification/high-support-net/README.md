# R044 high-support net price: minimal computation archive

`young_branches.wl/log/json` records conditional scalar endpoint checks for `q=3,...,17`, run on 15 actual subkernels. Kernel IDs are `[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]`; actual PIDs are `[96776, 96777, 96778, 96780, 96781, 96782, 96783, 96784, 96785, 96786, 96787, 96788, 96789, 96790, 96791]`. The log's `False` outputs are empty scalar violation sets, while both branch checks are `True`. The exact entries, critical points and constants are in the JSON. The general integer-q proof is the source's two-case argument, not an extrapolation from this finite table.

`fraction_examples.py/json` provides six exact Fraction examples at `lambda=1`: for `q=3,4`, respectively, a positive-net case, the `P=v^(2q-1)` boundary case, and a negative-net case. These are scalar examples only. `precheck.wl/log` is excluded as duplicate material.

## What the estimates do and do not say

The high-support long-tail conclusion for `q>=4` is a **net-price** estimate after subtracting the paid cost and allocating one fixed share of `Z(K)`; it is not a gross-intersection estimate. The short-window `q=3` bound is the coarser gross estimate `C_K v^7 S(B)`. It is neither `v^(15/2) S(B)` nor paid from `Z(K)`; the source identifies the separate R043 missing-slot budget as the proposed absorber. The package does not verify that geometric budget or the non-overcounting needed to apply it.

The scalar scripts do not prove the existence of actual rectangles, helper geometry, bad-root ownership, or a global eighth-order estimate. The frozen R044 local net-price lemma has now passed independent review, and R043 actual slot geometry has been independently reviewed separately. The common positive budget needed to use the two together still requires R045; this archive supplies none of that geometry.

## Reproduction

From this directory:

```sh
wolframscript -file young_branches.wl
python3 fraction_examples.py
```

The Wolfram script launches 15 kernels and writes `young_branches.json` beside itself; use a disposable copy if preserving the archived output. The Fraction script prints JSON to standard output. No calculations were rerun for this archive.

## Whitelist and SHA-256

- `young_branches.wl` — 1229 bytes — SHA-256 `a63d1b8dbb1e77e31e149d5b6c76fb025d5ff6f9c5cc064f63953d4d6f46b32b`
- `young_branches.log` — 424 bytes — SHA-256 `199bf3d07f21ad2f616e63210b2f46ad8c80a156293447599bcbb0da3a645dd5`
- `young_branches.json` — 3026 bytes — SHA-256 `f12b37aaddd65917e8231e4657ba577e82e0e135869e899894675b01bdc5ef55`
- `fraction_examples.py` — 436 bytes — SHA-256 `472a13e39575c743d00bdd38cc7ea08e2c354abd2a4bdd7f12180e0f8f059956`
- `fraction_examples.json` — 1388 bytes — SHA-256 `c293c75314d63823d0dc8f898eda7357ce7526affa0ce1c7c5d0e21d786f9186`

Frozen source SHA-256: `ef3a6df6ee4581e990fd55fbe76c2b03f43b2e5bdbd0cb27f7eba18231cf0cc4`.
