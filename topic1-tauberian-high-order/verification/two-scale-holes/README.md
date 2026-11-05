# R054 two-scale hole: finite scalar computation archive

This archive preserves the four successful computation artifacts from `r054-computation` byte for byte. The frozen mathematical source is `scratch/r054-small-source-eighth/small_source_eighth.md`, SHA-256 `bc74ac7f20b63be1387bf1f465316fb61e00513b6567002246ac32095b4f00d1`.

## What the recorded calculations check

The Wolfram 15.0.0 ARM run recorded in `two_scales.log` (master PID 13682) checks the stated `n=2`, `u -> 0+` conditional-budget series algebra. The three displayed lambda-factor residuals are zero. Through order `u^10`, the ratio difference is `8192 u^10/27`; normalizing by the source scale `a v^8=1024u^10` gives the reported gain `8/27`. The coefficient is independent of symbolic safe parameter `c` through this order. The baseline coefficient is `15694/6561`, and the adjusted coefficient is `17638/6561`.

`fraction_series.py` and its actual `fraction_series.json` output independently use exact rational arithmetic for `c=0,1,7`. All three report coefficients 0 through 9 vanishing in the difference, normalized gain `8/27`, and the same old/new coefficients.

## Scope boundary

These are algebraic checks for the specified budget and finite series, not evidence that the geometric micro-hole host exists, that the `BH` protection or strict macro/micro witnesses hold, or that a finite truncation is valid. They do not establish an optimum or an eighth-order theorem. The source supplies the mathematical argument; these files record only the listed symbolic and rational checks.

## Reproduction

```sh
wolframscript -file two_scales.wl
python3 fraction_series.py
```

The Wolfram log's PID and outputs are from the recorded original run; the archive did not rerun either calculation. The scripts use relative/local inputs only and contain no private absolute paths.

## Files and SHA-256

- `two_scales.wl` — 1318 bytes — SHA-256 `7078cb0ec878a1fcbb2acff0e3b38015dc5c369e1ad651fd43f3c7c1fc847c4c`
- `two_scales.log` — 451 bytes — SHA-256 `85b5c051a1a2f9dc14bf8df9afcfa1830f63e1b05d89151ec76f0ae6c3d92cea`
- `fraction_series.py` — 1581 bytes — SHA-256 `08ffc8295dacb93e12c1c36bd07a7b290dc70eb028cf8dbfa41f84439ee3623b`
- `fraction_series.json` — 513 bytes — SHA-256 `624f1175e274fccca0f41e5899300be7cf35cdc68b0db28adacbe67860250c97`
