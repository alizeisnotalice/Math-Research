# R057 global eighth-order assembly: finite algebra checks

Frozen source: `scratch/r057-eighth-global-assembly/global_eighth_assembly.md`, SHA-256 `e3ed7614d8775cd57fe6792ab6235bfdb8ec0029ec94532564052ba67d78277c`. The two formal reviews are `scratch/r057-review-global/audit.md` (SHA-256 `76a80e950cfe11bbbaaa73137ddfe176558a5abc593572bf2606b46a3cdd9a4f`) and `scratch/r057-review-adversarial/audit.md` (SHA-256 `cd12441096e803a07674ecb7182148c43cc4ef55308a4cffb6e307abaa0737a9`). Public proof: `proofs/057-global-eighth-assembly.md`.

## Recorded Wolfram checks

`settlement.wl/log` records Mathematica 15.0.0 ARM, PID 21045. These are conditional algebra checks for fixed `n>=2`, with the normalization stated in the script. The log reports zero exact ratio-factorization residual; the generic coefficient residuals for `n>=3` are `{0,0,0}`; normal `n=2` coefficients are `{4/27,-176/243,15694/6561}` and inside `J_2=17638/6561`. It also checks the displayed degree bookkeeping, the fifth-order coefficient of `sigma`, a sixth-derivative product bound (`Reduce` violation set `False`), the normalized fresh-derivative identity, the eight-layer telescope, and the error-subtraction identity; the printed residuals are zero where stated. The log renders one `v^2` label as `vÂ²`; source code is unchanged and uses `v^2`.

The printed `n=2` J8 result is `{256, 1737988096/6561}` for the baseline and modified coefficient expression. `fraction_check.py/json` independently perform exact rational polynomial calculations for the `n=2` normal/inside J, baseline J8, and `q=0,1/2,1` sigma coefficients; its JSON reports `all_exact: true` and remainder bound `32/3`.

## Scope boundary

These scripts verify algebraic identities, finite coefficient extractions, selected rational examples, and the stated scalar Taylor bound. They do not prove the geometric hypotheses, actual-family cover or overlap properties, global CDF settlement, or the Tauberian theorem. Those conclusions depend on the mathematical source and its stated proof assumptions; the calculations are not a substitute for them.

## Reproduction

```sh
wolframscript -file settlement.wl
python3 fraction_check.py
```

The files below are the original code/output bytes; no calculation was rerun for this archive.

## Files and SHA-256

- `settlement.wl` — 2566 bytes — SHA-256 `505dd6488197f3dd1a7fb463ad18c5306c8136c20caf0ede3384c99699be0f1d`
- `settlement.log` — 968 bytes — SHA-256 `520489ad350b5d580b261b8b2d405da2b72af539163f263e4e93c3f2546be8f4`
- `fraction_check.py` — 1479 bytes — SHA-256 `0482fcf2efee7d08335ddf2623340294c566aa4333ace8a14350d6c616703841`
- `fraction_check.json` — 401 bytes — SHA-256 `793244cfb3221b56578b350e3bac2cfcc62947547acffd1b425ff408e7485263`

## Independent review calculations

`independent_check.py/json` preserve the non-author Sol exact Fraction run: 99 ordinary coefficients (n=2,…,12, j=0,…,8) have zero residual, with general bracket and sigma checks. `check_error_transport.wl/log` preserve the independent Astra Wolfram run (PID 23343): actual-mass error transport and nonlinear ratio residuals are zero; the delta coefficient >4 violation set is False. They support scalar transcription only, with the same geometric scope boundary above. Failed/unresolved preliminary runs are excluded; no old batch was rerun.

Reproduce: `python3 independent_check.py` and `wolframscript -file check_error_transport.wl`.

- `independent_check.py` — 2872 bytes — SHA-256 `6ea999890cc1958e35d31ec8579891ae53b8ca1b62b23a2df2ee8e3d35cf9a12`
- `independent_check.json` — 3549 bytes — SHA-256 `08eeff969891dfde0f7291bbb41490f03f33709cf1f3d9a995e9f42853644791`
- `check_error_transport.wl` — 1038 bytes — SHA-256 `9a58455eef1f1b9d66611c6323d1d9648bc96ba145bfa2b1a88907900db1488a`
- `check_error_transport.log` — 492 bytes — SHA-256 `68993ebacc68a7ed9d0327e2b79245fe1333d8e5c1aadbf7f0f514b23438300a`
