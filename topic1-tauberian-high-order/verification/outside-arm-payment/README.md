# R034 outside-arm payment: arithmetic archive

This package retains the symbolic power/exponent checks and the exact rational boundary examples cited in the independent review. Its scope is scalar algebra and finite arithmetic for the stated estimates. It does not establish the geometric tube/rectangle exclusion, source charging, or a general-dimensional seventh-order asymptotic.

## Symbolic checks

`power_bounds.wl/log` works on the stated domain `0<t<1`, `v=t^8`, `0<p1<=p2<=p3`, `r^3=p1 p2 p3`. The `Reduce` outputs labelled `failure` or `violation` are counterexample sets; `False` means the corresponding set is empty. The identities and inequalities checked include:

- ordering implies `p2 p3 >= r^2`;
- `r^3>t^57` forces `r>t^19`;
- for `z>t^38`, the truncated multiplier `Min[1,z/t^24]/t^16` is at least `t^-2`;
- the displayed `H` branch identity has residual zero;
- the stated `A <= max(1,C^2) v^2 H` counterexample set is empty under `0<v<1`, `C>0`, `0<z<=C^2`;
- for `q>=4`, `q+6q/(2q-1) >= 52/7`, with equality at `q=4` and strict inequality for `q>4`.

The log also gives the scalar example `v=t^3`, `p_i=t^7`: `A=t^21`, `H=t^12`, `R=1/t`, while `1/v=1/t^3`. Thus this payment factor cannot be uniformly replaced by `1/v` in that example. This is a scalar example, not a realized geometric counterexample.

## Exact Fraction boundary cases

`boundaries.py/json` has eleven exact rational cases. The first nine cover `A=v^(57/8)` equality, a strict-above-threshold sample, and truncation-saturated samples at `t=1/2,1/4,1/8`. The last two check the separate `v^7` scalar example at `t=1/2,1/4`. The review explicitly used these to cross-check threshold-side and saturation cases, so they are retained. They are finite checks and do not replace the preceding symbolic inequalities.

The source `power_bounds.log` has a display-encoding artifact in one superscripted label (`CÂ²`, `vÂ²`); the source WL expression is `cc^2`, `v^2`, and the recorded `Reduce` result for that counterexample set is `False`. This is not a claim that every printed status is `True`: empty counterexample sets are intentionally reported as `False`, while the q-factor line is a factorization.

## Reproduction

From this directory:

```sh
wolframscript -file power_bounds.wl
python3 boundaries.py
```

The Wolfram output is recorded in `power_bounds.log`; `boundaries.py` prints the recorded JSON object to standard output. No computations were rerun for this archive.

## Whitelist and SHA-256

- `power_bounds.wl` — 930 bytes — SHA-256 `7434cee958c30e4535b45ca2431f06eb5af0c815101c550ea1b7a08cf4593cf6`
- `power_bounds.log` — 444 bytes — SHA-256 `256f9ce0e559e1374a1a59bff9327bc7c4e0b022e55b35c64e4342ec2ccc466a`
- `boundaries.py` — 627 bytes — SHA-256 `00419a83fa0f9686b5e8fb228742143eb09c1bbbc967909e086f3733ff44e8b8`
- `boundaries.json` — 2045 bytes — SHA-256 `17f278a6f1c0924d9b8da06eb3d9b9d7559e453041387402fd41892e43f8d9b8`
