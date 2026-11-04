# R036 adaptive soft-slot computation archive

`soft_price.wl/log` contains the scalar Wolfram checks used with the adaptive geometric argument. On the stated domain it checks the normalized soft-weight bound by reducing the violation sets for `f > omega`, `f > 1`, `f > 27 Min[x,y]/4`, and the quarter inequality; the recorded `False` outputs mean those violation sets are empty. The completing-square factorization is also recorded.

The independently reviewed result is limited to the adaptive geometric capacity statement, scalar calibration, bounded gross branch, and substitution under the explicit local-price condition `(P)`. The actual all-dimensional price `(P)` is not proved in R036 and is left to R037 review. This package therefore does not establish an unconditional full-dimensional seventh-order estimate or the Tauberian next coefficient.

The log also contains a layer-integral identity, returned as `True`. It is an extra algebra check only: the reviewed proof does not use layer-cake or integrate fixed-threshold estimates to obtain adaptive geometry. The computation does not prove geometric helper allocation, slot uniqueness, or full capacity; those are supplied by the proof and its independent audit.

## Reproduction

From this directory:

```sh
wolframscript -file soft_price.wl
```

Compare the output with `soft_price.log`. No computation was rerun for this archive.

## Whitelist and SHA-256

- `soft_price.wl` — 722 bytes — SHA-256 `fa613b8d2cbcebb09842ba28b7b85525f1c284e8642ceda3e9b4d580ba600a64`
- `soft_price.log` — 358 bytes — SHA-256 `bb09cc8d821720a06b50c823b9ffac6e0540347a5a5b897f57067f41570ea25a`
