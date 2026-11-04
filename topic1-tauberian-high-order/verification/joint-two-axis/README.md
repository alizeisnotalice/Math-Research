# R045 joint two-axis price: scalar subwindow check

This three-file archive retains only the successful `normal_subwindow.wl/log` and this README. The computation checks a strict **subwindow only**:

- `nu=1/100`, `c=nu^2/1000=10^-7`;
- `1/2<=b<=1`, `1/6<=x,y<=1`, and `|z|<=1/12`.

Within that subwindow, the log reports empty violation sets (`False`) for the positive-part domain and for the endpoint coercivity checks at `b=1/2` and `b=1`; it also records the coefficient of `z^2` and the stated scalar bounds. These are the actual software checks, not a numerical test of geometry.

This does **not** certify the full-domain inequality (9) for every `b,x,y,z` or supply the full constant `c_nu`. In the source, the remaining `c_nu` is obtained by a separate analytic ratio/compactness argument (`analyticratio`); it is not inferred from this fixed subwindow. The frozen R045 local interface has now passed independent Sol and adversarial review. Its geometry is supplied by those analytic proofs, not this subwindow computation; the new eight-layer global settlement remains separate.

## Reproduction

```sh
wolframscript -file normal_subwindow.wl
```

Compare with `normal_subwindow.log`. No calculation was rerun for this archive. The two invalid first-run logs are excluded.

## Whitelist and SHA-256

- `normal_subwindow.wl` — 800 bytes — SHA-256 `78929c7e2e7a007ff44229e6a3c8a5dd39efe2bbeded78fa8d6870571fca1855`
- `normal_subwindow.log` — 418 bytes — SHA-256 `58af04f4930eb874ada1b28e71fd8ee516bba2d6eba1f0d4617e52ad7097e865`

Frozen source SHA-256: `8f0ad8b01bb443789b33f55188c7a04d000a3ceaee3fd61d0a83eb458c9bf8ea`.
