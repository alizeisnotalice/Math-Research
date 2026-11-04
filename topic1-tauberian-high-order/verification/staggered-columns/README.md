# R032 staggered-column calculation archive

This package is a calculation-side filter for the proposed three-dimensional staggered-column construction. It retains two complementary checks: exact symbolic Jacobian/two-box algebra and a nonseparable periodic rational Simpson sample. The second is retained because it is the only concrete exact sample in this package that exhibits nonzero local cell contributions with exact period-average cancellation; it helps screen the tempting periodic-column route. It is not a general no-go theorem.

## Exact symbolic checks

`jacobian.wl/log` checks the determinant identity `det = 1 - pair + cycles`, and the symbolic three-cycle has the form
\[
u^7(ade+bcj+u^2(adf+bcg)),
\]
so the order-three term starts at order `u^7`. It also records the mixed two-box intersection `8 b^2 s u^7`. In the explicit affine example the shift is `(1-2u^3, 1-2u^3, 1-2u)`; the determinant is positive for `0<u<1/8`, all overlap widths are positive, and the intersection product is `8u^7`. This is a two-box intersection with three active coordinates, not a strict-halo lower bound.

## Periodic sample and interpretation

`periodic_fraction.py/json` defines a nonseparable `3x3x2` periodic trilinear perturbation and computes exact tensor-Simpson averages of its cell Jacobians at four rational values of `u`. All diagonal derivatives are exactly one. In each sample the mean determinant is one, the mean pair term is zero, and the signed cycle average is zero, while the maximum absolute cell-cycle coefficient divided by `u^7` is positive and the `xz` pair coefficient ranges from `-2/9` to `2/9`. This illustrates local nonzero values canceled by period averaging. It is a diagnostic for this specified periodic family only; it does not rule out arbitrary nonperiodic phase arrays or establish that every n=3 column phase has no periodic variation.

The scripts check algebra and exact finite parameter-cell averages. They do not prove global injectivity for the periodic sample, construct a halo, or certify all-radius halo behavior. The affine overlap example likewise does not supply source-cost charging or the missing opposite-face geometry.

The independent Astra audit `r032-review/staggered_column_audit.md` reports PROOF-VALID for its stated frozen proof. The separate Sol final audit is now PROOF-VALID for the same restricted construction class. This calculation archive alone does not establish the geometric proof. No proof of an all-radius halo statement is encoded here.

## Reproduction

From this directory:

```sh
wolframscript -file jacobian.wl
python3 periodic_fraction.py
```

The Wolfram output is recorded in `jacobian.log`; the Python script prints the exact rational sample object as JSON. The retained `periodic_fraction.json` is the recorded output; no computation was rerun to prepare this archive.

## Whitelist and SHA-256

- `jacobian.wl` — 965 bytes — SHA-256 `8e7262effd37f9142779ede203ddb52ff507b63e60ed952dccea7a531c484eaf`
- `jacobian.log` — 444 bytes — SHA-256 `8ad9f858239424717dd751c9aed03336661d9a775502aeebe7874403c0d52f08`
- `periodic_fraction.py` — 2230 bytes — SHA-256 `2bf35620be19ed457067c19d4a354c93b45207322a5007a04de3237ac8475fea`
- `periodic_fraction.json` — 1684 bytes — SHA-256 `b8656ebb8bde913d4a82e75b39f42225356a1c19009f9f5987a68d067bc4c1cb`
