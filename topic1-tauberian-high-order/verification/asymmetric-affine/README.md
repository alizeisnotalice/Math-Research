# R047 asymmetric affine screen: computation archive

This package preserves the three successful Wolfram runs only. The invalid first `pair_symbolic` log is excluded.

## Scope of the screening result

The source screens only the R041 `B2` two-round bad-set feedback rule with its natural axis-step rectangular extension. It does not cover arbitrary affine holes, arbitrary witness choices, other finite input families, or the Tauberian constant `C_n`. Its conclusion is that this specified class gives no larger eighth-order lower coefficient than the R041 value; it is not an eighth-order upper bound for `C_n`.

## Pair moments and geometry-error separation

`pair_symbolic.wl/log` computes the Hessian/feedback gradient and completed-square identity for the scalar amplitude expansion, along with the exact two-arm `U0` moment and the second-feedback moment. `exact_u1.wl/log` performs the exact piecewise integration for the linearized `U1` pair model, then checks the `lambda^4`/`lambda^5` coefficients and exchange symmetry. These are exact algebraic integrals for the displayed two-dimensional model. They are separate from the source's geometric comparison between model slices and actual bad sets; the computations do not establish that comparison.

The uniform small-cap remainder in the source, `M1=Q4/f-Q5/f^2+O(e^6/f^3)` for `P/e,Q/e` in a fixed positive compact range, comes from the analytic split at `h=PQ/f`: the proof controls the `r<=h` cap and `h<=r<=Q` region with the full truncated kernel. In the small cap `r<=h`, the proof bounds the unweighted area difference by `O(e^3/f)` on a cap of width `O(e^2/f)`; after the full weight kernel `P+r=O(e)` and outer `1/f` factor, this yields `O(e^6/f^3)`. The large-cap estimate is handled separately. Uniformity follows from these analytic bounds for `P/e,Q/e` in the stated compact range, not from the finite symbolic examples.

## Axis-step check

`axis_steps.wl/log` checks the implicit-root/product expansion and gain residual only for `n=2,3`. The `n=3` defining-product residual has the recorded `-d1*d2*(d1+d2)*u^12` term, while the gain residual through `u^9` is zero. This finite-dimensional check is not the general-`n` root/remainder proof; that proof is supplied analytically in the source.

## Reproduction

```sh
wolframscript -file pair_symbolic.wl
wolframscript -file exact_u1.wl
wolframscript -file axis_steps.wl
```

Compare each run with its adjacent `.log`. The retained logs are the successful recorded outputs; no computation was rerun for this archive.

## Whitelist and SHA-256

- `pair_symbolic.wl` — 2035 bytes — SHA-256 `9368d37f8079cfe9385665be968910ae3dc4250fae47a7061a880ee15a7e9d2d`
- `pair_symbolic.log` — 832 bytes — SHA-256 `f53838b9039449194cd9392e79814222aab48e7eca4e9170d4a9998f4406dc0d`
- `exact_u1.wl` — 1113 bytes — SHA-256 `6c247e5296632e997aa1fdb47c04659edc033de8f2a859a10eb8a796eaec61fb`
- `exact_u1.log` — 600 bytes — SHA-256 `8812827054390db6a909de441917edd2a3065a8aa5ad4444ad3f914c8feae4d6`
- `axis_steps.wl` — 698 bytes — SHA-256 `c48474cb168813d60597d41b2be41f0b0d3cae9cbcecdf4848292cd806cc6b89`
- `axis_steps.log` — 486 bytes — SHA-256 `bd806208b6f7fa9fbe92dff927028c58f07415fd0a599d3eb78e62fb70aea1f7`

Frozen source SHA-256: `ce64039150e8fc980ecd71824c2247fd5837a9010178a6c1f43637305907b121`.
