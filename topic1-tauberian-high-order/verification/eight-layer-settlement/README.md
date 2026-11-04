# R046 conditional eight-layer settlement: algebra archive

`eight_layers.wl/log` checks a conditional eight-layer telescoping ledger for fixed `n>=2`, `tau=b v/c`, `c>0`, with layer inputs `E7=Q T+e7 v^7` and `E6=(k+e6)v^6`, where `e6,e7 -> 0`. It records zero telescoping residual, the endpoint order `A0 b^8 m v^9/(aa c^7)`, zero unnormalized residual through `v^7`, and the remaining normalized `v^8` coefficient `C8 + b (b e6+c e7)/c^2`. The sixth-layer difference is also recorded exactly.

This is algebra under the stated layer hypotheses. It does not prove that a new actual finite family supplies those layer bounds, nor prove a corresponding upper bound. The new-family supply and the upper-bound argument remain subject to proof and independent review. The R38 denominator certificate is reused and is intentionally not duplicated here.

## Reproduction

```sh
wolframscript -file eight_layers.wl
```

Compare output with `eight_layers.log`. No computation was rerun for this archive.

## Whitelist and SHA-256

- `eight_layers.wl` — 974 bytes — SHA-256 `21ce70e6224edbe3d718d5acb0b82f038730c7a69e531e30b31f3d3841ca89e4`
- `eight_layers.log` — 466 bytes — SHA-256 `aed4aec82c6214fcaebe8254fae66b3b487fb80e0c0de950049588eea2c65893`

The R046 conditional proof has now passed two independent full reviews. Its exact H45 condition is supplied by the separately reviewed R045 version; the combined global conclusion is recorded in rounds/046.md. This software output itself still verifies only the conditional algebra.
