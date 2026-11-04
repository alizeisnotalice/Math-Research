# R049 eighth-stability: reverse-ratio identity check

This package contains only `reverse_ratio.wl/log` and this README. The Mathematica run checks one algebraic identity under the printed conditions `aa>0`, `vv>0`, and `av-aa(1+qq) != 0`; `Together[lhs-rhs]` returns zero. The actual run reports PID `5942`; the supplied execution record is exit status 0.

This check certifies only the displayed reverse-ratio identity. It does not verify the geometric or asymptotic premises, establish a numerical limit, or prove a Tauberian conclusion. The R049 stability review was frozen before this supplementary run; the run is post-review arithmetic evidence. The unrelated `local_variation` files are intentionally not included.

## Reproduction

```sh
wolframscript -file reverse_ratio.wl
```

Compare with `reverse_ratio.log`. No computation was rerun for this archive.

## Whitelist and SHA-256

- `reverse_ratio.wl` — 469 bytes — SHA-256 `5d5331eb9aa56423fc80229bc010d848a96f12df9569ab4560a687396f212c77`
- `reverse_ratio.log` — 315 bytes — SHA-256 `ed35cb2c46ff2b48dc48499ed6d74196dabb193e3c4ae2270a30b4638bb108f4`
