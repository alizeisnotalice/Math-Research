# R055 whole-height complement: scalar clip increment check

The frozen source is `scratch/r055-eighth-complement/whole_height_complement.md`, SHA-256 `0684d01f3665a5a064b4786b6eb1d78ea5780ccdae64f96df3ad01647224a1a1`.

`clip_increment.wl` asks Wolfram Reduce over real `H,u,w >= 0` for the violation sets of
`0 <= min(H,u+w)-min(H,u) <= w`, and of equality to `w` under `H >= u+w`. In the preserved successful Wolfram 15.0.0 ARM output, both violation sets are `False` (empty).

This archive checks only those scalar identities. It does not certify the `O_T(v^3)` estimate for free base variables in a coordinate split, actual same-source CDF bookkeeping, the order of limits in `T` and `m`, or any Tauberian conclusion.

Reproduce with:

```sh
wolframscript -file clip_increment.wl
```

The source script and recorded output are copied unchanged; no computation was rerun.

## Files and SHA-256
- `clip_increment.wl` — 289 bytes — SHA-256 `274bba84c1a184907a6e118ad9193ab2b5c05e725f38302767d775da572462ea`
- `clip_increment.log` — 145 bytes — SHA-256 `6b95cb0c57a112a139e44d6632ac64854b24bd66271506e611ab693237a373d4`
