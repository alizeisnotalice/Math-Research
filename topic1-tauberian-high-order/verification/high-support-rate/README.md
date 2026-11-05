# R061 high-support tail exponent algebra

Frozen source: `scratch/r061-high-support-rate/high_support_rate.md`, SHA-256 `a037675512160f4717f49f37f466f0579fb386e8cdc67078fa1a4a34b8b9c663`. Independent review: `scratch/r061-review/audit.md`, SHA-256 `f2a74b5ba9bb0d07c508ef42d5ae568a451016f74a195b3917c9be99aea6cb02`.

The preserved Mathematica 15.0 ARM run (PID 32510) checks exponent algebra only: the balance equation gives `kappa=3/5`, both balanced exponents are `43/5`, the `q>=4` violation sets for lower bounds `43/5` and `158/15` are `False`, the `q=4` value is `158/15`, and the Young exponent margin is `7/5`. The script also prints the already-supplied `q=3` scalar Young certificate.

This calculation is not a geometric proof. It does not establish stability, rectangle widths, target exclusion, actual-family overlap, or the complete `C_n` remainder. The general argument and its hypotheses are in the frozen source.

Reproduce with:

```sh
wolframscript -file tail_exponents.wl
```

The code and log are byte-for-byte copies of the successful run; no calculation was rerun.

## Files and SHA-256
- `tail_exponents.wl` — 632 bytes — SHA-256 `ae000b254e74721a81d11b915b1ee10a30fdcc6321e0eb4d1a011c8d0538e2d9`
- `tail_exponents.log` — 441 bytes — SHA-256 `a7a17575d651c037041734694dfd00ae505f4287ef2d11c9303000bc76bfa4d2`
