# R037 general-feedback computation archive

This package contains the two successful Wolfram runs used for the local price computation: `local_price.wl/log` and the supplementary `check_extra.wl/log`. It excludes the earlier `local_price-first-unresolved.log`, which is discussed below and is not evidence of a successful symbolic solution.

## Successful computation scope

`local_price.wl/log` records the scaled polynomial price checks, clipping identity, low-price and large-maximum violation reductions, and a scale-ratio reformulation. The `False` results are empty violation sets. The equality-classification line returns the condition `b == 1 && x == 1/3 && y == 1/3`; the supporting derivative line returns the expression `q^a Log[q]`. These are the actual outputs, not a blanket claim that Mathematica returned `True` for every assertion.

The successful scale-ratio check assumes `q=r_j/r_i>=1` and `beta=q^(n-2)>=1`, and reduces the remaining algebraic violation `q beta < 1`. The derivative expression supports the analytic monotonicity argument for powers on `q>=1`; it does not solve the earlier raw real-power system by itself.

`check_extra.wl/log` checks a clipped convexity violation set and records the symbolic predicate `q^a >= 1` under `q>=1,a>=2`. That second output remains a predicate rather than a software-returned `True`; it is included as a supporting expression, not represented as an automated proof result.

## Explicitly unresolved earlier attempt

The excluded `local_price-first-unresolved.log` contains `Reduce::nsmet` and an unevaluated two-scale real-power `Reduce[...]` expression. The two-scale comparison is handled by ratio variables and a direct power-monotonicity argument in the reviewed proof; the unresolved first output is not counted as verification. It is deliberately not copied into this minimal success package.

## Scope boundary

The two R037 audits report PROOF-VALID for the local price and the new fixed reconstruction (N=7) interface at the frozen statement. The package does not establish an all-input global seventh-order Tauberian upper bound, the complete global ledger, or a global next coefficient. The recorded calculations are scalar algebraic supports; they do not replace the geometric proof.

## Reproduction

From this directory:

```sh
wolframscript -file local_price.wl
wolframscript -file check_extra.wl
```

Compare each run with its `.log`. No calculation was rerun to prepare this archive.

## Whitelist and SHA-256

- `local_price.wl` — 1265 bytes — SHA-256 `826310f6e9c2aa9fe1fae812bac39d347ecc78bf5557fe7cdb89e445f673a375`
- `local_price.log` — 589 bytes — SHA-256 `c8e716862b988981cc6f022bc6503d7eae29ab764ba372b112ccd382eb2aa4f0`
- `check_extra.wl` — 272 bytes — SHA-256 `9bce563dc3d578dfb5417fd4a42c238bed0aadd905d52b1c85dc04cd2548be12`
- `check_extra.log` — 59 bytes — SHA-256 `852e06a6a5b82d4be00d9bae030df4e978b73ee94d8420f334bf21bbddf5b859`
