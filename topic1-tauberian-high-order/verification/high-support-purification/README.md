# R050 high-support purification: exponent check

`exponent_checks.wl/log` verifies the three scalar exponent expressions used in the analytic proof:

- `beta_q + 2 - 19/2 = (q-3)(4q-5)/(2(q-1))`;
- `q + q/(q-1) + 5 - 19/2 = (q-3)(2q-3)/(2(q-1))`;
- `q + 6 - 9 = q-3`.

For real `q>=3`, the recorded violation sets for each expression being negative are all `False` (empty), and all three expressions equal zero at `q=3`. This is a scalar algebra check supporting the real-q analytic proof; it is not geometric evidence and does not replace the proof of the high-support estimate.

The independent R050 audit reports PROOF-VALID for its frozen conditional statement and scope. This package contains only the exponent computation; it does not assert any broader input class or add a new proof step.

## Reproduction

```sh
wolframscript -file exponent_checks.wl
```

Compare with `exponent_checks.log`. No computation was rerun for this archive.

## Whitelist and SHA-256

- `exponent_checks.wl` — 247 bytes — SHA-256 `98db29a22f3e7ff941c0766aa8d00e2128d06bc2c2d3c690c67ff7b75c281717`
- `exponent_checks.log` — 214 bytes — SHA-256 `e2b508e31d1fe8234a14679fa2538574e88a1fb3a14fa1a2e5f38d68e83da07e`
