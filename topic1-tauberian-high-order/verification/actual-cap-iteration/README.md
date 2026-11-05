# R052 actual-cap operator: scalar computation archive

This archive contains only the successful `operator.wl/log` and this README. The frozen source has SHA-256 `d655c1e160c88914774248c81688c6bdde3f27614387c4eb7f0b9c1a0a0a5296` and the completed independent Sol review is archived in `reviews/052-actual-cap-iteration.md`. This package records computation scope only.

## Checked scalar facts

For `h>=0`, `clip(x)=Min[h,Max[0,x]]`. The two `Reduce` queries test violations of monotonicity and the 1-Lipschitz bound; both outputs are `False`, meaning the scalar violation sets are empty.

The WL then checks the finite symbolic recursion with `lambda=cc v`, `mass=aa area v^5`, and `rho=rr area v^7+ee area v^9`: it prints `(v/area) mandatory` and `(v^2/area) e1` with their lowest nonzero order `v^9`, and separately prints the unscaled term `lambda^2 mass=aa area cc^2 v^7`. These are algebraic expansions only; they do not establish the geometric hypotheses feeding the recursion.

Two log labels render a superscript `²` as `Â²` (`vÂ²` and `lambdaÂ²`). The source WL uses `v^2` and `lambda^2`, and the printed mathematical expressions/leading orders are otherwise recorded as shown. The copied WL/log are unchanged; this archive does not repair or rerun the original output.

## Scope boundary

The calculation checks only the scalar clipping facts and this finite recurrence. It does not prove actual-cap Tonelli estimates, source-chart or coordinate-splitting properties, the R051 same-witness necessity interface, or physical allocation by the strong CDF. Those are proof-level inputs checked in the separate completed review.

## Reproduction

```sh
wolframscript -file operator.wl
```

Compare with `operator.log`. No calculation was rerun for this archive.

## Whitelist and SHA-256

- `operator.wl` — 531 bytes — SHA-256 `467b24066f4037fd936a88ff5a935c3701fb0fd964beace6dd0efb0fd9fe5207`
- `operator.log` — 308 bytes — SHA-256 `24d5fcaa296b1874d73cf1f0c34743f2f72b390ddd7eb6e11ebc4313da99d7dd`
