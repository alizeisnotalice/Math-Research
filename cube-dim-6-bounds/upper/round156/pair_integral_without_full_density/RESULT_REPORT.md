# C156 actual result: whole true joint-law upper budget

This is an append-only report. The frozen producer, proof, inputs, release, saved fractions and root records remain unchanged. No bounds, roots, optimizer or saved numerical checker were rerun to prepare this report. Only JSON reading and SHA-256 provenance checks were performed here.

## Actual execution and provenance

Terra ran the single authorized batch in `Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round156/pair_integral_without_full_density`. PID 89280 exited 0, without retry. The operator receipt records wall time 56.724943374982104 seconds, child user CPU 43.502114 seconds and child system CPU 0.419483 seconds. The exact result SHA-256 is `f9a0689da14fb34d48c3ffb65b6abaef06193906ff71cc3cab77420074bee21f`.

The release SHA-256 is `93ce0baf5531e12335c5887c4d778edf1c39175256536bfc4d7704eb4145052f`; it binds main approval and the independent concrete preflight review, whose SHA-256 is `58bfc6f585076240a0915c957389bf29eb312b1c44081ce0d13518096a5491a5`. Those were preflight approvals; the actual saved result still requires its separate independent saved-data audit. The appended manifest pins the complete actual directory and its external source files. No claim of a completed saved-data audit is made here.

## Full-support mechanism and finite family

The round130 identity is used with its complete positive weights and factor `kappa=46080/pi^3`: sum the 21 cone weights alpha and all 496 unordered pairs in each cone. For each pair, its two coordinates are the positive pushforward of six independent Exp(1) variables. A nonsingular two-column seed has density supremum at most the reciprocal determinant; convolution with the remaining probability law preserves this bound. Two disjoint nonsingular two-column seeds admit the Cauchy convolution bound `1/(4 sqrt(D1 D2))`. Zero columns, repeated directions and singular seed choices contribute through the unchanged positive remaining-law factors, or the nonsingular seed is skipped. No signed density term is used.

The kernel is positive, decreasing and logconcave on `[0,infinity)` as proved in the frozen PROOF. With `p=k(0)`, its normalized form satisfies `k(a+b)/p <= [k(a)/p][k(b)/p]`. The two fixed anchors `0` and `sqrt(h)`, at the exact `h=1/128`, give exponential majorants. The finite family combines empty, two-column and four-column seeds with the remaining positive one-dimensional factors, taking a minimum of valid caps separately for every actual pair. It integrates the entire Exp probability law, without truncation, boundary loss, probability renormalization or density generation.

The actual batch covered 10416 raw pairs, evaluated 229314 caps, formed 156240 exact column determinants and 468720 four-seed pairing products. It made 147 new Cauchy square-root cache evaluations; 73524 choices were pruned by exact range comparisons. There was one new closed I0 evaluation and one positive kernel anchor evaluation, with exp=3, erfc=2 and erf=0. Old 343456 functionals, old six-dimensional grids, native kernels, Monte Carlo and optimizers were not run.

## Result and exact comparison scope

For each family, the saved exact rational interval encloses its positive *budget expression*. Its upper endpoint is a candidate valid upper bound for the true full beta2=mu2, pending the independent saved-data audit. Its lower endpoint is only a lower bound for that budget expression. It is not a true mu2 lower bound and does not exclude other upper-bound mechanisms or other parameters.

The complete exact fraction to be consumed downstream is the value of `result.json["budgets"]["finite_best"]["valid_true_beta2_upper"]`, with the result SHA-256 above. The entire fraction is retained verbatim in that immutable source; downstream code must read it as a Fraction. It must not reconstruct it from this report's decimal display.

Terra's saved `operator_display_summary.json` reports the display of finite_best as `1.6654804884606890991806385E-7`; the flat budget display is `0.0000016880413139347619776344893`. These strings are display only. The exact comparisons saved in result.json are:

| Budget | BU < 7/10^8 | BU < 8/10^8 | BU < saved146 legal beta2 cap |
|---|---|---|---|
| flat | false | false | false |
| finite_best | false | false | true |

The budget-expression lower endpoint lies above both 7/10^8 and 8/10^8 for both registered families. Thus these two fixed families cannot reach those thresholds at the frozen parameters. This is not a lower bound on actual beta2, M6 or C6, and it is not a full-family exclusion.

The finite_best upper endpoint improves the old legal beta2 cap if independently confirmed. This permits a new normalized common-parameter model, rather than continuation of round146's optimistic lower-endpoint coefficients. C160 must first round the complete upper fraction upward to a dyadic cap, verify the exact changed feature identity, and receive separate phase gates. No new M upper, global g proof or C6 bound is registered by C156.
