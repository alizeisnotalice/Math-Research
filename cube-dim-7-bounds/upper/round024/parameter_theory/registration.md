# R024 proposed finite scalar contract — not automatically released

Independent clone baseb135609efca22428ac8f6a9c9485d114a4a03846. Theory owner writes only upper/round024/parameter_theory, executes no mathematical numerical work/MATLAB/commit/push, and does not change any sealed round. Main must explicitly release the numerical owner after locking this proof/spec; this document alone does not start it.

Fixed functions n7,A=(1+2/a)^7,B=(1+1/a)^7(1+7loga),F=(sqrtA+sqrtB)^2,R=0.5*((a+1)/(a+2))^(5/2)*(a-7loga)/sqrt(1+7loga). Use512-bit Arb, exact rational inputs and positive square roots. Lock parameters.json and code hashes before executing. No MC, kernel, olda36tau7recomputation, added grid points, derivative search or library optimizer.

Prescribed sequence:

1. Evaluate exactly one new fixed expression F(36), with the closed form in proof(2), and retain its raw outward dyadic interval. The explicit tau_*(36)formula is proved analytically; it need not be another separately evaluated scalar row.
2. Evaluate R(35)-1 and R(37)-1. Require strict negative/positive intervals respectively. Failure stops the root-localization stage with honest results; do not change the starting bracket.
3. Exactly80midpoint sign evaluations. Keep endpoints and midpoint as exact Fractions, require a strictly separated sign at each step, update the corresponding endpoint and record input, signed interval and updated bracket. A ball containing0is unresolved, not equality; stop with no guess, retry or automatic precision expansion.
4. After80successful steps take the exact dyadic midpoint a_c of the final bracket and evaluate exactly one F(a_c). Record its exact numerator/denominator, the bracket width2^-79, the raw value interval, and its rational upper floor(U*10^12)+1 over10^12. Apply the same strict grid ceiling to F(36).

At most84scientific function evaluations are in scope: two Fvalues and82R-sign values. Cached log/sqrt intermediates at repeated exact a values may be reused, but actual library call counts must be reported. An algebraic cross-check F=A+B+2sqrt(AB)may be used internally without creating new parameter rows. Do not automatically add an interval evaluation of Fover the root bracket or claim its minimum-value enclosure from the candidate interval.

The accepted R023 result/manifest supplies its full winner interval and grid bound, read-only. Compare new upper endpoints to old lower endpoints only for a claimed strict numerical improvement; do not rerun the old row or compare to a display decimal. Nonseparation records an unresolved strict comparison and does not invalidate the theorem or authorize more parameters.

Acceptance: source consistency, all positive/finite square-root inputs, both initial signs, all80signs, exact bracket nesting/halving, exact candidate, strict outward grid ceiling and honest old-interval comparison. On failure preserve partial data and stop. A checker may verify the same fixed inputs, signs and intervals; no new science rows. Report actual small scalar CPU/wall/RSS without a saturation claim. Final language must say unique optimum of this parameter family, never optimum/lower bound of the actual C7.
