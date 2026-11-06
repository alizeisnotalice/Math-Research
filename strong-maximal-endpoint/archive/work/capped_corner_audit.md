# Capped-corner continuous witness: independent audit

Date: 2026-11-06. This audits the current exact-rational implementation in
`cert_capped_corner.py`. It does not assert an optimal value of `C_2`.

## Verdict

No blocking geometry, branch-selection, or logarithm-bound error was found.
The program supplies a lower bound for the **actual continuous** strong
maximal halo of its stated step function. Its comparison is only with the
family of constant functions on one axis-parallel rectangle of amplitude at
most 16.

## Geometry

For the right strip, fix `y` in `(j,j+1)` and use `[a,x] x I_y`, where
`x>4` and the vertical interval contains `y`. Since its mass is independent
of x, the average condition is `x<a+mass(I_y)/(tau*|I_y|)`. Expanding the
three choices `I_y=[l,r]`, `[l,y]`, and `[y,r]` gives exactly the constant
and Möbius functions appended in `strip`, after subtracting 4 from the
right boundary.

For the northeast corner, `[a,x] x [b,y]` has fixed mass `M_ab`; the
condition is `y<b+M_ab/(tau*(x-a))`, exactly the function in `corner`.
The interior, four strips, and four corners are disjoint apart from null
boundaries. Reflections and transposition supply their proper orientations.

The current `lower_integral` selects a maximizing branch at each midpoint
using only `Fraction` arithmetic. It then explicitly checks that the
chosen Möbius denominator does not vanish between its segment endpoints.
Because every branch is monotone on such a segment, its endpoint minimum is
a lower height for that branch and hence for the pointwise maximum. This is
valid even without proving that the selected branch is globally maximizing
on the whole segment. No floating-point value participates in the witness
area computation.

## Logarithmic bounds and comparator

The atanh sum has terms k=0,...,63 at u=1/3. The positive omitted tail is
bounded by `2*u^129/(129*(1-u^2))`, so `244+968*log(2)_upper` is an Orlicz
cost upper bound.

For a constant amplitude A on any rectangle, direct product-halo
integration at threshold 1 gives
`(1/A+4 log(A))/(1+log(A))`. Writing s=log(A), its derivative has the sign
of `4-exp(-s)*(s+2)`, positive for s>=0. Thus A=16 controls all 1<=A<=16,
and using the lower log bound in its denominator produces a valid upper
comparison ratio.

The reported strict gap therefore has a valid interpretation inside this
comparison family: 2.970596100510... exceeds 2.956286971652....
