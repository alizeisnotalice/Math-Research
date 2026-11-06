# r028 feedback Fourier band audit — 2026-11-06

Read-only audit of the final “Fixed-`mu` saddlepoint Fourier lemma: proof
record” in `work/strong-maximal-global-adversary.md`.

## Result

Accepted as a *pointwise* moderate-lower-tail / one-step feedback lemma.  It
does not settle the controlled-occupation issue, so it does not by itself
prove a global feedback or endpoint theorem.

For `s=k mu+x`, with `-C sqrt(k log N) <= x <= -c sqrt(k)`, the negative
tilt has mean `s/k`; untilting gives

`f'(s)/f(s) = lambda + O(1/k)`,

with `lambda=-x/(mu k)+O(x^2/k^2)>0`.  Thus the centered one-step
difference has leading term `mu f'(s)/2`, and its signs agree with the
claimed feedback: `a>0`, `b>0`, and the optimizing displacement is
`-a/b = x/(2k)(1+o(1))`.  The leading gain is
`p x^2 f(s)/(8 k^2)`; retaining a conservative `p/12` is justified.

The Fourier statements have the right centering and orders.  The odd
Gaussian contribution to the first derivative vanishes, giving
`g'(0)=O(k^-3/2)`, while the second derivative gives
`-g(0)/(k v_lambda)+O(k^-2)`.  The exponential factor in the global
derivative bound is sufficient for the Taylor interval used here.

## Rate wording to retain

The coarse displayed Taylor remainder `O((1+log N) f/k)` alone yields a
relative error `O(log N/sqrt(k))` at the bottom of the band.  This still
tends to zero because `k >= (log N)^4`, and is enough for a qualitative
strict gain.  The sharper stated

`rho = lambda + 1/(k lambda) = O(sqrt(log N/k))`

requires retaining the stronger consequence of the Fourier derivative
estimates, namely `O((lambda^2+1/k)f)` for the relevant Taylor remainder.
That stronger estimate is compatible with the preceding bounds, but should
be cited explicitly rather than inferred from the coarse remainder.
