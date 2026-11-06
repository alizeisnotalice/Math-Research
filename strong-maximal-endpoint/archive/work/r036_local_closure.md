# r036: nonlinear Doob correction — current closure audit

Date: 2026-11-06.  The Doob/Jensen module in `work/r036_doob.md` correctly
controls a nonlinear *Jensen defect* under its stated accessible-state
weighted shift hypotheses.  It does not, by itself, close the all-action
Bellman supersolution.

## What follows from the supplied estimates

For `W=v exp(w)`, `w=U/v`,

`W'=exp(w)[(1-w)v'+U']`.

On the natural tube, the stated `U'/d=o(1)` makes `W` decreasing there, and
the displayed derivative estimates give the claimed local `W''` expansion.
The Doob kernel identity is exact:

`K_c W(r)=K_c v(r) Ehat_{c,r}[exp(w(r-Z_c))]`,

where `Ehat` has the density from r036.  Under the verified weighted
Lipschitz and boundedness hypotheses, its Jensen excess is at most the
reference-summable source

`C f_t(r)(log t)^4/t^(3/2)`.

This part is valid conditionally on the global shifted-density estimates
invoked in r036.

## Precise missing all-action term

Taking logarithms in the preceding identity gives

`log(K_cW/v_t)=log(K_cv/v_t)+K_cU/(K_cv)-U_t/v_t+J_c`,

with `0<=J_c` the Jensen defect.  If the proposed Poisson potential obeys
`U_t=P U_{t-1}+g_t`, the first-order expansion of this expression necessarily
contains

`(K_c-P)U_{t-1}`

in addition to the reference action gain `(K_c-P)v_{t-1}`.  Thus a
supersolution still needs a uniform inequality of the form

`g_t(r) >= (K_c-P)v_{t-1}(r)+(K_c-P)U_{t-1}(r)
            + (controlled quadratic/Jensen errors)`             (C)

for every action `c`, or an equivalent exact cancellation.

The present weighted Lipschitz estimate only gives, after comparing the
two action kernels,

`|(K_c-P)U_{t-1}(r)| \lesssim
 f_t(r)(log t)^(3/2)/sqrt(t)`

without an additional action-derivative cancellation.  On a natural tube,
the intended quadratic gain is of order

`f_t(r) log(t)/t`.

The available bound is larger by a factor of order `sqrt(t log t)`, so it
cannot be absorbed by the proposed `f_t(log t)^4/t^(3/2)` Jensen source or
by the reference-summable errors from r035.  The global root for `v` does
not remove this term: it localizes the maximizing action for `v`, whereas
the perturbation changes the action objective.

## Conclusion

The candidate has not yet proved a nonlinear `W` supersolution.  A new
calculation is needed: either an exact parameter expansion/cancellation for
`K_cU-PU` matching the quadratic source, or a different potential whose
action transport has a favorable sign.  Until (C) is established, no
all-action closure, no sharp payoff bound, and no endpoint coefficient may
be claimed.  This is a specific gap, not an objection to the Doob/Jensen
or reference-Poisson summability modules themselves.

## Addendum: the moment cancellation repairs the linear term

The preceding conclusion was based only on a Lipschitz transport bound. The
new exact action derivative, assuming the stated operator identity, is much
sharper. With `b=m(a)-a`, Taylor expansion at `r-a` gives

`P U(r)-U(r-a)=-b U'(r-a)+.5 E_a(Z-a)^2 U''(r-a)+R_3`.

Hence the proposed global derivative bounds and `t<=N` imply

`|partial_c K_cU|_{c=a} <= C f_t(r)(log t)^2/t`.

The sign of the first term is `-b` (immaterial for the bound). In fact,
`b|U'|<=C f_t(logt)^2/sqrt(Nt)<=C f_t(logt)^2/t`, and the second and third
terms have that scale or better. If `B_U=p_aU'+O(bU'+U'')`, then
`B_U=o(d_t)` on the tube. Therefore the corrected *linear* candidate
`I=v+U` has derivative `A_v+A_U`, curvature
`B_v+B_U<=-(p_a/2)d_t`, and Young adds only
`C f_t(logt)^4/t^2` to the quadratic reference source. This is
reference-summable. The weighted shift estimates also make
`K_cU-K_aU=o(f_t)` outside a fixed action neighbourhood; combined with a
written quantitative reference exterior gap, this is a viable local
all-action closure for `I`.

For the exponential candidate the stated curvature scale needs correction.
Put `w=U/v`. Although the linear term cancels in
`e^w(1-w)-1=O(w^2)`, the available scales
`U=O(f_t logt)` and `f_t/v=O(sqrt(logt/t))` yield
`w=O((logt)^(3/2)/sqrt(t))`. Thus its curvature correction is generally
`O(f_t(logt)^3/t)`, not `O(f_t(logt)^2/t)`. This is still harmless after
Young (`O(f_t(logt)^6/t^2)` is summable), but a nonlinear proof must use
that larger error, or prove a sharper cancellation, and explicitly carry
the perspective/Doob action comparison. The shorthand `e^w g-Jensen`
does not itself establish it.

## Local MAX lemma (conditional interfaces explicit)

Let `P=K_a`, `v_t=P v_{t-1}`, and
`U_t=PU_{t-1}+g_t+C_0t b^2/N`.  On a natural tube, assume quantitatively:

1. `v'=-d`, `v''=O(f(|eta|+t^-1))`, an integrable weighted third derivative,
   and the strict mean gap `m(c)<=mu-kappa_0` outside a fixed neighbourhood
   of the unique maximum `mu`, while `m(a)=mu+O(b^2)`.
2. The three global weighted shift estimates
   `|U'|<=f log^(3/2)t/sqrt(t)e^(kz)`,
   `|U''|<=f log^2t/t e^(kz)`,
   `|U'''|<=f log^K t/t^(3/2)e^(kz)`, with `k<1/4`.
3. The r035 reference-action gain and its stated errors.

Then, for large enough fixed `t0`, a nonnegative global extension of

`g_t=(1+2eps)(p_a/8)f_t((r-t ell)/t)^2
      +B f_t log^6t/t^(3/2)`

gives locally `sup_c K_cW_{t-1}(r)<=W_t(r)` for
`W_t=v_t exp(U_t/v_t)`, if `B` is fixed sufficiently large.  The compact
action computation has `A_U=O(f log^2t/t)`, `B_U=o(d)`, and the exponential
remainder costs only `Ceps f log^6t/t^2`.  The strict mean gap gives
`K_cv-K_av<=-kappa d` outside the compact action set, while the weighted
estimates give `|K_c(W-v)-K_a(W-v)|=o(f)` uniformly there.  At `a`, Doob and
Jensen yield

`P W_{t-1}<=v_t exp(PU_{t-1}/v_t)+C f log^4t/t^(3/2)`,

and positivity gives exactly
`W_t-v_t exp(PU_{t-1}/v_t)>=g_t+C_0b^2/N`.
Choose `B` first and then `t0` so
`B sum_{j>=t0}log^6(j)/j^(3/2)<=1`; this prevents a bootstrap of potential
constants.  The affine term settles the per-step `b^2/N` error and is
`O(b^2)` at terminal time.  This is strictly a local lemma: a global outer
transition and proofs of the three weighted interfaces are not supplied.
