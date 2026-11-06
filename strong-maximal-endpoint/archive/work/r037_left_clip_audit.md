# r037: left coarse clip — local audit

Date: 2026-11-06.  Here `L` denotes the affine/Gamma left screen and
`a=L-v` (called `a_t` in r035), while the nonlinear inner branch is
`I=v exp(U/v)` and the actual candidate is `W=L wedge I` (the right screen
is irrelevant to this note).

## 1. The clip has an exact, favourable algebraic interface

For every positive linear action operator,

`K_cW=K_cI-K_c[(I-L)_+] <= K_cI`.                                  (1)

Thus no monotonicity of `I` is needed to get an upper comparison with the
inner branch.  Conversely `K_cW<=K_cL`, so if the left screen is a global
supersolution then `K_cW<=L_t`.  If both branch inequalities hold at the
current point, positivity gives the exact minimum rule

`K_c(L wedge I) <= (K_cL) wedge (K_cI) <= L_t wedge I_t`.

At a crossing there is no boundary distribution or differentiability term:
the leakage in (1) has the favourable sign.  This is the correct way to
treat the clip boundary; it does not require the raw inner branch to be
globally decreasing.

## 2. The relevant deficit is not exactly `a-U`

Put

`E=v(e^(U/v)-1)-U >= 0`,  `C_true=L-I=a-U-E`.

The linear ratio `q=U/a` is nevertheless a valid *one-sided* locator:
`I<=L` implies `U<=a`, hence the actually active inner branch necessarily
has `q<=1`.  The converse fails in a strip adjacent to `q=1`, since
`E>0`.  A proof must not identify `{q<=1}` with the active set.  This is
not a problem for (1), but it matters if derivatives of the deficit are
used.

The lower saddle estimates from r035 give, uniformly for
`A sqrt(log(t)/t)<=lambda<=C t^(-1/3)`,

`a/f = 1/[lambda(1+lambda)]-1/(N lambda^2)+o(1/lambda)`,
`a'/f = 1/(1+lambda)-1/(N lambda)+o(1)`,
`a''/f = lambda/(1+lambda)-1/N+o(lambda)`.                         (2)

The `-R-Put/N` pieces in the exact identity for `a` must be retained to
obtain (2); replacing `a` by an unspecified `F` is insufficient near the
coarse transition.

With the bridge expansion for `U`, `q=kappa t lambda^3(1+o(1))`; hence
actual activity forces `lambda<=C t^(-1/3)`.  Since
`d lambda/dr=-(t v_lambda)^(-1)`, differentiating `U=qa` yields

`(a-U)' = (1-q)a' -q'a
        = f[(1-q)(1+O(lambda))+c_1 lambda+o(lambda)]`,

where `c_1>0`; the crucial positive term is `-q'a`.  A further derivative
gives

`(a-U)'' = lambda f[(1-q)(1+O(lambda))+c_2 lambda+o(lambda)]`,

with `c_2>0` (the `-2q'a'` term supplies the boundary buffer; `-q''a` is
smaller by one power of lambda).  These are the needed discrete
boundary-layer expansions.  In particular, the apparent degeneration at
`q=1` is not real.

## 3. Exponential correction and choice of the lower cutoff

One must transfer the preceding estimates from `a-U` to `C_true`.  On this
range the tilted density bound implies

`a/v <= C f/lambda <= C t^(-1/2)lambda^(-1)
 exp(-c t lambda^2)`.

After choosing the fixed lower band constant `A` large enough, this is
`o(lambda)` uniformly.  In the active set `U<=a`, so `w=U/v=o(lambda)`.
Differentiating `E=v(e^w-1)-U` with the same saddle bounds gives

`|E'|=o(lambda f)`, `|E''|=o(lambda^2 f)`.

Consequently the positive derivative/curvature buffers above persist for
the true nonlinear deficit even in the thin layer where the clip switches.
This is exactly why the lower tube has to begin at a sufficiently large
multiple of `sqrt(log(t)/t)`.  It also shows why a proof that works only
with `C=a-U` is incomplete unless it explicitly accounts for `E`.

## 4. Consequence for the local MAX calculation

The usual compact-action Taylor calculation applied to the decreasing,
convex true cost has its leading gain reduced to

`C f lambda^2[(1-q)+O(lambda)]`.

It is bounded by the reference quadratic source after any fixed inflation
of that source; the quadratic/third-order remainders are smaller uniformly
in the stated range.  For paths whose previous state lies on the other
side of the clip, (1) supplies an additional negative leakage rather than
a positive crossing term.  Thus the proposed negative-intermediate local
MAX extension is viable, conditional on writing the uniform saddle
remainders in (2) and the analogous bridge derivatives for `U`.

This does **not** establish a global supersolution: it says nothing about
the positive large-deviation region, nor about continuing an arbitrary
controlled trajectory after it exits this left band.
