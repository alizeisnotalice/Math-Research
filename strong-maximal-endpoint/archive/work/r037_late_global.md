# r037: late-time global coarse splice audit

Date: 2026-11-06.  This audits the proposed late-time proposition, not an
initialization of the Bellman comparison.

Let `a=a_N` be the r036 matched parameter, `delta=mu-m(a)`, and let
`I_t=v_t^a exp(U_t/v_t^a)+Lambda t b^2/N`.  Let
`O_t=O_t^- wedge O_t^+` be the two r034 global barriers and
`W_t=O_t wedge I_t`.  The following argument is valid for every fixed
`gamma in (0,1)` after choosing fixed constants in the stated order.

## Outside the natural tube

The exact terminal decomposition used in r034 gives, on the left,

`O_t^- - v_t^a <= t delta/N + P_a(S_t<r)+B_t^Gamma`.

Here the omitted reference `B` term has the favourable sign.  At
`r<=t m(a)-A sqrt(t log t)`, uniform Cramer/Chernoff and the Gamma piece give

`O_t^- - v_t^a <= t delta/N+C t^(-K)+C exp(-ct)`.

Both estimates are uniform in `N>=t` because the matched `a_N` remains in a
fixed compact neighbourhood of `mu`.  On the right, for the same (or a
larger) fixed `A`,

`O_t^+ <= C t^(-K)`

whenever `r>=t m(a)+A sqrt(t logt)`.  Indeed
`t delta=O(t log N/N)=o(sqrt(t logt))` uniformly for `t>=N^gamma`, so the
latter region is, after reducing the fixed constant in the radius, to the
right of `t mu+A' sqrt(t logt)`, where r034's positive-Chernoff screen
applies.

Choose `Lambda` so `Lambda b^2-delta>=b^2`; this is possible because
`delta=O(b^2)`.  Since `I_t>=v_t^a+Lambda t b^2/N`, it is enough to make
the displayed tails `o(t b^2/N)`.  With `b^2 asymp logN/N`, at the worst
time `t=N^gamma` this is exactly

`N^(-gamma K)=o(N^(gamma-2) logN)`,

which follows from `K>2/gamma-1`.  (Wolfram locally verified this limit.)
Thus, for all large `N` and all `N^gamma<=t<=N`, `I_t>=O_t^-` on the left
exterior and `I_t>=O_t^+` on the right exterior.  Since `O_t<=O_t^-` and
`O_t<=O_t^+`, respectively,

`W_t=O_t` throughout the exterior.

## Pointwise late-time Bellman comparison

At an exterior current state, positivity gives

`T W_(t-1)<=T O_(t-1)<=O_t=W_t`.

At a current state in the matched natural tube, r036 supplies
`T I_(t-1)<=I_t`; independently `T W_(t-1)<=T O_(t-1)<=O_t` and
`T W_(t-1)<=T I_(t-1)<=I_t`.  Hence

`T W_(t-1)<=O_t wedge I_t=W_t`.

No hypothesis about where the *previous/future* state lies is needed; it is
only the pointwise domination `W_(t-1)<=O_(t-1),I_(t-1)` and positivity of
each `K_c` that is used.  This also avoids the old crossing/occupation
gap.

## Scope and the remaining gap

This proves `T W_(t-1)<=W_t` for every state and every
`ceil(N^gamma)<t<=N`, for sufficiently large N.  The strict inequality in
the polynomial tail comparison has enough margin to absorb all fixed
constants and the `t` versus `t-1` tube-centre displacement.

It does not prove an endpoint upper bound by itself.  At the new initial
layer `t0=ceil(N^gamma)`, one only knows `U_true,t0<=O_t0`; the clipped
value `W_t0=O_t0 wedge I_t0` can lie below the true Bellman value.  A proof
must supply a dominating initialization/transition there without losing the
candidate terminal correction.  The old fixed-4t shift cannot simply be
substituted with a shift of order `N^gamma`.
