# r029 Poisson-corrector route — preliminary audit (2026-11-06)

## Status

No algebraic counterexample found.  The corrector is a plausible way to
replace a controlled-occupation lower bound, but it is not a completed
endpoint result until the four listed proof interfaces are formalized.

## Checks that pass

Let `p_c=(1+exp(c))^-1` and `m(c)=1+c p_c`.  Direct differentiation gives

`m'(c)=p_c-c p_c(1-p_c)=-p_c(c-m(c))`.

Thus the first parameter derivative of `K_c U` has its `U'` coefficient
of size `p_c(c-m(c))`; at `c0` it is indeed controlled by
`|c0-m0|`, rather than by a constant.  A local MATLAB calculation at
`c=mu-10^{-3,-4,-5}` gave `|m'(c)|/|c-m(c)|` converging to `p_mu`, as
required.

For the Duhamel sum `U_t=sum_{j<=t}P^{t-j}g_j`, the proposed split is
consistent.  If `j<=t/2`, a common tilt at mean `r/t` cancels the moment
factors; the global tilted prefix derivative estimate and `phi_j=O(1/j)`
give `f_t(r) exp(C|u|/sqrt(t))/(j t^{ell/2})`, whose sum costs only
`log N`.  If `j>t/2`, differentiating the smooth compact-band source gives
`|g_j^(ell)|<=C f_j/j^(1+ell/2)`, and convolution gives the same required
bound after summing `O(t)` indices.  The all-shift factor is integrable
against the exponential action tail once `t` is large.

Since `psi=theta(2-theta)>=theta`, the source dominates the action-error
factor on the transition portion of the cutoff.  With
`t>=L^K`, `K>2`, and `t<=N/(C1 L)`, the displayed relative error
`O(|c0-m0|L+L/sqrt(t))` tends to zero.

## Interfaces still needed before a theorem

1. Fix the exact `t`/`t-1` convention.  Prove the actual greedy gain,
   initially expressed through the `t-1`-future density, dominates the
   source written with `f_t`; this needs a uniform density-ratio lemma.
2. State the Fourier/derivative bounds uniformly for the triangular kernel
   `c0=c0(N)` and both signs of the local saddlepoint tilt, not only the
   earlier fixed-`mu` lower-tail lemma.
3. Write the exact supersolution induction
   `K_{pi_t}(V_{t-1}-U_{t-1}) <= V_t-U_t`, including the `g_t` index.
4. Prove the conditional tilted local-limit bridge uniformly under the
   continuous conditioning `S_N=q`, for
   `j in [L^K,N/(C1 L)]`.  Its cutoff second-moment lower bound is the
   indispensable justification for `U_N(q)>=const f_N(q) log N`.

These are proof obligations, not evidence against the route.  In
particular, no uncontrolled path-occupation assertion is needed if the
supersolution and bridge are established.
