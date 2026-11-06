# r035: coarse-cut monotonicity and the exact annular residual

Notation: `P=K_mu`, `v_t=P^t phi_N`, and let `U_t` denote the centred
quadratic Poisson correction from r032, so that the candidate inner branch
is `I_t=v_t+U_t`. Let
\[
 O_t=O_t^-\wedge O_t^+,qquad W_t=O_t\wedge I_t . \tag{1}
\]
Both outer branches are global MAX supersolutions and nonincreasing.

## 1. A useful exact pasting/residual identity

The candidate (1) is globally defined (the exact r032 formula defines `U`
at every positive state), but it is not automatically a supersolution. At a
state with `W_t=I_t`, the *actual* residual is
\[
\begin{aligned}
 {\cal T}W_{t-1}(r)-I_t(r)
 \leq \sup_c\{&(K_c-P)v_{t-1}(r)+(K_c-P)U_{t-1}(r)-g_t(r)\\
 &-K_c[(I_{t-1}-O_{t-1})_+](r)\}. \tag{2}
\end{aligned}
\]
This follows from `W=I-(I-O)_+`. The final term has the favorable sign.
Formula (2), rather than a reference occupation estimate, is the exact
annular quantity needing a bound. In the strict interior where the old inner
comparison holds, the first line is the local Taylor residual; in an annulus
the negative leakage term may help but cannot simply be omitted if one seeks
sharp constants.

The generic sufficient monotonicity condition is also simple:

**Cut monotonicity lemma.** If `O` is nonincreasing and `I` is continuous,
locally absolutely continuous, and `I'(r)<=0` at almost every point with
`I(r)<O(r)`, then `O wedge I` is nonincreasing.  This is proved component by
component on the open set `{I<O}`; at endpoints continuity supplies the
ordering. Therefore it is enough to show that the coarse screen cuts off all
positive-derivative portions of `I`; one need not prove global monotonicity
of raw `I`.

## 2. Why the lower coarse screen is plausibly exactly the needed cut

For `r<rho_N`, the Gamma correction in `O^-` vanishes. Writing
\(\psi_N(u)=1-e^{-u}-u/N\) on this range gives the exact identity
\[
 a_t(r):=O_t^-(r)-v_t(r)
   =P^t\psi_N(r)=F_t(r)-R_t(r)-{\operatorname{Put}_t(r)\over N}. \tag{3}
\]
Here `Put_t=E(r-S_t)_+`.  Thus `I_t<=O_t^-` is exactly `U_t<=a_t`; it is not
a loose comparison.

Consider the lower saddle range `eta=-lambda<0`, with
`lambda sqrt(t)->infinity` but `lambda=o(1)`, and with `r` still a fixed
fraction of `t` away from zero. The common-tilt local limit expansion, after
integrating its uniformly controlled shifted-density ratio, gives
\[
 {a_t\over f_t}={1\over\lambda(1+\lambda)}-{1\over N\lambda^2}
 +o\!\left({1\over\lambda}\right),\qquad
 {R_t\over f_t}={1\over1+\lambda}+o(1). \tag{4}
\]
The first formula is simply the tilted integral of
`1-exp(-u)-u/N`.  The centred Poisson identity and the tilted posterior
moment expansion give
\[
 {U_t\over f_t}={p_\mu\mu^2\over8}t\lambda^2
 +O(t\lambda^3+\log t),\quad
 {U_t'\over f_t}=\lambda{U_t\over f_t}
 -{p_\mu\mu\over4}\lambda
 +o(\lambda)+O(\log t/t). \tag{5}
\]
The derivative in (5) is important: it has the negative derivative of the
quadratic factor, not merely `lambda U/f`.

Consequences, conditional only on making the usual uniform saddle remainders
in (4)--(5) explicit:

* If lambda is bounded below by a fixed positive constant, `U/a` is of order
  `t`, so raw `I` is discarded by `O^-` for all large `t`.
* In the intermediate lower range, `U<=a` implies `t lambda^3=O(1)`.  Insert
  `U/f<=a/f` into (5):
\[
 {U_t'\over f_t}\leq {1\over1+\lambda}
 -{p_\mu\mu\over4}\lambda+o(\lambda)
 <{R_t\over f_t}\leq {-v_t'\over f_t}. \tag{6}
\]
  Hence `I_t'=v_t'+U_t'<0` wherever the inner branch is selected in this
  range.  Notice that the margin comes from the derivative of the quadratic
  source; dropping it would give no strict inequality.
* In the central range `lambda=O(sqrt(log(t)/t))`, r033's derivative bounds
  already give `U_t'/(-v_t')=o(1)`.  The corresponding positive-saddle range
  is analogous: its possibly positive term is `O(|eta|f_t)`, whereas
  `-v_t'` is at least a fixed multiple of `f_t` on the natural tube.

This provides a concrete route to proving the selected-branch hypothesis of
the cut monotonicity lemma. It is **not yet a complete global proof**: one
must supply uniform tilted-density/posterior derivative bounds through the
transition `r approx rho_N` and show that the positive side is either cut by
`O^+` or has the stated derivative bound.

## 3. Actual computation, limited evidence

`r035_ratio_probe.wl` evaluates the exact finite mixture formula for (3) and
the exact r032 Poisson formula, at `N=256`, `t=16,32,64`, on grids below
`rho_N`. Its sampled maximum successive increments of `U/a` were
`-1.68e-4`, `-1.08e-4`, and `-7.11e-5`, respectively; the ratio falls from
large values near zero to below `5e-4` near the sampled right endpoint.
This is numerical evidence for the predicted decreasing crossing, not a
certificate and not a replacement for (4)--(6). WolframScript was used and
exited; MATLAB was not used.

## 4. Scope

The natural candidate is therefore (1), not a monotone majorant. A monotone
majorant `sup_{x>=r} W_raw(x)` would be an *upward* modification and has no
automatic MAX supersolution property; it should not be introduced unless its
uplift is separately controlled. The present calculation makes the annular
residual (2) explicit and identifies a plausible analytic route to eliminate
the raw-inner monotonicity objection. It does not establish the all-action
inequality (2), and hence does not prove the sharp Bellman upper bound.
