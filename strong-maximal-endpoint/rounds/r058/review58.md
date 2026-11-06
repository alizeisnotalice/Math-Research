# Independent audit of r058

## Verdict: accept (scoped coefficient and fixed-family claims)

I checked the displayed proof against `r057/common_window57.md` and the
relevant r054/r056 identities, rather than treating the stored green checks as
a proof.  Within the stated full-positive-quadrant, one-receiver, two-equal-
block, cap-two, target-one, regular stationary setting, the conclusions in
`proof58.md` and `budget58.md` follow.  They do **not** establish a uniform
Taylor/IFT radius, a nonzero lower bound for `K`, a multi-receiver result, or
an endpoint theorem for the strong maximal constant; the round explicitly
keeps all of those outside its scope.

### Mixed response coefficient

The sign and centering in
\[
 K=-B\bigl(k_-,s_{t,-}+(s_+-\beta)k_-\bigr),
 \qquad \beta=B(k_+,s_+)/B(1,k_+),
\]
are consistent with the full lower-coarea-plus-soft response form.  In
particular, differentiating the contrast response needs the mixed score
derivative `s_{t,-}`; a Fisher bound for `s` alone would not suffice.  The
algebraic checks were independently rerun from `code/analytic58.wl` and its
output exactly matches `rounds/r058/analytic58.json`.

The central estimate is supported by the r054 compact `C^1` tilt bridge:
the exact score formula is `s=bF(c,d)+O_{C^1}(1)`, with
`F_d(mu,mu)=0` (equivalently `u'(mu)=0`) and `d_t=O(1/b)`.  On the chosen
moderate band this gives `s_t=O(sqrt(L/N))`; the identity
`tk+cs=b m(c)-1` simultaneously gives `k-1=O(sqrt(L/N))`.

The tail argument is adequately finite-weight.  The shifted-Gamma formulas
give the required uniform-in-level pushforward maxima only for
`1,k,k^2,|s_t|,k|s_t|`; combined with fixed-order inverse-Gamma moments and
Chernoff/Hölder on the bad block, this is enough for every weight actually
occurring in (1).  It is important, and correctly stated, that this is not a
claim for arbitrary polynomial weights.  Enlarging the fixed band constant
after fixing the finite list of weights yields the claimed negligible tails.

Using the r057 Ward normalization and profile estimate gives
`r \gtrsim N/L`, `B00 \asymp 1/r=O(L/N)`, while the inherited multiplier
bound gives `r<=2N`.  The bootstrap `beta=O(L)` is then valid, and the
moderate integrand is `O(L^2/N)`.  Hence
\[
 |K|=O(L^3/N^2).
\]
Finally, r056/r057 give `Gaa=Nm''(c)/r+E[HZ_+]<=-p_mu/4` for large `N`, and
`Gd=2/r>=1/N`; thus the advertised coefficient bound
`|K^2/(8GaaGd)|=O(L^6/N^3)` follows.  This is only a coefficient estimate,
as the text says.

### Fixed-arity Orlicz family

The cell accounting uses the actual `Phi_N`.  There are `M` cells at
`1-delta` and `M-1` at `1+delta`; after the linear contributions cancel, the
exact modular is
\[
 J_N=M^2-\delta+(M-1)D_N(\delta),
 \quad D_N=\Phi_N(1+\delta)-(1+\delta).
\]
The common interval average is `1-delta/M^2`, whereas deletion of any one
fiber supplies a prefix on which every remaining fiber averages one.  Thus
the point regret is exactly `delta/M^2` in the stated construction.

The derivative lower bound for
`log(1+log(1+s))`, followed by the elementary normalization inequality,
correctly yields (4)--(6).  In particular a regret at least
`kappa (log N)^2/N` forces a normalized log-cost of order at least
`(log N)^2` (for all sufficiently large `N`, uniformly over the permitted
`M>=2`).  I reran `code/budget58.wl`; its 29 exact checks and serialized
output match `rounds/r058/budget58.json`.  The numerical rows remain only
diagnostics, exactly as labeled.
