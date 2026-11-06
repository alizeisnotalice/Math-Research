# r056 independent review

## `shared56.md`

Accepted in its final stated scope.  The exact regret identity and the
layer-cake lower bound are correct and tie-invariant.  The final version
correctly distinguishes two notions which must not be conflated:

* for the unrestricted continuous common-action set in the two-slab example,
  the intermediate action covers both fibers up to height `1/3`, so the
  actual layer deficit is `1/2` only on `(1/3,1/2)` and its integral is
  `1/12`;
* the endpoint payoff matrix has regret `1/8`, but it cannot by itself
  certify the continuous layer coverage statement.

The finite endpoint matrix does exactly certify each *grid* outer interval's
`A_J-S_J`; it does not verify the universal condition (4) for all fractional
outer intervals without a separate outer-interval arrangement/continuity
argument.  The document now states this finite-grid limitation, retains the
actual Orlicz denominator, and makes no general endpoint claim.

## Pending `common56.md`

For the proposed common-stationary theorem, the algebraic route is sound
provided the final text explicitly records the following:

\[
 T_{++}=N m''(c)P/\lambda+E[H Z_+],
 \qquad Z_+=e_1+e_2+2s_1s_2,quad EZ_+=0,
\]
and the r050 stationary Hessian is
\[
 (\log P)''=-{\lambda\over P}(T_{++}+B_{++})
 +{\lambda\over P}{B_{0+}^2\over B_{00}}.
\]
Positive semidefiniteness gives `B_{0+}^2/B_{00}<=B_{++}`.  The r053
soft/coarea proof applies to the sum score as well as the difference score,
because `E(s1+s2)^2=2I_b`; hence its scaled response is `o(N)`.  The r054
profile estimate and `||Z_+||_2=O(log^2 N)` make the Ward residual scaled
`o(N)`.  At a stationary point this proves common curvature
`N p_mu+o(N)`.

The resulting statement is only for regular stationary points in the r054
window.  If the value is `C^2` on the entire interval, positivity at every
zero of `P'` implies at most one stationary point; it does not prove that a
stationary point exists.  For fixed `H`, the action second derivative is
`T_{++}` and is strictly negative after using `P/lambda>=1/(2N)`, so the
fixed-profile mean has a local maximum.  None of this extends to arbitrary
input fibers or `C_N`.

## Final `common56.md`

Accepted.  The final file includes precisely the required off-stationary
Schur/gradient terms and does **not** claim a whole-window curvature bound.
At `M_+=0`, the r050 formula reduces to
\[
 (\log P)''_+=-{\lambda\over P}T_{++}
 -{\lambda\over P}\left(B_{++}-{B_{0+}^2\over B_{00}}\right),
\]
with the bracket nonnegative.  Its Ward decomposition, zero mean of
`Z_+`, profile estimate, and sum-score response bound give
`Np_mu+o(N)` uniformly at regular stationary points.  The uniqueness
statement correctly remains conditional on `C^2` regularity over a connected
window and asserts no existence.

The fixed-input consequence is also correctly normalized: `lambda/P<=2N`
implies `NP/lambda>=1/2`; since the normalized curvature is negative,
this is enough to deduce a strict local maximum of the stationary source's
mean under common receiver action.  It is explicitly not an arbitrary-input
or strict-target theorem.

## Final angular IFT corollary

Accepted.  In the final `common56.md`, the bounded KKT profile is convolved
against the zero-extended block densities, so derivatives are legitimately
placed on those densities; `b>=6` supplies the stated local `C^4`
regularity despite the profile's KKT jump.  Exchange symmetry makes
`F(v,a)` even in `v`, and the fixed-input stationary calculation gives
`F_a(0,c)=0`, `F_aa(0,c)<0`.  The implicit-function expansion is therefore
\[
 a(v)=c-{F_{avv}(0,c)\over2F_{aa}(0,c)}v^2+o(v^2),\qquad
 F(v,a(v))-F(v,c)=-{F_{avv}(0,c)^2\over8F_{aa}(0,c)}v^4+o(v^4).
\]
The factor and sign are correct.  The file correctly makes no claim that
the mixed coefficient is nonzero, nor any uniform angular/N-scale bound.

## Exact counterexample to a converse layer certificate

The proposed two-fiber construction is valid, with all inner intervals
required to contain `x_0=0`.  Over the two equal outer halves take
\[
 g_1={\bf1}_{[0,1]},\qquad
 g_2={\bf1}_{[0,1]}+3{\bf1}_{[1,2]}.
\]
Their adaptive values are `1` and `2`, hence `A_J=3/2`.  A common interval
`[0,r]` has mean
\[
 1\;(r\le1),\qquad(3-1/r)/2\;(1<r\le2),\qquad5/(2r)\;(r>2),
\]
so `S_J=5/4` and the common-action regret is `1/4`.  Wolfram 15 actually
maximized this piecewise formula and returned `5/4`.

Nevertheless `D_J(t)=0` for every height: `[0,1]` covers both active fibers
for `t<1`, while `[0,2]` covers the only adaptive-active second fiber for
`1<=t<2`.  Thus the layer deficit remains a sufficient lower certificate
only; no converse or characterization of positive common-action regret is
possible from this heightwise coverage statistic.  The exact `Phi_2` cost
of the input is `2+3(1+log 3)=5+log 27`, as required for a fully specified
finite step example.
