# r075 review: power-policy weighted pushforward is sharp

## Verdict

The proposed power policy is correct.  It gives an injective, smooth,
radially monotone critical-window policy for which the *remaining* Jacobian
loss after the exact rectangle-volume factor is exactly `1/q`.  Thus `q=1/N`
or `q=1/N^2` produces unavoidable bare losses `N` or `N^2`.  Its `q=0`
limit supplies the stronger singular-pushforward test.  These are geometric
sharpness statements only: the probes below are not asserted to be legal
physical source costs satisfying all efficiency constraints.

## Exact power-policy calculation

Fix `r_0>0`, write `r=|x|`, `u=log(r/r_0)`, and for `0<q<1` put

\[
 \lambda(u)=e^{-(1-q)u},\qquad
 T(x)=\lambda(u)x,\qquad |T(x)|=t=r_0e^{qu}.              \tag{1}
\]

Choose the action by

\[
 e^{-c(u)}=(1+e^{-c_r})\lambda(u)-1.                     \tag{2}
\]

Restrict `u` to the interval on which the right side gives
`c(u) in [c_r-rho,c_r+rho]`; it is a nontrivial interval of length
`Theta(rho)`.  The policy is smooth and increasing in `u`.  This is fully
consistent with a shrinking critical action window.

Since `x dot grad log lambda=-(1-q)`, the exact determinant is

\[
 \det DT=q\lambda^N.                                      \tag{3}
\]

The radial map `u mapsto log(t/r_0)=qu` is strictly increasing, so the map
is one-to-one on a cone-annulus.  The area formula therefore gives, for
every nonnegative `g` supported in `T(A)`,

\[
 \boxed{\ \int_A\lambda(x)^N g(Tx)\,dx
       ={1\over q}\int_{T(A)}g(Y)\,dY.\ }                \tag{4}
\]

This directly audits the relevant cancellation: the `lambda^N` on the left
is precisely the rectangle-volume factor that re-anchoring can cancel, but
the slope factor `q` in (3) remains.  Consequently no uniform geometric
pushforward bound follows from volume matching.  Taking `q=N^{-1}` or
`q=N^{-2}` gives exact losses `N` or `N^2`, respectively, while every action
in (2) stays in the same critical window after the indicated restriction.

## Collapse limit and Gaussian tube test

At `q=0`, `lambda=e^{-u}=r_0/r`, hence `T(x)` lies on the sphere of radius
`r_0`.  In polar coordinates, with `domega` denoting the angular measure,

\[
 dx=r_0^N e^{Nu}\,du\,domega,\qquad
 \lambda^Ndx=r_0^N\,du\,domega.                           \tag{5}
\]

Thus a cone-annular set with radial log-width `Delta u` has nonzero weighted
receiver mass proportional to `Delta u`, even though its image has zero
Lebesgue measure.  To see the failure through nonsingular test functions,
use

\[
 g_\epsilon(Y)=\exp\!\left[-{\log(|Y|/r_0)^2\over2\epsilon^2}\right]. \tag{6}
\]

After removing the common angular factor and `r_0^N`, its target integral is

\[
 \int_{-\infty}^{\infty}e^{-v^2/(2\epsilon^2)}e^{Nv}\,dv
 =\sqrt{2\pi}\,\epsilon\exp\{(N\epsilon)^2/2\}.        \tag{7}
\]

For fixed `N`, this tends to zero as `epsilon` tends to zero, whereas the
left side of the weighted pushforward on any fixed `Delta u` interval equals
its positive mass because `g_epsilon(Tx)=1`.  A compactly supported radial
cutoff can be added to (6) without changing this conclusion.  This confirms
the claimed singular limit.  The Gaussian is a geometric probe, not a claim
about an admissible source modular or an actual-policy cost.

## Scope checks

* The policy `c(u)` is monotone, but the meaningful no-collapse condition is
  positivity of the image-radius slope `d log t/du=q`, not monotonicity of
  `c` itself.  The `q=0` policy is a direct counterexample to conflating the
  two.
* Formula (4) concerns only a weighted receiver pushforward.  It does not
  eliminate the separate receiver-dependent source pullback
  `f(Y+(k/lambda)X)`, which exists even for a constant policy.
* Hence the calculation rules out a bare re-anchoring/Fubini improvement at
  sharp precision.  It neither constructs a full `C_N` counterexample nor
  disproves a theorem using additional action-indexed source packing.
