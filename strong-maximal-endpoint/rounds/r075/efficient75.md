# r075 review: low-cost successful policy still has exponential re-anchor loss

## Verdict

The proposed construction is valid.  It supplies a bounded, compact source
which is strictly successful and locally `j`-efficient on a nonempty open
receiver set under a measurable critical-window vector policy.  Its
re-anchoring map is an affine injection—there is no folding—yet the weighted
pushforward left after exact rectangle-volume cancellation has density
`mu^N exp(o(N))`.  Thus bare re-anchoring transport is exponentially too
expensive even under local success and low local cost.  This remains only a
geometric/pushforward obstruction, not a global upper-bound counterexample.

## Fixed source and an open locally good policy patch

Use the r074 compact source at the receiver `x_0=(1,...,1)`:

\[
 F_N(y)=\tau h_{*,3/2,c_a}(S(y))\mathbf1_{W_N}(S(y))
        \mathbf1_{R_{c_a}(1)}(y),\qquad
 c_a=c_r+{(\log N)^{1/4}\over\sqrt N}.                  \tag{1}
\]

Here `tau=1+O(N^{-3})`, `S=-sum_i log|y_i|`, and r074 proves

\[
 A_{c_a}F_N(1)>1,\qquad j_{c_a}\Phi_N(F_N)(1)
 \le J_{\rm ref}N^{o(1)}.                                 \tag{2}
\]

Let `r_a=e^{-c_a}` and define, coordinatewise near `x_0`,

\[
 c_i(x)=c_a+\log x_i.                                     \tag{3}
\]

Then `e^{-c_i(x)}x_i=r_a`, so the actual rectangle in coordinate `i` is

\[
 R_{c_i(x)}(x_i)=[-r_a,x_i].                               \tag{4}
\]

At `x_0` this is exactly the r074 geometry.  The source is compact away
from the axes and its sum support is contained in the positive interior
window `W_N`.  Hence moving the upper faces in (4), the actions in (3), and
the conditional normalizers changes both the received mean and the full
vector local cost continuously.  One explicit representation is of the form

\[
 {1\over\prod_i x_i}\int_{R_{\mathbf c(x)}(x)}
 {\Phi_N(F_N(y))\over
 H_{\mathbf c(x)}(\sum_i\log(x_i/|y_i|))}\,dy,             \tag{5}
\]

on the compact support; here H_vector(s)=Gamma(N)^{-1} sum_J (s−sum_{i in J}c_i)_+^{N−1}, including the empty J. Its denominator stays positive there.  Dominated
convergence applies (moving face boundaries have zero measure).

At x_0 the full local cost equals the wide-window cost, because this fixed source is supported on W_N in total increment. Using the full local cost in the continuity argument avoids changing the window with the receiver. Since the first inequality in (2) is strict and the second is `N^{o(1)}`,
for every fixed `epsilon>0` and all sufficiently large `N` there is a
nonempty open set

\[
 E_N\subset(1-\min\{N^{-2},\rho/4\},
             1+\min\{N^{-2},\rho/4\})^N,\qquad
 \rho=\sqrt{\log N/N},                                   \tag{6}
\]

on which every receiver is strictly successful and

\[
 j_{\mathbf c(x)}\Phi_N(F_N)(x)\le J_{\rm ref}N^\epsilon.\tag{7}
\]

This is an existence-of-an-open-neighborhood statement; it does not claim
that the entire displayed cube has a uniform continuity modulus.  Since
`log x_i=O(N^{-2})` on the containing cube and
`c_a-c_r=o(rho)`, all actions in (3) lie in the critical window for large
`N`.  No maximizing property is asserted or needed: this is a valid
successful policy.

## Exact affine re-anchoring and the surviving weight

For the fixed reference `c_r`, set `r_r=e^{-c_r}`.  The re-anchor multiplier
and receiver map are

\[
 \lambda_i(x)={1+e^{-c_i(x)}\over1+r_r}
 ={x_i+r_a\over x_i(1+r_r)},\qquad
 T_i(x)=\lambda_i(x)x_i={x_i+r_a\over1+r_r}.              \tag{8}
\]

Thus `T` is a globally injective affine map with

\[
 \det DT=(1+r_r)^{-N}.                                    \tag{9}
\]

There is no folding or multiplicity in this example.  Nevertheless, after
the local rectangle-volume cancellation contributes the weight
`prod_i lambda_i`, the pushforward is, for every nonnegative probe `g`,

\[
 \begin{aligned}
 \int_{E_N}\prod_i\lambda_i(x)\,g(Tx)\,dx
 &=\int_{T(E_N)}D_N(Y)g(Y)\,dY,\\
 D_N(Tx)&={\prod_i\lambda_i(x)\over|\det DT|}
 =\prod_{i=1}^N\left(1+{r_a\over x_i}\right).           \tag{10}
 \end{aligned}
\]

The direction of the determinant in (10) is essential: `dx=|det
DT|^{-1}dY`.  Uniformly for `x in E_N`,

\[
 \log D_N(Tx)=N\log\mu+o(N),                              \tag{11}
\]

because `x_i=1+O(N^{-2})`, `c_a=mu+o(1)`, and
`1+e^{-mu}=mu`.  Hence the density in (10) is
`mu^N exp(o(N))`, exponentially large.  This verifies the claimed loss even
though the source in (1) is actual locally low-cost and the re-anchor map is
as regular as possible.

More generally, if a coordinatewise repair gave image-radius slopes
`q_i lambda_i`, the same cancellation leaves the product `prod_i q_i^{-1}`.
A mere positive lower bound on every `q_i` is therefore insufficient for a
subpolynomial full-vector loss: it permits fixed `q_i=q<1` and hence
`q^{-N}`.  One needs a substantially stronger collective condition, not
only nonfolding in each coordinate.

## Scope

Equation (10) holds for arbitrary probes `g`; it is not an assertion that a
particular `g` is a source modular.  The set `E_N` may be extremely small,
and this construction does not create enough receivers or control their
shared global source budget to refute a global strong-maximal upper bound.
It proves only that actual local efficiency, strict success, a critical
action window, and even affine injective re-anchoring do not by themselves
justify a bare weighted fixed-reference transport step.
