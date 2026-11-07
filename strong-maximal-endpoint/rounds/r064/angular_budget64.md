# Exact angular-budget push for deep noncentered rectangles

## Outcome

For every fixed receiver radial coordinate, averaging over the receiver's
within-block angles gives an **exact action-independent modular contraction**.
The factor from the uniform-Dirichlet push exactly cancels the ratio of the
source and receiver hyperbolic Jacobians.  This is stronger than a
`1+O(N^2/L)` freezing argument.

It does not by itself prove a global physical upper bound.  The resulting
radial budget is controlled only after averaging over receiver angles at a
*fixed* radial receiver.  Integrating it over all receiver radii produces
an exact upper-cone multiplicity factor \(\prod_jv_j\) at a source point.
Controlling that multiplicity is a separate packing/Carleson problem.

## 1. Uniform Dirichlet transport at fixed increment

Let \(\Delta_{b-1}=\{\xi_i\ge0:\sum_i\xi_i=1\}\), with uniform
probability measure \(\sigma_{b-1}\).  Fix `V>0`, a receiver angle
\(\xi\), and a nonnegative increment vector `z` with \(\sum_i z_i=s\).
Set
\[
 T_{V,z}(\xi)=\Theta=\frac{V\xi+z}{V+s}.
\]
Its image is the shifted simplex
\[
 \Delta_z=\{\eta\in\Delta_{b-1}:\eta_i\ge z_i/(V+s)\},
\]
and its pushforward has the exact density
\[
 \frac{d(T_{V,z})_\#\sigma_{b-1}}{d\sigma_{b-1}}(\eta)
 =\left(\frac{V+s}{V}\right)^{b-1}\mathbf1_{\Delta_z}(\eta).
 \tag{1}
\]
Indeed the affine map scales `(b-1)`-dimensional Lebesgue measure by
\((V/(V+s))^{b-1}\), and the normalization of uniform Dirichlet measure
is the same on both simplices.  Thus every nonnegative `G` satisfies
\[
 \int_{\Delta_{b-1}}G(T_{V,z}(\xi))d\sigma(\xi)
 =\left(\frac{V+s}{V}\right)^{b-1}
 \int_{\Delta_z}G(\eta)d\sigma(\eta).                    \tag{2}
\]

## 2. Exact conditional angular modular inequality

Use the block notation of `angular63.md`.  At receiver radial coordinate
\(V_j=L\alpha_j+u_j>0\), receiver angle \(\xi_j\), and common action `c`,
let \(S_j\) be the block increment sum and define the conditional radial
terminal
\[
 h_{u,\xi,c}(s)
 =E_c\left[F(u+s,\Theta)\mid S=s,\xi\right].             \tag{3}
\]
Here `F` is the action-independent raw terminal and `Theta` is the exact
transport in (8) of `angular63.md`.  The terminal in (3) generally depends
on `c` and `xi`; only the raw `F` does not.

Let \(d_j=b_j-1\), and write \(d\sigma(\xi)=\prod_jd\sigma_j(\xi_j)\).
Conditional Jensen, followed by (2) for each block and then averaging the
conditional increment vector `z`, gives for every `s>=0`
\[
 \begin{aligned}
 \int \Phi_N(h_{u,\xi,c}(s))d\sigma(\xi)
 &\le E_c\left[\int\Phi_N(F(u+s,\Theta))d\sigma(\xi)
 \mid S=s\right]\\
 &\le \prod_j\left(\frac{V_j+s_j}{V_j}\right)^{d_j}
 \int\Phi_N(F(u+s,\eta))d\sigma(\eta).                  \tag{4}
 \end{aligned}
\]
The last integral is over the full source-angle product simplex; replacing
the shifted simplices by it is the only inequality after Jensen.

Multiply (4) by the exact receiver radial weight
\(e^{-\sum u_j}\prod_jV_j^{d_j}\), then by \(e^{-\sum s_j}\), and
integrate `s`.  The Dirichlet factor cancels exactly:
\[
 \boxed{
 e^{-\sum u_j}\prod_jV_j^{d_j}
 \int_{\mathbb R_+^m}e^{-\sum s_j}
 \int\Phi_N(h_{u,\xi,c}(s))d\sigma(\xi)ds
 \le \mathcal O_L^+(F;u), }
 \tag{5}
\]
where the physical source-tail modular is
\[
 \mathcal O_L^+(F;u)
 =\int_{s\ge0}e^{-\sum_j(u_j+s_j)}
 \prod_j(V_j+s_j)^{d_j}
 \int\Phi_N(F(u+s,\eta))d\sigma(\eta)ds.                 \tag{6}
\]
This is the requested `e^{-s}ds` radial modular: the `e^{-s}` factor is
essential.  Formula (5) has no action-dependent constant and no deep-limit
error.  Restoring the common factor \(2^Ne^{-L}/\prod_j\Gamma(b_j)\) turns its
right side into the corresponding part of the exact Lebesgue--Orlicz cost
from the physical Jacobian.

## 3. What (5) does and does not permit

For a fixed `u`, suppose a scalar radial budget theorem says that every
bounded terminal `h` with
\(E_c h(S)>1\) has
\[
 \int e^{-\sum s_j}\Phi_N(h(s))ds\ge\mathfrak q_N^{-1}.
 \tag{7}
\]
Apply it only to the successful angular fibres
\[
 E_u(c)=\{\xi:E_c h_{u,\xi,c}(S)>1\}.
\]
Then (5) implies the exact fibrewise estimate
\[
 e^{-\sum u_j}\prod_jV_j^{d_j}\,\sigma(E_u(c))\,
 \mathfrak q_N^{-1}\le\mathcal O_L^+(F;u).                \tag{8}
\]
Thus a positive angular fraction of successful rectangles at the same
radial receiver has to pay actual source Orlicz cost above that receiver.
This is compatible with the `r026` conditional-density philosophy, but is
stronger in that the within-block angular Jacobians have not been discarded.

There are two crucial limitations.

1. Pointwise success at one angle `xi` has zero \(\sigma\)-measure.  It
does not give a lower bound for the integrated left side of (5).  Any use
of (8) needs an angular thickening, a continuity margin, or a different
way to charge isolated successful rectangles.

2. The modular budget in (5) is integrated over `xi`; it is not a uniform
budget for each terminal \(h_{u,\xi,c}\).  Hence a radial curvature or
sharing theorem cannot be applied fibre by fibre unless an additional
selection/efficiency argument supplies individual radial budgets.

## 4. Global radial integration: the exact obstruction

It is tempting to integrate (5) in all receiver radial coordinates
\(u_j>-L\alpha_j\).  The right side can be evaluated exactly.  Let
\[
 A(w)=\int\Phi_N(F(w,\eta))d\sigma(\eta),\qquad
 v_j=L\alpha_j+w_j.
\]
Fubini and the substitution \(w=u+s\) give
\[
 \begin{aligned}
 \int_{u>-L\alpha}\mathcal O_L^+(F;u)du
 &=\int_{w>-L\alpha}e^{-\sum w_j}\prod_jv_j^{d_j}A(w)
 \left|\{u:-L\alpha<u\le w\}\right|dw\\
 &=\int_{w>-L\alpha}e^{-\sum w_j}\prod_jv_j^{d_j+1}A(w)dw.
 \tag{9}
 \end{aligned}
\]
The full physical Orlicz cost contains only \(\prod_jv_j^{d_j}\), not the
extra \(\prod_jv_j\).  This is an exact upper-cone multiplicity: a source
point at radial coordinate `w` is charged by every receiver `u<=w`.
Therefore global radial integration does not contract the physical budget.

For example, a nonzero bounded source supported on 0<w_j<1 has the ratio of the integral in (9) to its physical source modular between ∏_j Lα_j and ∏_j(Lα_j+1). It therefore diverges like L^m. This only demonstrates failure of the unrestricted charging step: taking source amplitude below one produces no receiving halo and is not a counterexample to any endpoint constant.

For a restricted receiver set `U`, the extra factor is instead the
Carleson multiplicity
\[
 M_U(w)=|\{u\in U:-L\alpha<u\le w\}|.                    \tag{10}
\]
Any physical upper-bound proof using (5) must control `M_U` or allocate
source modular mass to receiver cones without reusing it.  This is exactly
where a genuine halo packing argument, rather than a radialization error
estimate, is required.

## 5. Relation to action sharing

For common `c`, the raw terminal in (3) remains fixed as `c` changes, but
both the conditional radial terminal and its successful-angle set change
through the kernel.  Formula (5) provides an action-independent **average
budget**, not the fixed-terminal hypothesis of the common-curvature lemma.
The angular transport identities in `angular63.md` quantify the missing
terms.  A possible route is to combine (8) with an efficiency theorem that
forces a positive angular measure of near-critical fibres whose conditional
terminals vary little with action, together with a Carleson bound for (10).
Neither property follows from (5) alone.

Thus the exact Dirichlet cancellation is real progress toward the physical
noncentered-rectangle representation, but the remaining barriers are
angular thickness/selection and radial cone packing, not a hidden
Jacobian error.

## Appendix: orthant-labelled raw terminals

Sections 1--5 used one common raw terminal on every orthant, hence described
a reflection-invariant base case.  No such symmetry is needed.  Let
\(\sigma\in\{\pm1\}^N\) be the receiver orthant label and let the physical
input be specified by an arbitrary bounded family
\(F_\sigma(w,\eta)\).  In (1)--(6), replace the factor `2^N` and a single
angular integral by the sum over all labels \(\sigma\); this is the exact
Lebesgue--Orlicz decomposition for a non-even input. Equivalently one may use the uniform probability 2^(−N) on the labels and keep the physical multiplier 2^N e^(−L)/∏Γ(b_j); an unnormalized label sum uses multiplier e^(−L)/∏Γ(b_j).

Retain the latent Bernoulli signs in (6).  Conditional on the increment
vector `z`, a coordinate with \(J_i=1\) crosses to the opposite source
orthant.  Thus the exact physical rectangle average is
\[
 \mathcal A_{L,c}F_\sigma(u,\xi)
 =E_c\left[F_{\sigma(-1)^J}(u+S,\Theta)\right].           \tag{11}
\]
For fixed `s`, define the conditional radial terminal using this joint
law, rather than only the marginal law of `z`:
\[
 h_{u,\sigma,\xi,c}(s)
 =E_c\left[F_{\sigma(-1)^J}(u+s,\Theta)\mid S=s,\xi\right].
 \tag{12}
\]

Conditional Jensen and the Dirichlet push (1)--(2) apply for each fixed
`J`.  After summing receiver labels, the map
\(\sigma\mapsto\sigma(-1)^J\) is a permutation of the orthants.  Therefore
the exact cancellation in (5) continues to hold with
\[
 \int\Phi_N(h_{u,\xi,c})d\sigma(\xi)
 \quad\hbox{replaced by}\quad
 \sum_\sigma\int\Phi_N(h_{u,\sigma,\xi,c})d\sigma(\xi),
\]
and with the right side replaced by the corresponding sum of the actual
source Orlicz tails of \(F_{\sigma'}\).  No reflection invariance, sign
averaging approximation, or label loss is involved.

For action policies depending measurably on (u,sigma,xi), assume every block size b_j≥2 and use the labelled
joint likelihood envelope of `envelope64.md`.  Its factor multiplies the
right side as specified in Section 6 of envelope64.md, while the source-label permutation above
still preserves the summed source budget.  The angular-thickness and radial
cone-packing limitations of Sections 3--5 remain unchanged.

For arbitrary nonnegative inputs with finite physical Φ_N modular, the same budget inequalities follow by replacing every F_σ by min(F_σ,M) and using monotone convergence. The measurable policy may be held fixed during this limit. This extends input amplitudes; it does not enlarge the fixed-anchor/common-within-block rectangle family or resolve receiver-cone allocation.
