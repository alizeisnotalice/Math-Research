# Audit — conditional angular-envelope normalizer

## Verdict

**PASS, with the sharpened asymptotic.** Uniformly at
the critical block level and in the stated near-\(\mu\) parameter window,
\[
 Z_b(s)\asymp\sqrt{b\log b},\qquad
 \log Z_b(s)=\tfrac12\log b+\tfrac12\log\log b+O(1),
 \qquad Z_b(s)=\int_{\sum z_i=s}\sup_c g_c^s(z)\,dz.
\]
The result proves an unavoidable polynomial normalizer for unrestricted
conditional angular scalarization.  It does **not** bound \(Q_N\), \(C_N\),
or construct a near-efficient physical counterexample.

## Upper bound

Away from the downward jump hyperplanes \(c=z_i\), differentiation of the
conditional density gives
\[
 \partial_c\log g_c^s=bp_c-\partial_c\log f_{b,c}(s).
\]
The exact Gamma--Bernoulli score identity is
\[
 bp_c-\partial_c\log f_{b,c}(s)
 =(b-1)E_c[K/G\mid G+cK=s]\ge0.
\]
As \(c\) increases, all jump contributions to \(g_c^s\) are downward.
Therefore pointwise positive variation, followed by integration over the
simplex, gives
\[
 Z_b(s)\le1+\int_a^d
 \{bp_c-\partial_c\log f_{b,c}(s)\}\,dc
 =1+b\log{1+e^{-a}\over1+e^{-d}}
 -\log{f_{b,d}(s)\over f_{b,a}(s)}.
\]
At \(s=b\mu-d_b+O(\log b)\) and
\(c-\mu=O(\sqrt{\log b/b})\), the uniform tilted local theorem gives
\(f_{b,c}(s)\asymp(\log b)/b\), with endpoint log-ratio \(o(1)\).
The first term is \(O(\sqrt{b\log b})\). Hence
\[
 Z_b(s)\le C\sqrt{b\log b},\qquad
 \log Z_b(s)\le\tfrac12\log b+\tfrac12\log\log b+O(1).
\]
The sign and the use of positive variation are essential; differentiating
through the moving jumps without retaining their negative direction would
be invalid.

## Lower bound

For \(c<d\), direct integration gives the stated one-step affinity
\[
 A(c,d)=
 {1+(\sqrt2-1)e^{-c}+(2-\sqrt2)e^{-d}
  \over\sqrt{(1+e^{-c})(1+e^{-d})}}.
\]
Its one-sided logarithmic derivative at the diagonal is strictly negative,
and compactness of the near-\(\mu\) window yields
\(A(c,d)\le e^{-\kappa(d-c)}\) for a fixed \(\kappa>0\).

The normalized geometric-mean increment density is an exponential mixture
with shifts \(0,c,d\):
\[
 q_h(z)={e^{-z}\{1+\alpha1_{z>c}+\beta1_{z>d}\}\over T(c,d)},\quad
 \alpha=\sqrt2-1,\quad\beta=2-\sqrt2.
\]
For \(c\le d\), its mean is
\[
 M(c,d)=1+{\alpha ce^{-c}+\beta de^{-d}\over T(c,d)}.
\]
At \((\mu,\mu)\), both one-sided partial derivatives vanish:
\(\partial_cM=\alpha m'(\mu)=0\) and
\(\partial_dM=\beta m'(\mu)=0\). Thus throughout the window
\[
 M(c,d)=\mu+O(\rho^2),\qquad \operatorname{Var}q_h=\mu+O(\rho).
\]
The uniform exponential-plus-three-shift Fourier moderate-deviation proof
therefore gives
\[
 f_{h,b}(s)={\varphi(z_b)\over\sqrt{\mu b}}
 \{1+O((\log b)^{3/2}/\sqrt b)\},
\]
the same asymptotic as \(f_{b,c}(s)\) and \(f_{b,d}(s)\). Hence
\[
 \int\sqrt{g_c^s g_d^s}
 =A(c,d)^b(1+o(1))
 \le C e^{-\kappa b|c-d|}. \tag{*}
\]

Choose a grid with spacing \(D/b\), where \(D\) is fixed so large that
\(2C e^{-\kappa D}/(1-e^{-\kappa D})<1/2\). It has
\(K\asymp b\rho\asymp\sqrt{b\log b}\) points. Under each grid law, the
maximum-likelihood error is bounded by the sum of the affinities in (*),
hence is at most \(1/2+o(1)\). Therefore
\[
 Z_b(s)\ge\tfrac12K,
\]
which matches the upper order.

## Scope

This only rules out a cost-free replacement of arbitrary conditional
angular kernels by a common scalar terminal.  Turning its
\(\sqrt b\)-scale normalizer into any statement about physical endpoint
efficiency would require a separate link from conditional scalar budgets
to actual angular/radial packing and score accounting.

## Global-in-s policy envelope

The proposed all-level extension is valid.  For a parameter interval
\([a,d]\Subset(0,\infty)\), the positive conditional-density flux is
\[
 a_b(c,s):=bp_c-\partial_c\log f_{b,c}(s)
 =bp_c\,{f_{b-1,c}(s-c)\over f_{b,c}(s)}.
\]
For \(s<c\) it is zero. For \(s\ge c\), write
\(R_{b-1,c}(t)=e^t f_{b-1,c}(t)\). This is nondecreasing (equivalently,
use the positive shifted-Gamma polynomial representation). For
\(0\le z\le c\),
\[
 q_c(z)f_{b-1,c}(s-z)
 ={e^{-s}\over1+e^{-c}}R_{b-1,c}(s-z)
 \ge p_c f_{b-1,c}(s-c).
\]
The convolution formula therefore gives
\[
 f_{b,c}(s)\ge c\,p_c f_{b-1,c}(s-c),\qquad
 0\le a_b(c,s)\le b/c.
\]
(The degenerate \(b=1\) conditional simplex is handled separately and has
normalizer one.) Positive variation now yields, uniformly for **every**
\(s\),
\[
 Z_{b,s}([a,d])
 \le1+\int_a^d a_b(c,s)\,dc
 \le1+b\log(d/a). \tag{G}
\]

For an arbitrary measurable angular policy with each block action in
\([a_j,d_j]\), Jensen first moves the Orlicz function inside the
conditional increment law. Pointwise domination by the supremum
conditional density, then the fixed-increment Dirichlet push, gives the
exact actual source-tail modular bound with multiplicative factor
\[
 \prod_j\{1+b_j\log(d_j/a_j)\}.
\]
This is policy-uniform and requires no critical-sum localization, provided
every block has \(b_j\ge2\). It is an
envelope inequality, not an action-independent terminal or common-curvature
reduction. In particular it supplies neither a \(Q_N/C_N\) estimate nor a
solution to the receiver-cone multiplicity/Carleson obstruction.

## Joined physical-budget statement

The physical restoration in angular_budget64 is now exact.  With
\(v_j=L\alpha_j+w_j\), the full hyperbolic Jacobian contributes the common
factor
\[
 {2^N e^{-L}\over\prod_j\Gamma(b_j)}.
\]
Thus multiplying the radial source-tail expression in its (6), and the
policy-envelope version of its (5), by this factor restores the
corresponding exact Lebesgue--Orlicz pieces.  No additional
\(\alpha_j\)-factor or asymptotic correction is missing: the powers are
already \((L\alpha_j+w_j)^{b_j-1}\).

Combining the all-level envelope (G), Jensen, and the fixed-increment
Dirichlet push gives, for any measurable within-block policy satisfying
\(c_j(\xi)\in[a_j,d_j]\),
\[
 \text{policy receiver modular at fixed }u
 \ \le\
 \prod_j\{1+b_j\log(d_j/a_j)\}\,
 \mathcal O_L^+(F;u),
\]
with the same physical common factor on both sides.  This is an exact
source-tail budget comparison up to the displayed likelihood-envelope
factor, valid without critical-sum localization.

Finally, the example supported on \(0<w_j<1\) is correctly only a method
obstacle.  The ratio of unrestricted cone charging (angular_budget64 (9))
to physical source cost lies between
\(\prod_jL\alpha_j\) and \(\prod_j(L\alpha_j+1)\), while a source of height
below one has no level-one halo.  It proves that unrestricted radial
charging has multiplicity \(L^m\); it is not an endpoint counterexample.

## Orthant-labelled extension

**PASS for arbitrary (not necessarily even) physical input, with
\(b_j\ge2\) explicitly retained.**  Keep the latent labels \(J_i\) rather
than marginalizing them.  Their joint conditional density is
\[
 \widetilde g_c^s(z,J)=
 {e^{-s}\over(1+e^{-c})^b f_{b,c}(s)}
 \prod_{i:J_i=1}{\bf1}_{\{z_i>c\}}.
\]
It is normalized: summing \(J\) produces
\(e^{-s}2^{K_c(z)}/(1+e^{-c})^b\), whose simplex integral is
\(f_{b,c}(s)\).  For a fixed \((z,J)\), the density has the same
nonnegative absolutely-continuous logarithmic derivative \(a_b(c,s)\) as
before.  If \(J\ne0\), it is supported only below
\(\min_{J_i=1}z_i\), so its sole support-boundary jump as \(c\) increases
is downward; if \(J=0\), there is no such jump.  Consequently the joint
positive-variation argument is exact and yields the same normalizer
\[
 \widetilde Z_{b,s}=1+\int_a^d a_b(c,s)\,dc
 =1+b\log{1+e^{-a}\over1+e^{-d}}
   -\log{f_{b,d}(s)\over f_{b,a}(s)}
 \le1+b\log(d/a).
\]

The restriction \(b\ge2\) is substantive in this labelled version.  For
\(b=1\), when \(s\in(a,d)\), the \(J=0\) conditional mass jumps upward
from \(1/2\) to \(1\) at \(c=s\); the joint envelope is \(3/2\), not the
scalar-simplex value \(1\).  Thus the displayed block product must not be
claimed for arbitrary non-even inputs with a singleton block unless it is
replaced by an explicit \(b=1\) factor.  The intended macroscopic-block
application, with all \(b_j\ge2\), is unaffected.

For the physical signs, the exact source label is
\(\sigma'=\sigma(-1)^J\).  For every fixed \(J\), this is a permutation of
the receiver orthants.  Therefore summing over \(\sigma\) after conditional
Jensen and the fixed-increment Dirichlet push returns precisely the sum of
the actual source Orlicz tails \(F_{\sigma'}\).  The argument remains valid
when the policy is merely measurable in \((\sigma,\xi)\), since the joint
envelope is uniform in that choice.  This closes the former reflection/
evenness gap; it does not alter the separate angular-thickness and
cone-packing limitations.

Final integration check: the marginal b=1 sentence is explicitly distinguished from the joint b=1 exception, and the policy theorem with arbitrary orthant-labelled inputs assumes every b_j≥2. Section 6 now contains the displayed joint density, exact normalizer and policy reduction. Nonnegative unbounded inputs follow by monotone truncation at fixed measurable policy; this does not enlarge the fixed-anchor/shared-block operator family.
