# Conditional likelihood envelope for action-adaptive block angles

## Result

Let \(Z_c=E+cJ\), with \(E\sim\operatorname{Exp}(1)\) and
\(J\sim\operatorname{Bernoulli}((1+e^c)^{-1})\), and condition a block
of `b` independent increments on its sum `S=s`.  If
\(g_c^s\) denotes the conditional density on the simplex
\(\{z_i\ge0:\sum z_i=s\}\), then for b≥2, s>0 and 0<a<d<∞,
\[
 \mathcal Z_{b,s}([a,d])
 :=\int_{\sum z_i=s}\sup_{c\in[a,d]}g_c^s(z)\,dz_{-b}
 \le1+b\log\frac{1+e^{-a}}{1+e^{-d}}
 -\log\frac{f_{b,d}(s)}{f_{b,a}(s)}.                       \tag{1}
\]
At the critical lower-tail radius
\[
 s=b\mu-d_b+O(\log b),\qquad
 |c-\mu|\le K\sqrt{\log b/b},                              \tag{2}
\]
this gives
\[
 \mathcal Z_{b,s}=O(\sqrt{b\log b}),
 \qquad
 \log\mathcal Z_{b,s}\le\tfrac12\log b+	frac12\log\log b+O(1).
 \tag{3}
\]

Thus a measurable policy which chooses the block action as a function of the receiver angle can evade an action-independent conditional budget only through
this explicit polynomial likelihood-envelope factor, rather than through
the exponentially crude `2^b` factor.  This is a quantitative bridge to
the exact angular modular inequality of `angular_budget64.md`; it is not
yet a general upper bound for the strong maximal constant.

## 1. Exact conditional density and its one-sided action variation

The one-increment density is
\[
 q_c(z)=\frac{e^{-z}(1+\mathbf1_{\{z>c\}})}{1+e^{-c}},
 \qquad z>0.                                                \tag{4}
\]
On the sum simplex, write \(K_c(z)=\#\{i:z_i>c\}\).  The product density
is exactly
\[
 \prod_{i=1}^bq_c(z_i)=e^{-s}\frac{2^{K_c(z)}}{(1+e^{-c})^b}.
 \tag{5}
\]
Hence, relative to `dz_{-b}` on the simplex,
\[
 g_c^s(z)=\frac{e^{-s}2^{K_c(z)}}{(1+e^{-c})^bf_{b,c}(s)}. \tag{6}
\]

Away from the finitely many crossing values \(c=z_i\), the derivative of
the log density is independent of `z`:
\[
 \partial_c\log g_c^s(z)=bp_c-\partial_c\log f_{b,c}(s).
 \tag{7}
\]
Use the exact Gamma--Binomial conditional score identity
\[
 \partial_c\log f_{b,c}(s)
 =bp_c-(b-1)E_c\left[\frac K G\ \middle|\ S=s\right],   \tag{8}
\]
where `K` is the Bernoulli count and \(G=s-cK\) is the Gamma remainder.
Consequently the right side of (7) is nonnegative.  As `c` increases
through a crossing, \(K_c(z)\) drops by one, so \(g_c^s(z)\) has a
downward jump by a factor two.  Thus all upward variation in `c` is the
absolutely continuous variation in (7).

In the next integral, ∂_c g denotes its classical absolutely continuous derivative away from the moving hyperplanes, not its full distributional derivative (which includes negative atoms). For a positive function with only downward jumps,
\[
 \sup_{c\in[a,d]}g_c^s(z)
 \le g_a^s(z)+\int_a^d\partial_cg_c^s(z)\,dc.
\]
Integrating on the simplex and using \(\int g_c^s=1\) gives
\[
 \mathcal Z_{b,s}([a,d])
 \le1+\int_a^d\{bp_c-\partial_c\log f_{b,c}(s)\}\,dc,
\]
which is exactly (1).

## 2. Critical-window size

In the window (2), the ordinary compact-tilt score estimate from the
Gamma--Binomial kernel gives uniformly
\[
 \partial_c\log f_{b,c}(s)=O(\log b).                     \tag{9}
\]
The interval length is `O(sqrt(log b/b))`, whereas
\(bp_c=O(b)\).  Integrating (1) therefore gives
\[
 \mathcal Z_{b,s}([a,d])
 \le1+C\sqrt{b\log b}+O\!\left(\frac{(\log b)^{3/2}}{\sqrt b}\right),
\]
which proves (3).  This argument uses the conditional score sign (8), not
a pointwise comparison of the different conditional densities.

For reference, the envelope can be evaluated explicitly when `b=2` and
\(s>2d\).  Formula (6) has value \(1/[2(s-c)]\) on each outer interval
and \(1/(s-c)\) on the middle interval \([c,s-c]\).  Maximizing pointwise
in `c` and integrating gives
\[
 \boxed{
 \mathcal Z_{2,s}([a,d])
 =1-\frac{d-a}{s-d}+2\log\frac{s-a}{s-d}.}                \tag{10}
\]
This finite-dimensional identity checks both the downward-jump mechanism
and the orientation of the density ratio in (1).

### Uniform-in-s corollary

The same envelope has a global estimate.  Define
\[
 A_b(c,s)=bp_c-\partial_c\log f_{b,c}(s).
\]
The moving-jump flux identity is
\[
 A_b(c,s)=bp_c\frac{f_{b-1,c}(s-c)}{f_{b,c}(s)}\ge0,        \tag{10a}
\]
with a zero numerator when `s<c`.  It follows by differentiating all `b`
factors in the convolution: the denominator derivative contributes
\(p_cq_c\), while the jump moving through `z=c` contributes
\(-p_c\delta_c\).

For `s>=c`, the likelihood ratio
\(R_{b-1,c}(t)=e^tf_{b-1,c}(t)\) is nondecreasing.  Restricting the
convolution for \(f_{b,c}(s)\) to `0<=z<=c` gives
\[
 q_c(z)=\frac{e^{-z}}{1+e^{-c}},\qquad
 f_{b-1,c}(s-z)\ge e^{-(c-z)}f_{b-1,c}(s-c),
\]
and hence
\[
 f_{b,c}(s)\ge c p_c f_{b-1,c}(s-c),
 \qquad 0\le A_b(c,s)\le b/c.                             \tag{10b}
\]
For `s<c`, \(A_b=0\).  Inserting this into the variation proof of (1)
yields, for all `s>0` and `0<a<d`,
\[
 \boxed{\mathcal Z_{b,s}([a,d])\le1+b\log(d/a).}          \tag{10c}
\]
For b=1 the marginal, unlabelled angle is trivial and Z_{1,s}=1. This statement does not apply to the joint sign-labelled normalizer of Section 6: that normalizer is 3/2 when a<s<d. Its policy bound below requires every block size to be at least two.

## 3. Matching lower envelope mass in the critical window

The upper estimate has the correct order.  Fix a nontrivial critical window
\[
 [a_b,d_b]=\left[\mu-K\sqrt{\frac{\log b}{b}},
 \mu+K\sqrt{\frac{\log b}{b}}\right],\qquad K>0,          \tag{11}
\]
and let `s` satisfy (2).  Then
\[
 \mathcal Z_{b,s}([a_b,d_b])\asymp\sqrt{b\log b}.          \tag{12}
\]

To prove the lower bound, take `c<d` in this window and form the one-step
geometric-mean density.  With
\[
 \alpha=\sqrt2-1,\qquad \beta=2-\sqrt2,\qquad
 C_{c,d}=1+\alpha e^{-c}+\beta e^{-d},
\]
one has
\[
 \sqrt{q_c(z)q_d(z)}=
 \frac{C_{c,d}}{\sqrt{(1+e^{-c})(1+e^{-d})}}\,h_{c,d}(z), \tag{13}
\]
where \(h_{c,d}\) is the probability law of
\[
 E+W,\qquad
 P(W=0)=C_{c,d}^{-1},\quad
 P(W=c)=\alpha e^{-c}C_{c,d}^{-1},\quad
 P(W=d)=\beta e^{-d}C_{c,d}^{-1}.                           \tag{14}
\]
Thus its mean is exactly
\[
 M(c,d)=1+\frac{\alpha ce^{-c}+\beta de^{-d}}{C_{c,d}}.   \tag{15}
\]
At `c=d=mu`, both partial derivatives of `M` vanish: they are `alpha`
and `beta` times \(m'(\mu)=0\).  Hence throughout (11),
\[
 M(c,d)=\mu+O(\log b/b),\qquad
 \operatorname{Var}(E+W)=\mu+O(\sqrt{\log b/b}).           \tag{16}
\]
The uniform three-shift Fourier moderate-deviation proof therefore gives
\[
 \frac{f_{h_{c,d}^{*b}}(s)}{\sqrt{f_{b,c}(s)f_{b,d}(s)}}=1+o(1),
 \tag{17}
\]
uniformly in (11) and (2).  More precisely, all three laws have density φ(z_b)/√(μb) times [1+O((log b)^(3/2)/√b)] at this radius. Their tilted characteristic functions share the exponential factor (1+t²/(1−θ)²)^(-1/2); the remaining discrete three-point factor has modulus at most one and uniformly bounded cumulants. Thus the uniform Fourier proof of r063/refined63.md applies, while (16) moves the standardized critical exponent by O((log b)^(3/2)/√b).

The conditional Hellinger affinity consequently is
\[
 \begin{aligned}
 A_s(c,d)&:=\int\sqrt{g_c^s(z)g_d^s(z)}\,dz_{-b}\\
 &=\left\{\frac{C_{c,d}}
 {\sqrt{(1+e^{-c})(1+e^{-d})}}\right\}^{b}
 \frac{f_{h_{c,d}^{*b}}(s)}{\sqrt{f_{b,c}(s)f_{b,d}(s)}}.
 \end{aligned}                                             \tag{18}
\]
On a fixed compact neighborhood of `mu`, direct differentiation of the
braced factor at the diagonal gives a strictly negative one-sided slope.
Consequently there are fixed \(C,\kappa>0\) such that
\[
 A_s(c,d)\le C e^{-\kappa b|c-d|}                           \tag{19}
\]
for all sufficiently large `b` and all `c,d` in (11).  The important point
is the linear `b|c-d|` exponent, caused by the moving jump in (4), rather
than a Gaussian `b(c-d)^2` exponent.

Choose a grid \(c_k=a_b+kD/b\) in (11), where `D` is a sufficiently large
fixed constant.  Its cardinality is
\[
 K_b\asymp\sqrt{b\log b}.                                  \tag{20}
\]
For probability densities \(g_k=g_{c_k}^s\), Cauchy--Schwarz gives
\[
 K_b=\int\sum_kg_k
 \le\left(\int\max_kg_k\right)^{1/2}
 \left(\sum_{k,\ell}\int\sqrt{g_kg_\ell}\right)^{1/2}.  \tag{21}
\]
By (19), the double sum in (21) is at most a fixed constant times `K_b`.
Therefore \(\int\max_kg_k\ge cK_b\).  Since the full supremum over the
interval dominates the grid maximum, this proves the lower half of (12).

Combining (3) and (12),
\[
 \log\mathcal Z_{b,s}([a_b,d_b])
 =\tfrac12\log b+\tfrac12\log\log b+O(1).                 \tag{22}
\]

## 4. Tangent angular model and the meaning of the factor

The envelope mass is not merely a proof artifact.  Fix `s`, an interior
receiver angle \(\xi_0\), and a large radial coordinate `V`.  In tangent
coordinates
\[
 X=V(\xi-\xi_0),\qquad Y=(V+s)(\eta-\xi_0),
\]
the exact affine transport \(\eta=(V\xi+z)/(V+s)\) becomes
\[
 Y=X+z-s\xi_0.                                               \tag{23}
\]
Thus the conditional angular action is, up to the harmless fixed
translation \(s\xi_0\), convolution by \(g_c^s\).

For nonnegative \(\phi\in L^1\) on the tangent hyperplane, define
\[
 \mathcal M_s\phi(X)=\sup_{c\in[a,d]}(g_c^s*\phi)(X).
\]
Then
\[
 \|\mathcal M_s\|_{L^1\to L^1}=\mathcal Z_{b,s}([a,d]).  \tag{24}
\]
The upper bound follows from Fubini and
\(\sup_c\int g_c^s(z)\phi(X+z)dz\le
\int(\sup_cg_c^s(z))\phi(X+z)dz\).  For the reverse inequality, use a
positive approximate identity.  First restrict to a finite action grid,
where its convolution values converge in `L1` to the corresponding
densities and the integral of their maximum converges to the grid envelope;
then take a countable dense grid and monotone convergence.

In particular, let \(\phi\) be the indicator of an `epsilon`-cube in the
tangent hyperplane.  It has \(0\le\phi\le1\),
\(\Phi_N(\phi)=\phi\), and \(0\le\mathcal M_s\phi\le1\), so
\[
 \frac{\int\Phi_N(\mathcal M_s\phi)}{\int\Phi_N(\phi)}
 \longrightarrow\mathcal Z_{b,s}([a,d])\qquad(\epsilon\downarrow0).
 \tag{25}
\]
This is a low-efficiency packet: its source mass and its received average
both vanish with the packet volume.  It does not contradict any critical
average-one or near-efficient physical theorem.  A finite-depth physical
realization additionally requires a narrow radial band and a deep angular
packet limit; (23)--(25) do not by themselves provide that realization.

## 5. Connection with the angular modular inequality

In `angular_budget64.md`, after conditioning a block on `S=s`, an angular
source policy is tested against the conditional increment law.  If that
policy may depend on the action `c`, then for every nonnegative measurable
test `G` on the sum simplex,
\[
 \sup_{c\in[a,d]}E_c[G(Z)\mid S=s]
 \le\int G(z)\sup_{c\in[a,d]}g_c^s(z)dz_{-b}.             \tag{11}
\]
The dominating measure on the right has total mass \(\mathcal Z_{b,s}\).
Therefore the price for replacing an action-dependent conditional angular
choice by a single action-independent envelope is at most the polynomial
factor in (3), block by block.

Combining (11) with the exact uniform-Dirichlet push of
`angular_budget64.md` leaves its source/receiver Jacobian cancellation
intact; only the conditional angular measure acquires the factor
\(\mathcal Z_{b,s}\).  For several blocks the factors multiply.  Their
logarithms contribute at most
\[
 \frac12\sum_j\log b_j+\frac12\sum_j\log\log b_j+O(m).   \tag{12}
\]
This is much smaller than an exponential angular-policy loss, but it can
still be comparable with the half-logarithmic gains and losses in balanced
block arguments.  A future physical proof must decide whether this factor
is absorbable by an efficiency or packing estimate; (1)--(12) alone do not
settle that question.


## 6. Joint signs and arbitrary non-even source inputs

Assume b≥2, s>0 and 0<a<d. Retain the latent Bernoulli vector
J∈{0,1}^b. Relative to counting measure in J and dz_{-b}, the exact joint
conditional density is
\[
 \widetilde g_c^s(z,J)=
 \frac{e^{-s}}{(1+e^{-c})^b f_{b,c}(s)}
 \prod_{i:J_i=1}1_{\{z_i>c\}}.                              \tag{26}
\]
Summing J recovers (6), so this is normalized. For fixed (z,J), its
positive trajectory has log derivative A_b(c,s)≥0 until the first crossing
c=min_{J_i=1}z_i, where it is killed and remains zero. If J=0 it is never
killed. Because f_{b,c}(s) is continuous in c for b≥2, there is no positive
normalizer jump. Thus pointwise supremum equals the initial density plus
its entire absolutely continuous positive variation. Integration gives the
exact joint normalizer
\[
 \widetilde Z_{b,s}([a,d])
 =1+\int_a^d A_b(c,s)dc
 =1+b\log\frac{1+e^{-a}}{1+e^{-d}}
   -\log\frac{f_{b,d}(s)}{f_{b,a}(s)}
 \le1+b\log(d/a).                                         \tag{27}
\]
In the critical full window (11), the exact formula and the score bound
also give \(\widetilde Z_{b,s}=1+2Kp_\mu\sqrt{b\log b}+o(1)\).
The marginal envelope can be smaller, but Section 3 shows that even after
forgetting signs it still has the same logarithmic growth order.

For the physical input let σ∈{±1}^N be the receiver orthant and allow
arbitrary source terminals F_σ(w,η). The physical average evaluates
F_{σ(-1)^J}(u+S,Θ). Conditional Jensen, followed by (26)'s joint envelope,
removes any measurable choice c_j=c_j(u,σ,ξ) in [a_j,d_j]. For fixed J,
σ↦σ(-1)^J is a bijection, preserving the uniform orthant probability.
The Dirichlet push still applies for every fixed z,J. Therefore, if all
b_j≥2, the exact inequality (5) of angular_budget64.md holds for arbitrary
non-even inputs and these policies with the right side multiplied by
\[
 \prod_j[1+b_j\log(d_j/a_j)].                              \tag{28}
\]
The orthant average is included on both sides; using an unnormalized
orthant sum instead changes both sides by the same factor 2^N. No actual
source modular mass is discarded. This remains an angle-integrated,
fixed-radial-receiver bound, not a per-fibre critical budget or a global
halo-packing theorem.

For b=1, the marginal angle is still trivial, but the joint signs have an
exception: when a<s<d, J=0 conditional mass rises from 1/2 to 1 at c=s,
while J=1 is killed. Hence the joint envelope mass is 3/2. The moving-jump
normalizer argument in (27) cannot omit this upward jump. Singleton blocks
must either be excluded as above or be treated with their separate 3/2
factor.
