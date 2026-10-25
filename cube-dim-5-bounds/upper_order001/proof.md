# Maximum-order conditioning: exact interface and certified single-point failure

Let β=3/4, and let S_j, j=1,…,16, be independent positive stable variables
with Laplace transform exp(−λ^β). Let v_j be the unnormalized body vectors
in {1}×{−1,1}⁴. Their frame matrix is Σ_j v_jv_jᵀ=16I. Put A=Σ_j S_jv_jv_jᵀ.
The Gaussian-mixture and time normalization are those explicitly derived
in ../upper_ridge001/proof.md; in particular

    κ=3Γ(7/4)/(2^(7/2)π^(5/2)).

Write F and f for the stable CDF and density. This law has a continuous
density, as also follows from its Hankel inversion below. Ties for the
maximum therefore have probability zero. If the unique maximum is in
direction v and has value z, the other 15 variables are independent with
law S conditioned on S≤z. Their conditional moments are

    m(z)=E[S | S≤z],
    n(z)=E[S⁻¹ | S≤z],       a(z)=1/n(z).

The probability element for this maximum winner is f(z)F(z)^15 dz.
It integrates to 1/16, and there are 16 disjoint winner events. No entropy
change or tilt is made. The factor 1/16 appearing in an integrated bin
weight must not be applied again to the sum of the 16 events.

## Conditional Jensen and time integration

Conditionally on this event, put B_c(z)=c(16I−vvᵀ)+zvvᵀ. Log-determinant
concavity and deterministic weighted least squares give respectively

    E[log det A | v,z] ≤ log det B_m(z),
    E[xᵀA⁻¹x | v,z] ≤ xᵀB_a(z)⁻¹x.

The second direction follows from expectation of a minimum ≤ minimum of
the expected quadratic cost. It is not matrix-inverse Jensen. All moments
needed here are finite for fixed z>0: the first moment is at most z, and
the inverse moment is at most E S⁻¹/F(z). An independent five-column
subframe controls the negative log determinant; the positive part is
bounded by the trace. Thus ordinary exponential Jensen is legitimate.

At stable time t, each variable scales by t^(4/3), while the identity of
the maximum is unchanged. Tonelli and the exact time Gamma integral give

    G(x) ≥ H_order(x)
      = κ Σ_v ∫₀∞ f(z)F(z)^15 φ_v(z) dz,
    φ_v(z)=(16m(z))⁻²(11m(z)+5z)⁻¹/²
             [q_v/(16a(z))+p_v/(11a(z)+5z)]⁻⁷/⁴,
    p_v=(v·x)²/5,       q_v=|x|²−p_v.

This is a sum over a disjoint partition of the original probability
space. It is valid for every x, with nonnegative extended integrals. The
present certificate concerns only x=(1,1,1,1,0). At that point the exact
groups are (p,q,count)=(16/5,4/5,2), (4/5,16/5,8), and (0,4,6). The local
code recomputes these counts from the 16 integer body vectors.

## Exact CDF and truncated-moment series

Define

    t_k(z)=π⁻¹(−1)^(k+1) Γ(βk) sin(πβk) z⁻βk/k!.

The stable tail is T(z)=1−F(z)=Σ_{k≥1}t_k(z), as proved by the Hankel
integral and absolutely dominated expansion in ../upper_ridge002/proof.md.
The truncated moments have the exact expressions

    M(z)=E[S 1_{S≤z}]=Σ_{k≥1} [βk/(1−βk)] z t_k(z),
    N(z)=E[S⁻¹ 1_{S≤z}]=n₀−Σ_{k≥1} [βk/(1+βk)] t_k(z)/z,
    n₀=Γ(1/β)/β,       m=M/F,       a=F/N.

There is no undetermined constant in the first expression. To see this,
the Laplace transform of M(z) is

    ∫₀∞ exp(−λz)M(z)dz = βλ^(β−2)exp(−λ^β).

Separate the leading transform βλ^(β−2), whose inverse is
βz^(1−β)/Γ(2−β). The remainder
βλ^(β−2)(exp(−λ^β)−1) has order λ^(2β−2)=λ^−1/2 near zero, so its small
Hankel arc vanishes. Expanding on the cut is absolutely justified by
exp(−zr) r^(β−2)(exp(r^β)−1), integrable at both zero and infinity.
The inverse terms are

    β(−1)^(k−1) z^(1−βk)/[(k−1)! Γ(2−βk)].

The reflection and Gamma recurrence identities turn these into the M
series above, including its leading k=1 term. For N, integrate the density
tail termwise on [z,∞) against s⁻¹, which is absolutely justified for
every fixed z>0. Subtract this tail moment from n₀. The identity for n₀
follows directly by integrating exp(−λ^β) over λ>0.

For β=3/4, k divisible by 4 contributes zero; each remaining residue
class is alternating. The tail-term absolute ratio at k+4 versus k is
at most z⁻³/(k+4). For the M multiplier, once k≥2,

    |[β(k+4)/(1−β(k+4))]/[βk/(1−βk)]|
      = (k+4)(βk−1)/(k(βk+2)) < 1.

For the N tail multiplier, its ratio is less than (k+4)/k. Thus the
eventual absolute ratios of the three series are bounded respectively by
z⁻³/(k+4), z⁻³/(k+4), and z⁻³/k. The calculation uses z≥1 and 64 terms
per residue. At the first omitted indices 257,258,259 and thereafter,
all these ratios are strictly below one. Each remainder is enclosed
between zero and its first omitted term with its actual sign. Arb at
256 bits propagates these intervals through F,M,N,m,a.

## Stieltjes bins and full complement control

Both m(z) and a(z) are increasing: truncating an increasing variable at
a larger upper threshold increases its conditional mean, while the
conditional mean of S⁻¹ decreases. Also m≥a, a≤a₀=1/n₀, and m≤z.
For the registered bins b_j=(65/64)^j, j=0,…,1024, the exact weight per
winner is

    ∫_[b_j,b_(j+1)] f(z)F(z)^15 dz
      = [F(b_(j+1))^16−F(b_j)^16]/16.

Each certified weight is positive. Use m and z at the upper endpoint to
lower-bound the determinant factor and use a and z at the lower endpoint
to upper-bound the inverse quadratic form. This gives a lower bound for
φ throughout each bin. Reversing these endpoint choices gives an upper
bound. No derivative of the stable CDF and no sampled minimum is needed.

For z<1, B_a(z)≤16zI and |x|²=4 imply Q≥1/(4z). Since m≥a and N(z)≤n₀,

    φ_v(z)≤ C a(z)⁻² z^(5/4)
            ≤ C n₀² z^(5/4)F(z)⁻²,
    C=4^(7/4)/(256√5).

Multiply by fF^15, use z^(5/4)≤1, integrate fF^13, and sum all 16 winners.
The full omitted low contribution to H_order is therefore at most

    κ·16 C n₀² F(1)^14/14.

Let U=b_1024. For z≥U, m(z)≥m(U), a(z)≤a₀ and every q_v>0 at this point.
Consequently,

    φ_v(z)≤(16m(U))⁻²(5U)⁻¹/²[q_v/(16a₀)]⁻⁷/⁴.

Multiplication by [1−F(U)^16]/16 bounds the whole upper tail per winner.
The lower calculation discards both omitted intervals by nonnegativity.

## Result and structural diagnosis

The complete formula satisfies the outward-rounded strict bounds

    0.0000420630244543289 < H_order(11110)
      < 0.0000431012868806624 < 0.000052 < 0.0001.

The lower contributions from the three projection groups are respectively
approximately 14.4108e−6, 17.3802e−6, and 10.2720e−6, with multiplicities
already included. The omitted low and high intervals contribute at most
1.17851e−7 and 3.59695e−15. Finer quadrature or recovering the omitted
positive tail cannot overcome this formula's certified ceiling.

This construction removes KL but still replaces all 15 residual variables
by separate first and inverse moments. Its large-z behavior makes this
loss explicit: M(z) is asymptotic to
βz^(1−β)/[(1−β)Γ(1−β)], while a(z) tends to a₀. Thus m grows like z^1/4,
whereas the harmonic parameter stays bounded. At 11110 every q_v≥4/5,
so φ_v(z) decays like z⁻¹; the maximum-density factor decays like z⁻⁷/⁴.
Rare large residual variables are not retained as separate geometric
directions by this single-maximum layer. This identifies the specific
moment replacement and its tail behavior; it does not prove that all
order-conditioned bounds or top-two refinements must fail.

The upper endpoint bounds H_order only, never actual G. No full spatial
cover, new C5 upper bound, top-two layer, or parameter grid is claimed.
