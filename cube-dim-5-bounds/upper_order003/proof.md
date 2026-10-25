# Top-three conditional interface and the exact scope of the pilot

Let β=3/4, let the 16 body columns v_j∈{1}×{−1,1}⁴ satisfy VVᵀ=16I,
and let the original stable variables be iid with Laplace transform
exp(−λ^β). The Green normalization is

    κ=3Γ(7/4)/(2^(7/2)π^(5/2)).

## Disjoint events and conditional matrix

There are C(16,3)=560 unordered triples of directions. The event that a
given triple contains the three largest stable variables is disjoint from
all other such events, apart from null ties. Retain its three labeled
values z1,z2,z3, and put t=min(z1,z2,z3). The remaining13 variables are iid
conditioned on S≤t. If m(t)=E[S|S≤t] and a(t)=1/E[S⁻¹|S≤t], the conditional
matrix is

    B_c=c(16I−UUᵀ)+U diag(z1,z2,z3)Uᵀ,
    U=[v1 v2 v3].

The event density for this one unordered triple is

    f(z1)f(z2)f(z3)F(t)^13 dz1 dz2 dz3

on the **full positive cube**. There is no additional division by6.
The 13 remaining directions span R⁵: UUᵀ has operator norm at most its
trace15, hence16I−UUᵀ≥I. Finite first/inverse moments justify conditional
log-determinant Jensen and the deterministic least-squares inverse bound.
The exact time Gamma integral therefore gives

    G(x) ≥ H_top3(x)
      =κ Σ_{unordered triples} ∫∫∫ f(z1)f(z2)f(z3)F(t)^13 K_U(z) dz,
    K_U=(det B_m)^−1/2 (xᵀB_a⁻¹x)^−7/4.

This is a rigorous analytic lower-bound interface. No entropy penalty
appears. The numerical integral reported below is not a rigorous value
or enclosure of that interface.

## Gram-minor formulas, including degenerate cases

Write A=16c and h_i=z_i−c. For the physical endpoint choices used in this
interface, c is m(t) or a(t), so h_i≥0. For J⊆{1,2,3}, let

    g_J=det Gram(v_j:j∈J),     g_∅=1,
    q_J=det Gram(v_j:j∈J, x), q_∅=|x|².

The determinant lemma and a rank-one update by xxᵀ imply

    D_c=Σ_J A^(3−|J|) g_J ∏_(j∈J) h_j,
    N_c=Σ_J A^(3−|J|) q_J ∏_(j∈J) h_j,
    det B_c=A²D_c,        xᵀB_c⁻¹x=N_c/(A D_c).

All Gram minors are nonnegative. These identities do not invert the
retained Gram matrix. If some retained vectors or the augmented family
are dependent, the corresponding minors are simply zero; D_c remains
strictly positive because its empty-subset term is A³. Thus no artificial
singularity is introduced by Gram degeneration.

The exact integer enumeration in this task found that all560 actual
retained triples have rank3. Their full Gram determinants are48 for160
triples,64 for240 triples,and112 for160 triples. The augmented minors
may vanish and are retained as exact zeros in the implementation.

At x=11110, the coefficients are completely specified by the three
pair Gram determinants, the full Gram determinant, and the seven
nonempty augmented Gram minors. Canonicalizing this 11-component tuple
under all six label permutations gives28 classes. The integer
enumeration calculates every minor and multiplicity; these multiplicities
sum to560. Equal coefficient tuples imply equal kernels, without an
unproved symmetry reduction.

## Exact bin weights for a possible subsequent certificate

The following weights are proved here but no strict box integration has
yet been run. Let I_i=[b_i,b_(i+1)], let Δ_i=F(b_(i+1))−F(b_i), and write
Δ_i^(k)=F(b_(i+1))^k−F(b_i)^k. This notation is a difference of powers.

For one assignment of labeled variables to bins:

* If i>j>k, the weight is Δ_i Δ_j Δ_k^(14)/14.
* If i=j>k, the weight is Δ_i² Δ_k^(14)/14.
* If i>j=k, writing H=F(b_(j+1)), the weight is
  Δ_i·2[HΔ_j^(14)/14−Δ_j^(15)/15].
* If i=j=k, the weight is
  3[H²Δ_i^(14)/14−2HΔ_i^(15)/15+Δ_i^(16)/16].

The final expression follows by integrating the minimum CDF coordinate:
3∫_L^H(H-u)²u^13du. The factors2 and3 are internal to the repeated-bin
boxes. A sum over distinct label assignments has6,3,3,1 assignments in
the four cases above. These factors must not be applied twice.

## Numerical pilot: order-statistic normalization

For16 iid continuous variables, the CDF of the third-largest variable
has distribution Beta(14,3). Given this smallest retained CDF value u,
the other two CDF coordinates are independent uniforms on[u,1]; sorting
them gives the ordered top-three values. The ordered-value density is
3360 f(s)f(r)f(t)F(t)^13 on s≥r≥t, where3360=560×6.

Consequently one can sample this normalized ordered distribution and
evaluate

    H_top3=κ E_order [(1/3360) Σ_classes count·Σ_six_labels K_class].

The division by3360 belongs to this normalized expectation representation.
It is consistent with the full-cube event sum above and is not an extra
factor on that original sum.

The registered computation used one scrambled3-dimensional Sobol batch
of16384 points, seed20261026. It used a fixed4096-interval geometric
CDF/moment table from1/4 to2^40. Each table entry came from the rigorously
enclosed stable series in upper_order002; values were then interpolated
numerically. No point required extrapolation. This distinction matters:
rigorous table values do not make numerical interpolation or a Sobol
integral rigorous. There is no proved quadrature error or unsampled-tail
bound in this batch.

The numerical mean was7.816013872981878e−5. Eight consecutive block means
ranged from7.80906339793226e−5 to7.827768301216427e−5. These block extrema
are descriptive diagnostics, not confidence intervals or certified
bounds. The largest individual sampled contribution was3.80247e−4,
further illustrating why sampled values cannot be used as a global
integrand bound. Retained stable values ranged from0.6875 to7.5528e7,
and the smallest retained value ranged from0.6875 to232.85.

The estimate is about61.1% of the target128e−6. It does not presently
justify an expensive strict top-three integration as the leading route.
It also does not prove a top-three formula obstruction or an upper bound
on actual G. The rigorous top-three interface remains available; a future
strict certification must separately bound every box and all omitted tails.
