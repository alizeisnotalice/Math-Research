# Optimal null-space columns and Stieltjes correction

Fix a top-two conditional configuration at x=11110. The two retained
values are deterministic z_j, the other14 variables are independent with
law S|S≤t, t the smaller retained value. Put m=E S, n=E S⁻¹, a=1/n, and
let n_j=n at random indices and n_j=1/z_j at retained indices. Let
B_a=Σ_j n_j⁻¹v_jv_jᵀ. All coefficients below depend only on these fixed
conditional parameters and x, never on an unretained realization.

## Independently selected feasible columns

For each random index l choose k_l∈ker V with k_(l,l)=1. For deterministic
w0 satisfying Vw0=x and conditional functions h_l with E h_l=0, set

    w=w0+Σ_random_l k_l h_l(S_l).

This is feasible for every realization. Different l may use different
null-space vectors; no common projection or mutual orthogonality is
required. The weighted least-squares cost expansion has no cross terms
between h_l,h_k when l≠k. If its denominator is S_l or S_k, the other
centered factor has zero expectation; otherwise both do. Its lth
increment is

    2w0_l E[h_l(S)/S]+E[(1/S+K_l)h_l(S)²],
    K_l=Σ_(j≠l)k_(j,l)² n_j.

Retained indices have no h correction and retain their deterministic
inverse coefficient1/z_j throughout this calculation.

## Exact best column and Schur guard

Minimizing K_l subject to V_(-l)k_(-l)=-v_l is a deterministic weighted
least-squares problem. Its positive definite Gram matrix is

    B_(-l)=B_a−a v_lv_lᵀ,

and its minimum is K_l*=v_lᵀB_(-l)⁻¹v_l. Explicit minimizing coefficients
are k_(j,l)=−n_j⁻¹v_jᵀB_(-l)⁻¹v_l for j≠l, with k_(l,l)=1. Because the
other15 body columns span R⁵, B_(-l) is positive definite. Sherman–Morrison
gives

    h_l=v_lᵀB_a⁻¹v_l,
    K_l*=h_l/(1−a h_l).

There is a uniform Schur denominator guard in this conditional setting.
The retained values satisfy z_j≥t≥a, hence B_a≥16aI. Therefore

    0<a h_l≤5/16,       1−a h_l≥11/16.

The normalized fixed-projection column P e_l/P_ll, with
P=I−VᵀV/16, is feasible and has

    K_l^fixed=Σ_(j≠l)P_jl²n_j/P_ll².

Consequently K_l*≤K_l^fixed. Both are strictly positive in the actual
body frame. The code evaluates the fixed value as an original positive
sum, and checks the Schur guard and K ordering at every sampled box.

## Optimal centered function for a fixed column

For K≥0 define the positive conditional Stieltjes moments

    A0=E[1/(1+KS)],       A1=E[S/(1+KS)].

The unique weighted-quadratic minimizer, up to null sets, is

    h*(s)=−w0_l(1−λs)/(1+Ks),       λ=A0/A1.

It is centered by construction. Since0<S≤t and E S⁻¹<∞, the weighted
Hilbert norm and linear functional are finite, so the minimization is
legitimate rather than a formal manipulation of divergent expressions.
The identity A0+KA1=1 yields the effective inverse coefficient

    n_eff,l=A0/A1,       a_eff,l=A1/A0.

At K=0 the formula gives a_eff=m. At K>0 all denominators and both
integrals are positive. Two opposite-monotonicity covariance inequalities
give a≤a_eff≤m. The function a_eff(K) is nonincreasing: after normalizing
the weight1/(1+KS), its derivative is

    −Cov_weighted(S,S/(1+KS))≤0.

Thus the optimal column yields a weakly larger a_eff than the fixed
column and a weakly smaller inverse quadratic. These comparisons hold
pointwise in every conditional configuration, independently of numerical
quadrature. The arbitrary-function class also contains the prior linear
feature, so it weakly improves that analytic bound.

After minimizing over the deterministic feasible w0,

    E xᵀA⁻¹x≤xᵀ[Σ_random a_eff,l v_lv_lᵀ+Σ_retained z_jv_jv_jᵀ]⁻¹x.

The effective matrix remains between B_a and B_m. The old determinant
Jensen bound, time Gamma constant, and full-support top-two tail bounds
are therefore unchanged. In particular q_eff≥1/(20t) at11110, and the
previous complete omitted-domain bounds still apply. All optimized
functions and columns use only conditional parameters; their displayed
dependence on their own S_l is the permitted pathwise feasible witness.

## Fixed numerical inner quadrature and limitations

The pilot reuses exactly the256 outer top-two boxes. Its two fixed inner
levels are32 and64 positive CDF-increment bins fromδ=1/4 to t, geometric
in s. Each interval is represented by its geometric midpoint. The common
low cell[0,δ] is represented by its exact probability F(δ) and exact
conditional mean m(δ). A0 and A1 are both positive weighted sums; A1 is
not evaluated by subtracting A0 from1. The conditional normalizer F(t)
cancels in their ratio.

This positive inner quadrature is numerical. The low-cell representation
is also an approximation, not an exact integral. Its error has a useful
analytic bound: since1/(1+Ks) is convex with second derivative at most
2K² on s≥0,

    0≤I0_low−F(δ)/(1+K m(δ))≤F(δ)K² m(δ)(δ−m(δ)).

The corresponding I1 error is the negative of this difference divided
by K for K>0. This fact is recorded for future certification but is not
propagated into a full pilot integral bound. Differences between32 and64
levels are convergence diagnostics only, not certified errors or
confidence intervals. The same applies to the outer midpoint integration.
