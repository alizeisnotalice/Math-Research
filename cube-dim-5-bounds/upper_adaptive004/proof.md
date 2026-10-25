# Feasible random least-squares correction under order conditioning

This proof precedes the registered numerical pilot. Let V be the full
5×16 body frame, VVᵀ=16I. Fix the retained directions and stable values
z_j as conditional parameters. Every unretained variable is independent
with law S conditioned on S≤t, where t is the smallest retained value.
Set m=E S, n=E S⁻¹, a=1/n, d=mn−1≥0. For retained j, n_j=1/z_j; for
unretained j, n_j=n. All expectations below are conditional on these
fixed parameters. There is no measure change or entropy penalty.

## Pathwise feasibility and legal dependence

Define P=I−VᵀV/16. Then VP=0, P²=P, and P_ll=11/16. For a deterministic
vector w0 satisfying Vw0=x, let δ_l=S_l−m at unretained indices and
δ_l=0 at retained indices. For deterministic b, define

    w=w0+P diag(b)δ.

Both w0 and b may depend on x and all fixed conditional parameters, but
not on the unretained realization. The vector w itself depends on the
realization through the displayed δ, and satisfies Vw=x for **every**
realization. Thus the deterministic weighted least-squares identity gives

    xᵀA⁻¹x=min_(Vu=x) Σ_j u_j²/S_j ≤ Σ_j w_j²/S_j

sample by sample. This use of random feasible coefficients is legitimate;
there is no filtration or causality requirement in this algebraic minimum.
The later optimization of b is over the conditional deterministic class
specified above, so independence and centered-moment calculations remain
valid. Finite bounded positive moments and finite inverse moments justify
the expected expansion and all finite sums.

## Diagonal expected quadratic form

For unretained l, E[(S_l−m)/S_l]=1−mn=−d. Every cross term with a different
denominator index retains an independent centered factor. Therefore

    E cost =Σ_j n_j w0_j²+2Σ_l c_l b_l+Σ_l H_l b_l²,
    c_l=−P_ll d w0_l,
    H_l=P_ll² m d+Var(S) Σ_(j≠l)P_jl² n_j.

For l≠k, each term in E[(S_l−m)(S_k−m)/S_j] is zero: if j=l or k, the
other centered factor has mean zero; otherwise both do. Hence there are
no off-diagonal b_l b_k terms. For the diagonal,
E[(S−m)²/S]=m²n−m=md, which proves the formula for H_l. The retained
indices are never optimization variables b_l; their deterministic costs
1/z_j remain in the sum defining H_l.

If H_l=0, then P_ll²md=0 and d=0, since P_ll>0 and m>0. Thus c_l=0.
Choose b_l=0 and define its correction to be zero, avoiding0/0. In the
nondegenerate case, minimizing b_l gives the exact cost reduction
P_ll²d²w0_l²/H_l.

## Certified variance upper bound and effective weights

Since0≤S≤t, Var(S)≤m(t−m). Define the conditional deterministic upper bound

    Hbar_l=P_ll²md+m(t−m) Σ_(j≠l)P_jl² n_j ≥ H_l.

If Hbar_l>0, choose b_l=−c_l/Hbar_l. Its true expected reduction is
c_l²(2/Hbar_l−H_l/Hbar_l²), at least c_l²/Hbar_l. If Hbar_l=0, then again
c_l=0 and the correction is defined as zero. Thus

    E xᵀA⁻¹x ≤ Σ_j n_eff,j w0_j²,
    n_eff,l=n−P_ll²d²/Hbar_l   (unretained l),
    n_eff,j=1/z_j             (retained j).

For d>0, Hbar_l≥P_ll²md implies that the correction is at most d/m.
Consequently

    1/m ≤ n_eff,l ≤ n,
    a ≤ a_eff,l=1/n_eff,l ≤ m.

These inequalities also hold in the d=0 case. Minimizing the deterministic
upper quadratic over Vw0=x is now legal and yields

    E xᵀA⁻¹x ≤ xᵀB_eff⁻¹x,
    B_eff=Σ_unretained a_eff,l v_lv_lᵀ+Σ_retained z_jv_jv_jᵀ.

The old harmonic matrix satisfies B_a≤B_eff≤B_m. Thus the improved inverse
bound is no larger than the old one, pointwise in every retained
configuration and for every x. Keeping the same log-determinant Jensen
bound produces a weakly larger valid Green lower formula, without
assuming monotonicity between different order layers.

The frame identity gives the optional expression

    Σ_(j≠l)P_jl²n_j
      =[55n+Σ_retained(v_j·v_l)²(1/z_j−n)]/256.

The implementation instead sums the original nonnegative terms, so it
does not subtract nearly equal quantities. Retained coefficients remain
exactly1/z_j and do not receive an adaptive correction.

## Time scaling and full-support tails for the top-two pilot

Under a common scale λ of all stable variables, m,t scale by λ, n by
λ⁻¹, d is invariant, Hbar scales by λ, and n_eff scales by λ⁻¹. Hence
B_eff scales by λ. The existing time Gamma integral and κ are unchanged.
For the top-two conditional layer at11110 the kernel is therefore

    K_ad=(det B_m)^−1/2(xᵀB_eff⁻¹x)^−7/4.

It is defined on the full positive retained-value square. Moreover
a_eff,l≤m(t)≤t. If s≥t are the retained values and v is the direction of
s, then B_eff≤t(16I−vvᵀ)+s vvᵀ. Because11110 has squared distance at
least4/5 from every body line,

    xᵀB_eff⁻¹x ≥1/(20t).

The unchanged determinant satisfies det B_m≥(16m(t))³·16st. Hence exactly
the same uniform bound as in upper_order002 holds:

    K_ad(s,t)≤C0 m(t)^−3/2 t^(5/4)s^−1/2,
    C0=20^(7/4)/(4·16^(3/2)).

The low-minimum and high-maximum integrals in upper_order002/proof.md
therefore remain valid, proving integrability and full-support tail
control for this new formula at11110. In particular, for the same
registered256 boxes their combined Green omission is at most the already
certified low bound7.226e−52 plus high bound6.763e−7. This controls the
omitted domain, but does not turn a midpoint calculation on the retained
boxes into a strict integral certificate.

## Pilot scope

Only the previously used256-box top-two partition will be evaluated.
The first pilot computes full5×5 matrices and numerical solves at one
geometric midpoint per box, with exact-series probability weights. It
will compare baseline and adaptive kernels through the same matrix path.
The actual48-element stabilizer of11110 has11 orbits on the120 unordered
pairs; that exact action, rather than an assumed Gram-only equivalence,
will reduce the adaptive computation. No global spatial grid, new order
layer, or strict quadrature claim is authorized by this pilot.
