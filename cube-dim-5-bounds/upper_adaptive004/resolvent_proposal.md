# Uncomputed extension: optimal centered one-variable corrections

This is an analytic extension only; no new quadrature or numerical
candidate was run. Keep the same conditional setup and fixed null-space
projection P. Instead of b_l(S_l−m), allow an arbitrary deterministic
conditional function h_l(S_l) with E h_l=0, and use

    w=w0+P h(S).

The same independence argument eliminates all cross terms between
different l. For a fixed random index l, put p=P_ll and
B_l=Σ_(j≠l)P_jl²n_j. Its cost increment is

    2w0_l p E[h(S)/S]+E[(p²/S+B_l)h(S)²],    E h=0.

This is a one-variable weighted quadratic minimization, whose optimizer
is explicitly

    h(s)=−w0_l p(1−λ_l s)/(p²+B_l s),
    λ_l=R0_l/R1_l,
    R0_l=E[1/(p²+B_l S)],
    R1_l=E[S/(p²+B_l S)].

The choice of λ enforces E h=0. Completing the square, or substituting
the optimizer and using p²R0+B_lR1=1, gives the effective inverse
coefficient

    n_eff,l=R0_l/R1_l,       a_eff,l=R1_l/R0_l.

It is positive and satisfies a≤a_eff,l≤m. Indeed the two relevant
covariances have opposite monotonicities:
E[S/(p²+B_lS)]≤m E[1/(p²+B_lS)], while
E[S⁻¹]E[S/(p²+B_lS)]≥E[1/(p²+B_lS)]. The B_l=0 limit gives a_eff=m.
The exact optimum weakly improves every linear centered feature, and
the same matrix order and top-two tail bounds remain valid.

This formulation needs bounded conditional resolvent integrals R0,R1,
not a new inverse-moment series or a conditional variance estimate.
The functions are positive and monotone in s, so the existing certified
CDF bins could provide positive interval sums for a future test. All
optimized functions depend only on the fixed conditional parameters and
their own S_l; no dependence on another unretained realization is added.

A future pilot would still require separate registration and a careful
bound for the conditional low-s omission. No quantitative gain beyond
the completed linear-feature pilot is asserted here.
