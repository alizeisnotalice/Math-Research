# A distinct next mechanism: feasible random least-squares coefficients

This is an analytic proposal, not a computed improvement. It targets the
inverse-quadratic loss without adding another order layer or changing the
probability measure. The formulas below hold conditionally on any fixed
set of retained stable values; the other variables are independent with
the same upper truncation t. Their mean and inverse mean are m,n, and
d=mn−1≥0. Retained variable j has deterministic value z_j and inverse
mean n_j=1/z_j; each random variable has n_j=n.

For the full16-column frame define

    P=I−VᵀV/16,       VP=0,       P_ll=11/16.

Take any deterministic feasible w0 with Vw0=x. For each unretained index
l let δ_l=S_l−m, and let δ_l=0 at retained indices. The random coefficients

    w=w0+P diag(b)δ

are feasible sample by sample, since Vw=x. Hence
xᵀA⁻¹x≤Σ_j w_j²/S_j for every realization. Expanding the expected cost,

    E cost = Σ_j n_j w0_j² + 2Σ_l c_l b_l + Σ_l H_l b_l²,
    c_l=−P_ll d w0_l,
    H_l=P_ll² m d + Var(S) Σ_(j≠l) P_jl² n_j.

There are no off-diagonal b_l b_k terms: if l≠k, at least one independent
centered factor retains zero mean even when the denominator is S_l or
S_k. Likewise only j=l contributes to the linear term. The identity
E[(S−m)²/S]=m²n−m=md gives the displayed diagonal entry.

Eliminating b exactly would decrease the effective inverse coefficient
at each random index from n to

    n_eff,l=n−P_ll²d²/H_l.

The remaining deterministic minimization over w0 then gives the valid
improved bound

    E xᵀA⁻¹x ≤ xᵀ[V diag(1/n_eff,j)Vᵀ]⁻¹x,

with n_eff,j=1/z_j at retained indices. The correction is at most d/m,
because H_l≥P_ll²md. Thus n_eff,l≥1/m>0, and the effective direct weights
lie between the old harmonic mean a=1/n and the arithmetic mean m.

An exact conditional second moment is not required for an initial bound.
Since0≤S≤t, Var(S)≤m(t−m). Replace H_l by the explicit upper bound

    Hbar_l=P_ll²md+m(t−m) Σ_(j≠l)P_jl²n_j.

Choosing b_l=−c_l/Hbar_l yields an actual expected reduction at least
c_l²/Hbar_l. Therefore n−P_ll²d²/Hbar_l is already a rigorous improved
inverse coefficient using only the currently certified m,n,t. For d=0,
set the correction to zero rather than dividing0/0.

The frame gives a convenient formula for the remaining sum:

    Σ_(j≠l)P_jl²n_j
      =[55n+Σ_(retained j)(v_j·v_l)²(1/z_j−n)]/256.

Although the retained correction is nonpositive, the original sum of
squares proves this whole quantity nonnegative. A numerical implementation
should preserve that certified nonnegative structure. All existing
log-determinant and time-integration steps then remain applicable with
the improved inverse matrix. The resulting claim is uniform in x for
each retained configuration.

Suggested next finite test, requiring a separate registration: reuse the
same single point and fixed top-three quadrature design, compare the old
inverse bound with this feasible-random-coefficient bound, and record the
change in the full formula. Only if the numerical gain is competitive
should a strict positive-box certificate be built. No such test has been
run in the present directory.
