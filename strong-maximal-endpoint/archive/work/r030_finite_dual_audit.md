# r030: finite real receiver budget — exact dual and certification interface

Date: 2026-11-06.  This note concerns the stated finite-receiver convex
problem, not a claim that a current floating-point grid is certified.

Let `b=1+eta`, let

`D=[a_1,infty) x [a_2,infty)`, `dnu=exp(-(t_1+t_2))dt`,

and consider

`P_A = inf { int_D Phi_n(H)dnu : 0<=H<=A,
                 int_D R_j H dnu >= b, j=1,...,M }`.

Here `R_j>=0` must belong to `L^1(nu)`.  In the intended shifted
Gamma--Bernoulli application, the choice `a_i=min_j u_{ji}` together with
support in `[0,infty)` gives `int_D R_j dnu=1`.  This normalization is
essential: the bare statement `A>1+eta` is Slater only after it (or the
more general condition `A int R_j dnu>b` for every `j`) is verified.

## Exact dual

For `lambda in R_+^M`, put `q_lambda(t)=sum_j lambda_j R_j(t)` and

`m_A(q)=min_{0<=h<=A} {Phi_n(h)-q h}`.

Then the Lagrange dual is

`D_A = sup_{lambda>=0} [ b sum_j lambda_j + int_D m_A(q_lambda)dnu ]`.

Weak-duality gap for any primal feasible `H` and any `lambda>=0` is exactly

`P(H)-d(lambda)
 = int_D[Phi_n(H)-q_lambda H-m_A(q_lambda)]dnu
   + sum_j lambda_j( int_D R_jH dnu-b ) >=0`.

For `n>=2`, the kink and cap selection are as follows.  Write

`alpha_A=Phi_n'(A)`,
`Phi_n'(h)=(1+log h)^(n-2)(n+log h)` for `h>1`.

The subgradient at one is `[1,n]`; hence a pointwise minimizer of
`Phi_n(h)-qh` is

* `h=0` for `q<1`;
* any `h in [0,1]` for `q=1`;
* `h=1` for `1<q<=n`;
* `h=(Phi_n')^{-1}(q)` for `n<q<alpha_A`;
* `h=A` for `q>=alpha_A`.

At the two boundary values the indicated closed selections are valid.  In
particular, treating the kink as a smooth inverse loses the full `[1,n]`
plateau and is incorrect.  The corresponding exact value is

`m_A(q)=0` for `q<=1`, `1-q` for `1<=q<=n`,
`Phi_n(h_q)-qh_q` in the soft interval, and `Phi_n(A)-qA` above the cap.

## Strong duality and attainment

Assume `R_j in L^1(nu)` and `A int R_j dnu>b` for all `j`.  Choose a
constant `a'` strictly between `b` and `A` in the normalized case (or
choose any constant with `a'int R_j>b`).  It is a strict feasible point.
The box `0<=H<=A` is weak-* compact in `L^infty(nu)` because `nu(D)<infty`;
the constraints are weak-* continuous and the integral of the bounded
convex integrand is weak-* lower semicontinuous.  Thus the primal attains.
Finite-constraint Slater separation gives zero gap and a multiplier.
Moreover the dual maximum attains: choosing the strict constant `a'` in
the Lagrangian gives

`d(lambda)<=Phi_n(a')nu(D)-sum_j lambda_j(a'int R_j dnu-b)`,

so its continuous concave objective is coercive on `R_+^M`.

This is an infinite-dimensional statement only under the listed topology
and integrability/strict-feasibility assumptions; a finite-grid Farkas
certificate alone does not establish it.

## Continuous certificates from finite positive quadrature

Positive node weights by themselves are *not* a continuous certificate.
They become one after adding rigorous cell and tail enclosures.  A practical
form is a compact rectangle

`B=[a_1,U_1]x[a_2,U_2]`,

partitioned into cells `C_k` with exact (or outward enclosed) masses
`w_k=nu(C_k)>0`.  For a tabulated, cellwise constant candidate `H_k`, use
interval arithmetic or derivative bounds to obtain

`R_j(t) >= Rlow_jk` and `R_j(t)<=Rup_jk` on each cell.

If the candidate is zero off `B`, then

`I_j=int R_jH dnu >= sum_k w_k Rlow_jk H_k`;

therefore `sum_k w_k Rlow_jk H_k >= b` certifies continuous feasibility.
Its cost is exactly `sum_k w_k Phi_n(H_k)` (and is an upper bound on the
continuous optimum).  More generally, if a positive quadrature produces
numbers `Ihat_j,Jhat` and verified errors
`|I_j-Ihat_j|<=e_j`, `|J-Jhat|<=e_Phi`, require
`Ihat_j-e_j>=b` and report `J<=Jhat+e_Phi`.

For a dual lower certificate, enclose
`q_lambda(C_k)` from above by `qbar_k`.  Since `m_A` is nonincreasing,

`int_{C_k}m_A(q_lambda)dnu >= w_k m_A(qbar_k)`.

On `D\B`, use `m_A(q)>=-Aq`; if
`tau_j>=int_{D\B}R_j dnu` is a rigorous Gamma--binomial tail bound, then

`d(lambda)>=b sum lambda_j
 +sum_k w_k m_A(qbar_k)-A sum_j lambda_j tau_j`.

This is a certified lower bound on the *continuous* primal optimum.  The
difference between the certified primal upper and dual lower bounds is a
continuous optimality-gap certificate.  Equivalent quadrature-error form:
if `dhat` approximates the dual integral with verified error `e_D`, then
`d(lambda)>=dhat-e_D`.

## Monotone-profile tail option

If the actual witness additionally must be coordinatewise nondecreasing,
do **not** set it to zero outside `B`: that destroys monotonicity.  Instead
set `H=A` on `D\B` and zero below the lower corner.  If the values inside
`B` are coordinatewise nondecreasing and lie in `[0,A]`, this extension is
still coordinatewise nondecreasing.  Then:

* use the same box lower bounds alone to certify averages;
* add the exact/bounded cost tail `Phi_n(A)nu(D\B)` to the primal upper;
* retain the dual tail lower bound `-A sum lambda_j tau_j` above.

For this monotone extension the source tail must be included in cost; it
cannot be silently discarded.  The required `tau_j` is computable by
directed Gamma--binomial tail integration/interval arithmetic, not a
Monte-Carlo estimate.

## Numerical cross-check

Local MATLAB evaluated the `n=4,A=2` rule: the left/right slopes at one are
`1` and `4`, `Phi'(2)=13.4540673605`; sample minimizers for
`q=(0.8,1,2,4,6,14.454...)` were respectively
`(0,1,1,1,1.216312...,2)`, agreeing with the four-regime rule.  This is a
sanity check only, not a proof or a certificate for any receiver family.
