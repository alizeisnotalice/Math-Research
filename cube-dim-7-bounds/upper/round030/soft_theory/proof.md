# Round 030: soft eligibility kernels, a complete interface and a mixing obstruction

## Scope and outcome

The operator is unchanged: Lebesgue measure, centered axis-parallel cubes with a common radius in all coordinates, every positive radius, and arbitrary L1 input. This note concerns upper-bound selection mechanisms only. It proves no lower bound for the operator constant C_n.

There is a complete soft-kernel version of the finite selection/first-rejection argument. Its finite-buffer load is a joint supremum of explicit shell integrals. In the infinite-buffer limit, all nonincreasing kernels considered below have the **same** threshold mixing measure in coverage and load. Cauchy--Schwarz then proves that their resulting optimized budget cannot improve the existing hard-threshold family U_n. This includes finite mixtures, a power rolloff starting immediately at distance 1, and exponential rolloff.

The obstruction is to this particular joint payment interface, not to soft selection in general and not to C_n. A finite-buffer joint supremum involving thresholds below 1+2/a remains a genuine possible mechanism beyond the obstruction.

Source checkout: f65a5bf92911ebed4cf803071bd7693978e5e582. The required source files are `paper/main.tex` and `revision_v3/asymptotic/proof.md`; hashes are in `verification.json`. The local sparse checkout inherits an unrelated missing source blob and shows `D cube-lowdim-bounds/AGENTS.md`; this note does not modify or stage that item.

## 1. Kernels and the first-rejection coverage proof

Let μ be a probability measure on (1,infinity) and put

    κ(s) = ∫ 1_{s<q} μ(dq),       s >= 0.

Thus 0<=κ<=1, κ is nonincreasing, κ(s)=1 for s<=1, and κ tends to zero at infinity. Atoms and the strict inequality are allowed. Assume

    ∫ q^n μ(dq)<infinity,
    ∫_{1<q<2} log(1/(q-1)) μ(dq)<infinity.                 (M)

Fix a>2, τ>0, and use exactly the paper's increments I_j, volumes v_j and weights w_j=v_j/(2ar_j)^n. Replace eligibility by

    Σ_{i<j} w_i κ(||y_j-y_i||∞/(a r_i)) <= τ.             (S)

Choose a largest remaining eligible radius. Eligibility can only be lost as selected nonnegative terms are added, so selected radii are nonincreasing.

For an unselected Q_eta, stop just before the first selected radius smaller than r_eta, or at termination if no such selection occurs. The **prefix** preceding this event already violates (S), and all of its radii are at least r_eta. For x in Q_eta and every such prefix term,

    ||x-y_j||∞/(a r_j) < ||y_eta-y_j||∞/(a r_j)+1/a.

Therefore the transported kernel κ_cov(s)=κ(max(s-1/a,0)) has prefix sum >τ at x. Its normalized integral is, by Tonelli,

    (2ar_j)^(-n) ∫ κ_cov(||x-y_j||∞/(ar_j)) dx
      = ∫ (q+1/a)^n μ(dq) =: A_{n,a}(μ).                 (C)

Markov's inequality and disjointness of increments give

    |union of all input cubes| <= (1+A_{n,a}(μ)/τ) Σ_j v_j.

This is the same κ and μ as in (S); no separate coverage kernel is chosen. The argument never assumes that a rejected cube is still eligible immediately before termination.

## 2. The exact finite-buffer load surrogate

Let W(x)=Σ_j w_j 1_{aQ_j}(x). Pick the last among the smallest-radius contributors Q_k; normalize x=0, r_k=1. All other contributors precede k. The eligibility sum controls their κ-weighted part by τ, leaving

    W(0) <= τ+a^(-n)+∫ Σ_{j<k, 0 in aQ_j}
                        w_j 1_{||y_k-y_j||∞ >= aq r_j} μ(dq).       (L1)

For q>=2 this residual is zero, since ||y_k-y_j||∞<a(1+r_j)<=2ar_j. For 1<q<2 put

    t(q)=1/(q-1),  α=a/(a+1),  c_a(q)=(aq-1)/(a+1),
    m(z)=max(1,||z||∞),  u=y_k/(a+1),  z=y/(a+1).

For each actual residual increment, exactly the paper's triangle inequalities give

    1<=r_j<t(q),  ||z||∞<r_j,
    ||u-z||∞>c_a(q)r_j>=c_a(q)m(z).

Consequently, with H defined in the paper,

    D_{n,a}(μ)=a^(-n)+α^(-n) sup_{||u||∞<=α}
                 ∫_{1<q<2} H(u;c_a(q),t(q)) μ(dq),      (J)
    W(x)<=τ+D_{n,a}(μ).

The supremum stays outside the μ integral in (J). It is not legitimate to select a different u, different kernel, or different μ for coverage and load. The paper's final integration against |f|, finite covering and inner regularity now give

    C_n <= (1+A_{n,a}(μ)/τ)(τ+D_{n,a}(μ)),
    C_n <= (sqrt(A_{n,a}(μ))+sqrt(D_{n,a}(μ)))^2.          (F)

Formula (J) is computable even when c<=1. By coordinate reflection take 0<=u_i<=α. For rho>=1 define

    L_i(rho)=min(2rho, 2c rho, (1+c)rho-u_i),
    N_i(rho)=1_{|rho-u_i|<c rho}+1_{rho+u_i<c rho}.

All lengths above are positive because u_i<1 and c>0. Then

    H(u;c,t)=1-2^(-n) product_i L_i(1)
      +2^(-n) ∫_1^t rho^(-n)
        [2n(2rho)^(n-1)-Σ_i N_i(rho) product_{ell!=i} L_ell(rho)] d rho.   (SH)

Strict-boundary exceptions occur at finitely many rho and do not affect the integral. Formula (SH) follows directly by cube-face parameterization; it also holds for n=1 with empty products equal to one.

## 3. Common-corner regime and the limit for all kernels satisfying (M)

If q>1+2/a then c_a(q)>1, so the paper's coordinate monotonicity applies simultaneously to every such layer. The common maximizer is u=α(1,...,1). For q<2 its value is precisely J_a(t(q)) from the paper. In particular, for μ supported in (1+2/a,2),

    D_{n,a}(μ) = ∫ [a^(-n)+α^(-n) J_a(t(q))] μ(dq).       (CC)

This is equality for the surrogate (J), not an assertion that actual discrete configurations attain the surrogate.

For thresholds outside this common-corner range use the always-valid bound

    0<=H(u;c,t)<=1+n log t.

Define an explicit finite-a upper surrogate by taking J_a(t(q)) on 1+2/a<q<2, 1+n log t(q) on 1<q<=1+2/a, and zero on q>=2. Add a^(-n) and multiply its integral by α^(-n). It bounds (J). For every fixed q>1, the second range eventually disappears. Its bound is integrable by (M); for a>=2, α^(-n)<=(3/2)^n. Dominated convergence yields

    A_n(μ)=∫q^n μ(dq),
    B_n(μ)=∫ b_n(q) μ(dq),
    b_n(q)=B_n(1/(q-1)) for 1<q<2, and b_n(q)=0 for q>=2,
    C_n <= F_n(μ):=(sqrt(A_n(μ))+sqrt(B_n(μ)))^2.          (LIMIT)

This step handles kernels which start decreasing immediately above s=1; falsely requiring a positive plateau beyond 1 would miss them. It uses actual finite-a selections before passing to the parameter limit and therefore does not introduce an infinite-buffer selection algorithm.

## 4. Strict mixing obstruction

Write a_n(q)=q^n. Cauchy--Schwarz gives

    sqrt[(∫a_n dμ)(∫b_n dμ)] >= ∫sqrt(a_n b_n) dμ,
    F_n(μ) >= ∫(sqrt(a_n(q))+sqrt(b_n(q)))^2 μ(dq).        (CS)

For 1<q<2 the integrand is exactly the old hard-budget F_n(t(q))>=U_n. For q>=2 it equals q^n>=2^n>=U_n (the last comparison follows by t down to 1). Hence

    F_n(μ)>=U_n.                                        (OBSTRUCTION)

If the ratio a_n(q)/b_n(q) varies on positive μ mass, the first inequality in (CS) is strict; mass on q>=2 together with mass below 2 also makes it strict. The argument is an exact analytic inequality for the upper-budget family, not a numerical search and not a lower bound on C_n.

For finite a in the common-corner range, the same argument compares with the component finite-a hard budgets. The existing finite-a gap theorem then gives each component budget > its infinite-buffer hard budget >=U_n. Thus that finite-a soft family cannot improve U_n either.

## 5. Three checked classes

1. Finite step kernels: μ=Σ p_l delta_{q_l}, q_l>1, p_l>=0, Σp_l=1. Both costs are sums with these same p_l. The mixing obstruction applies without a search over step locations or weights.

2. Power rolloff touching 1: for 0<h<=1 and β>0,

       κ(s)=1                         for s<=1,
            1-((s-1)/h)^β            for 1<s<1+h,
            0                         for s>=1+h.

   Here q=1+hV with V density βv^(β-1) on (0,1), and

       A_n=Σ_{k=0}^n binom(n,k) β h^k/(β+k),
       B_n=∫_0^1 βv^(β-1) B_n(1/(hv)) dv.

   The k=0 term is one. The log moment equals log(1/h)+1/β, so (M) holds. A linear ramp is β=1. All β, not just checked numerical values, obey (OBSTRUCTION).

3. Exponential rolloff touching 1: κ(s)=1 for s<=1 and exp(-(s-1)/h) for s>1, h>0. Its threshold law has q=1+V, V exponential of mean h. Thus

       A_n=Σ_{k=0}^n binom(n,k) k! h^k,
       B_n=∫_0^1 exp(-v/h)/h B_n(1/v) dv.

   All moments and the required logarithmic moment are finite. Thresholds v>=1 pay coverage and contribute zero residual. Again (OBSTRUCTION) applies for every h.

The displayed one-dimensional quadratures are complete formulas, not invitations to splice the coverage of one class with another class's load.

## 6. Genuine remaining gap

For finite a, thresholds 1<q<=1+2/a have c<=1. The face count N_i then can increase with u_i; the common-corner proof is unavailable. Formula (SH) makes the problem explicit. It is possible in principle that maximizing positions differ by layer, so

    sup_u ∫H_q(u) μ(dq) < ∫sup_u H_q(u) μ(dq).

The mixing obstruction does not rule out exploiting this strict joint gain at finite a. A proof or counterexample requires evaluating the **joint** expression (J), including its a-dependent coverage cost. It cannot be inferred from independent hard-layer minima. The present note claims no such improvement.

There is an exact counterexample to extending the common-corner lemma to this range. In dimension one set a=10, q=11/10, hence α=c=10/11 and t=10. Formula (SH) gives

    H(0;10/11,10)=1/11+log(10),
    H(10/11;10/11,10)=1/2+(1/2)log(10).

Their difference is (1/2)log(10)-9/22>0; indeed log(10)>1 gives a difference >1/11. Thus the old corner is not even a maximizer for this layer. This is an auxiliary-integral counterexample, not an operator lower bound. It identifies why the finite-a noncommon regime cannot be silently included in (CC).

Further departures include changing the increments or using multilevel payments whose final coverage/load pair is not an affine mixture of the same hard pairs. These would need a new selection-to-L1 proof. Merely optimizing the three families above in the limit cannot produce a new mechanism, a coefficient below 1, or a seven-dimensional budget below U_7.

## 7. Verification and reproduction

Run `verify.py` with the existing local math-research-tools Python environment. It checks the power and exponential coverage formulas symbolically for n=1,...,8, the exact two-point Cauchy--Schwarz identity, and representative n=7 component/mixed budgets at 40 decimal digits. It also evaluates the two continuous soft formulas without optimizing them. These computations corroborate algebra and quadrature only; the all-n claims are proved above. `verification.json` records the actual results, versions and source hashes. No MATLAB session or large computation was started.
