# 11O: tensor Chebyshev/convex-chord integral certificate

Fixed n3 body+face frame, tau=1/3, eps=9/5000. Old11H is closed incomplete. This new representation computes true H_eps=integral_eps^infty p_t; no cost/min/tau search. Four representative patches actually passed, not yet a full I certificate. Source identities and whole-face short-time lemma require independent acceptance. Current U3/U4 unchanged by this pilot.

## Representation and all degeneracies

The locked96cone proof gives48 antipodal pairs. For each, rays r_j, positive k_j=psi(r_j), b_j=r_j dot(1,y,z)/k_j and weight w=2/prod(k_j). The exact point formula is H_eps=(8pi^3)^-1 sum w f_eps[b0,b1,b2], f_eps(b)=log(eps^2+b^2)/2. It follows by integrating the three resolvents; repeated knots are continuous Hermite limits f[a,a,b]=(f(b)-f(a)-f'(a)(b-a))/(b-a)^2 and f[a,a,a]=f''(a)/2. Exact Fraction equality chooses these cases. Zero projections are ordinary points because eps>0; all signed terms remain. Code independently constructs pairs from locked cb11_cones.json, rather than importing a mutable point driver. Identical absolute f arguments share a cache, using f evenness exactly; no approximation is introduced. Formula is the same already used by minimum's regularized evaluator, disclosed shared theory/arithmetic, not an independent derivation or audit.

## Uniform analyticity and node approximation

Each ray exactly satisfies k>=3(|r_y|+|r_z|). Cone linearity and the triangle inequality imply psi(xi)>=3(|xi_y|+|xi_z|) globally. For complex y,z with max imaginary part <=s<3eps, absolute Fourier integration gives |H_eps|<=M=15/[1088pi^3(eps-s/3)^2], using p_t(0)=15/(544pi^3t^3). This is a whole polydomain estimate, not a sum of natural interval cone ranges.

Every patch has halfwidth h=1/256. Degree32 tensor interpolation uses33 roots in each direction and rho=5/2. Its ellipse height s=h(rho-rho^-1)/2; eps-s/3>0 exactly. At384bit, M<2373.630728. The one-variable degree32 roots interpolation error is <=4M rho^-32/(rho-1): coefficients <=2M rho^-k and roots alias all Chebyshev modes to a signed mode or zero, so the tail difference is <=2 times the coefficient tail. Cardinal DCT expression gives Lebesgue norm<=2*32+1. Applying I-IyIz=(I-Iy)+Iy(I-Iz) yields

    Epoly=(2*32+2)*4M rho^-32/(rho-1)
          <7.706293505929e-8.

This bound holds for the actual exact-root interpolant, including both coordinates and all true G singular lines. It is not a smoothness assumption on G itself.

Positive roots are floored to exact dyadics with128 fractional bits after Arb proves a unique floor; negative nodes are exact reflections and the middle root exactly0. Physical coordinate displacement <h*2^-128. Cauchy radius R=eps at real points gives each coordinate derivative <=D=15/[1088pi^3(eps-R/3)^2 R]. Integrating along real coordinate segments proves total value displacement <=2Dh2^-128, actual upper <3.938477182864e-36. Every point enclosure is inflated by this amount before any matrix multiplication, so node perturbation and transform amplification are fully propagated.

## Exact polynomial mean and whole-patch chord

DCT matrix A maps true-node values to tensor Chebyshev coefficients. Exact integer recurrence builds T_k(2u-1)=sum_l p_kl u^l. Degree32 Bernstein transform is B_jk=sum_{l<=min(j,k)}p_kl binom(j,l)/binom(32,l), all Fraction. Thus C=B A V A^T B^T encloses the coefficients of the single actual exact-root interpolant P, despite correlated node balls. Since Bernstein bases are positive and sum1, their min/max enclose P on the entire closed square. Taking outward dyadic endpoints gives

    L=max(0,min_coeff_lower-Epoly), U=max_coeff_upper+Epoly.

True H>=0 is supplied by the heat semigroup. Exact Chebyshev means on the normalized patch are mu0=1, mu_odd=0, mu_even=-1/(k^2-1). Therefore mean(P)=mu A V A^T mu^T, computed without numerical quadrature. Mean(H)<=m=min(U,meanP_upper+Epoly). Actual checks verify m>=L; if this were to fail the program stops rather than raising the mean silently.

For q=3/2, convexity on [L,U] bounds H^q by its chord. The average obeys

    mean H^q <= L^(3/2) + S(m-L),
    S=(U+sqrt(UL)+L)/(sqrt(U)+sqrt(L)).

For U>L this is exactly (U^(3/2)-L^(3/2))/(U-L), without small denominator cancellation. For U=L>0, m=L and the value reduces to L^(3/2). U=0 would give zero; the current positive heat kernel gives U>0, checked explicitly. L=0 is valid. Monotonicity in m permits its upper bound. The final area-weighted upper endpoint is an exact dyadic, also saved as a Fraction for lossless summation.

## Triangle multiplicity and final short-time budget

Indices0<=j<=i<128 denote closed squares. Coordinate-swap symmetry pairs off-diagonal squares, weight2; diagonal squares are whole squares and weight1. Exact weighted area sum is128^2/128^2=1. They cover the full positive quadrant up to shared edges of Lebesgue measure0, with no missing triangle half or doubled diagonal. This integral step specifically uses constant Lebesgue density and the process's reflection/permutation invariance.

H_eps is finite and analytic across all true G lines. There are no removed hats. The difference T_eps=G-H_eps is nonnegative with whole-face Lq bound from main11M:

    Etime=2^-4/3*(3/4)*[(256/3)/pi^3]^(1/3)
          *[(16/3)/pi]^(2/3)*eps^(4/3).

The analysis uses the same density's sup and exact plane marginal, then Tonelli/Minkowski and reflection; it includes all short-time lines and the spatial infinity tail. This agent reviewed the analysis, but independent review remains a separate dependency. If every8256 patch contributes a certified bound summing to B, then I<= (B^(2/3)+Etime)^(3/2). The code computes this with Arb and compares to exact29/50000. Pilot four patches do not supply this B. Minimum g, general L1/allr transfer, and final maximal constant remain other workstreams.

## Actual finite pilot and reproducibility

Four patches (127,57),(127,127),(127,0),(0,0),1089 nodes each: actual384bit run5.85890675s total. Each ~1.46s. Max node radius3.938477194062e-36; max Bernstein radius5.878496676302e-17. Range bounds finite and positive. Their weighted I_H upper contributions respectively4.689856449922e-8,5.789640498317e-8,7.353195370358e-8,4.726469891163e-8. No increase to512bit was needed. Author selfcheck (no new point evaluations) verifies exact Bernstein means, all1089 DCT inverse enclosures, record chord recomputation. It is not independent audit.

Code cb11_integral_cheb_check.py is immutable for the fixed full batch. Stage four writes four-result and a622KB cache with exact midpoint/radius dyadics, already including node displacement. Full parts read these four matrices only; all other nodes are actually evaluated, no full node cache. A part uses patch serial index i(i+1)/2+j modulo parts. Each process writes a distinct part file. No subagents are spawned. Parent coordinates up to12 local processes and actual memory/CPU; ideal runtime from pilot is17minutes on12 cores, not an asserted completion time.

Commands in this directory, using existing venv Python:

    python cb11_integral_cheb_check.py --stage part --part 0 --parts 12
    ... distinct part1 through11 ...
    python cb11_integral_cheb_selfcheck.py --parts 12
    python cb11_integral_cheb_check.py --stage aggregate --parts 12

Author selfcheck additionally checks exact expected patch sets in every part, frozen checker hashes, nonnegative L<=mean<=U, all chord bounds, and exact Fraction/packed endpoint agreement. Aggregate requires exactly all8256 unique patches and sums exact rational contributions. Do not call four stage again merely for a baseline. Portable MATLAB calibration will compare fixed recorded chord/budget values only; Python Arb is the actual rigorous arithmetic. Source/manifest hashes and execution receipts distinguish authored pilot, parent full run, and independent audit. No claims about weighted backgrounds; no ledger/state/other project modifications or commits.
