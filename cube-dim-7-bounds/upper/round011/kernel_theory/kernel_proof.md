# R011: deterministic complete-fan evaluation of the true truncated Green kernel

2026-10-27. Sol proof/specification only. This owner runs no node, divided difference, kernel, quadrature, MC or MATLAB computation. The original7D body64/unit-weight/common-time Cauchy model is fixed. R010's accepted positive ordered-cone decomposition is input:135cells,212top simplices on46rays. These counts are source claims to be hash/checker validated, not recomputed here. No full-domain g or new C7 claim follows from a point certificate.

## 1. Compile the structural certificate without changing scale

For catalog ray r let w_r be primitive integer, xi_r=Mw_r with detM=1, and q_r=psi(xi_r)>0 integer. For simplex S of7rays define

    alpha_S=|det W_S|/product_(r in S)q_r>0.

Validate all212unique simplex tuples and their source cell signs/q equalities. Store xi_r,q_r,alpha_S as exact integers/rationals. Neither the lambda denominator N_r nor the u-coordinate volume determinant belongs in alpha: both normalisations cancel as proved in R010.

Set A_sum=sum_S alpha_S, an exact positive rational. The untransformed ordered-cone contribution is

    F_epsilon(y)=sum_S alpha_S DD6[f_epsilon](b_r(y):r in S),
    b_r(y)=y·xi_r/q_r,
    f_epsilon(b)=log(epsilon-i b).                                  (1)

R010 proves the complete group formula

    H_epsilon(x)=(2pi)^(-7)sum_(Q in B7)F_epsilon(Q^T x).             (2)

No frequency cutoff, probability quadrature, shape-expectation error or short-time correction appears in this exact epsilon>0 formula. The only numerical error in a complete evaluation is its outward arithmetic error. The mathematical formula still requires the accepted complete fan sources; a ray list alone would not justify(2).

For completeness, the time integral has this sign because for distinct real knots the partial-fraction coefficient of1/(t-i b_j) is i^(-6)/product_(k!=j)(b_j-b_k), and i^(-6)=-1. Integrating from epsilon to infinity cancels the common upper-end logT terms, leaving the positive DD6log(epsilon-i b). Coincident knots follow by continuity at epsilon>0. Cone coordinates give Jacobian|detR| and factors(tq_i-i x·r_i)^(-1), so extracting product q_i gives exactly alpha. Also psi>=c|xi| for some c>0 ensures integral exp(-epsilon psi)/psi finite, justifying absolute Fourier/Fubini exchange; cone boundaries have zero Lebesgue measure. Thus signed contributions may cancel but the complete expression is the actual H.

## 2. Correct phase-orbit/stabiliser reduction

Let O_x={Q^T x:Q in B7}. Every distinct phase y has the same fibre size c_x=|Stab(x)|, so

    H_epsilon(x)=c_x(2pi)^(-7)sum_(y in O_x)F_epsilon(y).

If x has z zero coordinates and positive absolute-value multiplicities m_lambda,

    c_x=2^z z! product_(lambda>0)m_lambda!,
    |O_x|=2^(7-z)7!/[z!product m_lambda!].                         (3)

Nonzero coordinate signs are uniquely determined when matching a permutation; only zero coordinates allow free sign flips. This proves the formula even if the physical input contains negative entries.

For x!=0 the orbit splits into opposite pairs. Choose O_x^+ by requiring the first **nonzero** coordinate of y positive. Since F(-y)=conjugate(F(y)),

    H_epsilon(x)=2c_x(2pi)^(-7)sum_(y in O_x^+)Re F_epsilon(y).      (4)

This is reduction by the actual phase Q^T x, not by a symmetry of a single frequency cone.

Fixed value pilots:

- x=e1: orbit{+/-e_j:j1,...,7}, size14,c_x=2^6·6!=46080. Half-orbit has the7positive axes.
- x=(1,...,1): orbit all128sign vectors,c_x=7!=5040. Half-orbit has first coordinate+1 and the other6signs arbitrary.
- Generic input with7distinct nonzero absolute values: orbit645120,c_x=1; half-orbit322560.

At x=0 the orbit has one point fixed by negation. **Do not apply the factor2 half-orbit rule.** Directly,

    H_epsilon(0)=645120 A_sum/[6(2pi)^7epsilon^6]
                =840 A_sum/[pi^7epsilon^6].                        (5)

This gives a rational positive coefficient independent of the DD implementation. Similarly p_t(0)=5040 A_sum/(pi^7t^7), and the spectral unit-ball volume is128A_sum. The known determinant512basis bound implies the rational sanity inequality5040A_sum<=1/512. It is a check, not an independent proof of fan coverage.

## 3. Exact Hermite divided differences and their sign

For real rational knots b0,...,b6, sort them exactly and let D[i,0]=f(b_i). For order k=1,...,6,

    if b_i=b_(i+k): D[i,k]=f^(k)(b_i)/k!,
    otherwise: D[i,k]=(D[i+1,k-1]-D[i,k-1])/(b_(i+k)-b_i).          (6)

Sorted endpoint equality means the entire subblock is equal. The other denominator is an exactly nonzero rational, converted directly to a ball; never obtain it by subtracting two rounded knot balls or replace a small difference by0.

The correct confluent values are

    f_epsilon^(k)(b)/k!=−i^k/[k(epsilon-i b)^k], k>=1.              (7)

The minus sign is outside i^k. Real epsilon>0 places every logarithm in the right half-plane. At all b=0, DD6=1/(6epsilon^6), verifying the dimension7 positive sign against integral_epsilon^infinity t^(-7)dt.

Because the knots and recurrence coefficients are real,

    Re DD6[f_epsilon]=DD6[h_epsilon],
    h_epsilon(b)=0.5log(epsilon²+b²).                              (8)

A real-Arb implementation of(8) is exactly equivalent to the complex formula and can avoid unnecessary imaginary arithmetic. Its Hermite values are the real parts of(7); for rational epsilon,b and k>=1 these are exact Gaussian-rational real values. At b0 the odd derivatives vanish and the even coefficients are(-1)^(m+1)/(2m epsilon^(2m)). The raw cone or phase contribution may be negative; keep every signed term.

## 4. Common integer knots remove Fraction sorting from the full loop

For physical input x choose a common exact denominator d_x, x_j=n_j/d_x. Put L=lcm_r q_r, B=d_xL and zeta_r=(L/q_r)xi_r, all exact integers. For a phase y_j=eta_j n_(pi(j))/d_x,

    b_r(y)=a_r/B,
    a_r=sum_j eta_j n_(pi(j))zeta_(r,j) integer.                   (9)

All46knots have the **same denominator**, so sorting/equality/differences require only integer operations. There is no per-knot approximate comparison or millions of unnecessary Fraction gcd reductions.

Let beta=epsilon B and ftilde(a)=log(beta-i a). Constant logB disappears from a positive-order divided difference, and scaling gives

    DD_b^6 f_epsilon=B^6 DD_a^6 ftilde.                            (10)

The real version uses0.5log(beta²+a²). All recurrence derivatives are now with respect to a, so use−i^k/[k(beta-i a)^k] and multiply by B^6 **once at the top**. Do not insert inconsistent B powers into cached suborders. The all-zero test is B^6/[6(epsilon B)^6]=1/(6epsilon^6).

Precompile alpha_S B^6 per physical input. Per phase, cache the function values for unique integer a_r, exact Hermite values on demand, and DD submultisets. A DD cache key is the full sorted integer tuple with multiplicities, plus the fixed epsilon/precision context. Symmetry of divided differences permits exact coalescing of identical top knot multisets, summing their positive alpha coefficients exactly. Do not merge distinct phases merely because their46unordered values look alike; the weighted simplex/multiset incidence must also agree exactly.

Each uncached7-knot table has28entries (21recursive differences). Across212simplices there are at most5936such entries before sharing, but the actual cache counts must be recorded. A full generic run requests68,382,720top DDs before such exact caching; do not describe all requests as distinct recomputations if the cache reduces them.

## 5. Exact phase ordering for the generic dyadic point

The physical point is fixed by the seven normalized binary64 hex strings in `testcases.json`, converted first to exact dyadic rationals. Its nonzero absolute values are distinct. Let pi run over lexicographic permutations of indices0..6. Let eta1=+1, with eta2,...,eta7 lexicographic in{-1,+1} (−1first). Then y_j=eta_j x_(pi(j)) enumerates exactly322560opposite-pair representatives.

Use phase_id=64*permutation_rank+sign_rank. There are5040natural chunks, one permutation/64signs each. The first4chunks are the256-phase generic pilot and are included once in any later complete run; the remaining5036chunks are not started without explicit main release. A partial phase sum is not positive in general and is not an H value. Its midpoint or extrapolation cannot be compared to a g target.

The phase representation is also an exact signed-permutation Q: Q_(pi(j),j)=eta_j. It can be retained symbolically for derivative work; arbitrary phase sorting would lose that information.

## 6. Precision/error acceptance, with no blind clipping

The proposed deterministic value target is raw full-H interval width<=10^(-12), an absolute target, not a claimed achieved error. All source weights and knots are exact; only transcendental evaluations and ball arithmetic have rounding radii.

Start the fixed pilots at256bits. If a prescribed pilot fails only its width target, one transparent same-input promotion to512bits is allowed in the proposed contract; record both attempts and count the extra evaluations. No equality threshold, row deletion, negative-cone clipping or midpoint replacement is permitted. Failure at512bits is reported and stops the batch pending a new main decision. Different precision alone is not an independent arithmetic kernel.

For a complete symmetry-orbit pilot, require the actual raw value interval narrow enough and positive; compare its epsilon monotonicity and origin peak sanity without changing a raw interval to make a test pass. An analytic intersection with positivity may be valid mathematically but must not hide a width failure.

For the256-phase generic pilot, require partial width<=10^(-12)*256/322560 as a **planning gate**. Multiplying this by the remaining phase count estimates cost/width but does not certify unseen cancellation. Full generic acceptance still requires the actual complete orbit width<=10^(-12). If it fails after the fixed run, preserve its genuine enclosure and report unresolved width; do not automatically rerun68million DDs at another precision.

## 7. Derivative interface and a stabiliser caveat

For affine knots b_j(x) and a fixed real direction v let c_j=(Qxi_j·v)/q_j. Analyticity of divided differences gives

    D_v^k f[b0,...,b6]
      =k! sum_(alpha0+...+alpha6=k) product_j c_j^(alpha_j)
          f[b0 repeated1+alpha0,...,b6 repeated1+alpha6].          (11)

This follows by differentiating the Cauchy representation of DD and extends to confluent knots. For k1 it is the sum of seven order7DDs, each with one duplicated knot. Do not differentiate an equality branch as though coincident knots stay coincident in a perturbation; use(11) or an analytic AD implementation that preserves the multiknot function. Hermite data are available from(7) at all orders.

A value-phase stabiliser fold is **not automatically valid for arbitrary derivatives**. Several Q with identical Q^Tx can have different Q^Tv. Their derivative coefficients/tensors must be summed or stabiliser-averaged. The simple value multiplicity is safe for a direction fixed by the physical stabiliser, such as radial v=x. At a generic distinct-absolute-value point the stabiliser is trivial and the exact permutation/sign record recovers Q. Full generic group pairing remains valid for derivatives by conjugacy.

The sharp Euclidean floor established from the sealed R003 spectral catalog is

    psi(xi)>=c|xi|, c=32sqrt(2).

Its proof, primary computation provenance and cross-checks are in extensions.md. Consequently the global Fourier derivative bound, independently of evaluating(11), is

    ||D_v^kH_epsilon||infinity
      <=Gamma(6+k)/(120pi^4·c^(7+k)epsilon^(6+k)), |v|=1.

The complex spatial tube has |Im x|<c epsilon. This replaces the earlier safe c8 bound; for k1 the bound improves by exactly2^20. The sharpness check is a new scalar functional of an already complete catalog, not a new vertex enumeration. The primary thread actually evaluated it; this theory owner did not.

This supplies a deterministic future spatial covering/Taylor interface. It does not turn the present single-point computation into a full-face g certificate.

## 8. Scope

The default proposed numerical work is fixed scalar DD tests, origin checks, four complete symmetry-orbit H values, and a256-phase generic pilot. A complete generic H at epsilon1/1024 requires a separate main release after pilot width/cost review. Generic epsilon0G, even moments and alternative1Dtime integration are derived in the separate `extensions.md` and do not start automatically.

Every accepted complete evaluation is of the true original H_epsilon, with no shape-expectation approximation. The remaining g/I/all-L1 transfer uses the original full-G source and is not established by a point interval. No new global g or C7 upper bound is asserted in this specification.
