# R049: independent audit of the actual operator threshold

2026-10-28. **ACCEPT.** The locked proof closes all sufficient conditions of the accepted R045 soft inverse-column family at d=2^-24, for every integer n with L=log n>=2^52. Consequently

    t=n/[(1-2^-43)L],
    C_n<=F(A_n(t),B_n(t)-2^-40 n)
        <n log n-n log log n-2^-41 n.

This is an operator threshold, including same-history eligibility and general L1/all-positive-common-radii passage. It is not inferred merely from the scalar comparison, and is not claimed to be a minimal threshold. No blocking error or requested author correction was found. The coordinator retains final proof/integration responsibility.

## Locked sources and audit method

Threshold proof read in full: `/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round049/explicit_threshold/THRESHOLD_PROOF.md`, SHA-256 `b8c385a90cefbe3b234f5269be749ad13984243e01d9109658bb38e27304e33c` (224 lines).

Its eventual guards were checked against the entire actual frozen R045 mechanism `/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round045/soft_inverse_column/PROOF.md`, SHA-256 `add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9`. The independently proved local formula/AC lemma has hash `3fa91ce60a9ef05df270d1fa85a314de6a8573f8223fe1ea69f8433d6b6caa39`. All three actual hashes are recorded in `microcheck.json`.

This audit independently checked the numerical constants and every infinite tail, rather than accepting the author's 37 base-point PASS entries. The author's verification code was not executed. Only this independent audit directory was written; no upper source, Git state, old numerical baseline, production calculation, n=7 work, or MATLAB was touched.

## 1. Geometric guards hold on the whole tail

Write L0=2^52, d=2^-24, epsilon0=2^-20, sigma=2^-43 and h=epsilon0/(Ld). Taylor positivity gives n=e^L>=L^2/2 and n>=L^8/8! for all L>0. Thus t=n/x>=L/2, source radii rho>=t/2>8, target radii r<2, n>=17, nd>=8, and Ld>=1 throughout L>=L0. In particular t>1, the target and source bands lie inside ||z||<t, are disjoint, and the required source/target radius separation is strict.

The quantities h/d=2^28/L and zeta*t/h<=2d/(epsilon0 n^2) decrease to zero, and their base values are strictly below 1/8 and 1/2. Therefore epsilon=zeta*t<=h/2, h+epsilon<d/2, kappa>1, beta<2, and g=d-h+epsilon*q lies strictly between 0 and d. All source q are between 1/2 and 1; source R=q+d-h>1, g<=R/2 and rho>=4R follow with substantial margin.

For each fixed target, the baseline q0=r-d+h may extend beyond the original source band. The entire required extended interval [q0-2h,q0] is nevertheless inside (1/2,1), its R lies in [r-2h,r] strictly above 1, and its rho>4R. Here g'=epsilon<=d, and T=(|w|-(kappa-1)rho)_+<=g+2h<2d. These statements cover the coordinate formula/derivative domain on **all** q, not only the original source interval. The bounds strengthen as L grows: h decreases, the lower q/R margins increase, and the crude rho lower bound (L/2)(1-3d/4-h) increases. No eventual geometric guard in the frozen source remains unquantified.

## 2. Free/pinned baseline and target mass

At q0, R=r. The exact one-coordinate cancellation E J=1 uses a free uniform target coordinate; its second moment estimate remains nonasymptotic. In the author's derivation, E T^2=g^3/(3R), B<=g^2/(3R), and A<=4T/3 give

    E(J-1)^2 <= [2g^4/(9R^2)+8g^3/(27R)]/rho^2
                <=3g^3/rho^2.

The log mean is nonpositive by Jensen, not zero. The log inequalities give the claimed -6/12 mean/variance constants. One pinned target-face coordinate contributes a deterministic log bounded by 2d/rho and no variance. With n/rho<=2L, summing the n-1 free coordinates and that pin gives precisely

    |E S| <=(24d^3 L^2+4dL)/n,
    Var S <=48d^3 L^2/n.

Replacing n by L^8/8! gives M(L) and V(L), positive sums of L^-6 and L^-7, respectively L^-6; all decrease on the entire tail. The independent check verifies directly V(L0)/(1/8-M(L0))^2<1/8, as well as M(L0)<1/32. The probability estimate therefore holds uniformly over every target radius and face for all L>=L0. No covariance or missing pinned variance is assumed.

For Kw=1+Bin(n-1,2d/r), mean<=1+2nd and variance<=2nd. The direct Chebyshev bound is f(z)=2z/(2z-1)^2 with z=nd>=8. Its derivative is -2(2z+1)/(2z-1)^3<0; at z=8 it is 16/225<1/8. This independently establishes the manuscript's slightly coarser counting bound. The two target good events need not be independent.

For all target r, 1/r-1/t >=1/(1+3d/4)-2/L>1/2, so the excluded original corner fraction is at most (3/4)^(n-1)<1/8 for n>=17. The target annulus mass is n log[(1+3d/4)/(1+d/4)]>=nd/4; after the three bad-event bounds the retained fraction is at least 5/8. This exceeds the asserted a_*=lambda(A)>=nd/16. The upper mass bound a_*<=nd also holds.

## 3. Whole matching bridge and inverse-row law

The frozen local formula is valid on the extended domain of Section 1. It is locally Lipschitz/absolutely continuous across the true capture-onset kink, so its a.e. derivative estimate integrates across that kink. The q-dependence keeps epsilon=zeta*t fixed. Summing the derivative bound on a target with Kw<=4nd gives at most 1280Ldh=5/4096<1/8 over the whole extended interval. Together with the baseline event, this yields product J in [1/2,2] at every matched q. No continuum statement is inferred from sampling.

Individual J bounds are separately supplied. A(R,g,T)<=4T/3 extends to T>=g because A(T)=A(g)+T-g; with T<=2d, J>=1-4d/(3rho)>1/2. Also J<=1+d^2/(3Rrho)<2. Thus the product bound is not incorrectly used to infer coordinatewise bounds.

The complete inverse angular coefficient includes each pinned source face:

    U=(product J)/(2 ell_*) sum_i N_i/J_i,
    ell_*=1-g/(2R) in [3/4,1], 1<=N_i<=2.

It lies in [n/8,8n]. The source law used for clipping is capture/G normalized separately on each captured face, with free normalizers J>=1/2 and deterministic pin loss at most d. Under that law E[-log ell]<=8g^2/(9Rrho)<=d^2/rho, since actual R>1. Hence E[-log G]<=d+2Ld^2<=3Ld^2. Markov at 12Ld^2 loses at most 1/4; this is an inverse-row probability, not a raw-source probability. It proves the weighted angular lower bound 3n theta/32. Matching intervals have ordinary length at least h and at most 2h and q in (1/2,1); their source log length is between h and 4h. The full weighted row is therefore between 3k/32 and 32k without any remaining eventual assumption.

## 4. Cap removal and tau: global Taylor bounds

The corner-complement exponent obeys

    (n-1)(1/rho-1/t)/2=(1-1/n)xv/(2q)
                          >=xd/16>=Ld/32.

Source log width is at most d. Since b<=1, the entire unweighted corner mass is a valid row-removal bound. Division by k=nh exp(-12Ld^2) yields

    Cap/k <=2^20 Ld^2 exp[-L(d/32-12d^2)]
          <=2^20 Ld^2 exp(-Ld/64)
          <=24*2^44/(L^3 d^2).

All powers of d and the factor 64^4 are correct. The last expression is globally decreasing; the independent check at L0 gives the stronger bound <1/64. Thus raw row 3k/32 minus cap loss is even >5k/64, and the claimed retained row >=3k/64 is safe. The upper row 32k and pointwise column 4k are inherited directly under the same bounded kernel. No source-concentration hypothesis is needed.

For x in [L/2,L], the original cap bound gives B>=n[L-log L-4/L]>8n on L>=64 and B<=nL; the latter follows from x>e. Delta=d/65536<1, so B-Delta n>0. The chosen tau=sqrt[A(B-Delta n)] is legitimate, and

    tau/k <=2^20 d L^(3/2) exp(-4d^2L).

The source factor theta has been paid: the exponent is -16d^2L+12d^2L=-4d^2L. Applying e^u>=u^64/64! globally, this upper bound is a constant times L^(-125/2), decreasing on the full tail. At L0, u=64 and the bound is 2^74*64!/2^384. The independent exact check gives the stronger result <1/128, hence certainly <3/128. Neither this nor the cap argument relies on numerical exponential evaluations or a finite tail grid.

## 5. Same history, finite a, and all inputs/scales

With target mass nd/16, retained row 3k/64, max degree 32k and tau<=3k/128, the energy deficit is at least 3nd/65536. Keeping only delta=nd/65536 leaves a fixed factor-three safety margin. This is the accepted missing-measure lemma on the common actual history, not a separate independent continuum selection.

For all a>=8n^4, t<n implies a>2t. With beta<2,

    alpha beta-kappa-2/(a+1)
       =zeta-(beta+2)/(a+1)>n^-3-4/(8n^4)>0.

A target actual increment has radius <3 (from cr_j<alpha+r and c>1), while each source actual increment has radius >=rho>4. Thus every included source precedes the same target. The strict capture/centre guard yields that target's original eligibility relation. Bounded b<=1 means its weighted actual row is <=tau, with no inverse-theta charge and no condition on the actual measure's distribution.

On both supports lambda_a=alpha^-n lambda exactly; ambient energy scales by alpha^-2n, degrees and target mass by alpha^-n. Since alpha^-n>=1, the missing-mass bound does not deteriorate. Delta therefore works for every such a and every original finite history. Targets/sources lie in the original corner cap and inside t, as checked in Section 1.

The witness remains the a-dilate of each original selected small cube, for the same general f in L1. Compact strict-level sets are covered by finitely many such cubes. For fixed n,t,tau only the resulting scalar bounds are sent to a->infinity; no convergence of histories, kernels or selectors is required. Inner regularity returns the full level set. All positive common radii remain available; there is no bounded-scale or compact-support restriction on the operator.

## 6. Complete root cost and strict endpoint saving

The new parameter, not the old x=L-3log L, is used throughout. With sigma=32d^2 and Delta=d/65536,

    F(A,B-Delta n)
       <=T-Delta n+2sigma n+n R(L),
    R(L)=1/n+exp(-sigma L)+2sqrt(L)exp(-sigma L/2).

At the endpoint, 2sigma=2^-42=Delta/4, so obtaining the half-budget saving requires R(L)<2^-42. Taylor positivity of orders 2,16,16 bounds the three terms by positive constants times L^-2, L^-16 and L^(-31/2). Each decreases globally. At L0, the exact combined upper bound is

    2^-103+16!/2^144+2^27*16!/2^128.

The independent check gives the stronger inequality <2^-43. This proves, strictly, F(A,B-Delta n)<T-Delta n/2=T-2^-41 n on the entire tail. Both A and the positive square-root cross term were retained; the proof does not infer a full-cost ordering by comparing deficits alone.

## 7. Actual checks, completeness and limits

`microcheck.py` was run locally using Python and Fraction, PID 83221. Its 15 independent compressed gates all PASS. They check actual endpoint identities; whole-bridge base domain, rho and epsilon guards; n/nd conditions; direct baseline Chebyshev probability; the worst z=8 counting probability; target cap; bridge constant; stronger cap, tau and scalar base margins; and the full scalar endpoint margin. `microcheck.json` records the exact bounds, hashes and runtime. The script neither imports nor runs an author helper or receipt.

The finite exact checks corroborate base constants only. The proof of the infinite range uses the monotonic inequalities explicitly audited above: L^2/L and L^8, decreasing h and epsilon/h, decreasing M/V, the analytic counting function, fixed corner exponential, and the global Taylor power envelopes L^-3, L^(-125/2), L^-2, L^-16 and L^(-31/2). All remaining R045 interfaces are nonasymptotic once these conditions hold. In particular no unresolved O_d term or eventual start remains in the fixed endpoint theorem.

This acceptance activates the R049 operator domain L>=2^52. It leaves the frozen R045 general fixed-d theorem unchanged and does not assert this explicit start for other d, global optimality of d, public priority, or a lower bound for C_n. Existing old explicit bounds retain their own domains and may still be combined by taking a minimum.
