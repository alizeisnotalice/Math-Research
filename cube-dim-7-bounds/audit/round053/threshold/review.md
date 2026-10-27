# R053 fixed-tuple operator threshold: independent audit

2026-10-28. **ACCEPT.** With the already accepted R051 full mechanism, the locked R053 proof establishes, for every integer n with L=log n>=2^32,

    a=73/8, s=293/16, d=1023/2^20,
    sigma=s d^2, t=n/[(1-sigma)L],
    C_n<=F(A_n(t),B_n(t)-gamma n)
        <n log n-n log log n-2^-16 n.

Here gamma is the exact rational coefficient in Section 1 below. This is a complete original-operator threshold, not an acceptance of a candidate objective or a finite base-point experiment. No blocking error or requested author correction was found. The main coordinator retains final proof/integration responsibility.

The threshold source retains the historical wording “conditional on R051 full acceptance.” That prerequisite is now satisfied by the locked full review identified below. This audit does not rewrite the frozen source or its generation history. It does not certify a grid-best/global-optimality identity, a varying tuple, a uniform threshold for the continuous family, or a lower bound for C_n.

## Sources and actual audit work

All three files were read in full and their actual hashes recomputed:

| File | SHA-256 |
|---|---|
| parameters `upper/round053/explicit_candidate_threshold/THRESHOLD_PROOF.md` | `26adb4e31f0026637ea17d8f7943bc2933dce97b8655e121f2604b6ddeea2775` |
| geometry `upper/round051/flexible_budget/PROOF.md` | `aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6` |
| soft `audit/round051/flexible_budget/proof_review.md` | `66355a1d58d3b6a832e844645b9e61ca9e3c5a2f27bc0958947054aa9144f960` |

The full absolute paths and runtime data are recorded in `microcheck.json`. Only this auditor's `audit/round053/threshold/` directory was written. No author source, Git state, shared MATLAB, n=7 work, production search, R052 enumeration, old grid, or author 38-check program was run or changed.

## 1. Exact reconstruction of the fixed budget

Rebuilding directly from the R051 formulas, with epsJ=epsD=epsS=epsT=1/1024, epsA=1/64 and eps0=2^-20, yields:

    P0=3579/4096, qmin=1-3d/4,
    ulo=3579/8200,
    uhi=8796093022208/7674592741893,
    M*=max(4,2uhi/qmin)=4,
    m*=64449/134315936,
    b*=11148919/38310400,
    Gamma=718536680631/20582788938137600,
    gamma=735063024285513/21076775872652902400,
    sigma=306632997/17592186044416,
    c_raw=6468305130094850733787980147/
          370780099533215320770594052505600.

Every listed rational matches the threshold source exactly. All domain inequalities are strict where needed: a>3, s>2a, d<=1/16, d<1/(32a), sigma<1/2; the error fractions are positive with epsJ/epsA/epsD<1; b*>0. In particular Gamma>gamma>0 and gamma<1. The exact c_raw exceeds c=2^-16. The positive scalar charge is sigma/(1-sigma), not an omitted change of t or a replacement by 2sigma. The deficit gap Gamma-gamma=epsD Gamma is preserved independently of the scalar gap c_raw-c.

## 2. All geometric and extended-coordinate guards

Let L0=2^32. Global Taylor positivity gives n=e^L>=L^2/2 and n>=L^8/8!. With sigma<1/2, x lies in [L/2,L], x>1, t>=L/2 and t<n. Source radii exceed 8 and targets are below 2, so target/source supports are disjoint, lie inside the original radius-t envelope, and have the strict >4/<3 separation required for actual histories.

At the base h/d<1/8 and zeta*t/h<=2d/(eps0 n^2)<1/2; these upper bounds decrease with L. Thus epsilon=zeta*t<=h/2, h+epsilon<d/2, kappa>1 and beta<2. Source R=q+d-h>1 and g=d-h+epsilon*q lies in (0,d], with g<=R/2.

For a fixed target, q0=r-d+h can be just outside the source band. The entire extension [q0-2h,q0] still has

    q>1/2, q<1, 1<R<=r<2, rho>4r>=4R,
    g'=epsilon<=d, T=(|w|-(kappa-1)rho)_+<=g+2h<2d.

The extremal lower q is 1-3d/4-h and lower R is 1+d/4-2h; both improve as h decreases. The lower rho bound (L/2)(1-3d/4-h) increases. These prove every local-formula/derivative/pinned-face domain on the whole tail, including the extended baseline q. Neither the old tiny d restriction nor the old soft exponent is used. All remaining n, nd, Ld guards are already strict at L0 and strengthen for larger L.

## 3. New target probability and the exact mass coefficient

The free-coordinate baseline has E J=1 and log mean/variance bounded by -6g^3/rho^2 and 12g^3/rho^2. The pinned coordinate contributes a deterministic absolute log at most 2g/rho and zero variance. Summing on each face, uniformly over r, gives

    |E S| <=(24d^3 L^2+4dL)/n,
    Var S <=48d^3 L^2/n.

Using n>=L^8/8! gives M(L)=8![24d^3/L^6+4d/L^7] and V(L)=48*8!d^3/L^6. Both decrease globally. Independently, M(L0)<1/32 and V(L0)/(1/8-M(L0))^2<epsA/3. Chebyshev therefore bounds failure of the baseline log event by epsA/3 at every face/radius throughout L>=L0. No pinned coordinate is assigned a free uniform law, and no log mean is silently set to zero.

For Kw=1+Bin(n-1,2d/r), the direct counting bound 2z/(2z-1)^2, z=nd, is <=8/(9z) for z>=2. The algebraic difference is nonnegative there (14z^2-32z+8>=0). At the base, 8/[9(L0^2/2)d]<epsA/3; the denominator only grows.

The target corner-complement has fraction at most (3/4)^(n-1) per face: 1/r-1/t>=1/(1+3d/4)-2/L>1/2. Since n>=33 and (3/4)^32<epsA/3, its bad rate is also below epsA/3. A union bound suffices; the three events are not assumed independent. Thus the good set retains at least 1-epsA=63/64 of **each** full face, not the old 5/8.

The annular mass is n log[(1+3d/4)/(1+d/4)]. Integrating 1/r over its width d/2 bounds it below by nd/[2(1+3d/4)]. Consequently lambda(A)>=m*n with the exact m* above. There is no substitution of a coarse nd/4 mass or old retained fraction in the new budget.

## 4. Uniform matching product, individual precision, and degrees

The accepted local formula is locally Lipschitz/absolutely continuous at capture onset, with the derivative bound valid a.e. The geometric extension in Section 2 permits integration across its kinks. Summing on Kw<=4nd and integrating over length <=2h gives the exact matching log change 1280Ldh=5/4096. Baseline 1/8 plus that change is 517/4096; exp(-z)>=1-z and exp(z)<=1/(1-z) give the stated P0 product enclosure for every matched q, not just a sampled set.

Individual J is separately controlled. The piecewise closed form has A(T)<=4T/3 also for T>=g, and B<=g^2/(3R), so on T<=2d,

    |J-1|<=2d/rho<=4dL/n<=4d*8!/L^7.

The last bound is <epsJ at L0 and decreases as L^-7. This pointwise uniform estimate does not involve an n-coordinate failure union and is never raised to the n-th power to replace the independent product estimate.

Pinned source faces remain in the exact angular identity

    U=(product J)/(2ell_*) sum_i Ni/Ji,
    1-d/2<=ell_*<=1, 1<=Ni<=2.

Thus ulo*n<=U<=uhi*n with the reconstructed sharp rational coefficients. The inverse-row law has normalizers J>=1/2 (also available coarsely from the coordinate formula), free mean logloss <=d^2/rho and deterministic pin loss <=d. Since Ld>=1,

    E_inv[-logG]<=d+2Ld^2<=3Ld^2.

With theta=exp(-aLd^2), Markov under this **inverse captured law**, not raw source probability, retains at least 1-3/a=49/73. Matching log widths are in [h,2h/qmin]. Hence the pre-cap row bounds are (1-3/a)ulo*k and 2uhi*k/qmin. The full-column envelope gives column<=4k pointwise. The maximum-degree constant is exactly M*=4, rather than a copied old value 32.

## 5. Distinct source-cap and tau payments on the full tail

The O(nd exp(-Ld/32)) in the generic proof can be replaced by the explicit coefficient used here. Indeed the original corner-complement shell coefficient is at most (n/2)exp(-Ld/32), because

    (n-1)(1/rho-1/t)/2=(1-1/n)xv/(2q)
                           >=xd/16>=Ld/32.

Its source log-band width is <=d, so the full mass is <=nd exp(-Ld/32). Bounded b<=1 pays removal by this unweighted mass for every target. Dividing by k gives

    Cap/k<=2^20 Ld^2 exp(-eta L), eta=d/32-a d^2>0.

Taylor order four gives 24*2^20*d^2/(eta^4 L^3). The coefficient is correct, the power decreases on the whole tail, and its exact L0 value is <epsS. Thus final row >=[(1-3/a)ulo-epsS]k.

The chosen tau is real: on the x-range used here, B>8n, B<=nL and gamma<1. A<=n exp(-sigma L) implies

    tau/k<=2^20 d L^(3/2) exp(-nu L),
    nu=(s/2-a)d^2=d^2/32>0.

This exponent includes payment of the soft source factor. At L0, nu L0=1046529/8192>64 and the prefactor is 1023*2^48<2^58. Global Taylor order 64 yields the decreasing power L^(-125/2); its base value is <2^58*64!/2^384<epsT. The coarse prefactor bound is only used to compare the base value; subsequent L are covered by the actual Taylor power function, not by pretending that the prefactor remains constant. The source-cap exponent eta and selection exponent nu are different and are both checked.

Therefore final row minus tau is at least b*k, and target mass/degree yield an ambient deficit at least m*b*n/M*=Gamma n. The strict margin Gamma-gamma=epsD Gamma survives; epsD is not inadvertently reused as a source-cap or scalar error.

## 6. Actual finite histories and the general operator

For every original dilation mathfrak a>=8n^4, t<n implies mathfrak a>2t; beta<2 and zeta=n^-3 give

    alpha beta-kappa-2/(mathfrak a+1)
      =zeta-(beta+2)/(mathfrak a+1)
      >n^-3-4/(8n^4)>0.

An actual uncaptured target increment meeting A has radius <3, and every included source increment has radius >=rho>4. Thus each source precedes the **same actual target**. The strict centre-displacement guard puts that target centre in the original dilated source cube. Source factor b<=1 preserves the target's original eligibility budget tau; no theta^-1 factor and no distributional assumption on actual mu is introduced.

Both supports lie in D infinity, where lambda_mathfrak-a=alpha^-n lambda exactly. The single missing measure lambda_mathfrak-a-mu_mathfrak-a supplies the bipartite energy deficit; source/target sets are disjoint. Energy scales by alpha^-2n and degrees/target mass by alpha^-n. Since alpha^-n>=1, the lower deficit Gamma n only improves. Thus gamma n is a valid reduction of the full actual uncaptured load for **all** such dilations and finite histories.

The witness remains the original selected cube's mathfrak-a dilate for the same arbitrary L1 input. Strict compact-level-set finite covers and the original witness payment stay unchanged. Holding n,t,tau fixed, only the proved scalar inequalities are sent to mathfrak-a->infinity; no selector or history convergence is used. Inner regularity then gives the full level set, and every positive common radius remains available. This closes the original operator statement at L0, rather than a continuum/model-only objective.

## 7. Complete scalar cost and its strict slack

At the actual x=(1-sigma)L,

    B<=T+[sigma/(1-sigma)]n+1,
    F(A,B-gamma n)<=T-c_raw n+nR(L),
    R(L)=1/n+exp(-sigma L)+2sqrt(L)exp(-sigma L/2).

The exact rational c_raw-c is positive. At L0, sigma L0>2^16; global Taylor bounds give

    R(L0)<=2^-63+16!/2^256+2^17*16!/2^240<c_raw-c.

The three upper envelopes decrease as L^-2, L^-16 and L^(-31/2) on all L>=L0. Thus the strict inequality persists on the whole tail. A and the positive square-root cross term are both paid. This is why the final saving is the strict interior coefficient 2^-16, rather than gamma or c_raw itself.

## 8. Full-cost ordering keeps old own-domains

For L>=Lcmp=2^40, the R049 actual deficit is 2^-40 n; the R042 actual deficit is 2^-31 nL^-1/3<=2^-40 n. At their own t values, x belongs to [L/2,L], and the old root arguments are positive. Their complete costs are >=T-n(4/L+2^-40), whose loss at Lcmp is 5*2^-40<2^-16 and decreases. This lower-envelope direction is correct and concerns formula costs, not a lower bound for C_n.

For older listed fixed-t instances, actual deficits are <=2n sqrt(logL/L) on their separate domains. The resulting common loss 4/L+2sqrt(logL/L) decreases when logL>1. Since log Lcmp=40log2<40<49, its base is <2^-38+14/2^20<2^-16. It follows that the new complete cost is strictly smaller there without comparing only displayed negative terms.

These are scalar sortings at Lcmp. An operator-level comparison also requires each old bound's own domain: R049 L>=2^52, R042 L>=2^108, R036 log L>=128, R034 its positivity/Y^2 guards, and fixed R040 its own unquantified N. In particular Lcmp does not certify R036's or R042's operator hypotheses. The all-t unit-core asymptotic comparison remains only eventual, with no fabricated explicit start for its O term. These qualifications are all present in the source.

## 9. Independent computation and completeness

`microcheck.py` was actually run with Python/Fraction, PID 84467, and all 23 compressed independent certificates PASS. It reconstructs one fixed tuple directly from the generic formulas, matches the exact stated rational fields, and checks domain/degree/sign conditions, geometric/extended-q guards, direct Chebyshev and target counting/cap rates, individual J precision, source cap, tau, complete scalar slack and both cost-comparison gaps. It imports no author helper and does not traverse a parameter list. `microcheck.json` preserves exact rational bounds, source hashes and runtime identity.

Those finite certificates verify constants only. Infinite tails are covered analytically above by increasing geometric margins; decreasing M/V, counting probabilities and fixed corner exponent; the a.e./AC whole-q bridge; and the global Taylor powers L^-3, L^(-125/2), L^-2, L^-16 and L^(-31/2). Every eventual condition in the frozen generic mechanism has an explicit replacement here. No residual O term or unknown N remains for this complete fixed tuple at L>=2^32. No claim about a continuous parameter optimum, enumeration completeness, or other tuples' thresholds is supplied.
