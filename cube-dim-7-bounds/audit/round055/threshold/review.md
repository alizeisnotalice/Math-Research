# R055 simple fixed tuple: independent operator-threshold audit

2026-10-28. **ACCEPT.** Conditional only on the already satisfied full acceptance of R051, the frozen R055 source proves for every integer n with L=log n>=2^27:

    a=9, s=19, d=2^-10, sigma=19/2^20,
    t=n/[(1-sigma)L],
    C_n<=Phi_simple=F(A_n(t),B_n(t)-gamma n)
        <n log n-n log log n-2^-16 n,

where gamma=9775559871/281956214374400. No blocking mathematical error or requested author correction was found. The coordinator retains final integration/proof responsibility. This accepts a standalone tuple and explicit threshold, not numerical list membership, production history, a new ranking or global optimality.

## 1. Frozen source and exact budget

The R055 proof was read in full:
`/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round055/simple_tuple_threshold/PROOF.md`, actual SHA-256 `5d2c61adf0e341b2a11daabec61daedd6e009b65a0139bb468711175f0288e9f`.

The R051 full-family source and accepted R053 threshold remain unchanged, with actual hashes respectively `aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6` and `26adb4e31f0026637ea17d8f7943bc2933dce97b8655e121f2604b6ddeea2775`. Their exact paths and hashes are in `microcheck.json`. The new tuple was rebuilt directly from the generic formulas, without a search helper, old table or enumerated list. It satisfies a>3, s>2a, d<=1/16, d<1/(32a), sigma<1/2, every error domain and positive bracket. With epsJ=epsD=epsS=epsT=1/1024, epsA=1/64 and eps0=2^-20:

    ulo=3579/8200,
    M=max(4,2uhi/qmin)=4,
    mass=63/131168,
    bracket=151679/524800,
    Gamma=9555777/275347865600,
    gamma=(1023/1024)Gamma=9775559871/281956214374400,
    c_raw=gamma-sigma/(1-sigma)
         =4893063658542547/295647162275777740800 >2^-16.

The uhi denominator and M branch were independently recomputed. Gamma>gamma>0, gamma<1, and the fixed epsD deficit margin is preserved. The exact positive scalar price is sigma/(1-sigma)=19/1048557; it is not replaced by an omitted price or a comparison of gamma values alone.

## 2. Geometry and whole extended q interval

For every L>=L0=2^27, global Taylor positivity supplies n=e^L>=L^2/2 and n>=L^8/8!. Hence x in [L/2,L] is >1, t>=L/2 and t<n, n>=33, nd>=2, and Ld>=1. Target radii r=1+u are below 2 and source radii rho=tq exceed 8; both lie below t and their supports are disjoint.

At L0, h/d=2^-27<1/8, and zeta*t/h<=2d/(eps0 n^2)<1/2. Both upper bounds decrease with L. Thus epsilon=zeta*t<=h/2, h+epsilon<d/2, kappa>1 and beta<2. For sources, R=q+d-h>1 and g=d-h+epsilon*q lies in (0,d], with g<=R/2 and rho>4r>=4R.

For each target, the actual local bridge is [q0-2h,q0], q0=r-d+h. Although q0 may extend slightly past the original source band, on that entire interval q>1/2, q<1, R>=r-2h>1, rho>4r, g'=epsilon<=d and T<=g+2h<2d. The lower q/R/rho bounds improve as L grows. These check the exact-coordinate, pinned-face and derivative domains without using any residual eventual qualifier from R051. In particular there is no need to identify this start with R053's different tuple threshold.

## 3. Face probabilities, exact target mass and J precision

On each baseline target face R=r, the n-1 free coordinates have E J=1 exactly, log mean in [-6g^3/rho^2,0] and variance <=12g^3/rho^2. The one pinned coordinate contributes a deterministic absolute log <=2g/rho and zero variance. Summing yields the stated M(L),V(L). After n>=L^8/8! they are positive sums of L^-6,L^-7, respectively L^-6, so decrease on the whole tail. Independently at L0, M<1/32 and the direct Chebyshev bound V/(1/8-M)^2<epsA/3. This is uniform over target face and radius; log mean is not treated as zero.

For Kw=1+Bin(n-1,2d/r), the exact Chebyshev bound 2z/(2z-1)^2, z=nd>=2, is <=8/(9z). The latter is <epsA/3 already using n>=L^2/2 at L0 and decreases thereafter. The target corner fraction is <=(3/4)^(n-1): 1/(1+3d/4)-2/L>1/2 gives ell_t(r)<3/4. Since n>=33, (3/4)^32<epsA/3. A union bound, without independence of these events, therefore retains at least 1-epsA=63/64 of every face.

Full target mass is n log[(1+3d/4)/(1+d/4)], bounded below by nd/[2(1+3d/4)]. Thus the sharp mass coefficient 63/131168 is justified. No old 5/8 retention or nd/4 replacement is used in the new budget.

The accepted local derivative is an a.e. estimate with absolute continuity across onset kinks. All its domains were checked above, so integrating it over the entire bridge gives 1280Ldh=5/4096. Combined with baseline 1/8, the log product is bounded by 517/4096, giving P0<=product J<=1/P0 for every matched q.

Individual precision is separate: the exact piecewise formula gives |J-1|<=2d/rho<=4dL/n<=4d*8!/L^7. At L0 this is <epsJ and it decreases thereafter. This pointwise all-coordinate/all-q bound requires no probabilistic union over n coordinates and is not exponentiated n times to infer the product bound.

## 4. Inverse probability, source cap, tau and full deficit

The exact complete source-face coefficient retains the pin:

    U=(product J)/(2ell_*) sum_i Ni/Ji,
    1-d/2<=ell_*<=1, 1<=Ni<=2.

It gives ulo*n<=U<=uhi*n. In the normalized inverse capture/G law, free-coordinate normalizers have J>=1/2; mean logloss per free coordinate is <=d^2/rho, and the pin loss is <=d. Thus E_inv[-logG]<=d+2Ld^2<=3Ld^2. Markov at aLd^2 retains 1-3/a=2/3. This is the inverse-row law, not raw source probability. Matching source log width is in [h,2h/qmin], giving the claimed lower/upper rows and M=4 together with the pointwise column cap.

For this tuple, eta=d/32-a d^2=23/2^20 and nu=(s/2-a)d^2=2^-21 are distinct and positive. The explicit unweighted original corner mass is <=nd exp(-Ld/32), by (n-1)(1/rho-1/t)/2>=Ld/32 and source log width <=d. Since b<=1, it bounds removal in every weighted row. Consequently

    Cap/k<=L exp(-eta L).

Global Taylor order four produces a decreasing L^-3 envelope. At L0 its exact value is 24*2^27/2944^4<epsS, since eta L0=2944. No hidden O coefficient is left.

B>8n, B<=nL and gamma<1 justify the real positive choice tau=sqrt[A(B-gamma n)]. With A<=n exp(-sigma L),

    tau/k<=2^20 d L^(3/2) exp(-nu L).

At L0, nu L0=64 and the prefactor is 2^(101/2)<2^51. Global Taylor order 64 gives a decreasing L^(-125/2) envelope whose base is <2^51*64!/2^384<epsT. This exponent pays the source soft factor; it is not borrowed from the cap exponent. The prefactor comparison applies at the base, while the actual decreasing power covers the tail.

Therefore row-tau>=bracket*k, maximum degree<=4k and target mass>=mass*n. The ambient deficit is >=Gamma n>gamma n, preserving the specified epsD margin. No actual-mu concentration or source independence is assumed.

## 5. Finite a, one history and the original operator

For every original dilation mathfrak a>=8n^4, t<n implies mathfrak a>2t and

    alpha beta-kappa-2/(mathfrak a+1)
      =n^-3-(beta+2)/(mathfrak a+1)>n^-3-4/(8n^4)>0.

Actual target increments have item radius <3 and included actual source items have radius >=rho>4. Every source therefore precedes the same actual target j, and the strict centre guard implies j's original eligibility capture. Multiplying by b<=1 preserves its actual row budget tau, without a theta^-1 factor.

Both supports lie in the original limiting corner cap, where lambda_mathfrak-a=alpha^-n lambda exactly. Energy scales by alpha^-2n and degrees/target mass by alpha^-n; since alpha^-n>=1, the full missing-measure bound cannot shrink. Thus gamma n is a valid deficit for every such finite history/dilation, not just the continuum envelope.

The original compact strict-level finite cover, actual disjoint increments, and mathfrak-a-dilated common-radius witnesses for the same arbitrary L1 input are unchanged. At fixed n,t,tau only already-valid scalar inequalities are sent to mathfrak a->infinity. Inner regularity then handles the whole level set. No bounded-support input, bounded-scale operator or selector/history limit is required. This closes C_n<=Phi_simple on all L>=2^27.

## 6. Complete scalar saving and R053 full-cost ordering

The actual buffer sigma/(1-sigma) is paid before selecting c=2^-16. Expanding the complete root gives

    Phi_simple<=T-c_raw n+nR(L),
    R(L)=1/n+exp(-sigma L)+2sqrt(L)exp(-sigma L/2).

At L0, sigma L0=2432>2^11 and sigma L0/2=1216>2^10. Global Taylor bounds of orders 2,16,16 give

    R(L0)<=2^-53+16!/2^176+2^15*16!/2^160<c_raw-2^-16.

The cross-term base coefficient is safe because 2sqrt(L0)=2^(29/2)<2^15. The three global envelopes decrease as L^-2, L^-16 and L^(-31/2), so the strict complete-cost saving holds on the whole tail. Neither A nor the positive square-root term is omitted.

The sharper full-cost lower envelope for the simple tuple is correct:

    B_simple>=T-n log(1-sigma_s)-4n/L
              >=T+sigma_s n-4n/L,
    Phi_simple>=T-n[gamma_s-sigma_s+4/L].

It follows directly by retaining -log(1-sigma_s) in the original cap lower bound; it is a lower bound for this formula cost, not for C_n. Independently reconstructing the R053 rational budget and applying its established complete upper envelope at Lr=2^32 gives

    R_raw(Lr)+4/Lr<c_raw,raw-(gamma_s-sigma_s).

The exact rational gap is positive (recorded in the receipt). Both positive losses decrease with L. Hence the complete raw root is strictly below the complete simple root for every L>=2^32, using each one's actual t and gamma, not only comparing c_raw or displayed 2^-16. Both operator thresholds apply there. On 2^27<=L<2^32, the simple proof supplies the earlier accepted domain; no unsupported ordering is claimed there.

## Actual checks and limits

`microcheck.py` was run locally with Python/Fraction, PID 84659; all 20 independent compressed certificates PASS. It reconstructs the simple tuple and R053 comparison budget from generic formulas, checks exact fractions/domain/degree branch, direct baseline probability, target rates, extended-q/individual-J guards, both exponential rates and base bounds, full scalar slack and full-cost comparison. It imports no author helper and does not run the author's 35 checks or any search/list traversal. `microcheck.json` records the actual source hashes, exact bounds and runtime.

The finite certificates corroborate base constants; all infinite tails and original-operator bridges are justified analytically above. No eventual O term remains for this tuple on L>=2^27. Only this private audit directory was written; no author source, Git, production run, old grid, n=7 or MATLAB was changed or used. No claim of minimal threshold, list completeness, parameter optimality or public priority is included.
