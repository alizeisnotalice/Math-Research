# R089 independent explicit-threshold audit: ACCEPT L=log n >= 2^27

Decision: ACCEPT the explicit threshold extension for the exact R089
s=5 tuple, using the frozen independently accepted R088 arbitrary-tau
theorem. The R088 condition was pending when this audit was registered;
its main acceptance and full review were subsequently read and hash-
verified BEFORE the new audit run. DEPENDENCY_ACCEPTANCE.json records
that discharge. Acceptance of R088 alone was not treated as acceptance
of R089: all base guards and every infinite-tail extension are audited
below. Final ledger activation remains with the main agent.

For every integer n with natural L=log n>=2^27, the resulting bound is

    C_n <= H_tau=(1+A_n(t)/tau)(tau+B_n(t)-(3193/2000)n)
        <T_n-(9/10)n,
    tau=n exp[-(12532/25000)L],
    t=n/[(4987/10000)L], T_n=n log n-n log log n.        (1)

The scope is the original Lebesgue centered axis-parallel common-radius
cube maximal operator, arbitrary L1 inputs and every positive witness
radius. The bound is H_tau, not an optimized square-root F bound.
The failed 2^26 receipt remains failed and untouched. No minimum-
threshold or claim at 2^26 is made.

## 1. Frozen source and independent calculation identity

The accepted R088 generic proof has SHA256
75c56ce9cbd59efedc1a6863d9a125b14f30088ab85106cb98760d0507c35151.
The archived main acceptance hash is
e8fb2a35fcb8a74239a91adace420422237d9d74c9b88ce504174a130ead3b5b,
and its full independent review hash is
51dbe1a81f04287963d6b3ffa34cb1a460704caf9439c226e45a3dfa8644483f.
Those texts were read in full. Their earlier acceptance explicitly left
R089's threshold unreviewed; this report supplies that separate step.

Both fixed producer registrations, both scripts and both actual receipts
were read and locked in REGISTRATION.json. The 2^26 author result is
FAIL on exactly face_retention; the separately registered 2^27 result
is PASS on 56 gates. Neither author program was executed by this audit.
The old floating search, worker optima and inferred workload statistics
are not used as mathematical evidence.

The independent registration fixes ONE full base test at L0=2^27 and
ONE scalar historical face-failure check at 2^26, with no grid or
fallback. The different recipe is

    N(L)=L^12/12!<=exp L=n,
    8 positive atanh-log terms plus exact geometric tails,
    T8(y)=sum_(j=0)^8 y^j/j!<=exp y.

The final check.py actually ran once and returned PASS on 20 compressed
source/base/history groups. It parsed author statuses and hashes but
did not execute or reconstruct the author's 56-gate program. Its exact
receipt is verification.json; displays below are orientation only.

## 2. Exact tuple and all fixed primary inequalities

Use

    s=5, qmin=1/5, qmax=49999/50000,
    sigma=5013/10000, chi=201/200, d=12532/25000,
    x=(1-sigma)L=(4987/10000)L, h=L^-2,
    trim=1/10000000,
    rlo=1+trim, rhi=5*qmax-trim=4.9998999,
    eB=eN=eJ=eRad=eS=eT=eD=epsA=e=1/100000,
    Ret=999/1000, m=4007/2500, gamma=3193/2000,
    m_face=247/250, pin p=3/5,
    price=3479/5000, c=9/10,
    scalar_allowance=7/10000.

Let mu=4-6log(5/3), d_mu=(chi-mu)/4. The independent log8
enclosures supply fixed bounds mu_lo<=mu<=mu_hi and
dlo<=d_mu<=dhi, with

    mu_hi<0.935046257448,
    dlo>0.017488435638,
    (mu_hi+2dhi)/(1-e)<m_face<chi,
    chi-m_face=17/1000.

The R088 degree definitions remain exactly

    P0=1-2e, Ulo=P0(1-e)/(1+e), Uhi=1/Ulo,
    M=2(1+e)Uhi,
    b=2Ret(1-e)Ulo-2e, Gamma=m*b/M.

Independent arithmetic proves

    (1-e)log(rhi/rlo)>1.6092026009>m=1.6028>gamma,
    m*Ret>gamma,
    0<P0<1, M>2, b>0,
    (1-e)Gamma>1.6010050502>gamma=1.5965,
    -log(1-sigma)<0.695750566445<price,
    gamma-price-c=7/10000.                            (2)

All comparisons are exact Fraction comparisons. The relatively coarse
log8 lower bound for target mass is still adequate. The fixed new rates
are independently recomputed as

    d-chi(1-sigma)=173/2000000>0,
    sigma-d=1/50000>0,
    d>0.                                             (3)

The stale sqrt_L=8192 display field in the producer's second registration
is explicitly marked unused there. No square root of L is used in this
audit or in the new tau/H_tau formulas, so it creates no mathematical
guard at 2^27. Its source bytes are preserved.

## 3. Geometry and the full five-event probability budget

At L0, N12 is about 7.1347972e88. Use n>=N12 throughout; no value
of exp(2^27) is materialized. The new independent guards verify

    2h<trim, rlo-2h>1, 1-h>0, (s-1)qmin-h>0,
    h<.001, epsilon=1/(x*n^2)<1,
    epsilon*qmax<=h/s, x/n<1, x*n^2>1,
    t>=N12/x>2rhi,
    rho>=N12/(s*x)>max(20,4rhi), 1+rhi<6.

Thus the complete matching q_±=(r±h)/s is inside the source band,
R=sq-h lies in [r-2h,r] and is >1 on active sources. The exact identity
g-(1-1/s)R=epsilon*q-h/s<=0 gives g/R<=4/5 and ell>=3/5.
Inactive edge columns receive no R>=1 local assertion. R>0 on the
full band, beta<2, kappa>1, target/source bands are disjoint, and both
are below norm t. The accepted weak uniform local interfaces used are
ONLY |J-1|<9/rho, |d_q log J|<120/rho, J1<2/rho,
|d_q J1|<60/rho and J2<1/rho. No stronger unneeded local estimate
from the R088 author's discussion is required for this threshold.

Every target event is evaluated in the full ambient target-face law,
with n-1 independent free coordinates and the true deterministic pin:

* Baseline log: absolute mean <=81s^2*x^2/N12+18s*x/N12,
  variance <=324s^2*x^2/N12. The former is <e/2, and Chebyshev
  gives failure <=4 variance/e^2<e/8; the actual bound is about
  2.03452e-59.
* Both endpoint counts: each failure <=4/(e^2*N12)<e/8, about
  5.60633e-79. The true target pin contributes exactly one.
* Original target cap D: the free fraction is
  1+x/(2n)-1/(2r)<=1+x/(2N12)-1/(2rhi)<19/20.
  Since n>=N12>=512, failure <=(19/20)^511<e/8, about
  4.13774e-12. Negative pinned faces lie in D.
* First profile: n-free mean <=mu*x, pin <=10x/n, variance
  <=100x^2/n. The new guards 10/N12<dlo/2 and
  400/(dlo^2*N12)<e/8 pay the exact R088 target event
  sum J1<=(mu+d_mu)x. The failure bound is about 1.83306e-83.

There are FIVE failures: log, two counts, cap and first profile. Their
sum is <=5e/8<epsA, without any independence between the events.
Hence lambda(A)>=(1-e)n log(rhi/rlo)>m*n.

Individual J error 9/rho<e, product-log transport 240/L<e,
whole-matching count error 2h+1/N12<3e/4, first-profile transport
120h<dlo<=d_mu and radial width
1/[1-(h/rlo)^2]<=1+e are all required. The first four were explicitly
tested in check.py; the last follows immediately from h<=2^-54,
rlo>1, e=10^-5: (h/rlo)^2<2^-108<e/(1+e), since
2^108>100001.
This last rational inequality is far below its budget, and is not an
unjustified eventual statement. The affine ellpin' is negative and
the accepted count calibration therefore applies over the ENTIRE
matching interval, not only its baseline.

## 4. True source-face moments and the preserved 2^26 failure

For every captured signed source face before cap or soft conditioning,
the free joint density factors as product 1_capture/(2rho*J_i*ell_i).
The exact free means and raw seconds are J1_i/J_i, J2_i/J_i. The
pin is deterministic and <p=.6, with no variance contribution.
The new source-face guard in Section 2 gives mean <=p+m_face*x.
The pointwise J2 bound gives TRUE centered variance

    sum_(i!=j)[J2_i/J_i-(J1_i/J_i)^2]
       <=s*x/(1-e)<10x.                               (4)

No second-profile target event, mixture variance or actual-increment
product law is used. Consequently the uniform facewise failure is

    f(x)=10x/[(17/1000)x-3/5]^2.

At L0=2^27, its exact value is <0.000516955743<1-Ret=.001.
Only after this bound may source faces be mixed.

The separately locked 2^26 author receipt failed exactly this guard.
This audit independently recomputed ONLY its same scalar envelope and
obtained f(x26)>0.001033912575>.001, confirming the failure. It did
not overwrite either receipt or rerun the other 55 old checks. Failure
of this Chebyshev gate is not a proof that C_n lacks a bound at 2^26,
and 2^27 is not claimed to be the smallest threshold.

## 5. Weighted cap, new tau L^2 and full H_tau cross term L

Let

    nu=1/s-(1+1/s)e>0,
    eta=(1-qmax)/(2qmax)=1/99998>0,
    K(L)=(nu-1/N12(L))*eta*x/(1+e).

The accepted whole-matching count, arbitrary omitted source pin (-1),
pre-cap independence, bounded soft weight and exact matched dq/q give
the complete weighted deletion <=M exp[-K(L)]*k. The new independent
T8 check proves at L0

    M/T8(K0)<7.35809e-13<e.                            (5)

The target radial weight cancels dr/r, all pin signs are counted and
the complete source column is <=2k<M*k, exactly as in the accepted
R088 bridge. These are the same final kernel and weights <=1.

For the SAME global tau=n exp(-dL) fixed before selection,

    tau/k=rhi*L^2*exp[-(173/2000000)L].                 (6)

The independent upper prefix envelope is
rhi*L0^2/T8((173/2000000)L0)<1.09950e-11<e.
There is NO sqrt(L) in (6); the prefactor is L^2.

For Q=B_n-gamma*n>0 the entire cover cost is

    H_tau=Q+tau+A_n+A_n Q/tau.

A_n<=n exp(-sigma L), Q<=B_n<=nL yield the full normalized remainder

    1/n+exp(-dL)+exp(-sigma L)+L exp[-(sigma-d)L].       (7)

The cross prefactor is L, and the cross rate is sigma-d=1/50000;
it is not the row rate in (6) or the old square-root exponent.
The independent complete prefix bound at L0 is

    1/N12+1/T8(dL)+1/T8(sigma L)+L/T8((sigma-d)L)
       <2.00132e-15<7/10000.                           (8)

The exact logarithmic price and slack in (2) therefore give the strict
H_tau<T_n-.9n once the whole-tail extensions below are established.
No implication from H_tau to a smaller F expression is used.

## 6. Every base guard extends to the entire tail

This section covers ALL L>=2^27 analytically; it is not a dimension
sample. N12=L^12/12! grows, x is proportional to L, and the lower
source-radius bound N12/(s*x) grows as L^11.

* epsilon*qmax/h is bounded by a constant times L^-23; epsilon,
  x/N12 and individual-J errors decrease. All geometry/radius guards
  from Section 3 persist. In particular beta<2 and kappa>1 persist.
* Baseline mean is a positive sum of L^-10 and L^-11, and its
  variance/failure envelope is proportional to L^-10. Count and
  first-profile bad probabilities are constant times L^-12. Their
  pin budgets improve, and target-cap fraction x/N12 decreases.
  The fixed bound (19/20)^511 remains valid since N12 increases.
* 240/L, 120h, 2h+1/N12 and the radial-width error decrease.
  Thus all full-matching product/count/profile bounds persist.
* For a=17/1000 and p=3/5, x0>p/a and
  d[10x/(ax-p)^2]/dx=-10(ax+p)/(ax-p)^3<0. Hence the face
  retention bound persists for every larger x, not only the base.
* Writing x=a_x L, the weighted-cap exponent is
  eta*a_x/(1+e) * [nu L-12! L^-11]. Its derivative is positive.
  T8 has positive coefficients, so M/T8(K(L)) decreases.

The tau and scalar-prefix derivatives can be checked directly with
the NEW degree-eight prefix; no author polynomial is imported. For y>=4,

    y*T8'(y)-2*T8(y)>=-2-y+y^3/6>0.                   (9)

The lower polynomial is 14/3 at y=4 and has derivative -1+y^2/2>0
thereafter; all omitted coefficients are nonnegative. Therefore
L^2/T8(cL) decreases whenever cL>=4, proving (6)'s payment for the
entire tail. At this tuple (173/2000000)L0>4 is checked exactly.

For y>=2,

    y*T8'(y)-T8(y)>=-1+y^2/2>0.                       (10)

Thus L/T8((sigma-d)L) decreases; its base argument exceeds 2.
Both pure reciprocal prefixes decrease, and 1/N12 decreases. This
proves the ENTIRE remainder (8) stays below its allowance. The proof
uses the actual L cross term and L^2 row payment, not an obsolete
L^(3/2) square-root envelope.

All fixed mass, mean-rate, degree and price inequalities in (2) are
independent of L and remain unchanged. Every finite base condition
used by the accepted R088 generic proof now has a verified full-tail
extension.

## 7. Positive Q, same finite history and all input/radius scope

For L>=L0>=8, the positive exponential series gives
exp(L/2)>=(L/2)^2/2>=L>=x, so log t>=L/2. The original cap identity
then gives B_n>=nL/4>gamma*n, as checked at the base and increasing.
Also B_n<=1+n log t<=nL: x>=2 and n>=512 imply n log x>1.
Hence the Q and every scalar upper direction used in Section 5 are
valid on the entire tail, even though gamma=1.5965>1.

For each fixed n, this positive tau is fixed before any input,
witness family, a_alg or selection history. R088's original finite
theorem accepts it. It is used consistently in actual eligibility,
actual weighted row, ambient row subtraction, captured old-weight
payment and the entire cover factor.

For EVERY a_alg>=8n^4 the checked x>=2 and n>=512 ensure
a_alg>2t and ccap>1. An actual target meeting A has radius
r_j<alpha+rhi<6. Every source has r_i>=rho>20>6 and precedes
that SAME target. Since beta<2,

    n^-3-(beta+2)/(a_alg+1)
      >n^-3-4/(8n^4)=(n-1/2)/n^4>0.

The actual center-distance estimate with r_j<r_i therefore applies
uniformly, giving actual weighted target row<=arad*tau<=tau.
No old square-root-generated history is reused. Both weights remain
bounded by one and genuine item weights are unchanged.

The finite-cap containment margin is positive, giving exact envelope
scaling c_a*lambda on disjoint A and P. Energy/degree/mass scale
as c_a^2,c_a,c_a. After (5), (6) the accepted energy-difference proof
therefore gives a TOTAL finite-cap deficit >=Gamma*n>gamma*n for
every original finite history, uniformly in a_alg. Taking only the
scalar a_alg limit of the full finite inequality yields precisely H_tau.
The terminal self term and complete cover factor are retained until
this scalar limit. No selected families must converge.

R088's independently accepted strict-level compact finite witness
cover, arbitrary positive terminal normalization, true selected
witnesses and inner regularity recover every L1 input and every positive
radius. Their hypotheses do not acquire a new scale or input restriction
from this explicit threshold. Combining this frozen bridge with Sections
2--6 proves (1), subject to final main-agent activation of this review.

## 8. Full old-cost comparison and the actual contribution of this branch

The separate main-agent primary_comparison registration, source and
actual PID 33608 receipt were read, not executed again. Its source hash
is 61e956d8a6a15a1816c60ac025a8cc576cf4775b30cb9a04685ada56854a71ca
and receipt hash is
174d7bb4d22c382717241894128172311a6bdac8285dbb8f830049884ab6878a.
Three independent log24 comparisons passed.

The original I_n(t)<=2t/(n-1) and nonnegative B core give the old
R081 complete-cost lower bound

    Phi81>T_n-[c81+5/L]n,
    c81=13/10+log(2/5)<.384.

Its 5/L factor is appropriate to x_old=.4L. R082 makes this exact
old Phi81 valid at L>=2^24. Therefore on the TRUE common domain
L>=2^27, (1) and the separate exact comparison imply

    Phi81-H89>.5159n.

This compares the complete old upper-bound function with the complete
new H89, not only two display constants, and it is not a C_n lower
bound. Neither branch receives a new validity domain from the comparison.

The same receipt checks the accurate prices

    c88=3/2+log(11/20) > c89=3193/2000+log(4987/10000).

Consequently this review does not claim that H89 eventually dominates
H88 over its full domain. R089 supplies an EXPLICIT .9 threshold for a
different tuple, whereas R088's fixed instance retains its own stronger
priced coefficient and unknown eventual threshold. Accepted branches
continue to be taken in the minimum on their own actual domains.

## 9. Preservation, reproducibility and decision boundaries

The 2^26 failure, its registration and all source receipts are immutable
inputs to this audit. The 2^27 registration is a separate source identity.
No author 56-gate program, floating search, old baseline or dimension
sweep was executed. This report does not retrospectively modify R088's
eventual theorem or claim that its own earlier audit certified R089.

To reproduce, run check.py with --source-dir pointing to the archived
R089 tau_compute folder and --generic-proof, --main-acceptance,
--full-review pointing to the frozen R088 files. The code verifies their
hashes before its ONE independent N12/log8/T8 run. Its exact inputs,
recipe and registration identity are retained in the accompanying files.

No gap was found in the threshold extension. ACCEPT L=log n>=2^27
for the fixed R089 H_tau tuple, retaining the source 2^26 failure.
No smallest-threshold, global-optimum or F assertion is included.
The paper and previous accepted results
were not edited. Only this independent clone's owned audit/round089/
threshold files were written; no shared MATLAB, commit or push was used.
