# R091 fixed c=1.1 at log n>=2^33 — independent ACCEPT

For EVERY integer dimension n with REAL L=log n>=2^33=8589934592,
the original Lebesgue centered axis-parallel common-radius cube maximal
operator, for arbitrary L1 inputs and all positive witness radii, satisfies

    C_n<=H_tau=(1+A_n(t)/tau)(tau+B_n(t)-gamma*n)
       <n log n-n log log n-(11/10)n,
    tau=n exp(-dL), t=n/[(1-sigma)L].

The fixed rational tuple is recorded below. The conclusion is the
COMPLETE H_tau bound, not an optimized square-root F bound. L need
not be an integer. This is a separate explicit-threshold theorem;
neither generic R090 acceptance nor the producer's 58 base PASS
gates alone establishes it. No claim of minimal threshold is made.

The underlying R090 proof is frozen at SHA256
`12605f79fc17c83516fd4423c69fd58b1340a2a780588490591795ee3ec59631`.
Its formal main acceptance, REPORT and archived independent full review
were read and bound in this audit's receipt. That generic proof supplies
the mathematical estimates; this review checks the NEW tuple, all
numerical prerequisites, their full-tail propagation and their original
finite-history use.

This is independent temporary staging, not a Git clone. Original sources,
old rounds, paper and Git state were read only. All new files are inside
`audit/round091/threshold`. No shared MATLAB, commit, push, author58
replay, worker15 execution, new scan or second threshold was used.

## 1. Candidate provenance and separate bounded arithmetic

Read in full: exact_report.md, exact_registration.json, exact_verify.py,
the actual result at L0=2^33 and exact_manifest.json. The registered
summary hash was checked without executing or repeating its workers.
The fixed exact candidate was created by explicitly documented decimal
quantization of model fields. It is a NEW rational tuple, not a recovered
prior exact certificate. No conclusion about an old numerical optimum
is used here.

The primary and secondary values are

    s=42.57987173=4257987173/10^8,
    chi=12.90458926=645229463/50000000,
    sigma=.92808239=92808239/10^8,
    d=.92807606=46403803/50000000,
    qmax=.99999352=12499919/12500000,
    gamma=3.73319131=373319131/10^8,
    m=3.74228319=374228319/10^8,
    Ret=.99939262=49969631/50000000,
    e=10^-5, m_face=12.77456974=638728487/50000000,
    trim=10^-7, price=6581/2500=2.6324, c=11/10.

The registration's `delta` is the tau exponent d, NOT the geometric
trim. All these constants stay fixed on the dimension tail. Its source
receipt PID 33459 records 58 PASS gates at 2^33; the planned fallback
2^34 was unnecessary and was not evaluated here.

This audit separately registered one new finite arithmetic check before
execution. `check.py` actually ran with local Python Fraction arithmetic,
PID 34233. It uses the weaker N12=L^12/12!<=exp L, single positive
Taylor term T8(y)=y^8/8!<=exp y and BINARY RANGE-REDUCED log20
intervals. Range reduction is useful because log(s*qmax) has a large
argument; all reduced log arguments are in [1,2). Seven independent
grouped checks pass. It imports or executes no producer or worker code.

Display-only independent values include

    N12>3.3693*10^110, rho_lower>1.2808*10^100,
    D_s=4s^2+6s+1 approximately 7508.66113655,
    d_mu_lower>0.1300827432887,
    source mean upper<12.644551<m_face,
    complete matching log shift<1.258*10^-7<e,
    baseline failure<2.383*10^-75,
    each count failure<1.188*10^-100,
    target-cap/T8 probability<2.856*10^-10<e/8,
    first-profile failure<2.307*10^-102,
    source-face failure<.000347214<1-Ret=.00060738,
    weighted-cap/T8 payment<3.395*10^-9<e,
    free-tau/T8 payment<1.144*10^-13<e,
    full-H/T8 remainder<4.533*10^-24<.00079131,
    paid Gamma-gamma>.00637009127.

Exact rational comparisons and all source identities are preserved in
receipt.json. The producer's actual run and this independent run remain
separate, and neither finite calculation substitutes for the following
all-tail argument.

## 2. Strict fixed-tuple theorem conditions

Independent rigorous logs enclose mu(s)=s-1-(s+1)log(2s/(s+1))
and its actual fixed d_mu=(chi-mu)/4. They prove

    s>=4, 0<sigma<1, 1/s<qmax<1,
    0<mu<chi, chi(1-sigma)<d<sigma,
    0<gamma<log(s*qmax), c<gamma+log(1-sigma),
    0<e<min(1/12,1/[2(s+1)]).

The R090 definitions of P0,Ulo,Uhi,M,b,Gamma are used exactly with
THIS tuple. Fixed constants satisfy

    (1-e)log(rhi/rlo)>m>gamma, Ret*m>gamma,
    (1-e)Gamma>gamma,
    (mu+2d_mu)/(1-e)<m_face<chi,
    -log(1-sigma)<2.6324,
    gamma-price-c=79131/10^8=.00079131>0.

Here rlo=1+10^-7 and rhi=s*qmax-10^-7. If [mu_lo,mu_hi] is
the independent log interval, d_lo=(chi-mu_hi)/4>0 is used for
pin, probability and transport payments. The mean upper bound uses
the exact relation mu+2d_mu=(mu+chi)/2. No rounded decimal replaces
the true fixed constant in the mathematical good-target event.

## 3. Arbitrary-s local domain, central deletion and whole matching

Put x=(1-sigma)L, h=L^-2 and epsilon=1/(x n^2). At L0 the
weaker N12 check pays

    N12>=10001, x>=2, x/N12<1, x*N12^2>1,
    0<h<=.001, 2h<trim, rlo-2h>1,
    1-h>0, (s-1)/s-h>0,
    1/(x*N12^2)<min(1,h),
    qmax/(x*N12^2)<=h/s,
    N12/x>=2rhi,
    rho_lower=N12/(s*x)>4(s+1),
    1+rhi<s+1, 2s/rho_lower<1/2 and <=e.

Thus every active q match has q>=1/s, 1<R<=s,
|w_i|<=r<=min(s,R+2h), rho>4(s+1), a0>0 and
g-(1-1/s)R=epsilon*q-h/s<=0. The complete q_±=(r±h)/s
interval stays inside [1/s,qmax]. Both annuli lie below t and are
disjoint. Full-band R,a0,G remain positive; unmatched edge sources
have zero columns and receive no active R>=1 estimate.

These are the NEW arbitrary-s hypotheses, not an s<=5 domain. g/R
may exceed .8. R090's fresh lemma uses ell>=1/2, Q<=log2<1,
|J-1|<=2s/rho, |d log J|<=D_s/rho, J1,J2<=s/rho and
|J1'|<=5s/rho. It covers central deletion |w|>R with the essential
-m_w term and AC integration across onsets. No old 9/rho, 120/rho
or target-radius cutoff at 5 is silently used.

For every larger REAL L, N12/x grows as L^11, x/N12 decreases
as L^-11, x*N12^2 grows as L^25 and
[qmax/(x*N12^2)]/h decreases as L^-23. Each trim, epsilon and
individual-J error improves, while rho and t grow. Thus the complete
uniform local domain is paid on the ENTIRE tail, not by dimension
sampling or a single-point extrapolation.

## 4. All five target events with the new s^4 moments

The accepted local full-target identities give baseline log-sum bounds

    absolute mean<=4s^4*x^2/n+4s^2*x/n,
    variance<=16s^4*x^2/n.

At L0 their N12 envelopes pay mean<=e/2 and four times variance/e^2
<=e/8. The adverse powers are L^-10,L^-11 and decrease thereafter.
Each endpoint count failure <=4/(e^2 n)<=e/8 decreases as L^-12.
These use s approximately 42.58, not the bounded-s constants.

The positive-face target-cap base is 1-1/(4s), not .95. t>=2rhi
gives the required fraction, and negative pinned faces already lie in D.
For actual n>=N12>=10001,

    (1-1/(4s))^(n-1)
      <=exp[-10000/(4s)]<=1/T8(10000/(4s))<e/8.

This fixed conservative bound holds for every larger n. The producer's
S24 estimate is stronger but is not needed for independent acceptance.

The true first-profile pin is <=s^2*x/n. The base pays
s^2/N12<=d_lo/2, so total mean <=(mu+d_mu/2)x. Free variance
is <=s^4*x^2/n and the excluded endpoint is (mu+d_mu)x. Its bad
probability is therefore <=4s^4/(d_mu^2 n)<=4s^4/(d_lo^2 N12)
<=e/8, decreasing as L^-12. All five events—log, two counts,
D and J1—are separately paid. Their mutual independence is not
assumed, and their union leaves target mass >m*n.

The COMPLETE matching q interval transports the log product by
2D_s*x*h and J1 by 10s*x*h. The base pays 2D_s*x*h<=e and
10s*h<=d_lo, improving as L^-1 and L^-2. Count calibration
2h+1/n<=3e/4 and the radial width bound also improve. Thus P0,
all individual J_i denominators, the whole-interval count ratio,
Ulo/Uhi and sum J1<=(mu+2d_mu)x are valid throughout matching.
No baseline-only or q-grid substitution is made.

## 5. True source variance, pin 1 and Ret on the whole tail

For each actual pre-cap SOURCE inverse face, free means and raw
second moments are J1_i/J_i and J2_i/J_i in the accepted PRODUCT
law. The source pin is deterministic and <=log2<1. The transported
mean is <1+m_face*x and the TRUE centered variance is bounded by

    sum_(i!=j)[J2_i/J_i-(J1_i/J_i)^2]
       <=s^2*x/(1-e)<2s^2*x.

No 10x or other s<=5 variance factor, pin variance, mixture variance
or actual-increment law replaces this expression. For the fixed
a_face=chi-m_face=.13001952>0, put u=a_face*x-1. The exact
facewise Chebyshev envelope is

    f(x)=2s^2*x/u^2.

The new base has u>0 and f(x0)<.000347214<.00060738=1-Ret.
Direct differentiation gives

    f'(x)=-2s^2*(a_face*x+1)/(a_face*x-1)^3<0.

Therefore this uniform retained fraction persists over every larger
L, target, matching q, omitted coordinate, face and sign. The pin 1
is retained in the gap. Mixing occurs only after the common facewise
bound, and soft retention is applied before cap conditioning.

## 6. Weighted cap and NEW free-tau payment

The s-dependent error condition yields
nu=1/s-(1+1/s)e>1/(2s)>0. The whole-matching count, full captured
lower strips and arbitrary positive source pin's -1 cost give the
accepted WEIGHTED loss bound

    caploss/k<=M exp[-(nu-1/n)*eta*x/(1+e)],
    eta=(1-qmax)/(2qmax)>0.

It uses bsoft<=theta/G pointwise, all source face coefficients and
the actual matched dq/q interval; theta cancels against k exactly.
An unweighted cap estimate would be inappropriate and is not used.

For n>=N12, the exponent is at least

    K12(L)=C[nu*L-12!*L^-11],
    C=(1-sigma)*eta/(1+e)>0.

The base has nu-1/N12>0 and K12>=4. Its derivative is
C[nu+11*12!*L^-12]>0, so M/T8(K12(L)) decreases. Its base is
<e, paying the complete source-cap loss on the entire tail. The
radial error decreases as h shrinks; complete ambient row/column
bounds and all target signs remain those of accepted R090.

For THIS fixed tau exponent d, the exact positive row-payment rate
is

    cT=d-chi(1-sigma)=.0000088423891314>0.

The SAME selection threshold tau=n exp(-dL) obeys

    tau/k=rhi L^2 exp(-cT L)<=rhi L^2/T8(cT L).

The new base has cT*L0>=4 and payment<e. The T8 envelope is a
constant times L^-6, so it decreases for every larger L. This is
L^2, not the old square-root L^(5/2) payment.

The producer's S24 alternative is independently valid: for y>=4,
yS24'-2S24>=-2-y+y^3/6>0. The lower polynomial is positive at
4 and has derivative -1+y^2/2>0 there and above. Thus its entire
L^2/S24(cT L) envelope decreases too. No hidden asymptotic claim
is needed after the new base domain is checked.

## 7. Entire H cost, positive Q and every finite actual history

For all L>=L0, exp(L/2)>=L^2/8>=x, with the last inequality
already valid for L>=16. Thus log t>=L/2,
B>=nL/4>gamma*n and Q=B-gamma*n>0. x>=2,n>=10001 also imply
B<=1+n log t<=nL. A<=exp x holds directly.

The complete cost is EXACTLY

    H_tau=Q+tau+A+A Q/tau.

Its normalized positive remainder is bounded by

    1/n+exp(-dL)+exp(-sigma L)+L exp[-(sigma-d)L].

Here sigma-d=.00000633>0. The independent T8 base pays the
complete upper envelope

    1/N12+1/T8(dL)+1/T8(sigma L)+L/T8((sigma-d)L)
       <79131/10^8.

The four all-tail powers are L^-12,L^-8,L^-8 and L^-7. Thus the
FULL H remainder stays strictly below the fixed allowance. For the
source S24 cross envelope, yS24'-S24>=-1+y^2/2>0 for y>=2,
and that base domain is covered. The cross term is proportional to
L, not sqrt L, and no optimization back to F is performed.

For every actual n on the domain, n>=N12>=10001 and t<n. For
every original a_alg>=8n^4, ccap>1. A genuine target item meeting
A has radius<alpha+rhi<s+1; an included source has radius>=rho
>4(s+1)>s+1. It precedes that SAME target in the original history.
The displacement is established for ACTUAL n algebraically:

    n^-3-(beta+2)/(a_alg+1)
      >n^-3-4/(8n^4)=(n-1/2)/n^4>0.

The n lower envelope is used only for n>=2 and beta<2, not to
infer that a rational expression sampled at N transfers to n. All
included pairs genuinely capture the target in its own eligibility
test. The same tau is fixed before any f, finite a, witnesses or
history; bounded weights give actual row<=arad*tau<=tau.

The finite envelope equals c_a*lambda on disjoint A,P, energy,
degree and mass scale by c_a^2,c_a,c_a, and the two missing-mass
terms are paid from the TOTAL finite-cap deficit. The accepted
R090 identity therefore gives Delta_a>=Gamma*n>gamma*n for every
original finite history. gamma>1 is harmless because Q positivity
was proved directly. The entire cover factor and terminal self-term
are retained in taking only the uniform NUMERIC a limit.

Finally, the general-L1 extension uses actual positive witnesses,
open shrunk covers of compact subsets, selected true dilates,
normalization at every positive terminal radius and inner regularity.
It needs no compact support, radius truncation, scale lattice,
actual-increment independence or convergence of histories. These
mathematical steps are supplied by the accepted frozen R090 proof;
all their numerical prerequisites are now explicitly valid.

## 8. Exact conclusion and independent acceptance object

The fixed price inequality -log(1-sigma)<6581/2500 and the identity

    gamma-6581/2500-11/10=79131/10^8

pay the ENTIRE H_tau remainder. Therefore the claimed strict bound
C_n<T_n-1.1n holds for every integer n with log n>=2^33.

ACCEPT this new explicit-threshold theorem. No missing target event,
local-size condition, inverse-face or finite-history numerical guard
was found. The proof uses actual s-dependent constants throughout;
generic R090 acceptance has not been substituted for the explicit
threshold argument.

All original identities, the new rational-candidate provenance, the
producer's PID 33459 and this audit's PID 34233 remain separate.
Neither 2^34 nor any worker or nearby tuple was evaluated by this
audit. Source hashes—including R090 PRIMARY_ACCEPTANCE—are recorded
and rechecked at handoff. No old result, paper or source is changed.
