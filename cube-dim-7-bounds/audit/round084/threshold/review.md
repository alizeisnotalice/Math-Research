# R084 explicit threshold conditions — independent CONDITIONAL ACCEPT

The registered L0=2^50=1125899906842624 covers the explicit numerical
conditions of the frozen R083 s=9/2 proof, and these conditions persist
on the entire dimension tail. **Threshold conditions: ACCEPT. Final
operator acceptance is conditional on the separate independent full-chain
acceptance of frozen R083 proof
`6d2dff6fedf7cde05d4823830bca674235bd84be4608fc35ea129aa44d649c0e`.**

Under that mathematical input, the explicit statement is, for every
integer n with L=log n>=2^50,

    C_n<=F(A_n(t),B_n(t)-(1501/1000)n)
       <n log n-n log log n-(11/20)n,
    t=n/[(387/1000)L].

This is a new threshold identity, not a rewrite of R083's original
eventual theorem. The full local, inverse-law and original-operator
arguments are being audited separately and are not inferred here from
55 base gates. This review audits their explicit envelope coverage and
full-tail propagation. No claim of minimal threshold is made.

The original files were read only. This directory is independent
temporary staging, not a clone. No producer program, old 55-gate replay,
parameter search, shared MATLAB, commit or push was used. All new files
are confined to `audit/round084/threshold`.

## 1. Locked source and independently registered new calculation

Read in full: R084 registration.json, README.md, verify.py and result.json,
and the frozen R083 full proof. The producer's actual PID 29861 records
55 PASS gates at the single pre-registered L0. The receipt's script,
registration and frozen-proof hashes agree with the unchanged source
files. Its scope correctly separates base arithmetic from proof.

A separate REGISTRATION.json in this audit directory was saved before
the new bounded run. The independent check used N12=L^12/12!<=exp L,
single positive Taylor term T8(y)=y^8/8!<=exp y, 20-term atanh logs
with rigorous geometric tails, and exact polynomial coefficients.
Actual PID 30060 returned PASS on eight grouped checks. No producer
module was imported or executed, and no second threshold was tried.

Useful independent display-only bounds at the fixed base are

    N12>8.6628*10^171,
    rho_lower>4.4181*10^156,
    baseline mean<2.175*10^-140,
    baseline failure<3.480*10^-129,
    each count failure<4.618*10^-162,
    first-profile failure<7.397*10^-162,
    source-face bad probability<5.100078*10^-6<e,
    weighted-cap/T8 payment<2.672*10^-73<e,
    tau/T8 payment<1.798*10^-48<e,
    complete-root/T8 remainder<1.346*10^-104<1/1000,
    paid Gamma-gamma>0.00180462097.

All decision inequalities were exact rational comparisons. The original
receipt uses the stronger N16/S16 recipe; these new checks deliberately
use different weaker envelopes and do not relabel its run identity.

## 2. Tau-prefix domain: both the true and conservative quadratics work

For S=S16, the true coefficient expansion is

    yS'-(5/2)S
      >=-5/2-3y/2-y^2/4+y^3/12+y^4/16.

The R084 README instead uses -y^2/2, explicitly as a weaker lower
bound. This is valid, because y^2>=0. Neither this weaker form nor
the accurate -y^2/4 makes the R082 quartic repair incorrect.

Write u=y-4>=0. For the exact low-degree polynomial the shifted
coefficients are

    53/6, 33/2, 27/4, 13/12, 1/16.

For R084's deliberately weaker polynomial they are

    29/6, 29/2, 13/2, 13/12, 1/16.

All coefficients in both lists are positive. Every omitted term of
degree >=5 in yS'-(5/2)S also has positive coefficient. Consequently
both proofs establish positivity for every y>=4, not merely at y=4.
Thus L^(5/2)/S16(cT L) decreases for cT L>=4.

Here cT=379/2000000>0 and cT*2^50>=4 is independently checked.
The scalar half-power prefix has yS'-(1/2)S>0 for y>=2, since already
its constant and linear terms sum to -1/2+y/2>0 there. Its source
domain sigma*L0/2>=2 is also checked. The independent T8 argument
below provides an additional full-tail proof without any prefix
derivative reliance.

## 3. New s=9/2 constants and explicit active domain

The operative inputs are

    s=9/2, qmin=2/9, qmax=9999/10000,
    sigma=613/1000, chi=1583/2000, x=.387L,
    h=L^-2, delta=10^-6, rlo=1+delta,
    rhi=4499549/1000000, e=epsA=10^-5.

The threshold is checked against the R083 constants, not s=4 values:
|J-1|<=7/rho, log mean 49n/rho^2+14/rho, log variance
196n/rho^2, and new J1/J2 bounds. The new base check covers

    N12>=506250, x>=2, x/N12<1, x*N12^2>1,
    2h<delta, rlo-2h>1, h<=.001,
    1-h>0, (s-1)/s-h>0,
    epsilon<1, qmax/(x*N12^2)<=h/s,
    N12/x>=50 and >rhi,
    rho_lower=N12/(s*x)>4rhi and >=14,
    1+rhi<6<4rhi.

The h<=.001 guard is also implied directly by 2h<10^-6. The rho>=14
guard pays 7/rho<=1/2; the much stronger 7/rho<=e is checked as
well. Thus there is no missing coarse local-domain size condition.
The actual full matching q_±=(r±h)/s lies inside the source band, and
active R lies in [r-2h,r], exceeds one and is <s. The relation
g-(7/9)R=epsilon*q-h/s<=0 is paid everywhere in matching, not just
at the baseline. Inactive edge columns vanish; no active R>=1 estimate
is extended to them. Full-band R,a0,G stay positive.

After replacing n by N12, rho_lower grows as L^11, x/N12 decreases
as L^-11, and x*N12^2 grows as L^25. The epsilon relative-width
ratio [qmax/(x*N12^2)]/h decreases as L^-23. All other size and trim
guards improve as h decreases. They therefore hold for every L>=L0.

## 4. Five target events with e=10^-5 and all whole-matching guards

At the baseline, rho>=n/(s*x). The R083 log envelopes become

    Mlog=49*s^2*x^2/n+14*s*x/n,
    Vlog=196*s^2*x^2/n.

Using N12, the new base has Mlog<=e/2 and 4Vlog/e^2<=e/8.
Their adverse upper envelopes are positive sums of L^-10 and L^-11,
so they decrease over the tail. Each count failure bound
4/(e^2 n)<=e/8 decreases as L^-12. The target-cap free fraction
1+x/(2N12)-1/(2rhi)<=.9 decreases. Negative pinned faces already
lie in D, t>=50 is covered, and n>=300 implies the per-positive-face
failure <=.9^299<=e/8.

The new first-profile mean uses

    mu=7/2-(11/2)log(18/11)<.7914,

verified by an independent LOWER log enclosure. Its deterministic
target pin is <=9/(8rho)<=81*x/(16n). The base N12>=506250 pays
81/(16N12)<=10^-5, so total mean <=.79141x. The free variance is
<=6561*x^2/(256n). The remaining gap to .79143x is x/50000,
giving exactly 64072265625/n<=e/8. This bound decreases as L^-12.

These are exactly five events: log, two counts, D and first profile.
No independence between events is needed. Their total failure is
<=5e/8<e, so the source's fixed strict mass inequality gives
lambda(A)>(1503/1000)n. The independent log20 check confirms that
fixed mass gap and the deficit/price constants.

The NEW whole-matching log shift is 180xh=(3483/50)/L. The base
pays it by e, and it decreases as L^-1. Individual 7/rho<=e,
the count-calibration bound 2h+1/n<=3e/4, and first-profile transport
80h<=10^-5 are all paid at the base and decrease. The whole-matching
profile is therefore <=.79144x, and the exact rounding

    .79144/(1-e)<.79147

holds. The pointwise second-profile bound gives true centered source
variance <=[27/(16(1-e))]x<2x. This checks the new constants without
claiming a full inverse-law proof from arithmetic; that law remains a
separate R083 full-chain input.

## 5. Narrow mean gap and source-face retention on the full tail

The actual source-face mean envelope in R083 is .6+.79147x, whereas
chi=.7915. Thus the correct Chebyshev distance is

    u(x)=(3/100000)x-3/5,
    f(x)=2x/u(x)^2.

The new base has u>0 and f(x0)<5.100078*10^-6<e=10^-5, which pays
the .99999 retained fraction. This is substantially tighter than the
old s=4 retention comparison, and it is evaluated with the NEW gap.
Let a=3/100000 and p=3/5. Direct differentiation gives

    f'(x)=-2(a*x+p)/(a*x-p)^3<0 when a*x>p.

Since x increases, the distance stays positive and the same uniform
facewise retention bound holds for every larger L. The subtraction
.6 is never dropped, and no asymptotic normal or mixture-variance
argument is substituted. Face probabilities are mixed only after
the R083 law supplies uniform facewise control.

## 6. Count 2/9, weighted cap, radial degrees and tau

The whole-matching count uses a0/R>=1/s=2/9, giving

    nu=1/s-(1+1/s)e=22221/100000,
    eta=(1-qmax)/(2qmax)=1/19998.

The actual R083 weighted-cap exponent is
(nu-1/n)*eta*x/(1+e). Its N12 lower envelope has the exact form

    K12(L)=C[nu*L-12!*L^-11],
    C=(387/1000)*eta/(1+e)>0.

The base pays nu-1/N12>0 and K12>=4. Its derivative is
C[nu+11*12!*L^-12]>0. Therefore M/T8(K12(L)) decreases, and its
already tiny base value is <e on the whole tail. This retains the
arbitrary positive source pin's -1 cost. It is the WEIGHTED captured
cap deletion, with theta cancellation after actual matched dq/q
integration. No unweighted cap rate or mass is used.

The radial-width error 1/[1-(h/rlo)^2]<=1+e is checked and decreases.
Subject to R083's deterministic complete row/column proof, raw rows
and both degrees keep exactly their stated constants; no source or
target face sign is discarded by this threshold audit.

The tau rate is newly verified as

    sigma/2-chi*(1-sigma)=379/2000000>0.

The full true tau payment is bounded by

    rhi*L^(5/2)*exp(-cT L)
      <=rhi*L^(5/2)/T8(cT L).

The base uses sqrt(L0)=2^25 EXACTLY and yields <e. This bound is a
fixed constant times L^(-11/2), so it decreases throughout the tail.
The source's own S16 payment also propagates by the correct quartic
argument in Section 2. No constant square-root base value is treated
as a uniform bound on sqrt L.

## 7. Positive scalar roots and entire cover cost

For all L>=L0, directly

    exp(L/2)>=(L/2)^2/2=L^2/8>=x=.387L.

The final inequality already holds for L>=16. Thus log t>=L/2,
B>=nL/4>1.501n, and the SAME global original-selection threshold
tau=sqrt[A(B-1.501n)] is positive. This proves root positivity
directly; gamma>1 requires no alteration to the threshold.

A<=exp x and B<=nL apply on this explicit domain. For the latter,
the original B<=1+n log t and x>=2,n>=2 give n log x>=n/2>=1.
The complete positive-root remainder is bounded by

    1/N12(L)+1/T8(sigma L)
      +2sqrt L/T8(sigma L/2).

At L0 the last numerator is exactly 2^26. Its base is <1/1000.
The three tail envelopes have powers L^-12,L^-8 and L^(-15/2), so
the full remainder stays strictly below 1/1000. The source S16
half-power domain sigma*L0/2>=2 is independently covered as well.

The fixed exact inequalities log(1000/387)<19/20 and
1.501-.95-.55=.001 then yield the entire F<T-.55n conclusion,
including all positive A/cross/self terms, if the R083 original-cover
bridge is independently accepted. Fixed paid Gamma>1.501 and mass
>1.503 are checked with independent exact log/rational calculations.

## 8. Conditions for EVERY finite actual n and original history

For each actual integer n with log n>=L0, n>=N12(L)>=506250 and
t<n. Every a_alg>=8n^4 therefore has a_alg>2t and ccap>1. The
radius guards give target item radius<1+rhi<6 and included source
increment radius>=rho>4rhi>6. This is the enlarged target/source
separation required by R083, not the old 5 or 3 bound.

The finite-history displacement is proved for ACTUAL n algebraically:

    n^-3-(beta+2)/(a_alg+1)
      >n^-3-4/(8n^4)
      =(n-1/2)/n^4>0.

This follows from beta<2 and n>=2. It is not an assertion that a
rational expression evaluated at N12 transfers automatically to n.
The size guards are enough for the R083 same-history inclusion
argument, exact finite-envelope equality on A and P, c_a scaling,
genuine deficit, complete cover factor and terminal self-weight.
Those mathematical implications still belong to the R083 independent
full audit; their numerical prerequisites now hold uniformly here.

The general-L1/all-positive-radius extension introduces no extra
threshold requirement once that original finite-history bridge is
accepted: its arbitrary finite positive witnesses, compact subcovers,
normalization and inner-regularity steps are dimension independent.

## 9. Conclusion and integration boundary

ACCEPT the numerical threshold and full-tail conditions for the fixed
R083 branch at log n>=2^50. No unresolved numerical condition or
tail gap was found. The final C_n theorem is explicitly conditional
on separate full-chain acceptance of the frozen R083 proof. This
audit is not a substitute for that local/law/operator review.

The R084 README's deliberately weaker quadratic lower bound is valid;
both it and the accurate R082 quartic correction are explicitly
verified above. No original document or receipt is rewritten.

R083's original unquantified-eventual identity and the producer's
R084 single-threshold execution remain distinct. All source hashes
are recorded and rechecked at handoff. No smaller threshold, global
optimum, old-domain extension or change to another accepted branch
is claimed.
