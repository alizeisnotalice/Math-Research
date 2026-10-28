# R087 fixed tuple at log n>=2^80 — independent ACCEPT

For every integer dimension n with the REAL number L=log n>=2^80,
the original centered axis-parallel Lebesgue cube maximal constant obeys

    C_n<=F(A_n(t),B_n(t)-(619/400)n)
       <n log n-n log log n-(111/200)n,
    t=n/[(37073/100000)L].

Thus the displayed coefficient is exactly .555. No restriction that
L itself be an integer is imposed. No minimal-threshold or global
optimality claim is made.

The mathematical framework is the independently and main-agent accepted
R085 fixed-tuple theorem, frozen PROOF.md SHA256
`ded160a902cdd543a7be9ebbb9c1f4d076c33f866d8aa96e9481f92d69bb41fa`.
Its accepted full review and independently rederived local lemma were
read in full. This review verifies the FIXED R087 tuple against that
framework and supplies an explicit all-tail threshold, not an inference
from a base PASS or a new derivation of the already accepted generic
family theorem.

This is an independent temporary staging directory, NOT a clone. All
original materials were read only; only this owned
`audit/round087/threshold` directory was written. No original or old-round
file, paper, Git state or shared MATLAB was modified. No commit/push,
parameter sweep, other threshold or author baseline execution occurred.

## 1. Frozen inputs, separate registration and actual new computation

Read in full: the two R087 registrations, primary.py, threshold.py,
README.md, both actual arithmetic receipts and the supplied manifest;
R085 PROOF.md, audit/full_review.md and audit/local_uniform.md. The
historical producer PIDs are 31183 (eight primary checks) and 31252
(56 secondary/threshold checks). Their primary receipt, registrations,
script and frozen-proof hash links agree with the current unchanged
files. Neither program was executed or imported by this audit.

The authoritative fixed choices are

    s=47/10, qmin=10/47, qmax=99999/100000,
    chi=42431/50000, sigma=62927/100000,
    gamma=619/400, c=111/200,
    trim=delta=e=1/10000000, Ret=999999/1000000,
    m=154753/100000, m_face=424309/500000, p=3/5,
    price=99229/100000, scalar slack=21/100000,
    L0=2^80, sqrt(L0)=2^40, h=L^-2.

All of them remain fixed as n grows. The definition
d_mu=(chi-mu(s))/4 is retained as a real constant, rather than
replaced by a rounded display.

This audit saved its own REGISTRATION.json before one new bounded run.
`check.py` actually ran with local Python Fraction arithmetic, PID 31660,
using weaker N12(L)=L^12/12!<=exp L, 32 positive atanh-log terms
plus exact geometric tails, and the single positive exponential-series
term T8(y)=y^8/8!<=exp y. Eight independently grouped checks passed.
These differ from the author 96-term/N16/S16 recipe; neither the
primary eight nor the threshold 56 checks was replayed.

Independent display-only results include

    N12(L0)>2.0345*10^280, rho_lower>9.6583*10^255,
    d_mu_lower>2.00807052158*10^-6,
    transported mean/(1-e)<0.848616068721<0.848618,
    baseline failure<2.827*10^-215,
    each count failure<1.967*10^-266,
    target-cap failure<=1.293358*10^-9<e/8,
    first-profile failure<4.876*10^-267,
    true source-face bad probability<5.578053*10^-12<1-Ret,
    weighted-cap/T8 payment<3.020*10^-137<e,
    tau/T8 payment<3.093*10^-91<e,
    complete-root/T8 remainder<2.024*10^-172<21/100000,
    paid Gamma-gamma>2.65954366*10^-5.

Every decision inequality was evaluated as an exact Fraction. The
receipts record source hashes and actual execution identities. The
following analysis, rather than these finite displays, supplies the
entire-tail and operator conclusion.

## 2. Strict parameter theorem and unrounded d_mu

For s=47/10 define Q=log(2s/(s+1)) and
mu=s-1-(s+1)Q. Independent log32 intervals give

    mu_lower<=mu<=mu_upper<chi,
    4<=s<=5, 0<sigma<1, 1/s<qmax<1,
    cT=sigma/2-chi(1-sigma)=130537/5000000000>0,
    0<gamma<log(s*qmax),
    c<gamma+log(1-sigma).

Thus the primary constraints of R085 hold with strict positive slack.
The trim and secondary choices satisfy

    rlo=1+delta<rhi=s*qmax-delta=4.6999529,
    (1-e)log(rhi/rlo)>m>gamma,
    Ret*m>gamma, 0<e<1/12,
    (1-e)Gamma>gamma,
    -log(1-sigma)<price,
    gamma-price-c=21/100000.

P0,Ulo,Uhi,M,b and Gamma are exactly the R085 quantities, with this
fixed e, Ret and m; no older round's deficit constants are imported.

Let [mu_lo,mu_hi] be the independently calculated interval. Define

    d_lo=(chi-mu_hi)/4, d_hi=(chi-mu_lo)/4.

Then 0<d_lo<=d_mu<=d_hi. Lower d_lo is used wherever a positive
pin, probability or transport budget must be paid. For the source
face mean there is the exact identity

    mu+2d_mu=(mu+chi)/2.

Its upper bound (mu_hi+chi)/2 preserves correlation between mu and
d_mu and gives

    (mu+2d_mu)/(1-e)<m_face<chi.

No numerical interval endpoint is silently substituted for the actual
constant defining the good-target event. All R085 fixed-tuple hypotheses
are therefore met by one fixed tuple, independently of n.

## 3. Uniform local domain and full matching for all larger L

At L0 the weaker N12 check proves

    N12>=400, x=(1-sigma)L>=2,
    x/N12<1, x*N12^2>1,
    2h<delta, h<=.001, rlo-2h>1,
    1-h>0, (s-1)/s-h>0,
    1/(x*N12^2)<1, qmax/(x*N12^2)<=h/s,
    N12/x>=2rhi,
    rho_lower=N12/(s*x)>20 and >4rhi,
    1<rhi<5, 9/rho_lower<=e and <=1/2.

These explicitly cover the accepted uniform local lemma's h<=.001,
rho>=20 and fixed-target derivative domain. At each target radius,
q_±=(r±h)/s stays inside [1/s,qmax]; active R lies in [r-2h,r],
exceeds one and is <5. The identity

    g-(1-1/s)R=epsilon*q-h/s<=0

gives the actual relative gap uniformly on every matching q. Thus
g/R<=1-1/s<=4/5 and ell>=3/5. The independent accepted local
formulas and AC derivative bounds, including central deletion |w|>R,
apply across the complete interval. Inactive edge sources have zero
columns and receive no inappropriate R>=1 assertion; full-band
R,a0,G remain positive. a0<qmax<1<rlo puts the true target pin
outside the endpoint count threshold automatically.

These statements persist over EVERY real L>=L0. N12/x grows as
L^11, x/N12 falls as L^-11, and x*N12^2 grows as L^25. The
relative epsilon-width envelope [qmax/(x*N12^2)]/h falls as L^-23.
h and all trim errors decrease, while source radii and t grow. No
discrete L grid or sampled tail is used.

## 4. All five full-target events and whole-matching transport

The accepted local value bound is 9/rho, with baseline free mean
and variance -81/rho^2<=E log J<=0 and Var(log J)<=324/rho^2,
and deterministic pin contribution <=18/rho. Consequently

    Mlog=81s^2x^2/n+18s*x/n,
    Vlog=324s^2x^2/n

bound the absolute log-sum mean and variance. At L0, the N12
envelopes pay Mlog<=e/2 and 4Vlog/e^2<=e/8. They decrease as
positive sums of L^-10,L^-11 thereafter. Each endpoint count
failure <=4/(e^2 n)<=e/8 has an L^-12 envelope.

Negative target pinned faces belong to D because t>rhi. On positive
faces t>=2rhi gives free complement fraction

    1+1/(2t)-1/(2r)<=1-1/(4rhi)<19/20.

Since n>=N12>=400, its failure is at most (19/20)^399<=e/8.
This bound is per signed face, with no hidden one-half factor.

Full-target Fubini gives free first-profile mean <=mu*x. The actual
pin contributes <=2/rho<=10x/n. The new base pays 10/N12<=d_lo/2,
so the total mean is <=(mu+d_mu/2)x. Free variance is
<=100x^2/n. The remaining d_mu*x/2 gap to the exact good-target
threshold (mu+d_mu)x is paid by

    400/(d_mu^2 n)<=400/(d_lo^2 N12)<=e/8.

This decreases as L^-12 because d_lo and d_mu are FIXED. These
are precisely five target conditions: baseline log, two endpoint
counts, D and first profile. The five events need not be mutually
independent. Their total failure <=5e/8<e gives the strict mass
bound lambda(A)>m*n with the fixed logarithmic-band inequality.
These are full-target ambient laws, never actual-increment laws.

For the complete matching interval the accepted local derivatives
give product-log drift <=240xh<=240/L and first-profile drift
<=120xh. The base checks 240/L<=e, individual 9/rho<=e,
2h+1/n<=3e/4 and 120h<=d_lo. All these improve on the tail.
Thus product J_i lies in [P0,1/P0], the endpoint-count ratio lies
in [1-e,1+e] throughout matching, and summing every captured sign
and omitted coordinate gives the stated Ulo*n<=U<=Uhi*n.
The complete-matching first-profile sum is <=(mu+2d_mu)x.
There is no baseline-only count or profile substitution.

## 5. TRUE source-face Chebyshev and fixed retention on the full tail

On each actual pre-cap inverse SOURCE face the free coordinates have
the accepted product density and exact first/raw-second moments
J1_i/J_i and J2_i/J_i. The source pin is deterministic with loss
<p=.6. The mean envelope is

    E V<=p+m_face*x,
    m_face=424309/500000=.848618,

and the TRUE centered variance, after subtracting squared means,
is bounded by s*x/(1-e)<10x. The pin adds no variance and no
variance of the face mixture is inserted. The source rate gap is
EXACTLY

    a=chi-m_face=1/500000=2*10^-6.

It follows uniformly on every retained target, every matching q,
every source face and sign that

    P_inv(V>chi*x | face)<=f(x)=10x/(a*x-p)^2

once a*x>p. The new base pays a*x0-p>0 and
f(x0)<5.578053*10^-12<10^-6=1-Ret. Direct differentiation gives

    f'(x)=-10(a*x+p)/(a*x-p)^3<0.

Hence this same retention holds on the ENTIRE tail. The 10x, the
small exact rate gap and the subtraction p=.6 are all retained;
the R083 2x or other old constants are not reused. Only after uniform
facewise control are faces mixed. min(1,G/theta)=1 on V<=chi*x
then gives the raw soft-row lower coefficient. No Gaussian model,
post-cap independence or actual-increment probability law is used.

## 6. Weighted cap, complete degrees and the SAME global tau

The full-matching count gives

    nu=1/s-(1+1/s)e>0,
    eta=(1-qmax)/(2qmax)=1/199998>0.

Every counted free coordinate captures its whole lower strip of
length 1-q. An arbitrary positive source pin omits at most one
counted coordinate; negative pinned faces cannot meet D-complement.
The accepted pre-cap inverse-product argument gives

    caploss/k<=M exp[-(nu-1/n)*eta*x/(1+e)].

Pointwise bsoft<=theta/G, the complete face sum and actual matched
dq/q integration are essential: theta cancels against k exactly.
This is the weighted captured deletion, not an unweighted source-band
mass. No requirement eta>chi is imposed.

For n>=N12 its lower exponent is

    K12(L)=C[nu*L-12!*L^-11],
    C=(1-sigma)*eta/(1+e)>0.

The base pays nu-1/N12>0 and K12>=4. Its derivative is
C[nu+11*12!*L^-12]>0, so M/T8(K12(L)) decreases. The independent
base is already <e; the complete cap payment stays paid for all
larger L.

The radial upper envelope 1/[1-(h/rlo)^2]<=1+e is checked and
decreases; its lower side is at least 1 and hence >=1-e. Thus raw
rows retain 2Ret(1-e)Ulo*k as their lower bound and M*k as upper
bound. The complete target-face column proof includes every pin
sign, uses r>=R and bsoft*G<=theta, and gives <=2k<M*k. Inactive
columns vanish. These deterministic degree arguments are covered
by the accepted R085 chain and have no further threshold condition.

Directly for all L>=L0,

    exp(L/2)>=(L/2)^2/2=L^2/8>=x=(1-sigma)L.

The last inequality holds already for L>=16. Therefore log t>=L/2,
B>=nL/4>gamma*n. The same single numeric selection threshold
tau=sqrt[A(B-gamma*n)] is positive even though gamma>1.
A<=exp x and B<=nL hold; B<=1+n log t<=nL also follows from
x>=2,n>=400 and n log x>=n/2>1.

Consequently the TRUE tau payment is

    tau/k<=rhi*L^(5/2)*exp(-cT L)
       <=rhi*L^(5/2)/T8(cT L).

At the base sqrt L0=2^40 exactly, and this is <e. The entire
envelope is a constant times L^(-11/2), so it decreases for all
larger L. The same tau is used in the original algorithm, the
actual-row eligibility inequality and the ambient deficit; no second
unpaid threshold is introduced.

The author's S16 alternative is valid as well. For y>=4,
yS16'-(5/2)S16 exceeds the lower polynomial with coefficients
-5/2,-3/2,-1/2,1/12,1/16. At y=4+u its coefficients are
29/6,29/2,13/2,13/12,1/16, all positive. The exact quadratic
coefficient is -1/4, and -1/2 is a valid weaker lower bound. The
base cT*L0>=4 is independently verified. This agrees with, rather
than contradicts, the earlier accurate quartic repair.

## 7. Finite history and genuine deficit for EVERY actual dimension

At each actual integer n on this domain, n>=N12(L)>=400, x>=2
and t<n. For EVERY a_alg>=8n^4, a_alg>2t and ccap>1. A genuine
target item meeting A has radius<1+rhi<6, while any included source
increment has radius>=rho>20>6. It precedes that SAME target in
the original largest-radius history. These are the correct uniform
family dimensions, not an inherited smaller-s separation argument.

The displacement margin is proved for ACTUAL n algebraically:

    n^-3-(beta+2)/(a_alg+1)
      >n^-3-4/(8n^4)
      =(n-1/2)/n^4>0.

The base substitution N12 is used only to guarantee n>=400 and
beta<2, not to infer monotonicity of this rational expression.
Thus each included true source center captures the true target in
its own original eligibility test. Integrating only genuine included
increment pieces costs no more than original item weights, and both
bounded kernel weights give actual row<=arad*tau<=tau.

On A and P the finite envelope equals c_a*lambda exactly, because
the original D cap lies in the finite cap with the positive algebraic
margin in R085. Energy, degree and target mass scale by c_a^2,c_a,c_a.
A and P are disjoint, so the two missing-mass terms cost at most
c_a*M*k times the TOTAL full finite-cap deficit. The accepted
energy identity therefore yields, for EVERY original finite history,

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
       >=[E-tau*lambda(A)]/(M*k)
       >=lambda(A)*b/M>=Gamma*n>gamma*n.

The c_a>=1 direction is correct. Paying cap and the SAME tau
leaves row-tau>=b*k, and the strict fixed paid-Gamma inequality is
independently verified. No gamma<1 restriction, deficit<=n bound,
new selector or actual-increment iid assumption is involved.

For arbitrary L1 input and every strict positive level, use |f| and
its actual positive-radius witnesses. Open shrunk cubes give finite
subcovers of compact subsets; selected dilates remain true witnesses.
Normalization works for every positive terminal radius. Lebesgue
inner regularity recovers the whole strict level set. No compact
support, truncated scale range, radius lattice, independent history
or selector limit is needed. The uniform finite-a proof is the
accepted original-operator bridge, now with all its numeric side
conditions explicitly paid.

## 8. Entire cover/root cost and the strict .555 conclusion

The full original finite cover factor, captured tau and terminal
self-weight are retained. Only in the family of scalar inequalities
already uniform over finite histories is a_alg allowed to tend to
infinity. Its exact result is

    C_n<=(1+A/tau)(tau+B-gamma*n)=F(A,B-gamma*n).

The complete positive remainder is

    R(L)=1/n+exp(-sigma L)+2sqrt L exp(-sigma L/2).

The independent base check pays

    1/N12(L0)+1/T8(sigma L0)
      +2^41/T8(sigma L0/2)<21/100000.

This includes the positive A and root cross-term and the terminal
scalar self contribution. The corresponding all-tail envelopes have
powers L^-12,L^-8 and L^(-15/2), so R(L)<21/100000 on every
larger L. The base square root is used only at L0. Equivalently,
the author's S16 half-power envelope decreases on y>=2, whose base
domain sigma*L0/2>=2 is also checked.

Finally

    -log(1-sigma)<99229/100000,
    gamma-99229/100000-111/200=21/100000.

These strict exact inequalities and the full R(L) payment prove
F(A,B-gamma*n)<T_n-(111/200)n. The claim concerns the complete
cover function, not a raw deficit, a model objective, the supremum
coefficient c_* or an unpriced display.

## 9. Decision and identity preservation

ACCEPT the fixed R087 branch for every natural dimension n with
log n>=2^80, with the full statement at the start of this report.
All five event bounds, uniform local conditions, actual inverse-source
retention, weighted cap, same tau, full root and actual finite-history
requirements have explicit base and entire-tail coverage. No numerical
or mathematical gap was found within the accepted R085 framework.

The author 8+56 run identities, the accepted generic-family theorem
and this new explicit-threshold audit remain separate. Original hashes
are recorded and rechecked before handoff. No other round, paper,
accepted branch domain or global-optimality claim is changed.
