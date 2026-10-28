# R081 proposed eventual theorem: a new s=4 covering geometry

Status: complete proposed eventual theorem, pending independent review.
No explicit dimension threshold is claimed. This proof rederives the
local profiles on the larger relative-gap domain; it does not use the
s=3/2 profile constants or a numerical model optimum.

There exists a finite L_* such that, for every integer n with
L=log n>=L_*, the original Lebesgue centered common-radius cube maximal
operator, for every L1 input and all positive witness radii, satisfies

    C_n <= Phi81=F(A_n(t),B_n(t)-(13/10)n)
        < T-(3/8)n,
    t=n/[(2/5)L], T=n log n-n log log n,
    F(A,B)=(sqrt A+sqrt B)^2.                           (1)

The assertion is eventual in log n, with a finite but unquantified
threshold determined by the finitely many explicit inequalities below.
It is neither a theorem at L=4096 nor an assertion about optimal s.

## 1. Fixed new geometry, weights and budgets

Fix

    s=4, qmin=1/4, qmax=99/100,
    sigma=3/5, chi=7/10, x=(2/5)L, t=n/x,
    beta=1+1/t, zeta=n^-3, epsilon=zeta*t=1/(x*n^2),
    kappa=beta-zeta,
    delta=1/100000, rlo=1+delta, rhi=99/25-delta,
    h=L^-2, theta=exp(-chi*x), k=n*h*theta/rhi.

Set ALL the fixed errors eB,eN,eJ,eRad,eS,eT,eD and epsA equal to
e=1/10000. Take Ret=999/1000, mL=137/100, gamma=13/10,
price=23/25 and c=3/8. Put

    P0=1-2e,
    Ulo=P0(1-e)/(1+e), Uhi=(1+e)/[P0(1-e)],
    M=2(1+e)Uhi,
    b=2 Ret(1-e)Ulo-2e,
    Gamma=mL*b/M.

The independent fixed-constant Fraction run verifies

    (1-e)log(rhi/rlo)>mL,
    (1-e)Gamma>13/10,
    log(5/2)<23/25,
    gamma-price-c=1/200.                               (2)

For orientation only, (1-e)Gamma is about 1.366988. There is NO
assumption Gamma<1 or gamma<1; positivity of B_n-gamma*n will be
proved directly in Section 8.

Use the original envelope and its pre-cap version

    D={||z||_infinity<t: some z_i<1-beta*m(z)},
    m(z)=max(1,||z||_infinity),
    d lambda=2^-n*m(z)^-n*1_D dz,

with lambda0 the same density on full annuli. Let A0 have target
radius r in [rlo,rhi], and let Pfull have source radius rho=tq with
q in [qmin,qmax]. Write P=D intersect Pfull. On the full source band
define

    R=4q-h, a0=(1-epsilon)q, g=R-a0,
    ell_R(z_i)=min(1,(R+kappa*rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i), bsoft=min(1,theta/G).

Extend bsoft by zero off Pfull, and let arad=r/rhi on A0, zero off.
After the retained target set A is defined below, use exactly

    K(w,z)=1_A(w)1_P(z)1_{|r-4q|<h}
       *1_{||z-w||_infinity<kappa*rho}*arad(w)*bsoft(z). (3)

Both weights are bounded by one; they do not change actual item weights
or the original largest-eligible-radius selection algorithm.

## 2. Uniform eventual active domain

For a target radius r, the COMPLETE matching interval and baseline are

    q_-=(r-h)/4, q_+=(r+h)/4, qb=q_+, R(qb)=r.

When 2h<delta this interval lies in [1/4,99/100], and all active
sources satisfy

    1<R in [r-2h,r]<4,
    a0>0, rho>4rhi, 0<epsilon<1,
    g-3R/4=epsilon*q-h/4<=0, 0<g<=3R/4.               (4)

These statements hold uniformly once epsilon*qmax<=h/4 and t*qmin
is large. Indeed epsilon*qmax/h is a constant times L/n^2 and tends
to zero. R>1 follows from rlo-2h>1; g>0 follows from 3q-h>0.
Also beta<2, kappa>1 and t>rhi eventually. Both bands are inside
norm<t, and are disjoint. Inactive edge sources may have R<=1 but
have zero final columns; none of the following local estimates is used
there. On the full source band R>0 and G>0 eventually.

All constants below are uniform on (4), |w_i|<=r<4 and q>=1/4.
In particular ell_R(z_i)>=1-g/(2R)>=5/8 on [-rho,rho].

## 3. New coordinate formulas and new derivative constants

Define, before source-cap or soft conditioning,

    J_j(w;q)=(2rho)^-1 integral_capture phi(z)^j/ell_R(z) dz,
    phi=-log ell_R, j=0,1,2, J_0=J.

Let Q(R,u)=-log(1-u/(2R)), m=(|w|-a0)_+, v=(g-m)_+ and
F(R,u)=2R*Q(R,u)-u. Direct integration of the two source boundary
strips gives, also when |w|>R,

    J=1+[F(R,g)+F(R,v)-m]/(2rho),
    J_1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J_2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).                 (5)

For example removing a segment of length m from the near source tail
leaves boundary length v=(g-m)_+; extra central deletion when m>g
has phi=0 and is already included by -m. This proves the formulas
without assuming |w|<=R throughout matching.

For 0<=u<=3R/4,

    0<=Q<=log(8/5)<1/2,
    |Q_R|<=3/(5R), Q_u<=4/(5R),
    0<=F<=R<=4, |F_R|<=11/5, F_u<=3/5.              (6)

The logarithm bound follows from e^(1/2)>1+1/2+1/8=13/8>8/5.
Here R'=4, g'=3+epsilon<=4 and a0'=1-epsilon in [0,1].
The positive parts satisfy |m'|<=1, |v'|<=5 a.e. They are Lipschitz,
so (5) is AC; genuine onset kinks do not invalidate integration of
the a.e. derivative estimates.

The numerator A=F(R,g)+F(R,v)-m in J has |A|<=12 and
|A'|<=24. Since rho'/rho=1/q<=4, this yields

    |J-1|<=6/rho, |J'|<=36/rho,
    |partial_q log J|<=72/rho<100/rho                 (7)

once 6/rho<=1/2. Similarly (5)--(6) give the NEW bounds

    0<=J_1<=1/rho, |J_1'|<=29/rho<40/rho,
    0<=J_2<=1/(3rho).                                 (8)

For the J_1 derivative, the terms from R', the two Q terms and rho'
are bounded respectively by 1/rho, 24/rho and 4/rho. The estimates
are deliberately coarse and do not import the s=3/2 constants.

## 4. New baseline target laws and a five-event good set

At qb, R=r. On each FULL TARGET face the n-1 free coordinates are
independent Uniform[-r,r], and one coordinate is pinned at +/-r.
Fubini gives E_W J=1: averaging the capture indicator over the full
target interval gives exactly ell_R(z), which cancels its inverse.
Since |J-1|<=6/rho_b<=1/2, elementary log bounds give

    -36/rho_b^2 <= E_W log J <=0,
    Var_W(log J)<=144/rho_b^2.

The actual target pin satisfies |log J|<=12/rho_b. Consequently the
baseline sum has absolute mean at most 36n/rho_b^2+12/rho_b and
variance at most 144n/rho_b^2. These tend uniformly to zero because
rho_b>=t/4 and n/rho_b^2<=16x^2/n.

For the first profile, the same Fubini cancellation gives exactly

    E_W J_1=E_fullsource phi
       =(2R/rho)[1-a+a log a], a=1-g/(2R)>=5/8.

The function H(a)=1-a+a log a has derivative log a<=0 on [5/8,1].
Since R/q<=4,

    n E_W J_1 <=8x H(5/8)=mu4*x,
    mu4=3+5log(5/8)<13/20.                            (9)

The strict inequality is independently checked by a positive atanh
lower enclosure log(8/5)>47/100; it is not taken from the model.
The true target pin contributes at most 1/rho_b. The variance of
the free sum is <=n/rho_b^2<=16x^2/n. For large L the total mean
is <=.655x, so Chebyshev gives

    P(sum_i J_1(w_i;qb)>(33/50)x)<=640000/n.           (10)

Here .655=.65+.005 and .66-.655=.005. This probability tends to
zero uniformly on target radius and face.

Let A be A0 intersected with exactly these FIVE conditions:

1. |sum log J(w_i;qb)|<=e.
2. |Nsum(q_-)-E_-|<=e*n/4.
3. |Nsum(q_+)-E_+|<=e*n/4.
4. w belongs to the original D.
5. sum J_1(w_i;qb)<=(33/50)x.

As before Nsum=sum_i[1+1_{|w_i|<a0(q)}] and
E_+/-=n+(n-1)a0(q_+/-)/r. Eventually a0<r, so the actual target
pin contributes exactly 1, while the remaining excess indicators are
independent Bernoulli(a0/r). Each count failure is <=4/(e^2*n).
The baseline-log failure tends to zero by its mean/variance bounds.
For the target cap, negative pinned faces lie entirely in D once
t>rhi; on positive faces its complement has free fraction
1+1/(2t)-1/(2r)<=9/10 for large L. Its failure is therefore at most
(9/10)^(n-1). The five failures are each <=epsA/8 for all sufficiently
large L, uniformly in radius/face. No independence between them is used.

Thus

    lambda(A)>=(1-epsA)n log(rhi/rlo)>mL*n.             (11)

All these probability laws are ambient full-target laws used to select
A; the actual increment measure is not assigned any such distribution.

## 5. Whole-matching products, counts and genuine source-face retention

For fixed w the matching width is 2h/s=h/2. From (7),

    sum_i |partial_q log J_i| <=400x,
    |sum log J_i(q)-sum log J_i(qb)|<=200xh=80/L.

For L large enough that 80/L<=e, the baseline event therefore gives
P0<=product_i J_i<=1/P0, with P0=1-2e. Indeed exp(-2e)>=1-2e
and exp(2e)<=1/(1-2e). Also |J_i-1|<=e eventually, uniformly.

The source pin inverse factor is ellpin=(1+a0/R)/2. Its derivative
is -(1-epsilon)h/(2R^2)<0, independently of s. Nsum is nondecreasing,
so Nsum/(2n ellpin) is nondecreasing. At q_+, R=r and its mean ratio
is at most 1. At q_- its mean ratio is at least 1-2h-1/n: subtract
a0(1/R-1/r)+a0/(nr) from 1+a0/R and use a0<=1,R>=1.
Count fluctuations change the ratio by at most e/4. For large L,
2h+1/n<=3e/4, so over the full matching interval

    1-e<=Nsum/(2n ellpin)<=1+e,
    U=[product_i J_i/(2ellpin)]sum_i N_i/J_i,
    Ulo*n<=U<=Uhi*n.                                   (12)

The expression sums all captured signed source faces; none is omitted.
The new first-profile AC bound (8) likewise gives a whole-matching
transport at most 80xh, since n/rho<=4x. Once 80h<=1/100,

    sum_i J_1(w_i;q)<=(67/100)x.                        (13)

Fix ANY captured signed SOURCE face j before cap deletion and soft
clipping. Its loss pin is deterministic and <=log(8/5)<1/2. The
free coordinate inverse laws are independent, with means J_1i/J_i
and raw second moments J_2i/J_i. From (8), (13), and J_i>=1-e,

    E(V|j)<=1/2+(671/1000)x,
    Var(V|j)<=sum_(i!=j) J_2i/J_i
             <=4x/[3(1-e)]<2x,
    V=-log G.                                         (14)

The fixed-constant script verifies .67/(1-e)<.671. No second-profile
target event, lower variance or Gaussian approximation is needed.
For u(x)=(29/1000)x-1/2>0, Chebyshev on each face gives

    P_inv(V>chi*x|j)<=2x/u(x)^2 ->0.                   (15)

Hence for all sufficiently large L the retained fraction is at least
Ret=.999, uniformly over targets, complete matching, faces and signs.
Only then mix faces. The ratio bsoft/(theta/G)=min(1,G/theta) equals
one on V<=chi*x, so the fullsource shell coefficient is between
Ret*theta*Ulo*n and theta*Uhi*n. Independence is not asserted after
cap or soft conditioning.

## 6. General-s weighted source cap and complete column bound

From a0/R>=1/4 and (12), the counted coordinates satisfy

    #{i: |w_i|<a0}>=nu*n,
    nu=1/4-(5/4)e>0.

At source rho=tq, the cap complement is exactly all
z_i>=-rho+(1-q). A negative source pin cannot lie there. For every
counted FREE coordinate its entire lower strip [-rho,-rho+(1-q))
is captured: w_i<a0 pays one endpoint, and rho>4rhi pays the other.
Its pre-cap inverse-law probability is >=(1-q)/[2rho(1+e)]. An
arbitrary positive source pin deletes at most one counted coordinate.
Independence and 1-y<=exp(-y) yield, on EVERY captured source face,

    P_inv(D^c|j)<=exp[-((nu-1/n)*eta/(1+e))*x],
    eta=(1-qmax)/(2qmax)=1/198.                         (16)

Multiplying bsoft<=theta/G pointwise, summing all face coefficients
and integrating only over matched dq/q cancels theta against k. The
exact radial relation is unchanged by s:

    (r/rhi)log((r+h)/(r-h))
       in [2(1-e)h/rhi,2(1+e)h/rhi]

for all sufficiently small h. Thus the complete source-cap loss obeys

    Lcap(w)/k<=M exp[-((nu-1/n)*eta/(1+e))*x] ->0.       (17)

Indeed nu-1/n>=nu/2 eventually and x=(2/5)L; the remaining exponent
has a fixed strictly positive coefficient. It is therefore <=eS=e
for all sufficiently large L. No condition chi<eta is needed.

Before cap deletion, raw rows lie between 2Ret(1-e)Ulo*k and M*k.
The column argument must also be checked at the enlarged target radius.
For a fixed source and matched target radius r>=R, let ell_i(r) be its
target-coordinate capture fraction. One has ell_i(r)>1/2 because
kappa*rho-|z_i|>=a0>0; also ell_i(r)<=ell_i(R), since
min(1,1/2+(kappa*rho-|z_i|)/(2r)) decreases in r. If ell_i(r)<1 at
most one target pin sign is captured, otherwise at most two. Summing
the complete target faces gives coefficient <=n product_i ell_i(r)
<=nG(z). Multiplication by arad cancels dr/r, the ordinary matching
target width is at most 2h, and bsoft*G<=theta. The complete column
is consequently <=2k<M*k. Inactive columns are zero. This uses no
old s=3/2 constant or probabilistic target law.

## 7. Existence of one common eventual threshold

All geometry guards in Section 2, individual-J, log-product, count,
profile-transport and radial-width guards are eventually valid by the
displayed uniform bounds. The bad target probabilities tend uniformly
to zero, and the bound (15) tends to zero as L tends to infinity.
The source-cap payment tends to zero exponentially by (17).

These are finitely many inequalities with fixed positive error budgets.
Therefore their thresholds have a finite maximum. This assertion uses
ordinary limits of explicit envelopes, not convergence of selectors or
histories. For example every occurrence of L^a/n tends to zero using
e^L>=L^k/k! with an integer k>a, and every L^a exp(-bL), b>0,
tends to zero by the same positive-series argument. No parameter grid
or unquantified local uniformity is being used.

The resulting threshold need not be small. In particular the present
coarse log-product constant requires L>=800000, and the Chebyshev
retention threshold may be larger. This proof does not transfer the
R080 threshold to the new geometry.

## 8. Tau payment, positive roots and the original finite-history bridge

The original cap identity with Icap<=log t implies
B_n>=(n/2)log t. Since e^(L/2)>=x for all sufficiently large L,
log t>=L/2 and

    B_n>=nL/4>(13/10)n.

Thus gamma>1 causes no problem. Define the SAME global selection
threshold tau=sqrt[A_n(B_n-gamma*n)]>0. The original scalar bounds
A_n<=exp x=n exp(-sigma*L), B_n<=nL give

    tau/k<=rhi*L^(5/2)*exp(-cT*L) ->0,
    cT=sigma/2-chi(1-sigma)=1/50>0.                   (18)

For sufficiently large L it is <=eT=e. After source-cap deletion and
this same tau are paid, row minus tau is at least b*k. Both ambient
degrees are at most M*k and target mass is >mL*n.

For every original algorithmic a_alg>=8n^4 use the usual normalization
by the last minimum positive-contribution radius and reflect the last
center into the positive orthant. With alpha=a_alg/(a_alg+1) and
ccap=(a_alg*beta-1)/(a_alg+1), a_alg>2t gives ccap>1. For an actual
target item j whose increment meets A,

    ccap*r_j<=||u-w||<alpha+rhi<5,

so r_j<5. Actual source increments have r_i>=rho>4rhi>5, and precede
that SAME target j. The displacement margin is still

    alpha*beta-kappa-2/(a_alg+1)
       =zeta-(beta+2)/(a_alg+1)>0,

by beta<2 and a_alg>=8n^4. Included pairs thus capture the true target
center in its own original eligibility test: their actual scaled center
distance is less than kappa*rho+(r_i+r_j)/(a_alg+1), which is at most
[kappa+2/(a_alg+1)]r_i<alpha*beta*r_i because r_j<r_i.
The two weights in (3)
are <=1, so the actual weighted target row is <=arad*tau<=tau.

On A and P the true finite corner envelope is exactly
lambda_a=c_a*lambda, c_a=alpha^-n>=1, because

    alpha-ccap*m-(1-beta*m)=[(beta+1)m-1]/(a_alg+1)>0.

Let E=integral K d lambda d lambda. The true uncaptured-increment
measure mu_a is dominated by lambda_a. On this support energy, degree
and target mass scale by c_a^2,c_a,c_a. A and P are disjoint, so in
the energy-difference identity both missing-mass terms cost at most
c_a*M*k times their missing masses, whose sum is bounded by the TOTAL
finite-cap deficit Delta_a. The actual energy is <=tau*c_a*lambda(A).
Hence for EVERY finite original history

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/[c_a*M*k]
       >=[E-tau*lambda(A)]/[M*k]
       >=lambda(A)*b/M>=Gamma*n>gamma*n.               (19)

No step requires gamma<1. General L1 inputs and all positive scales
are recovered exactly as in the original cover: use |f|, finite open
shrunk-witness coverings of compact subsets of a strict level set,
selected dilates as true witnesses, and inner regularity. Every positive
terminal radius is allowed. There is no compact-support, truncated-scale
or actual-increment iid restriction, and no selector/history limit.

## 9. Complete cover cost and the claimed display

Retain the original captured budget, terminal self-weight and ENTIRE
cover factor at every finite a_alg:

    C_n <=[1+(beta+1/a_alg)^n/tau]
       *[tau+a_alg^-n+alpha^-n J_a_alg(t)-gamma*n].

In this uniformly valid family of numeric bounds alone let a_alg tend
to infinity, obtaining

    C_n<=(1+A_n/tau)(tau+B_n-gamma*n)=Phi81.            (20)

The scalar bounds give

    B_n<=T-n log(1-sigma)+1=T+n log(5/2)+1,
    [A_n+2sqrt(A_n B_n)]/n
       <=exp(-sigma*L)+2sqrt(L)exp(-sigma*L/2).

The right-hand remainder plus 1/n tends to zero. By (2) it is <1/200
for all sufficiently large L, giving Phi81<T-(3/8)n. This is the
complete root cost with its positive cross-term paid, not a deficit-
only or model coefficient.

For comparison to R080's complete function, the original integral
I_n(t)<=2t/(n-1) gives at its own sigma=.11 and gamma=.39029

    Phi80>=T-n[.39029-.11+4/L]
          =T-n[28029/100000+4/L].

Since 3/8>28029/100000+4/4096, the new complete Phi81 is strictly
lower than Phi80 on their actual common domain. This does not quantify
L_* or claim a bound for the new geometry at 4096. All accepted branches
still retain their own validity domains and are taken in the minimum.

## 10. Actual computation and remaining acceptance responsibility

verify_constants.py is a bounded Fraction computation of the FIXED
constant gaps, not a parameter sweep or finite-threshold certificate.
It independently verifies the log bounds in (2), mu4<.65, the face
mean/variance rounding, positive weighted-cap and tau rates, the explicit
derivative-constant arithmetic, Gamma>gamma and the complete old-cost
comparison. verification.json records the actual run and exact fractions.
The integral identities, uniform target probabilities, eventual guards
and operator bridge are the mathematical proof above and remain subject
to independent whole-chain review by the main agent.

Only upper/round081/general_s in this independent clone was written.
No old result, paper, original checkout, Git commit/push, shared MATLAB
session, n=7 calculation or lower-bound task was changed or used.
