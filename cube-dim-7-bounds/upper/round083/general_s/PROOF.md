# R083 proposed eventual theorem: fixed s=9/2 and complete cost T-11n/20

Status: complete proposed eventual theorem, pending independent whole-chain
review. Its dimension threshold is finite but not quantified. The local
constants are rederived on the s=9/2 domain; no s=4 or s=3/2 local
constant is assumed. No accepted predecessor or paper is changed here.

There exists a finite L_* such that, for every integer n with L=log n
>=L_*, the original Lebesgue centered axis-parallel common-radius cube
maximal operator, for every L1 input and all positive witness radii, obeys

    C_n <= Phi83=F(A_n(t),B_n(t)-(1501/1000)n)
        < T-(11/20)n,
    t=n/[(387/1000)L], T=n log n-n log log n,
    F(A,B)=(sqrt A+sqrt B)^2.                            (1)

The proof supplies an eventual construction, not an explicit threshold,
an optimal parameter assertion, or a claim at L=4096.

## 1. Locked fixed tuple and scalar gaps

Fix exactly

    s=9/2, qmin=2/9, qmax=9999/10000,
    sigma=613/1000, chi=1583/2000,
    x=(387/1000)L, t=n/x,
    beta=1+1/t, zeta=n^-3, epsilon=zeta*t=1/(x*n^2),
    kappa=beta-zeta,
    delta=1/1000000, rlo=1+delta, rhi=4499549/1000000,
    h=L^-2, theta=exp(-chi*x), k=n*h*theta/rhi,
    eB=eN=eJ=eRad=eS=eT=eD=epsA=e=1/100000,
    Ret=99999/100000, mL=1503/1000, gamma=1501/1000,
    price=19/20, c=11/20.

Put

    P0=1-2e,
    Ulo=P0(1-e)/(1+e), Uhi=(1+e)/[P0(1-e)],
    M=2(1+e)Uhi,
    b=2 Ret(1-e)Ulo-2e,
    Gamma=mL*b/M.

The fixed arithmetic gaps, independently checked with 12-term positive
atanh-log enclosures and exact tails, are

    (1-e)log(rhi/rlo)>mL,
    (1-e)Gamma>gamma>1,
    log(1000/387)<19/20,
    gamma-price-c=1/1000.                              (2)

Here (1-e)Gamma is about 1.502804621; only the exact rational in the
receipt certifies its direction. No assumption gamma<1 or Gamma<1
is made. The producer's registered 64-term constant receipt was read
and locked, but its program was not run again.

## 2. New active geometry and exact kernel

Use the unchanged original corner envelope

    D={||z||_infinity<t: some z_i<1-beta*m(z)},
    m(z)=max(1,||z||_infinity),
    d lambda=2^-n*m(z)^-n*1_D dz.

Let lambda0 denote its full-annulus version before imposing D. Let A0
have target radius r in [rlo,rhi], and Pfull have source radius
rho=tq with q in [qmin,qmax]. Put P=D intersect Pfull. Define on the
FULL source band

    R=sq-h, a0=(1-epsilon)q, g=R-a0,
    ell_R(z_i)=min(1,(R+kappa*rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i), bsoft=min(1,theta/G).

Extend bsoft by zero off Pfull and let arad=r/rhi on A0, zero off.
For the retained target set A to be defined below, the FINAL kernel is

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
       *1_{||z-w||_infinity<kappa*rho}*arad(w)*bsoft(z). (3)

All fields are Borel and both weights lie in [0,1]. No actual item
weight or original largest-eligible-radius selection rule is changed.

For each r the COMPLETE matching interval and baseline are

    q_-=(r-h)/s, q_+=(r+h)/s, qb=q_+, R(qb)=r.

Once 2h<delta this interval lies in [qmin,qmax]. Uniformly on all
ACTIVE source matches, for all sufficiently large L,

    1<R in [r-2h,r]<9/2,
    rho>4rhi, 0<epsilon<1, a0>0,
    g-(7/9)R=epsilon*q-h/s<=0, 0<g<=7R/9.             (4)

Indeed epsilon*qmax/h is a fixed multiple of L/n^2 and tends to zero,
and t*qmin tends to infinity. Also beta<2, kappa>1 and t>rhi
eventually. Both bands lie inside norm<t and are disjoint. Inactive
edge sources can have R<=1 but have zero final columns; no local
R>=1 estimate is invoked there. R>0 and G>0 on all of Pfull.

The source-coordinate lower capture fraction is now

    ell_R(z_i)>=(R+a0)/(2R)>=11/18>1/2.                (5)

This enlarged relative-gap domain is the input to the new local bounds.

## 3. Fresh s=9/2 profiles, absolute continuity and derivative constants

For |w|<=r define the actual pre-cap captured coordinate profiles

    J_j(w;q)=(2rho)^-1 integral_capture phi(z)^j/ell_R(z) dz,
    phi=-log ell_R, j=0,1,2, J_0=J.

Write m=(|w|-a0)_+, v=(g-m)_+,
Q(R,u)=-log(1-u/(2R)), and F(R,u)=2R*Q(R,u)-u.
Integrating the two source boundary strips gives

    J=1+[F(R,g)+F(R,v)-m]/(2rho),
    J_1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J_2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).                 (6)

If m>g, extra central deletion has phi=0 and is already paid by -m.
Thus (6) also holds when |w|>R, as can happen inside matching. The
positive-part functions are Lipschitz, so the profiles are AC in q;
their onset kinks do not invalidate the a.e. derivative estimates.

For 0<=u<=7R/9 and 1<=R<=9/2,

    0<=Q<=log(18/11)<1/2,
    |Q_R|<=7/(11R), Q_u<=9/(11R),
    0<=F<=R<=9/2, |F_R|<=25/11, 0<=F_u<=7/11.        (7)

For a purely rational proof of the strict logarithm bound,
exp(1/2)>1+1/2+1/8+1/48=79/48>18/11. In particular it does not
use the insufficient quadratic prefix from the s=4 case.

At fixed n, epsilon is constant in q. Here R'=9/2,
g'=7/2+epsilon<=9/2, a0'=1-epsilon in [0,1],
|m'|<=1, |v'|<=11/2 and rho'/rho=1/q<=9/2. The numerator
A=F(R,g)+F(R,v)-m has

    |A|<=27/2,
    |A'|<=2*(25/11)*(9/2)+(7/11)*(9/2+11/2)+1
          =306/11.

Therefore, with 7/rho<=1/2,

    |J-1|<=27/(4rho)<7/rho,
    |J'|<=3897/(88rho)<45/rho,
    |partial_q log J|<=90/rho.                         (8)

The J_1 derivative is separately bounded by

    |J_1'|<=[9/8+(9/4)(72/11+81/11)+81/16]/rho
           =6597/(176rho)<40/rho,
    0<=J_1<=9/(8rho), 0<=J_2<=3/(8rho).               (9)

The three terms are respectively the R' contribution, the two Q
derivatives, and the rho' contribution. These constants have been
recomputed on (4); none is imported from a smaller s domain.

## 4. Full-target Fubini laws and exactly five good events

At qb, R=r. On each FULL target face the n-1 free coordinates are
independent Uniform[-r,r] and one true coordinate is pinned at +/-r.
Fubini gives E_W J=1: the full target capture probability is precisely
ell_R(z), canceling its inverse. By (8), for 7/rho_b<=1/2,

    -49/rho_b^2<=E_W log J<=0,
    Var_W(log J)<=196/rho_b^2.

This uses log(1+y)>=y-y^2 and |log(1+y)|<=2|y| for |y|<=1/2.
The deterministic target pin has |log J|<=14/rho_b. Hence the whole
baseline log sum has absolute mean <=49n/rho_b^2+14/rho_b and
variance <=196n/rho_b^2. Both tend uniformly to zero, since
rho_b>=t*qmin and n/rho_b^2<=s^2*x^2/n.

For the first profile, the same Fubini cancellation gives exactly

    E_W J_1=E_fullsource phi=(2R/rho)[1-a+a log a],
    a=1-g/(2R)>=11/18.

H(a)=1-a+a log a decreases on [11/18,1], and R/q<=s. Consequently

    n E_W J_1 <=2s*x*H(11/18)=mu*x,
    mu=7/2+(11/2)log(11/18)<7914/10000.               (10)

The strict upper bound is certified by an independent rational LOWER
enclosure of log(18/11). It is not inferred from the numerical model.
The true target pin contributes <=9/(8rho_b). Once n>=506250 this
is <=x/100000, because rho_b>=t/s. The total mean is therefore at
most (79141/100000)x. The free variance is bounded by

    n*(9/(8rho_b))^2<= (6561/256)*x^2/n.

Chebyshev at distance x/50000 gives the exact bound

    P(sum_i J_1(w_i;qb)>(79143/100000)x)
        <=64072265625/n ->0,                          (11)

uniformly in target radius and signed face.

Let A be A0 intersected with precisely these FIVE conditions:

1. |sum log J(w_i;qb)|<=e.
2. |Nsum(q_-)-E_-|<=e*n/4.
3. |Nsum(q_+)-E_+|<=e*n/4.
4. w belongs to D.
5. sum J_1(w_i;qb)<=(79143/100000)x.

Here Nsum=sum_i[1+1_{|w_i|<a0(q)}] and
E_+/-=n+(n-1)a0(q_+/-)/r. Since a0<r eventually, the true target
pin contributes exactly 1, while the other excess indicators are
independent Bernoulli(a0/r). Each count failure is <=4/(e^2*n).
The baseline-log failure tends to zero using its displayed mean and
variance. For the original target cap, negative pinned faces are in D
once t>rhi. A positive face's cap-complement free fraction is

    1+1/(2t)-1/(2r)<=8/9+1/100<9/10

once t>=50, since r<9/2. Thus its failure is at most (9/10)^(n-1).
Together with (11), the five failures are each <=epsA/8 uniformly
for all sufficiently large L. No independence between good events
is assumed. The target mass is consequently

    lambda(A)>=(1-epsA)n log(rhi/rlo)>mL*n.             (12)

These are ambient full-target laws used to select A, never a law for
actual uncaptured increments.

## 5. Whole-matching products, counts and the narrow true face gap

The matching width is 2h/s=4h/9. Since n/rho<=s*x, (8) gives

    |sum log J_i(q)-sum log J_i(qb)|<=180xh
        =(3483/50)/L ->0.

When this is <=e, the baseline event gives P0<=product_i J_i<=1/P0;
exp(-2e)>=1-2e and exp(2e)<=1/(1-2e). Also |J_i-1|<=e eventually.

Let ellpin=(1+a0/R)/2. Its derivative is
-(1-epsilon)h/(2R^2)<0 and Nsum is nondecreasing. Thus their ratio
Nsum/(2n ellpin) is nondecreasing. At q_+ its mean ratio is <=1,
and at q_- it is >=1-2h-1/n: a0<=1 and R,r>=1 bound the difference
a0(1/R-1/r)+a0/(nr). The endpoint fluctuations alter the ratio by
at most e/4. For 2h+1/n<=3e/4, this proves over the FULL interval

    1-e<=Nsum/(2n ellpin)<=1+e,
    U=[product_i J_i/(2ellpin)]sum_i N_i/J_i,
    Ulo*n<=U<=Uhi*n.                                  (13)

The expression sums every captured signed source face and every
omitted-coordinate term. By the freshly derived J_1 derivative (9),
its matching transport costs <=80xh. Once 80h<=1/100000,

    sum_i J_1(w_i;q)<=(79144/100000)x.                  (14)

Fix ANY captured signed source face j BEFORE cap or soft conditioning.
Its pin loss is deterministic and <=log(18/11)<1/2<3/5. Its n-1
free coordinates are independent in their actual normalized inverse-
capture laws, with means J_1i/J_i and raw second moments J_2i/J_i.
Indeed their joint density factors as the product of
1_capture(z_i)/(2rho*J_i*ell_R(z_i)) over i!=j.
From (9), (14), and |J_i-1|<=e,

    E(V|j)<=3/5+(79147/100000)x,
    Var(V|j)<=sum_(i!=j)J_2i/J_i
             <=[27/(16(1-e))]x<2x,
    V=-log G.                                        (15)

The fixed arithmetic verifies (79144/100000)/(1-e)<79147/100000.
The pin adds no variance. No lower variance or second-profile target
event is needed, and mixture variance is never substituted here.

The positive rate chi-79147/100000=3/100000 is small but strictly
positive. For u(x)=(3/100000)x-3/5>0, facewise Chebyshev gives

    P_inv(V>chi*x|j)<=2x/u(x)^2 ->0.                   (16)

For x> (3/5)/(3/100000), this upper bound decreases: its derivative
is -2(ax+p)/(ax-p)^3<0 with a=3/100000,p=3/5. Thus for all
sufficiently large L every face retains at least Ret=.99999. This
requires no Gaussian, compound-Poisson, model or finite-dimension
approximation. Only after uniform facewise control may the faces be
mixed. Since bsoft/(theta/G)=min(1,G/theta) is one on V<=chi*x,
the complete fullsource shell coefficient is between
Ret*theta*Ulo*n and theta*Uhi*n.

## 6. General-s weighted cap and exact complete column

From (13) and a0/R>=1/s=2/9,

    #{i: |w_i|<a0}>=nu*n,
    nu=1/s-(1+1/s)e=22221/100000>0.

At source rho=tq the original cap complement is exactly all
z_i>=-rho+(1-q). A negative source pin is excluded. For each counted
free coordinate the entire lower strip [-rho,-rho+(1-q)) is captured:
w_i<a0 pays one endpoint, rho>4rhi the other. Its pre-cap inverse
probability is >=(1-q)/[2rho(1+e)]. Any positive source pin omits
at most ONE counted coordinate. Independence, followed by 1-y<=e^-y,
therefore gives for every captured signed source face

    P_inv(D^c|j)<=exp[-((nu-1/n)*eta/(1+e))*x],
    eta=(1-qmax)/(2qmax)=1/19998>0.                    (17)

This is uniform over every retained target and the complete matching
interval. No signed empirical law or post-cap independence is assumed.

Pointwise bsoft<=theta/G. Multiply each face's inverse coefficient by
(17), sum the nonnegative face coefficients, and integrate over only
the actual matched dq/q. The exact radial identity

    (r/rhi)log((r+h)/(r-h))
       in [2(1-e)h/rhi,2(1+e)h/rhi]

holds uniformly for small h. Dividing by k cancels theta and gives

    Lcap(w)/k<=M exp[-((nu-1/n)*eta/(1+e))*x] ->0.       (18)

Eventually nu-1/n>=nu/2, so this is bounded by a fixed positive-rate
exponential in L and becomes <=eS=e. No condition chi<eta is used.
Before cap deletion raw rows lie between 2Ret(1-e)Ulo*k and M*k.

For completeness, the enlarged-radius column argument is deterministic.
For a fixed matched source and target radius r>=R, the target capture
fraction ell_i(r)=min(1,1/2+(kappa*rho-|z_i|)/(2r)) is >1/2 and
<=ell_i(R). If ell_i(r)<1 at most one target pin sign is captured;
otherwise at most two. Summing the COMPLETE target faces is at most
n product_i ell_i(r)<=nG(z). Multiplication by arad changes dr/r
to dr/rhi, the ordinary target matching width is <=2h, and
bsoft*G<=theta. The complete column is <=2k<M*k; inactive columns
are zero. The target/source restrictions only reduce upper degrees.

## 7. All eventual guards have one common finite threshold

The geometry and individual-J guards follow from (4), h=L^-2 and
t=n/x. The uniform baseline and first-profile bad probabilities are
bounded by explicit multiples of polynomial(L)/n; count failures
are O(1/n), target-cap failure is <=(.9)^(n-1). The log-product
transport is 69.66/L, the J_1 transport divided by x is <=80/L^2,
and the count and radial calibration errors tend to zero. The true
source-face failure (16) tends to zero as 1/x and decreases eventually.
The weighted source-cap payment (18) tends to zero exponentially.

All these estimates are uniform in the stated compact target-radius
band and complete matching domain. They are finitely many inequalities
with fixed positive budgets. Their thresholds have a finite maximum.
For rigor, polynomial(L)/e^L tends to zero by a sufficiently high
positive exponential-series term, and the same argument gives
L^a*e^(-bL)->0 for every fixed b>0. There is no numerical grid or
selector/history convergence argument behind this eventual threshold.

The threshold can be very large, particularly because (16) has only
a 3e-5 rate gap. It is not inherited from R080 or R081 and is not
claimed to be computed by the fixed-constant certificates.

## 8. Positive global tau and the original finite history at target radius <6

The original scalar cap identity with Icap<=log t gives
B_n>=(n/2)log t. For large L, e^(L/2)>=x, so log t>=L/2 and

    B_n>=nL/4>(1501/1000)n.

Thus the SAME global selection threshold

    tau=sqrt[A_n(B_n-gamma*n)]>0

is well defined despite gamma>1. The original scalar bounds
A_n<=exp x=n exp(-sigma*L), B_n<=nL yield

    tau/k<=rhi*L^(5/2)*exp(-cT*L) ->0,
    cT=sigma/2-chi(1-sigma)=379/2000000>0.             (19)

It is <=eT=e for all sufficiently large L. After cap and this SAME
tau are paid, final row minus tau is >=b*k. Both degrees remain
<=M*k and lambda(A)>mL*n.

For EVERY original algorithmic a_alg>=8n^4 normalize by the last
minimum positive-contribution radius and reflect its last center into
the positive orthant. Let alpha=a_alg/(a_alg+1),
ccap=(a_alg*beta-1)/(a_alg+1). Since a_alg>2t, ccap>1. An actual
target item j whose increment meets A satisfies

    ccap*r_j<=||u-w||<alpha+rhi<6,

so r_j<6. An actual source increment has r_i>=rho>4rhi>6; it
precedes this SAME j. The displacement margin remains

    alpha*beta-kappa-2/(a_alg+1)
       =zeta-(beta+2)/(a_alg+1)>0.

In fact its scaled actual center distance from the target is less than
kappa*rho+(r_i+r_j)/(a_alg+1), at most
[kappa+2/(a_alg+1)]r_i<alpha*beta*r_i since r_j<r_i. Thus it
captures the genuine target center in that original eligibility test.
Both weights are <=1, so actual weighted row is <=arad*tau<=tau.
No altered history, actual-increment law or unpaid row factor appears.

On A and P the finite corner envelope equals lambda_a=c_a*lambda,
c_a=alpha^-n>=1, exactly, since

    alpha-ccap*m-(1-beta*m)=[(beta+1)m-1]/(a_alg+1)>0.

The true uncaptured-increment measure mu_a is dominated by lambda_a.
Energy, degree and target mass on the kernel support scale by
c_a^2,c_a,c_a. For E=integral K d lambda d lambda, actual energy
is <=tau*c_a*lambda(A). In the energy-difference identity the two
missing-mass terms are bounded by c_a*M*k times missing masses on
A and P. These sets are disjoint, so their sum is bounded by the
TOTAL finite-cap deficit Delta_a. Consequently

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/[c_a*M*k]
       >=[E-tau*lambda(A)]/[M*k]
       >=lambda(A)*b/M>=Gamma*n>gamma*n.               (20)

This uses no gamma<1 restriction. It holds uniformly over all original
finite histories at each sufficiently large fixed n.

For arbitrary f in L1 use |f|. Every point of a strict maximal level
set has an actual positive witness radius, whose open shrunk cube
covers that point. Compact subsets have finite such coverings and
selected dilates remain true witnesses. The normalization applies to
every positive terminal radius; inner regularity gives the whole strict
level set. Compact-support, radius-truncation and selector-limit
assumptions are not introduced.

## 9. Complete cover/root cost and strict T-11n/20 conclusion

Retain the original captured budget, terminal self-weight and ENTIRE
cover factor at every finite a_alg:

    C_n <=[1+(beta+1/a_alg)^n/tau]
       *[tau+a_alg^-n+alpha^-n J_a_alg(t)-gamma*n].

Only in these numeric inequalities, already uniform over all original
finite families, let a_alg tend to infinity. The exact limit is

    C_n<=(1+A_n/tau)(tau+B_n-gamma*n)=Phi83.            (21)

The original scalar bounds retain the exact price:

    B_n<=T+n log(1000/387)+1,
    [A_n+2sqrt(A_n B_n)]/n
       <=exp(-sigma*L)+2sqrt(L)*exp(-sigma*L/2).

Together with 1/n this full positive-root remainder tends to zero,
and is <1/1000 for all sufficiently large L. Equations (2) and (21)
give Phi83<T-(11/20)n. The root cross-term is paid explicitly; the
claim is not an unpriced deficit or a numerical model coefficient.

## 10. Complete-function comparisons with R081 and R080

The comparison needs the accurate old positive price. For the original
scalar functions, with ell_t(rho)=1+1/(2t)-1/(2rho),

    I_n(t)=integral_1^t ell_t(rho)^(n-1) d rho/rho
       <=integral_0^(1-1/t) exp[-(n-1)v/2]/(v+1/t) dv
       <=2t/(n-1).

The nonnegative core of B_n and n/(n-1)<=2 imply for any t=n/x,

    B_n(t)>=n log t-2n/x.                              (22)

For R081, x_old=(2/5)L and gamma_old=13/10. A rational lower
log enclosure independently verifies log(5/2)>229/250=.916.
It follows, on its own true positive-root domain, that

    Phi81>=B_old-(13/10)n
       >=T-n[13/10-log(5/2)+5/L]
       > T-n[48/125+5/L].                              (23)

The term is 5/L, not 4/L: the old x is .4L and cannot be treated
as x>=L/2. Using only the weaker price sigma_old=.6 would give a
loss .7 and would not separate these complete functions.

For R080, x_old=(89/100)L and gamma_old=39029/100000. The lower
price -log(.89)>=.11 and 2/x_old<=4/L give

    Phi80>=T-n[28029/100000+4/L].                      (24)

The new independent fixed gates prove

    11/20>48/125+5/64,
    11/20>28029/100000+4/64.

The right sides decrease thereafter. Enlarge the unknown new L_* to
at least 64 if necessary. Then Phi83 is strictly smaller than Phi81
and Phi80 on each actual common validity domain, by (1), (23), (24).
These are lower bounds for old upper-bound FORMULAS, never lower
bounds for C_n. They neither activate an old branch outside its domain
nor quantify the new threshold. Other accepted branches retain their
own domains and are taken in the minimum; no global optimum is claimed.

## 11. Actual fixed computation, provenance and remaining proof status

The producer registration, verify.py and result.json were actually read;
REGISTRATION.json locks their SHA256 identities. Its 18 fixed-constant
PASS gates are historical arithmetic evidence, not an operator theorem.
The independently authored verify_constants.py uses 12 positive atanh
terms with exact tails, fresh derivative-constant arithmetic, target
probability constants and both complete old-cost comparisons. It does
not execute the producer's 64-term program or any parameter grid.

The independent script was actually run once and returned PASS on 19
fixed checks; verification.json preserves its exact rational values and
display-only decimals. This run did not evaluate a dimension threshold.

The new script is a bounded fixed-constant check only. It does not
compute L_* or replace the proofs of AC, full-target laws, source-face
independence, eventual uniformity, weighted-cap deletion, same-history
eligibility or arbitrary-L1/all-radii transfer supplied above. Those
steps remain subject to independent whole-chain review and main-agent
acceptance. All display decimals are separate from exact rational evidence.

Only the owned upper/round083/general_s directory was written. No prior
accepted result, original source, paper, shared MATLAB session, n=7 or
lower-bound task, Git commit or push was touched.
