# R080 proposed full theorem: weighted source cap at log n >= 4096

Status: complete proposed fixed-tuple proof, pending independent whole-chain
acceptance. The selected producer receipt is a conditional finite arithmetic
certificate; the proof below supplies the law, whole-tail and operator
transitions. No accepted ledger or paper is changed by this note.

The claimed conclusion is, for every integer n with L=log n>=4096,

    C_n <= Phi80=F(A_n(t),B_n(t)-(39029/100000)n)
        < T-(5473/20000)n,
    t=n/[(89/100)L], T=n log n-n log log n,
    F(A,B)=(sqrt A+sqrt B)^2.                           (1)

The operator is the original centered axis-parallel common-radius cube
maximal operator on Lebesgue measure, for every L1 input, every positive
level and every positive witness radius. The new mechanism is a weighted
source-cap deletion estimate. The old unweighted cap rate is negative at
this tuple and is not used.

## 1. Exact fixed tuple, fields and numerical identities

Fix

    s=3/2, qmin=2/3, qmax=247/250,
    sigma=11/100, chi=51/1000, x=(89/100)L, t=n/x,
    beta=1+1/t, zeta=n^-3, epsilon=zeta*t=1/(x*n^2),
    kappa=beta-zeta,
    delta=2^-19, rlo=1+delta, rhi=741/500-delta,
    h=1/(1536L), theta=exp(-chi*x), k=n*h*theta/rhi,
    eB=eN=eT=2^-19, eJ=eRad=eS=eD=2^-20,
    epsA=1/65536, Ret=999999/1000000,
    mL=196691/500000, gamma=39029/100000,
    price=5827/50000, c=5473/20000.

Define

    P0=1-eB-6/1536,
    Ulo=P0*(1-eN)/(1+eJ), Uhi=(1+eN)/[P0*(1-eJ)],
    M=2(1+eRad)Uhi,
    b=2 Ret(1-eRad)Ulo-eS-eT,
    Gamma=mL*b/M.

The exact values include

    x0=91136/25, h0=1/6291456, P0=522239/524288,
    Gamma=16492589442569454218876130874189746369 /
          42255181197909082324515317350400000000,
    0<gamma<(1-eD)Gamma<Gamma<1,
    gamma-price-c=1/10000.                             (2)

The independent compact script proves the new mass and price inequalities

    (1-epsA)log(rhi/rlo)>mL,
    -log(1-sigma)<price.                               (3)

It uses 12 positive atanh terms with an exact geometric remainder. The
actual producer best-case receipt used 64 terms, N64 and S64; that receipt
is locked and parsed but its program and 53 gates are not reexecuted here.

## 2. Complete matching geometry, ambient sets and final kernel

Use the original corner envelope

    D={||z||_infinity<t: some z_i<1-beta*m(z)},
    m(z)=max(1,||z||_infinity),
    d lambda(z)=2^-n*m(z)^-n*1_D(z) dz.

Write lambda0 for the same density on full annuli before imposing D.
At a full source shell rho>1 each signed maximal face has shell weight
d rho/(2rho) and n-1 independent Uniform[-rho,rho] free coordinates.

Let A0={w:r=||w||_infinity in [rlo,rhi]}, let Pfull be the FULL source
band q=||z||_infinity/t in [qmin,qmax], and P=D intersect Pfull.
The matching condition is |r-sq|<h. On all of Pfull set

    rho=tq, R=sq-h, a0=(1-epsilon)q, g=R-a0,
    ell_R(z_i)=min(1,(R+kappa*rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i), bsoft(z)=min(1,theta/G(z)).

All these fields are defined before cap restriction and outside matched
sources. Extend bsoft by zero off Pfull. Set arad(w)=r/rhi on A0, zero
off A0. Once the good-target set A below is specified, the FINAL kernel is

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
       *1_{||z-w||_infinity<kappa*rho}*arad(w)*bsoft(z). (4)

Both weights belong to [0,1], and all fields are Borel.

For an independent full-tail guard use Nstar(L)=L^32/32!<=e^L=n.
The compact run verifies at L0=4096 the active guards in its registration:
2h<delta, h<=2^-20, epsilon*qmax<=h/s, epsilon<1/2,
rho>=Nstar*qmin/x>4rhi, n/x>rhi, x/n<1, x*n^2>1,
and 2<=x<=L. Their whole-tail extensions are proved in Section 7.
Thus beta<2, kappa>1 and the target and source bands are disjoint.

For every target radius the COMPLETE matching interval and baseline are

    q_-=(r-h)/s, q_+=(r+h)/s, qb=q_+, R(qb)=r.

They lie in [qmin,qmax]. On every ACTIVE source,

    R in [r-2h,r], R>rlo-2h>1,
    g-R/3=epsilon*q-2h/3<=0, 0<g<=R/3.                (5)

Inactive edge sources may have R<1, but have zero final columns; no local
R>=1 profile estimate is applied there. On the full band R>0 and G>0.
In particular no old hcoef interval is inherited: the new hden=1536
domain is checked through these abstract active guards.

## 3. True target laws and all six retained events

For fixed w,q define the actual captured coordinate profiles

    J_j(w;q)=(2rho)^-1 integral_capture phi(z)^j/ell_R(z) dz,
    phi=-log ell_R, j=0,1,2, J_0=J.

The accepted relative-gap local lemmas bound, a.e. with AC in q,

    |partial_q log J|<3/rho, |J-1|<1/(3rho),
    0<=J_1<3/(50rho), |partial_q J_1|<1/(2rho),
    0<=J_2<1/(125rho), |partial_q J_2|<1/(10rho).        (6)

These are the same local identities used in the frozen R073 proof,
whose actual active hypotheses (5) hold here. In particular J_2 is the
second integrand profile, not J_1 squared or a whole-face raw moment.

On each full TARGET face at radius r, n-1 coordinates are independent
Uniform[-r,r] and one is the actual target pin. At qb, R=r. The local
Fubini identities give the baseline log mean and variance envelopes

    ML=(3/4)L^2/(qmin^2*Nstar)+(1/2)L/(qmin*Nstar),
    VL=(3/2)L^2/(qmin^2*Nstar).

Also, with Q0=log(6/5),

    mu1=1/2-(5/2)Q0<221/5000,
    mu2=1-5Q0-(5/2)Q0^2.

The exact finite-h baseline second-profile mean divided by x lies between

    (1-1/n)[mu2-h(mu2+1/25)]
    and mu2+1/(250*n*qmin).                            (7)

This uses the actual relation a-5/6=h/(3R)-epsilon*q/(2R), with
a=1-g/(2R). It is not inferred from only an upper relative-gap bound.

Let A be A0 intersected with all SIX conditions:

1. |sum_i log J(w_i;qb)|<=eB.
2. |Nsum(q_-)-E_-|<=eN*n/4.
3. |Nsum(q_+)-E_+|<=eN*n/4.
4. w belongs to D.
5. sum_i J_1(w_i;qb)<=(221/5000+1/100000)x.
6. (527/100000)x<=sum_i J_2(w_i;qb)<=(53/10000)x.

Here Nsum(q)=sum_i[1+1_{|w_i|<a0(q)}] and
E_+/-=n+(n-1)a0(q_+/-)/r. The true target pin contributes exactly 1
to each Nsum; its excess Bernoulli indicator is zero.

The new compact finite run verifies ML<=eB/2, 4VL/eB^2<=epsA/8;
each endpoint-count bad probability <=4/(eN^2*Nstar)<=epsA/8;
the target cap free fraction <=3/4, with (3/4)^74<=epsA/8;
the J_1 pin condition n>=9000 and bad bound 324000000/Nstar<=epsA/8;
and the second-profile means strictly between .00528x and .005295x,
with bad bound 5760000/Nstar<=epsA/8. These are the R073 elementary
Chebyshev arguments applied to the newly checked tuple. The last bound
uses mean distance at least x/200000 from the excluded endpoints and
true free-coordinate variance <=n/(125^2*rho_b^2).

The six failures total at most 6epsA/8<epsA; no independence of these
events is required. Consequently

    lambda(A)>=(1-epsA)n log(rhi/rlo)>mL*n.             (8)

All probability statements precede conditioning on A or D. No target
law is imposed on the actual uncaptured-increment measure.

## 4. Whole-matching inverse rows and Bernstein retention

The derivative ellpin'=-(1-epsilon)h/(2R^2)<0, the endpoint counts,
and the newly checked 2h+1/Nstar<=3eN/4 give on the WHOLE interval

    1-eN<=Nsum/(2n ellpin)<=1+eN.

The AC log bridge in (6) costs at most 6Lh=6/1536, so
P0<=product_i J_i<=1/P0. Individual J_i lie in [1-eJ,1+eJ] by
x/(3qmin*Nstar)<=eJ. Summing every captured signed SOURCE face gives

    U=[product_i J_i/(2ellpin)]sum_i N_i/J_i,
    Ulo*n<=U<=Uhi*n.                                   (9)

The other strong derivatives in (6) cost xh and xh/5. The new width
guards therefore yield throughout matching, for every w in A,

    sum J_1<=2211x/50000,
    5269x/1000000<=sum J_2<=5301x/1000000.              (10)

Fix any captured signed source face j BEFORE source-cap deletion and
soft clipping. Since kappa*rho=rho+a0 and |z_i|<=rho,
ell_R(z_i)>=1-g/(2R)>=5/6. Thus each loss phi_i lies in
[0,Q0], Q0<183/1000. The source pin has a deterministic loss in
the same interval. The n-1 FREE coordinates are independent in their
actual normalized inverse-capture laws with means J_1i/J_i and raw
second moments J_2i/J_i. For V=-log G this proves

    E(V|j)<=p+m*x, Var(V|j)<=v*x,
    p=183/1000, m=2211/[50000(1-eJ)], v=531/100000.      (11)

The true variance is sum_(i!=j)[J_2i/J_i-(J_1i/J_i)^2]; the new
guard (5301/1000000)/(1-eJ)<v pays its upper bound. The pin adds no
variance and no face-mixture variance is substituted for (11).

For the centered free losses Y_i=phi_i-Ephi_i, |Y_i|<=p. For k>=2,
E|Y_i|^k<=p^(k-2)E Y_i^2 and k!>=2*3^(k-2). Their independent
mgfs consequently give, for 0<=lambda<3/p,

    log E exp(lambda sum Y_i)
       <=lambda^2*v*x/[2(1-p*lambda/3)].

Write a=chi-m>0 and u(x)=a*x-p. Markov's inequality with the
admissible lambda=u/(v*x+p*u/3) gives

    P_inv(V>chi*x|j)<=exp[-H(x)],
    H(x)=u(x)^2/[2(v*x+p*u(x)/3)].                      (12)

The compact new base checks prove u(x0)>0 and S32(H(x0))>10^6,
where S32(y)=sum_(j=0)^32 y^j/j!. The base exponent is about
14.43058, only a display of an exact fraction. Moreover, for d=p/3,

    H'(x)=u[a*v*x+p*v+d*a*u]/[2(v*x+d*u)^2]>0.

Hence on the entire tail exp(H(x))>10^6 and the retained fraction is
strictly larger than Ret. Only after this uniform facewise estimate is
proved may the signed source faces be mixed. The soft ratio
bsoft/(theta/G)=min(1,G/theta) equals one on V<=chi*x, so the complete
fullsource shell row lies between Ret*theta*Ulo*n and theta*Uhi*n.

## 5. Weighted source-cap deletion and complete rows/columns

The new weighted-cap lemma is proved in the owned R079 weighted_cap
note. Its short derivation is included here to close the new mechanism.
By the whole-matching count in Section 4,

    mcount=#{i: |w_i|<a0}>=nu*n, nu=2/3-(5/3)eN,

because 2ellpin=1+a0/R>=5/3. At source radius rho=tq the cap
complement is exactly all z_i>=-rho+(1-q). A negative source pin
cannot lie there. For every counted FREE coordinate, w_i<a0 implies
the entire strip [-rho,-rho+(1-q)) is captured; rho>4r pays the other
capture endpoint. Its inverse-law probability is at least

    (1-q)/[2rho(1+eJ)],

since ell<=1 and J_i<=1+eJ. An arbitrary positive source pin omits
at most one counted coordinate. Independence in the pre-cap face law
and 1-y<=exp(-y) therefore prove, for every captured face,

    P_inv(D^c|j)<=exp[-((nu-1/n)*eta/(1+eJ))*x],
    eta=(1-qmax)/(2qmax)=3/494.                         (13)

The right side is uniform over the WHOLE matching interval and every
retained target. Negative source pins have probability zero. It uses
no signed target distribution or cap-conditioned independence.

Pointwise bsoft<=theta/G. Multiplying the inverse cap integral on each
source face by its full inverse coefficient and summing yields the
weighted cap shell loss <=theta*Uhi*n times (13). Integrate only over
the actual matched dq/q interval, and use

    (r/rhi)log((r+h)/(r-h))<=2(1+eRad)h/rhi.

After dividing by k=n*h*theta/rhi this gives

    Lcap(w)/k<=M*exp[-((nu-1/n)*eta/(1+eJ))*x]
               <=M/S32(K(L)),
    K(L)=(nu-1/Nstar(L))*eta*x/(1+eJ).                 (14)

The theta cancels exactly. No unweighted source-band log mass or factor
1536L is used. In particular the obsolete negative quantity
(eta-chi)(1-sigma) is neither a gate nor an exponent in this proof.

The independent compact base run obtains M/S32(K0)<eS, about
7.817854e-7. Section 7 proves this payment persists on the whole tail.

Before cap deletion the exact source log-width with arad gives raw rows
between 2Ret(1-eRad)Ulo*k and M*k. For a fixed matched source and
target radius r>=R, the complete sum of captured pinned TARGET faces
is at most n product_i ell_i(r)<=nG(z): if ell_i(r)<1 at most one
pin sign is captured, otherwise at most two, and ell_i(r)>1/2.
Multiplying by arad changes dr/r to dr/rhi, the ordinary target width
is <=2h, and bsoft*G<=theta. Thus the complete column is <=2k<Mk.
Inactive columns are zero; restrictions only reduce upper degrees.

## 6. The same global tau and its full-tail payment

For the original scalar cap identity, Icap<=log t gives
B_n>=(n/2)log t. Since e^(L/2)>=(L/2)^2/2>=L>=x on this domain,
log t>=L/2 and B_n>=nL/4>gamma*n. Also B_n<=nL and
A_n<=exp x=n exp(-sigma*L). Thus

    tau=sqrt[A_n(B_n-gamma*n)]>0,
    tau/k<=rhi*1536*L^(3/2)*exp(-cT*L),
    cT=sigma/2-chi(1-sigma)=961/100000>0.              (15)

The independent base check uses sqrt(4096)=64 exactly and proves

    rhi*1536*4096*64/S32(cT*4096)<eT,

about 3.538526e-8. This is a distinct, deliberately coarser certificate
than the producer's S64 payment. The full-tail justification is below.
After cap and this SAME tau are paid, final row minus tau is >=b*k.
Both ambient degrees remain <=M*k and lambda(A)>mL*n.

## 7. All-tail verification, rather than dimension sampling

Every base guard above was freshly evaluated for this exact tuple using
Nstar=L^32/32!, log12 and S32. For L>=4096, Nstar increases and
Nstar/L grows as L^31. After replacing n by Nstar, epsilon and the
individual-J errors are negative powers of L; the geometry/radius
lower guards increase. ML and VL are positive sums of L^-30 and
L^-31, and all bad probabilities are constants times L^-32. The
count error 2h+1/Nstar decreases. The target-cap fraction decreases
since x/Nstar decreases. The H2 lower mean envelope (7) increases
and its upper decreases, with positive factors already checked at the
base. Both profile-transport errors decrease. The fixed log-product
cost 6Lh=6/1536 is constant. The radial-width error decreases.

Bernstein retention persists by the explicit H'>0 proof in Section 4.
For the new weighted cap,

    K(L)=eta*(89/100)/(1+eJ) * [nu*L-32!*L^-31],

so K'(L)>0 once nu>0, and K0>4. S32 has positive coefficients,
so the upper payment M/S32(K(L)) decreases. There is no polynomial
prefactor to be paid in this new cap argument.

For tau and the scalar roots, set S=S32. For y>=4,

    y*S'(y)-(3/2)S(y)
      >=-3/2-y/2+y^2/4>=0,

because every remaining coefficient is positive. Therefore
L^(3/2)/S32(cT*L) decreases throughout the tail, with cT*4096>4.
The same argument implies L^(1/2)/S32(sigma*L/2) decreases, and
1/S32(sigma*L) decreases. Thus the tau payment and every scalar
remainder envelope checked at the base remain valid forever. No
finite grid is used to infer any infinite-tail statement.

## 8. Original finite history, true increments and arbitrary positive radii

Use the original largest-eligible-radius algorithm with the single numeric
threshold tau from (15). For EVERY algorithmic dilation a_alg>=8n^4,
normalize a positive load-support point by its last minimum positive-
contribution radius and reflect its last center into the positive orthant.
Write alpha=a_alg/(a_alg+1), ccap=(a_alg*beta-1)/(a_alg+1), and u for
the last normalized center. The true uncaptured-increment measure mu_a
is dominated by the FULL finite corner envelope lambda_a.

On A and P, lambda_a=c_a*lambda exactly, c_a=alpha^-n>=1, since
D is contained in the finite cap by

    alpha-ccap*m-(1-beta*m)
      =[(beta+1)m-1]/(a_alg+1)>0.

An actual target item j whose increment meets A satisfies
ccap*r_j<=||u-w||<alpha+rhi<3. Since a_alg>2t, ccap>1, so r_j<3.
An actual source increment has r_i>=rho>4, and hence precedes that
SAME target j. The displacement margin is

    alpha*beta-kappa-2/(a_alg+1)
       =zeta-(beta+2)/(a_alg+1)>0,

using beta<2 and a_alg>=8n^4. Included pairs thus really capture the
target center in its own original eligibility test. Integrating included
pieces costs no more than their genuine item weights. Both weights in
(4) are <=1, so the actual weighted target row is <=arad*tau<=tau
on mu_a-almost every target. No new selection rule or independent
history is introduced by the weighted-cap estimate.

Let E=integral K d lambda d lambda and Delta_a be the TOTAL full
finite-cap deficit lambda_a-mu_a. Energy, degree and target mass on
the support scale as c_a^2, c_a and c_a. The sets A and P are disjoint.
The energy-difference identity bounds its two missing-mass terms by
c_a*M*k times the sum of missing masses on A and P, hence by
c_a*M*k*Delta_a. On the other hand the actual energy is at most
tau*c_a*lambda(A). Thus

    Delta_a >=[c_a^2*E-tau*c_a*lambda(A)]/[c_a*M*k]
            >=[E-tau*lambda(A)]/[M*k]
            >=lambda(A)*b/M>=Gamma*n>gamma*n.          (16)

This holds uniformly for every original finite history and every such
a_alg. Actual increments and their original item weights are unchanged.

For arbitrary f in L1 use |f|. Each point in its strict maximal level set
has an actual positive witness radius; its open a_alg-shrunk cube covers
the point and finitely many cover any compact subset. Selected dilates
are still true witnesses. Normalization and (16) apply to every positive
terminal radius. Inner regularity then gives the full strict level set.
There is no compact-support, scale-truncation, selector-convergence or
actual-increment iid hypothesis.

## 9. Entire cover/root cost and strict display

The original captured budget tau, terminal self-weight and ENTIRE cover
factor give at every finite a_alg the accepted finite inequality

    C_n <=[1+(beta+1/a_alg)^n/tau]
       *[tau+a_alg^-n+alpha^-n*J_a_alg(t)-gamma*n].

Only in this family of numeric bounds, already uniform over all finite
histories, let a_alg tend to infinity. The exact limit is

    C_n <=(1+A_n/tau)(tau+B_n-gamma*n)=Phi80.            (17)

For the scalar remainder, B_n<=T-n log(1-sigma)+1 and

    [A_n+2sqrt(A_n B_n)]/n
       <=exp(-sigma*L)+2sqrt(L)*exp(-sigma*L/2).

The new independent base bound is

    1/Nstar(4096)+1/S32(sigma*4096)
      +128/S32(sigma*4096/2)<1/10000,                  (18)

with 128=2sqrt(4096), about 1.4931e-38. Every term decreases
on the tail by Section 7. Combining (2), (3), (17) and (18) proves
the strict second inequality in (1). The positive root cross-term is
explicitly paid rather than omitted.

## 10. Comparison with the complete old Phi69 and Phi78 functions

This comparison uses positive lower buffers for B, not only display
coefficients. The original scalar functions are

    ell_t(rho)=1+1/(2t)-1/(2rho),
    I_n(t)=integral_1^t ell_t(rho)^(n-1) d rho/rho,
    B_n(t)=1-[(1+1/t)/2]^n+n log t-(n/2)I_n(t).

For n>=2,t>1, the change v=1/rho-1/t gives

    I_n(t)<= integral_0^(1-1/t)
                exp[-(n-1)v/2]/(v+1/t) dv
           <=2t/(n-1).

Since the core 1-[(1+1/t)/2]^n is nonnegative, n/(n-1)<=2,
and x=(1-sigma_old)L>=L/2 for the old tuples,

    B_n(t_old)>=n log t_old-2n/x
       >=T-n log(1-sigma_old)-4n/L
       >=T+sigma_old*n-4n/L.                           (19)

Their true positive root domains are part of their accepted theorems;
F(A,B-gamma_old*n)>=B-gamma_old*n. Hence

    Phi69>=T-n[317/1000-83/1000+4/L]
           =T-n[234/1000+4/L],       L>=2^30,
    Phi78>=T-n[599/2000-25/256+4/L], L>=8192.           (20)

The compact run checks at these respective base values

    5473/20000>234/1000+4/2^30,
    5473/20000>599/2000-25/256+4/8192.

The right sides decrease thereafter. Therefore the new COMPLETE Phi80
is strictly smaller than Phi69 throughout the latter's true domain and
strictly smaller than Phi78 throughout the latter's true domain. The
numeric comparison margins are approximately .03965 and .071318.
Equations (19)--(20) lower-bound old upper-bound formulas, not C_n.

No older branch receives a new validity threshold from this comparison.
The full accepted envelope still takes the minimum over every accepted
branch on its own domain. There is no claim of global optimality among
continuous parameters, all possible kernels, or all covering methods.

## 11. Actual provenance, independent new guards and handoff

The frozen R073 proof hash is
50d48b85f02a14f50614f0eea5cedbd7ccd68e9b2f7e1564fe8b41fcd2af9c4b.
It binds the accepted relative J/AC, first-profile and second-profile
local lemmas, whose original hashes are c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1,
67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030,
and 4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8.
Their abstract active guards, not their old parameter thresholds, are
used and freshly checked here.

The R078 bounded-loss Bernstein proof and R079 weighted-cap proof were
authored in the same independent clone. The present note includes their
operative arguments in Sections 4--5, so no missing CP/BE theorem or
empirical-target interface is imported.

Producer files read in full: corrected R080 case.py, best-case receipt
case_988_1000.json, README, registration and correction_registration;
also the R079 basis registration/amendment. REGISTRATION.json locks
the four R080 producer identities before the independent compact run.
The initial batch's cost-slack error and subsequent corrected replay
remain in the producer archive and are not rewritten here. The stale
cap_rate field is an obsolete unweighted diagnostic, not an active gate.

The new verify.py actually ran ONCE, using ordinary Python Fraction
arithmetic with Nstar=L^32/32!, 12 atanh terms, S32 and exact base
square root 64, and returned PASS on 18 compact independent guards.
It covers the new active geometry, six target-event budgets, profile
transport/variance, Bernstein, target mass/degree, weighted cap, tau,
deficit, full root remainder and both COMPLETE old-cost comparisons.
It does not execute the producer 53-gate code or rerun a parameter grid.
verification.json preserves exact fractions and display-only decimals.

Reproduce against the separately archived producer directory:

    python3 verify.py --compute-dir ../compute

Only this owned upper/round080/weighted_threshold directory was written.
No MATLAB session, n=7 calculation, lower-bound work, Git commit/push,
paper update or accepted-ledger mutation was performed by this agent.
Independent whole-chain review and main-agent acceptance remain pending.
