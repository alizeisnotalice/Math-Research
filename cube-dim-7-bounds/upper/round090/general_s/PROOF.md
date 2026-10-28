# R090 arbitrary fixed s>=4: complete free-threshold family and limsup

Status: independently derived complete proof candidate for main-agent
review. Every s and all parameter gaps are fixed before n grows. Local
constants are rederived below and depend on s; no s<=5 constant or
threshold is imported. This independent staging directory is not a clone.
Only this directory is written. Earlier results and the paper are unchanged.

Define

    mu(s)=(s-1)+(s+1)log((s+1)/(2s)),
    C_free(s)=log[s/(1+mu(s))], T_n=n log n-n log log n.

For EVERY fixed real s>=4 and EVERY real c<C_free(s), there is one
fixed strict tuple and finite L_*(s,c) such that for every integer n with
L=log n>=L_*(s,c), the original Lebesgue centered axis-parallel
common-radius cube maximal operator, for arbitrary L1 inputs and all
positive witness radii, satisfies

    C_n <= H_tau=(1+A_n(t)/tau)(tau+B_n(t)-gamma*n)
        < T_n-c*n.                                    (1)

Here t=n/[(1-sigma)L] and tau=n exp(-dL), with fixed gamma,sigma,d
specified below. H_tau is the complete original cover cost at this same
algorithmic tau. A square-root F bound is not asserted at these parameters.

Consequently, for every epsilon>0, all sufficiently large n satisfy

    C_n<T_n-[-log(1-log 2)-epsilon]n,
    limsup_(n->infinity)(C_n-T_n)/n <= log(1-log 2).    (2)

No explicit dimension threshold, finite-domain coefficient exactly
-log(1-log 2), uniform threshold over s or epsilon, or global covering
optimality is claimed.

## 1. Fixed parameter order and arbitrary-positive-tau legality

For a_s=(s+1)/(2s), H(a)=1-a+a log a>0 for 0<a<1, so
mu(s)=2s H(a_s)>0. First fix s>=4 and primary parameters satisfying

    chi>mu(s), 0<sigma<1, chi(1-sigma)<d<sigma,
    1/s<qmax<1, 0<gamma<log(s*qmax),
    c<gamma+log(1-sigma).                              (3)

All are independent of n, f, witness families and histories. Existence
for every c<C_free(s) is proved in Section 9.

Choose fixed trim Delta>0 with
rlo=1+Delta<rhi=s*qmax-Delta and log(rhi/rlo)>gamma.
Choose fixed gamma<m<log(rhi/rlo), gamma/m<Ret<1 and
d_mu=(chi-mu(s))/4>0. Choose fixed

    0<e<min(1/12,1/[2(s+1)])

small enough that

    (1-e)log(rhi/rlo)>m,
    (mu(s)+2d_mu)/(1-e)<chi,
    (1-e)Gamma>gamma,                                 (4)

where

    P0=1-2e, Ulo=P0(1-e)/(1+e), Uhi=1/Ulo,
    M=2(1+e)Uhi,
    b=2Ret(1-e)Ulo-2e, Gamma=m*b/M.

These choices exist by continuity at e=0: Gamma tends to m Ret>gamma
and mu+2d_mu=(mu+chi)/2<chi. In particular b>0, M>2, and

    nu=1/s-(1+1/s)e>1/(2s)>0.                        (5)

Choose fixed m_face between (mu+2d_mu)/(1-e) and chi. The same e
pays target mass, product/count calibration, radial, cap and tau errors.

The original finite-cover theorem (R027 and revision_v3/cover/audit)
is for EVERY t>1,a>2t,tau>0. A largest-radius eligible cube is chosen
under the actual condition sum_previous w_i 1_(a beta Q_i)(center)
<=tau, with frozen w_i=|I_i|/(2a r_i)^n>=0. The first-rejection
prefix and captured-subsum arguments need only positivity of tau.
There is no square-root, size or finite-a restriction on tau. The
square root is an optional scalar optimization after that theorem.

Thus tau=n exp(-dL)>0 is legitimate. For each n fix this input-
independent tau BEFORE the actual selections, keep it the same for all
finite a, and use it in every eligibility test and actual row payment.
No old history is held fixed while its threshold is changed afterward.

## 2. Geometry with constants depending on the fixed s

For sufficiently large L put

    x=(1-sigma)L, t=n/x, beta=1+1/t,
    zeta=n^-3, epsilon=zeta*t=1/(x*n^2), kappa=beta-zeta,
    h=L^-2, theta=exp(-chi*x), k=n*h*theta/rhi.

Use the original corner envelope

    D={||z||_inf<t: some z_i<1-beta*m(z)},
    m(z)=max(1,||z||_inf),
    d lambda=2^-n*m(z)^-n*1_D dz.

Let lambda0 be its full-annulus version. A0 has target radius
r in [rlo,rhi]; Pfull has source radius rho=tq with q in [1/s,qmax].
Let P=D intersect Pfull. On ALL of Pfull define

    R=sq-h, a0=(1-epsilon)q, g=R-a0,
    ell_R(z_i)=min(1,(R+kappa*rho-|z_i|)/(2R)),
    G(z)=product ell_R(z_i), bsoft(z)=min(1,theta/G(z)).

Extend bsoft by zero off Pfull and let arad=r/rhi on A0, zero off.
The retained A below gives the Borel kernel

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
           *1_{||z-w||_inf<kappa*rho}*arad(w)*bsoft(z).

Both weights lie in [0,1]; actual increments, weights and the original
largest-eligible-radius rule are unchanged.

For target radius r the complete matching interval and baseline are
q_±=(r±h)/s, qb=q_+, R(qb)=r. Once 2h<Delta the whole interval
lies in [1/s,qmax]. On all active matches,

    1<R in [r-2h,r]<s,
    rho>4(s+1), 0<epsilon<1, a0>0,
    0<g, g-(1-1/s)R=epsilon*q-h/s<=0,                 (6)
    ell_R>=(1+a0/R)/2>=(s+1)/(2s)>1/2.

These are eventual: epsilon*qmax/h is a fixed multiple of L/n^2,
t/s tends to infinity, and h tends to zero. Also beta<2, kappa>1,
t>2rhi and h<=.001 hold eventually. On the full source band R>0
and G>0 eventually, but inactive edge columns can have R<=1 and
are zero; no R>=1 local lemma is applied there. Both bands lie below
t and are disjoint, since rho>4(s+1)>rhi.

## 3. Fresh local lemma for every fixed s>=4

Fix s,t,h,epsilon,w during q differentiation. Sufficient local guards
are 1/s<=q<1, 0<h<=.001, rho>4(s+1), 1<=R<=s,
0<=epsilon<1, epsilon*q<=h/s, 0<=g<=R and
|w|<=r<=min(s,R+2h).
They hold on the actual active domain (6). Note epsilon<=h<1.
The exact identity kappa*rho=rho+a0 is used throughout.

For the actual pre-cap inverse capture profiles put

    J_j=(2rho)^-1 integral_capture phi(z)^j/ell_R(z) dz,
    phi=-log ell_R, J=J_0,
    m_w=(|w|-a0)_+, v=(g-m_w)_+,
    Q(R,u)=-log(1-u/(2R)), F(R,u)=2R Q(R,u)-u.

The two source edge layers of width g are disjoint because rho>s>=g.
On each layer, the substitution Q=-log(1-u/(2R)) gives exact
integrals g+F(R,g), RQ(R,g)^2 and (2R/3)Q(R,g)^3 for
1/ell, phi/ell and phi^2/ell, respectively. By reflection take w>=0.
If w<=a0 there is no deletion. Otherwise source capture deletes its
left interval of length m_w. Since m_w<=s and 2rho-g>7s+8,
the deletion cannot reach the opposite edge. If m_w>g it deletes
the near edge plus central length m_w-g, where ell=1 and phi=0.
Consequently also in the central-deletion regime |w|>R,

    J=1+[F(R,g)+F(R,v)-m_w]/(2rho),
    J1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).                 (7)

The -m_w term in J is essential when v=0. Boundaries have zero
coordinate mass. These formulas are not restricted to the baseline.

For 0<=u<=R, elementary calculus gives

    0<=Q<=log 2<1,
    |Q_R|<=1/R, Q_u<=1/R,
    0<=F<=2R, -2<=F_R<=0, 0<=F_u<=1.                (8)

The strict logarithm bound follows from exp(1)>1+1=2. The sign
F_R<=0 follows from -log(1-z)<=z/(1-z); its lower bound follows
from Q>=0. Here z=u/(2R)<=1/2. No smaller-s numerical constant
has been imported.

For A=F(R,g)+F(R,v)-m_w, -s<=A<=4s, so |A|<=4s.
Almost everywhere R'=s, g'=s-1+epsilon<=s, |m_w'|<=1,
and v'=0,g',or g'-m_w', so 0<=v'<=s. Therefore

    |A'|<=6s+1, rho'/rho=1/q<=s,
    |J-1|<=2s/rho<1/2,
    |J'|<=(4s^2+6s+1)/(2rho),
    |d_q log J|<=(4s^2+6s+1)/rho=:D_s/rho.          (9)

The affine functions and positive parts are Lipschitz; Q,F are smooth
near the closed relative domain. Thus J,J1,J2 are AC. Since J>1/2,
log J is AC too. Uniform a.e. bounds integrate through every onset.

The values satisfy J1<=R/rho<=s/rho and
J2<=2R/(3rho)<=s/rho. For u=g or v, |d_q Q|<=2s/R.
Using the exact prefactor derivative (R/rho)'=h/(q*rho),

    |J1'|<=h/(q*rho)+4s/rho<=5s/rho.                (10)

All new constants depend only on the already fixed s. Equations
(7)--(10) hold throughout full matching, including central deletion.

## 4. Full-target laws and exactly five good events

At qb, R=r. On a full target face the n-1 free coordinates are
independent Uniform[-r,r] and one coordinate is its true deterministic
pin. Full-target Fubini gives E_W J=1 because the target capture
probability is ell_R(z), canceling its reciprocal. With (9),

    -4s^2/rho_b^2<=E_W log J<=0,
    Var_W(log J)<=16s^2/rho_b^2.

The true pin contributes absolute log loss <=4s/rho_b. Thus the
full log sum has absolute mean <=4s^2 n/rho_b^2+4s/rho_b and
variance <=16s^2 n/rho_b^2. These vanish uniformly over the fixed
target band and signed faces: n/rho_b^2<=s^2*x^2/n and
1/rho_b<=s*x/n.

The same Fubini cancellation for J1 yields exactly

    E_W J1=(2R/rho)H(a), a=1-g/(2R)>=a_s,
    H(a)=1-a+a log a.

Since H decreases on [a_s,1] and R/q<=s, the free sum mean is
<=mu(s)x. The actual pin adds <=s/rho_b<=s^2*x/n, eventually
<=d_mu*x/2. The free sum variance is <=n*s^2/rho_b^2
<=s^4*x^2/n. Chebyshev therefore gives

    P(sum J1(qb)>(mu+d_mu)x)<=4s^4/(d_mu^2 n) ->0.    (11)

Define A as A0 intersected with precisely FIVE conditions:
|sum log J(qb)|<=e; the two count deviations
|Nsum(q_±)-E_±|<=e*n/4; w in D; and
sum J1(qb)<=(mu+d_mu)x. Here
Nsum=sum_i[1+1_{|w_i|<a0(q)}], E_±=n+(n-1)a0(q_±)/r.
Eventually a0<r, so the true pin contributes exactly one. Free excess
indicators are independent Bernoulli(a0/r), giving each count failure
<=4/(e^2 n). The log failure tends to zero by its mean/variance.

For D, negative pinned target faces lie in D once t>rhi. On positive
faces, the cap-complement free fraction is

    1+1/(2t)-1/(2r)<=1-1/(4rhi)<1-1/(4s)<1

when t>=2rhi. Thus its failure is at most
(1-1/(4s))^(n-1), which tends to zero for this fixed s. Each of
the FIVE failures is eventually <=e/8 uniformly in radius and face.
Their mutual independence is unnecessary. Since lambda=lambda0 on D,

    lambda(A)>=(1-e)n log(rhi/rlo)>m*n.                (12)

These are ambient full-target laws only, not laws of actual increments.

## 5. Whole matching, true source variance and retention

The q width is 2h/s. Since n/rho<=s*x, (9) transports the
baseline log product by at most 2D_s*x*h. For fixed s this is
O(1/L) and eventually <=e. Together with the baseline event it
gives P0<=product J_i<=1/P0. Separately, pointwise (9) gives
every J_i in [1-e,1+e] eventually.

ellpin=(1+a0/R)/2 has derivative -(1-epsilon)h/(2R^2)<0,
and Nsum is nondecreasing. The q_+ mean ratio Nsum/(2n ellpin)
is <=1, while at q_- it is >=1-2h-1/n. The latter uses
a0<=1, R,r>=1 and R(q_-)=r-2h. Count fluctuations alter the
ratio by at most e/4. Once 2h+1/n<=3e/4, throughout matching

    1-e<=Nsum/(2n ellpin)<=1+e,
    U=(product J_i)/(2ellpin)*sum_i N_i/J_i,
    Ulo*n<=U<=Uhi*n.                                  (13)

All captured signed source faces and omitted-coordinate terms appear
in U. The first-profile transport from (10) is at most
10s*x*h. Once 10s*h<=d_mu,

    sum J1(w_i;q)<=(mu+2d_mu)x                       (14)

throughout complete matching.

Fix ANY captured signed source face BEFORE source cap and soft
conditioning. Its free normalized inverse-capture joint density is
the product of 1_capture/(2rho J_i ell_R), i!=j. Hence the free
coordinates are actually independent in this law. The pin loss is
deterministic and <=log2<p=1. With V=-log G,

    EV<=1+[(mu+2d_mu)/(1-e)]x<1+m_face*x,
    Var V=sum_(i!=j)[J2_i/J_i-(J1_i/J_i)^2]
         <=n*s/[rho(1-e)]<=2s^2*x.                   (15)

This is true centered variance. No pin variance, mixture variance,
whole-sum raw second moment or actual-increment law is substituted.
For fixed a_face=chi-m_face>0,

    P_inv(V>chi*x | face)
       <=2s^2*x/(a_face*x-1)^2 ->0                   (16)

once a_face*x>1. Thus every source face retains at least the fixed
Ret for all sufficiently large L, uniformly over retained targets,
matching q, faces and signs. Only after this facewise estimate may
faces be mixed. Since bsoft/(theta/G)=min(1,G/theta) is one on
V<=chi*x, the fullsource shell row is between
Ret*theta*Ulo*n and theta*Uhi*n.

## 6. Fixed-s weighted cap and entire column

By (13) and a0/R>=1/s, at least nu*n coordinates satisfy
|w_i|<a0. The positive nu in (5) was chosen with s dependence;
no uniform-in-s error is assumed. The source cap complement is exactly
all z_i>=-rho+(1-q). A negative source pin is excluded. Every
counted free coordinate captures its entire lower strip of length
1-q, using |w_i|<a0 and rho>4(s+1). Its pre-cap inverse-law
probability is >=(1-q)/[2rho(1+e)]. A positive pin omits at most
one counted coordinate. Independence before conditioning gives

    P_inv(D^c | face)<=exp[-(nu-1/n)*eta*x/(1+e)],
    eta=(1-qmax)/(2qmax)>0.                            (17)

Multiply bsoft<=theta/G by this bound, sum all face coefficients,
and integrate only over the actual matched dq/q interval. The radial
identity (r/rhi)log((r+h)/(r-h)) is uniformly within 1±e of
2h/rhi for small h. Division by k cancels theta. Therefore

    caploss/k<=M exp[-(nu-1/n)*eta*x/(1+e)] ->0.       (18)

For large n, nu-1/n>=nu/2. The resulting exponent is fixed positive
for this tuple, so cap loss is eventually <=e*k. There is no chi<eta
requirement and no independence asserted after conditioning. Raw rows
before cap deletion are >=2Ret(1-e)Ulo*k and <=M*k.

For the complete column, a fixed matched source has target fractions
ell_i(r)>1/2 for r>=R, decreasing in r, and <=ell_i(R).
The source radius rho>4(s+1)>r ensures the usual one-boundary
coordinate formula. At most one target pin sign is captured when
ell_i(r)<1, and at most two otherwise. Summing ALL target faces gives
<=n product ell_i(r)<=nG. arad changes dr/r to dr/rhi;
ordinary target matching width is <=2h and bsoft*G<=theta.
The entire column is <=2k<M*k, including every target face.
Restrictions to A,P only decrease upper degrees; inactive columns vanish.

## 7. Same free tau, fixed-s finite history and all positive scales

The exact free-threshold payment envelope is

    tau/k=rhi*L^2*exp[-(d-chi(1-sigma))L] ->0.         (19)

Its rate is positive by (3). This SAME global eligibility tau is
eventually <=e*k. After cap and tau are paid, row minus tau is
>=b*k. Both degrees are <=M*k and lambda(A)>m*n.

For every original a_alg>=8n^4, normalize by the last minimum
positive-contribution radius and reflect its center into the positive
orthant. Let alpha=a_alg/(a_alg+1) and
ccap=(a_alg*beta-1)/(a_alg+1). Eventually a_alg>2t and ccap>1.
An actual target item j meeting A satisfies

    ccap*r_j<=||u-w||<alpha+rhi<s+1,

so r_j<s+1. An actual source increment has
r_i>=rho>4(s+1)>s+1 and precedes this SAME j in the original
largest-radius history. This explicitly replaces the bounded-s target
and source constants by fixed-s bounds. The displacement margin is

    alpha*beta-kappa-2/(a_alg+1)
       =n^-3-(beta+2)/(a_alg+1)>0,

using beta<2 and a_alg>=8n^4. Since r_j<r_i, actual center
distance is <kappa*rho+(r_i+r_j)/(a_alg+1)
<alpha*beta*r_i. Included sources really capture the genuine target
center in its original eligibility test with this SAME tau. Both
kernel weights are <=1, so actual weighted row<=arad*tau<=tau.

On A and P the finite corner envelope equals c_a*lambda exactly,
c_a=alpha^-n>=1, since its cap margin is
[(beta+1)m(z)-1]/(a_alg+1)>0. Actual uncaptured-increment measure
is dominated by that finite envelope. Energy, degree and target mass
scale by c_a^2,c_a,c_a. A and P are disjoint; the two missing-mass
terms in the energy difference cost at most c_a*M*k times their sum,
at most the total finite-cap deficit Delta_a. Actual energy is
<=tau*c_a*lambda(A). Thus with E=integral K d lambda d lambda,

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
            >=[E-tau*lambda(A)]/(M*k)
            >=lambda(A)*b/M>Gamma*n>gamma*n.           (20)

This is uniform over all original finite histories at each sufficiently
large n and all these finite a. There is no gamma<1 constraint, actual-
increment probability law, altered selector or unpaid finite-a factor.

For general f in L1 use |f|. Each point of a strict maximal level set
has a genuine positive witness radius. Open shrunk witness cubes give
finite covers of compact subsets; selected dilates remain true witnesses.
Normalization allows every positive terminal radius. Inner regularity
yields the entire strict level set without compact-support, radius
truncation, scale lattice, history convergence or iid assumptions.

For every fixed s and tuple, all guards have finitely many explicit
envelopes: h->0, t->infinity, polynomial(L)/n->0, geometric target-
cap decay with base 1-1/(4s)<1, transports O_s(1/L), O_s(1/L^2),
(16)->0, and positive-rate cap/tau exponentials. Their fixed positive
budgets give one finite maximum threshold. No uniformity over varying
s or vanishing gaps is needed or asserted.

## 8. Entire original H_tau cost and sublinear remainder

Keep the captured budget, terminal self-weight and ENTIRE finite cover:

    C_n <=[1+(beta+1/a_alg)^n/tau]
          *[tau+a_alg^-n+alpha^-n J_a_alg(t)-gamma*n].

For each fixed n,t,tau let only a_alg tend to infinity in these numeric
bounds for the same original operator. The uniform finite-history deficit
(20) gives exactly the H_tau bound in (1). No selections must converge.

The original scalar functions are

    A_n(t)=(1+1/t)^n,
    B_n(t)=1-[(1+1/t)/2]^n+n log t
           -(n/2) integral_1^t
            [1+1/(2t)-1/(2rho)]^(n-1) d rho/rho.

B_n>=n log t/2>=nL/4>gamma*n eventually. Set Q=B_n-gamma*n>0.
The exact complete cost is

    H_tau=Q+tau+A_n+A_n Q/tau.

With P=-log(1-sigma), the original bounds give
A_n<=exp x=n exp(-sigma L), B_n<=T_n+P*n+1 and B_n<=nL
eventually. Therefore

    H_tau <= T_n+(P-gamma)n
       +n[exp(-dL)+exp(-sigma L)
          +L exp(-(sigma-d)L)+1/n].                   (21)

All terms in brackets tend to zero because 0<d<sigma. The positive
fixed allowance gamma-P-c pays this ENTIRE remainder and proves (1).
No cross term has been omitted. The algebraic fact F(A_n,Q)<=H_tau
does not imply C_n<=F(A_n,Q), so no square-root optimization is
made after the new-tau defect proof.

## 9. Strict finite choices for every c<C_free(s)

Fix s>=4 and c<C_free(s), and put D=C_free(s)-c>0. Choose
chi>mu so close that log(1+chi)<log(1+mu)+D/4. Choose
sigma>chi/(1+chi) so close that P<log(1+chi)+D/4.
Choose fixed d strictly between chi(1-sigma) and sigma. Choose
qmax close to one with qmax>1/s and log qmax>-D/8, and positive
gamma<log(s*qmax) with gamma>log(s*qmax)-D/8. Then

    gamma>log s-D/4,
    P<log(1+mu)+D/2,
    gamma-P-c>D/4>0.

Complete the trim, mass, Ret and s-dependent error choices of Section 1
before n grows. Every strict choice exists by continuity. Thus each
fixed s,c has a genuine fixed finite construction and its own threshold.
No boundary parameter or n-dependent s is inserted into the proof.

## 10. Finite s first, then the limiting coefficient

An exact rearrangement gives

    (1+mu(s))/s
      =1-(1+1/s)log2+(1+1/s)log(1+1/s).

Since 0<log2<1 and 0<=log(1+1/s)<=1/s, this positive ratio
tends to 1-log2>0. For example its difference from that limit has
absolute value <=(2+1/s)/s. Continuity of the logarithm therefore gives

    C_free(s) -> c_infty=-log(1-log2) as s->infinity.   (22)

Now fix epsilon>0. By (22), choose ONE finite s>=4 with
C_free(s)>c_infty-epsilon/2. Set c=c_infty-epsilon, which is
strictly smaller than C_free(s). Section 9 chooses ONE fixed strict
tuple for this s,c, and Sections 1--8 yield its own finite eventual
threshold. Consequently (C_n-T_n)/n< -c_infty+epsilon eventually.
The definition of limsup, with arbitrary positive epsilon, proves (2).

This quantifier order never substitutes s=infinity into a local lemma,
finite kernel, finite history or operator bound. It does not demand a
threshold common to all epsilon. The limit is a consequence of actual
finite-s theorems, not a formal limit of a model coefficient.

## 11. Actual finite mathematical auxiliary and provenance

REGISTRATION.json locks the original arbitrary-positive-tau sources
before one bounded auxiliary run. check.py ran with PID 32534 and PASS.
It uses exact Fraction polynomial bookkeeping for D_s=4s^2+6s+1,
the s-dependent count floor and variance factor, independent log24
enclosures proving 0<log2<1 and enclosing c_infty in
(1.18138,1.18140), and four direct quadrature cases at the single
fixed diagnostic s=10,q=3/5. There g/R>.8, outside the old s<=5
relative domain. Full capture, edge deletion and signed central deletion
agree with all three exact profile formulas to maximum error 2.22e-16.
Finite differences are diagnostics only. receipt.json records exact
values, PID, script/runtime identities and raw cases.

The computation does not prove a uniform local lemma or a dimension
threshold; the analytic proof supplies both the fixed-s reasoning and
the epsilon/limsup quantifiers. No old grid, baseline or author check
was replayed, and no shared MATLAB, commit or push was used. The
registered original source hashes were rechecked after the run and at
packaging. Main-agent review and acceptance retain final proof responsibility.
