# R088 free global threshold: fixed-s family and a 0.9n instance

Status: independently derived proof candidate for main-agent review.
This branch is restricted to fixed s in [4,5]. It changes the threshold
choice in the actual original algorithm, retains the complete H_tau cost,
and does not claim the square-root F bound at the new parameters.
The staging directory is independent, not a Git clone. Source identities
are registered; no previous source, result, paper or Git state is changed.

Let mu(s)=(s-1)+(s+1)log((s+1)/(2s)) and T_n=n log n-n log log n.
For every fixed s in [4,5] and every

    c < C_free(s):=log[s/(1+mu(s))],                    (1)

there is one fixed parameter tuple and finite L_* such that, for every
integer n with L=log n>=L_*, the original Lebesgue centered axis-parallel
common-radius cube maximal operator for arbitrary L1 inputs and all
positive witness radii satisfies

    C_n <= H_tau=(1+A_n(t)/tau)(tau+B_n(t)-gamma*n)
        < T_n-c*n.                                    (2)

Here the SAME global threshold used to build every finite selection is

    tau=n exp(-d L), x=(1-sigma)L, t=n/x,               (3)

with fixed d, sigma, gamma chosen below. No explicit L_* is certified.
There is also an explicit fixed s=9/2 tuple proving (2) with c=9/10
eventually. No s>5 extension or global optimum is included.

## 1. The original finite-a selection permits every tau>0

The original R027 theorem states its finite bound for fixed t>1,
a>2t and ANY tau>0. Its selection rule is: among remaining cubes
whose centers satisfy

    sum_(previous selected i) w_i 1_(b_alg Q_i)(center)<=tau,

select a largest radius, where b_alg=a(1+1/t),
I_i=Q_i minus previous selected small cubes, and
w_i=|I_i|/(2a r_i)^n are frozen nonnegative weights. Eligibility
can only be lost; the finite family terminates for every positive tau.
The first rejection prefix proves coverage with factor
1+(1+1/t+1/a)^n/tau. Nothing requires tau to be a square root,
bounded by one, related to a, or independent of n. It must be fixed
before that selection and independent of f, witnesses and history.

R027's captured previous terms cost at most this same tau by the actual
eligibility inequality. Its self term costs a^-n separately. The
square-root substitution appears only after the original finite theorem
and scalar a limit, as an optional numerical minimization. The detailed
revision_v3/cover/audit gives the same universal tau>0 quantifier.

Consequently (3), which is a positive input-independent real for each
fixed n, is legal. Fix n and this tau before taking any finite a or
any witness family. The same tau is used at all finite a, though the
selected families may vary. The new proof does not keep an old
square-root-generated history and then change its threshold afterward.

## 2. Strict fixed primary parameters and geometric errors

First fix s in [4,5] and

    chi>mu(s), 0<sigma<1, chi(1-sigma)<d<sigma,
    1/s<qmax<1, 0<gamma<log(s*qmax),
    c<gamma+log(1-sigma).                              (4)

The interval for d is nonempty exactly when sigma>chi/(1+chi).
All these constants are fixed before n grows. Choose a fixed trim
Delta>0 with rlo=1+Delta<rhi=s*qmax-Delta and
log(rhi/rlo)>gamma. Choose

    gamma<m<log(rhi/rlo), gamma/m<Ret<1,
    d_mu=(chi-mu(s))/4>0.

Choose fixed e in (0,1/12) so small that

    (1-e)log(rhi/rlo)>m,
    (mu(s)+2d_mu)/(1-e)<chi,
    (1-e)Gamma>gamma,

where

    P0=1-2e, Ulo=P0(1-e)/(1+e), Uhi=1/Ulo,
    M=2(1+e)Uhi,
    b=2Ret(1-e)Ulo-2e, Gamma=m*b/M.

Existence follows by continuity at e=0, where Gamma=m Ret>gamma.
Choose fixed m_face between (mu+2d_mu)/(1-e) and chi. The same
e pays all geometric, target-mass, cap and tau errors. In particular
nu=1/s-(1+1/s)e>1/10, b>0 and M>2.

For n large define beta=1+1/t, zeta=n^-3,
epsilon=zeta*t=1/(x*n^2), kappa=beta-zeta, h=L^-2,
theta=exp(-chi*x), k=n*h*theta/rhi. Use the original corner
envelope D, lambda and its full-annulus version lambda0:

    D={||z||_inf<t: some z_i<1-beta*max(1,||z||_inf)},
    d lambda=2^-n*max(1,||z||_inf)^-n*1_D dz.

A0 has target radius r in [rlo,rhi]; Pfull has source radius rho=tq
with q in [1/s,qmax], and P=D intersect Pfull. On the full source
band put R=sq-h, a0=(1-epsilon)q, g=R-a0,

    ell_R(z_i)=min(1,(R+kappa*rho-|z_i|)/(2R)),
    G=product ell_R, bsoft=min(1,theta/G).

Set arad=r/rhi on A0, zero elsewhere. After selecting A below use

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
           *1_{||z-w||_inf<kappa*rho}*arad(w)*bsoft(z).

Both weights are <=1; original item weights and eligibility rule
are unchanged. The matching interval q_±=(r±h)/s lies in the full
band once 2h<Delta, with baseline qb=q_+ and R(qb)=r. Uniformly
on matching, R in [r-2h,r] exceeds one, rho>20, r<5, and

    g-(1-1/s)R=epsilon*q-h/s<=0,
    ell_R>=(s+1)/(2s)>=3/5.

All these conditions, as well as beta<2, kappa>1, a0>0 and
t>2rhi, hold eventually. Inactive edge columns are zero; no local
R>=1 estimate is imposed there. R>0 and G>0 on Pfull eventually.
The target and source bands are disjoint and lie below norm t.

## 3. Local formulas, five events and face retention, independent of tau

The registered R085 uniform local lemma is used only for its geometric
and integral statements, not its old tau choice. It was independently
derived on s in [4,5] and is restated here with its proof ingredients.
For fixed t,h,epsilon,w in q differentiation, let

    m_w=(|w|-a0)_+, v=(g-m_w)_+,
    Q(R,u)=-log(1-u/(2R)), F(R,u)=2R Q(R,u)-u.

The pre-cap inverse-capture profiles satisfy exactly

    J=1+[F(R,g)+F(R,v)-m_w]/(2rho),
    J1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).

Each source boundary layer has width g and ell=1-u/(2R). Integrating
1/ell, phi/ell and phi^2/ell gives g+F, RQ^2 and (2R/3)Q^3.
Deletion has length m_w and cannot reach the opposite edge since
m_w<=5, rho>20 and g<=4. If |w|>R, deletion enters the central
ell=1, phi=0 region, v=0 and -m_w remains in J. Thus all formulas
hold throughout matching. Positive parts are Lipschitz, so the
profiles and log J are AC through onset kinks. On h<=.001,rho>=20,
the independent bounds imply the following safe uniform estimates:

    |J-1|<9/rho, |d_q log J|<120/rho,
    J1<2/rho, |d_q J1|<60/rho, J2<1/rho.             (5)

For example Q<.6, |Q_R|<=2/(3R), Q_u<=5/(6R),
-4/3<=F_R<=0 and F_u<=2/3. The affine and positive-part
derivatives give stronger |d_q log J|<=243/(4rho) and
|d_q J1|<9.002/rho. These are independent of sigma, d and tau.

At qb, full-target-coordinate Fubini cancels ell_R, yielding
E_W J=1 and E_W J1=(2R/rho)H(a), H(a)=1-a+a log a,
a>=a_s=(s+1)/(2s). The log baseline sum has absolute mean
<=81n/rho_b^2+18/rho_b and variance <=324n/rho_b^2.
Also its J1 sum mean is eventually <=(mu+d_mu/2)x and variance
<=100x^2/n. Thus the J1 bad probability is <=400/(d_mu^2 n).
The pin terms are deterministic; n/rho_b^2<=25x^2/n and
1/rho_b<=5x/n pay their errors uniformly.

Choose A by exactly five conditions: |sum log J(qb)|<=e;
the two endpoint deviations |Nsum(q_±)-E_±|<=e*n/4;
w in D; and sum J1(qb)<=(mu+d_mu)x. Here
Nsum=sum_i[1+1_{|w_i|<a0}], E_±=n+(n-1)a0(q_±)/r.
The target pin contributes exactly one; the free excess counts
are independent Bernoulli(a0/r), so each count failure is
<=4/(e^2 n). Negative pinned target faces lie in D; positive
faces have cap-complement probability <=(.95)^(n-1), since
rhi<5 and t>2rhi. The log failure vanishes by its mean/variance.
Each of the five failures is eventually <=e/8; their mutual
independence is unnecessary. Consequently

    lambda(A)>=(1-e)n log(rhi/rlo)>m*n.                (6)

Over matching, log-product transport is <=240xh<=240/L and
J1 transport <=120xh. These are eventually <=e and d_mu*x.
Individual J_i lie in [1-e,1+e], and product J_i in [P0,1/P0].
ellpin=(1+a0/R)/2 decreases with derivative
-(1-epsilon)h/(2R^2), while Nsum increases. The q_+ mean
ratio Nsum/(2n ellpin) is <=1 and the q_- mean ratio is
>=1-2h-1/n; fluctuations change it by <=e/4. Thus throughout
matching, once 2h+1/n<=3e/4,

    1-e<=Nsum/(2n ellpin)<=1+e,
    U=(product J_i)/(2ellpin)*sum_i N_i/J_i,
    Ulo*n<=U<=Uhi*n.                                  (7)

This includes all captured signed source faces and omitted coordinates.

For any captured source face before cap/soft conditioning, its
free-coordinate joint inverse law factors as
product_(i!=j) 1_capture/(2rho J_i ell_R). Its pin loss is
deterministic and <p=.6. True centered variance and mean satisfy

    EV<=p+m_face*x,
    Var V=sum_(i!=j)[J2_i/J_i-(J1_i/J_i)^2]<10x.

With fixed a_face=chi-m_face>0,
P_inv(V>chi*x | face)<=10x/(a_face*x-p)^2 ->0. Every
face therefore retains at least Ret eventually, uniformly on A
and matching. Only after this uniform facewise estimate may faces
be mixed. The soft ratio min(1,G/theta) equals one on V<=chi*x.
No Gaussian, mixture-variance or actual-increment law is assumed.

## 4. Weighted cap and degrees before paying the new tau

Equation (7) with a0/R>=1/s gives at least nu*n counted
coordinates. In source D-complement all z_i>=-rho+(1-q).
Negative source pins are excluded. Each counted free coordinate
captures its whole lower strip of length 1-q, with inverse
probability >=(1-q)/[2rho(1+e)]. A positive pin omits at most
one count. Independence before conditioning gives

    P_inv(D^c | face)<=exp[-(nu-1/n)*eta*x/(1+e)],
    eta=(1-qmax)/(2qmax)>0.

Multiply by bsoft<=theta/G, sum all full face coefficients, and
integrate only over the matched dq/q interval. The exact radial
factor (r/rhi)log((r+h)/(r-h)) is within 1±e of 2h/rhi.
After dividing by k, theta cancels and

    caploss/k<=M exp[-(nu-1/n)*eta*x/(1+e)] ->0.       (8)

The rate is fixed positive once nu-1/n>=nu/2. No chi<eta condition
is needed. The cap payment is eventually <=e*k. Before this deletion,
raw rows are >=2Ret(1-e)Ulo*k and <=M*k.

For a fixed source and matched r>=R, target capture fractions exceed
1/2 and decrease in r. At most one target pin sign is captured if
the fraction is <1, otherwise at most two. All target faces sum
to <=n product ell_i(r)<=nG. arad cancels dr/r, matching width
is <=2h, and bsoft*G<=theta, giving complete columns <=2k<M*k.
Restrictions only lower upper degrees; inactive columns vanish.

## 5. New tau payment and original-history bridge

For the new threshold (3), the exact envelope is

    tau/k=rhi*L^2*exp[-(d-chi(1-sigma))L] ->0.         (9)

The rate is positive by (4). Therefore this SAME global algorithmic
tau is eventually <=e*k. After (8) and (9), row minus tau is
>=b*k; both degrees are <=M*k, and target mass is >m*n.

For each original finite a_alg>=8n^4 normalize by the last minimum
positive-contribution radius and reflect the last center. Since
a_alg>2t eventually, ccap>1. An actual target meeting A has
ccap*r_j<alpha+rhi<6, so r_j<6. An actual source increment
has r_i>=rho>20>6 and precedes that SAME target in this new-tau
selection history. The displacement margin

    alpha*beta-kappa-2/(a_alg+1)
       =n^-3-(beta+2)/(a_alg+1)>0

and r_j<r_i imply actual center distance
<kappa*rho+(r_i+r_j)/(a_alg+1)<alpha*beta*r_i.
Every included source is eligible in that target's own original
test. Bounded weights give actual weighted row<=arad*tau<=tau.
This is exactly the tau from (3), not a post-selected threshold.

On A and P the finite envelope equals c_a*lambda, c_a=alpha^-n>=1,
by the cap margin [(beta+1)m(z)-1]/(a_alg+1)>0. Actual increment
measure is dominated by it. Energy, degrees and mass scale as
c_a^2,c_a,c_a. The two missing-mass terms on disjoint A and P
cost at most c_a*M*k times their sum, at most the total deficit
Delta_a. Actual energy is <=tau*c_a*lambda(A). Hence

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
            >=[E-tau*lambda(A)]/(M*k)
            >=lambda(A)*b/M>Gamma*n>gamma*n.           (10)

All estimates are uniform over all original finite histories at this
n and tau, and over a_alg>=8n^4. No actual-increment product law
is assumed. The original compact-subset cover by open shrunk true
witnesses, selected genuine dilates, arbitrary terminal normalization
and inner regularity recover general L1 inputs and every positive
radius. There is no compact-support, radius truncation, scale lattice,
or convergence of histories.

All guards use finitely many fixed budgets and explicit envelopes:
polynomial(L)/n, transports O(1/L), O(1/L^2), face failure O(1/x),
and positive-rate exponential cap/tau payments. Thus a common finite
eventual threshold exists for each fixed tuple. None is computed.

## 6. Full H_tau cost without square-root optimization

Retain the original captured budget, self term and entire finite cover:

    C_n <=[1+(beta+1/a_alg)^n/tau]
          *[tau+a_alg^-n+alpha^-n J_a_alg(t)-gamma*n].

For fixed n, t and the already chosen tau, take only the scalar
a_alg limit of these inequalities for the same operator. Uniform
finite-history deficit (10) gives precisely (2)'s H_tau bound.
No selected family must converge.

Eventually B_n>=n log t/2>=nL/4>gamma*n. Expanding the complete
cost, with Q=B_n-gamma*n>0, gives

    H_tau=Q+tau+A_n+A_n Q/tau.

The original scalar estimates A_n<=exp x=n exp(-sigma L),
B_n<=T_n+P*n+1, B_n<=nL, with P=-log(1-sigma), imply

    tau/n=exp(-dL), A_n/n<=exp(-sigma L),
    A_n Q/(tau*n)<=L exp[-(sigma-d)L].                (11)

Both rates in (11) are positive by (4). Consequently

    H_tau <=T_n+(P-gamma)n
             +n[exp(-dL)+exp(-sigma L)
                 +L exp(-(sigma-d)L)+1/n].            (12)

The bracket tends to zero. Its eventual size can be paid from the
fixed positive allowance gamma-P-c, proving H_tau<T_n-c*n.
This pays ALL cover and cross terms. Although algebraically
F(A_n,Q)<=H_tau, that direction does not turn C_n<=H_tau into
C_n<=F(A_n,Q). Such an F conclusion is not claimed; its optimizing
threshold need not satisfy the new defect-payment condition (9).

## 7. Every strict coefficient below log[s/(1+mu(s))]

Fix s in [4,5], c<C_free(s), D=C_free(s)-c>0. Choose in order:

1. chi>mu so close that log(1+chi)<log(1+mu)+D/4.
2. sigma>chi/(1+chi) so close that P<log(1+chi)+D/4.
3. d strictly between chi(1-sigma) and sigma, e.g. their midpoint.
4. qmax close enough to one that qmax>1/s, log qmax>-D/8.
5. positive gamma<log(s*qmax) with gamma>log(s*qmax)-D/8.

Then gamma>log s-D/4, P<log(1+mu)+D/2, so
gamma-P-c>D/4>0. Complete trim/mass/Ret/error choices of Section 2
before allowing n to grow. Sections 3--6 apply. No parameter drifts
with n and no boundary parameter is inserted into a finite construction.
This proves (1)--(2) for fixed s in [4,5]. No maximization or s>5
extension is needed for this branch.

## 8. One fully fixed s=9/2, c=9/10 instance

Use the exact tuple

    s=9/2, chi=4/5, sigma=9/20, d=89/200,
    qmax=9999/10000, Delta=1/1000000,
    rlo=1000001/1000000, rhi=4499549/1000000,
    m=1503/1000, Ret=99999/100000, e=1/100000,
    gamma=3/2, m_face=199/250, c=9/10.

d_mu=(chi-mu(9/2))/4 is fixed. Independent exact log24
enclosures verify

    0<mu(9/2)<chi,
    (1-e)log(rhi/rlo)>m>gamma,
    Ret>gamma/m,
    (mu+2d_mu)/(1-e)<m_face<chi,
    (1-e)Gamma>gamma,
    P=-log(11/20)<299/500,
    gamma-299/500-c=1/500>0.

The two new exponent rates are exactly

    d-chi(1-sigma)=1/200,
    sigma-d=1/200.

The source-face gap is chi-m_face=1/250. Thus its true-variance
Chebyshev envelope is 10x/((x/250)-3/5)^2 ->0. The cap rate
is positive. The new tau payment is rhi L^2 exp(-L/200)->0.
The entire normalized scalar remainder is bounded by
exp(-89L/200)+exp(-9L/20)+L exp(-L/200)+1/n, eventually
<1/500. Therefore, at one finite unquantified threshold,

    C_n <= H_{n exp(-89L/200)} < T_n-(9/10)n.

This is an actual fixed-tuple eventual statement, not the boundary
coefficient, a numerical model or an optimized F statement.

## 9. Actual bounded computation and provenance

REGISTRATION.json fixes all four source identities and the new exact
tuple before one bounded local Python run. check.py ran with PID
32025 and PASS on 13 exact arithmetic checks. receipt.json preserves
log24 rational bounds, exact rates, mass/face/deficit slack, and the
complete H_tau allowance. All source hashes were rechecked on completion.
No author baseline, old run or parameter scan was replayed; no dimension
threshold was evaluated. The Fraction evidence supports fixed arithmetic;
the analytic argument above supplies the eventual uniform theorem.

Original finite selection constraints were checked directly in R027
and revision_v3/cover/audit, not inferred from a later optimized formula.
The unchanged local geometry is locked to the R085 uniform proof and
the independent local_uniform.md source. Final acceptance of this new
candidate and any old-branch comparisons remains with the main agent.
