# R085 proposed parameter-family theorem and its analytic limiting coefficient

Status: complete proposed eventual parameter-family theorem, pending
independent whole-chain review. Every construction below uses fixed
parameters with strict positive slack. No dimension threshold, finite-
domain bound with the limiting coefficient, or global covering optimum
is claimed. Prior accepted files and the paper are unchanged.

For s in [4,5] define

    mu(s)=(s-1)+(s+1)log((s+1)/(2s)),
    C(s)=log[s/(1+2mu(s))],
    T_n=n log n-n log log n.

For EVERY fixed s in [4,5] and EVERY real c<C(s), there is a finite
L_*(s,c) such that for every integer n with log n>=L_*(s,c), the
original Lebesgue centered axis-parallel common-radius cube maximal
operator, for all L1 inputs and positive witness radii, satisfies

    C_n<=F(A_n(t),B_n(t)-gamma*n)<T_n-c*n,              (1)

for one FIXED tuple constructed below; t=n/[(1-sigma)log n].
Here F(A,B)=(sqrt A+sqrt B)^2 and gamma, sigma do not vary with n.

The maximum of C(s) on [4,5] is attained uniquely at

    s_*=1/(2 exp(-1/2)-1),
    c_*=-log(3-4 exp(-1/2))=0.5553395620244... .        (2)

Consequently, for every epsilon>0, all sufficiently large n satisfy

    C_n<T_n-(c_*-epsilon)n,
    limsup_(n->infinity) (C_n-T_n)/n <= -c_*.           (3)

The two occurrences of C (the maximal constant C_n and scalar function
C(s)) have their displayed distinct arguments throughout.

## 1. A fixed-tuple theorem with explicit parameter slack

First fix s in [4,5] and numbers satisfying

    chi>mu(s),
    0<sigma<1, sigma/2>chi(1-sigma),
    1/s<qmax<1,
    0<gamma<log(s*qmax),
    c<gamma+log(1-sigma).                              (4)

The construction will prove (1) for these fixed primary parameters.
mu(s)>0: if a_s=(s+1)/(2s), then 0<a_s<1 and

    mu(s)=2s H(a_s), H(a)=1-a+a log a>0,

since H(1)=0 and H'(a)=log a<0 for a<1. Thus chi>0.

Choose a FIXED trim delta>0 so small that

    rlo=1+delta<rhi=s*qmax-delta,
    log(rhi/rlo)>gamma.

Choose a mass constant m with gamma<m<log(rhi/rlo), and then choose
a fixed retained fraction Ret in (gamma/m,1). Put

    d_mu=(chi-mu(s))/4>0.

Now choose one FIXED error e in (0,1/12) sufficiently small that

    (1-e)log(rhi/rlo)>m,
    (mu(s)+2d_mu)/(1-e)<chi,
    (1-e)Gamma>gamma,                                 (5)

where

    P0=1-2e,
    Ulo=P0(1-e)/(1+e), Uhi=(1+e)/[P0(1-e)],
    M=2(1+e)Uhi,
    b=2Ret(1-e)Ulo-2e,
    Gamma=m*b/M.

Such an e exists. The first two inequalities are strict at e=0,
and (1-e)Gamma tends to m*Ret>gamma as e tends to zero. In
particular b>0 and M>2 for small positive e. The same e is used for
eB,eN,eJ,eRad,eS,eT,eD and epsA; it is not subsequently made
n-dependent. The count constant

    nu=1/s-(1+1/s)e

is positive, since e<1/12 and s<=5 imply nu>1/10.

Finally choose m_face strictly between (mu(s)+2d_mu)/(1-e) and chi.
All of these numbers are fixed before n is allowed to grow. The scalar
allowance

    kappa_scalar=gamma+log(1-sigma)-c>0                (6)

is also a fixed positive number.

## 2. Uniform geometry and the exact bounded kernel

For n with L=log n large, put

    qmin=1/s, x=(1-sigma)L, t=n/x,
    beta=1+1/t, zeta=n^-3, epsilon=zeta*t=1/(x*n^2),
    kappa=beta-zeta, h=L^-2,
    theta=exp(-chi*x), k=n*h*theta/rhi.

Use the original corner envelope

    D={||z||_infinity<t: some z_i<1-beta*m(z)},
    m(z)=max(1,||z||_infinity),
    d lambda=2^-n*m(z)^-n*1_D dz.

Write lambda0 for the full-annulus density before restricting to D.
Let A0 have target radius r in [rlo,rhi], and Pfull have source radius
rho=tq with q in [1/s,qmax]. Put P=D intersect Pfull and define on
ALL of Pfull

    R=sq-h, a0=(1-epsilon)q, g=R-a0,
    ell_R(z_i)=min(1,(R+kappa*rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i), bsoft=min(1,theta/G).

Extend bsoft by zero off Pfull, and let arad=r/rhi on A0, zero off.
The final Borel kernel, once A is specified, is

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
       *1_{||z-w||_infinity<kappa*rho}*arad(w)*bsoft(z). (7)

Both weights lie in [0,1]. Actual increments, their item weights and
the original largest-eligible-radius algorithm are unchanged.

For each target r the COMPLETE matching interval and baseline are

    q_-=(r-h)/s, q_+=(r+h)/s, qb=q_+, R(qb)=r.

For all sufficiently large L this interval lies in [qmin,qmax] and
every ACTIVE match satisfies

    1<R in [r-2h,r]<5, rho>20,
    0<epsilon<1, a0>0, beta<2, kappa>1, t>2rhi,
    0<g, g-(1-1/s)R=epsilon*q-h/s<=0.                (8)

Indeed 2h<delta, epsilon*qmax<=h/s, t*qmin>20 and the other
guards all hold eventually. In particular

    g/R<=1-1/s<=4/5,
    ell_R(z_i)>=(1+a0/R)/2>=(s+1)/(2s)>=3/5.          (9)

Inactive edge sources may have R<=1 but have zero final columns;
no local R>=1 lemma is applied there. R>0 and G>0 on Pfull for
large L. Both bands lie inside norm<t and are disjoint.

## 3. Fresh local bounds uniform on s in [4,5]

Fix a target coordinate |w|<=r and define its actual pre-cap inverse-
capture profiles

    J_j(w;q)=(2rho)^-1 integral_capture phi(z)^j/ell_R(z) dz,
    phi=-log ell_R, j=0,1,2, J_0=J.

Set m_w=(|w|-a0)_+, v=(g-m_w)_+,
Q(R,u)=-log(1-u/(2R)), F(R,u)=2R*Q(R,u)-u.
Direct integration on the two boundary strips gives

    J=1+[F(R,g)+F(R,v)-m_w]/(2rho),
    J_1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J_2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).                 (10)

The formulas also hold for |w|>R: central deletion beyond a boundary
strip has phi=0 and is included by -m_w. Positive-part onset kinks
are Lipschitz, so the profiles are AC and their a.e. derivatives
integrate over the whole matching interval.

For 0<=u<=4R/5, 1<=R<=5,

    Q<=log(5/3)<3/5,
    |Q_R|<=2/(3R), Q_u<=5/(6R),
    0<=F<=6, |F_R|<=38/15, 0<=F_u<=2/3.

The logarithm bound follows from exp(3/5)>1+3/5+(3/5)^2/2>5/3.
At fixed n, R'=s<=5, g'=s-1+epsilon<=5, |m_w'|<=1,
|v'|<=6, and rho'/rho=1/q<=s<=5. The numerator A in J has
|A|<=17 and |A'|<=101/3. It follows that

    |J-1|<=17/(2rho)<9/rho,
    |J'|<=178/(3rho)<60/rho,
    |partial_q log J|<=120/rho.                        (11)

The last estimate uses J>=1/2, ensured by rho>20. Independently,
the R', Q and rho' contributions to J_1' are bounded by
9/5, 95/2 and 9, divided by rho. Thus

    0<=J_1<=9/(5rho)<2/rho,
    |J_1'|<=583/(10rho)<60/rho,
    0<=J_2<=18/(25rho)<1/rho.                         (12)

These uniform constants are proved here on the full relative-gap
domain (9); the s=4 and s=9/2 constants are not assumed.

## 4. Full-target laws and a five-event retained set

At qb, R=r. On each FULL target face n-1 coordinates are independent
Uniform[-r,r] and the remaining coordinate is its true deterministic
pin. Averaging the capture indicator over the full target interval
gives exactly ell_R(z), so Fubini yields E_W J=1. Using (11) and
the elementary log bounds for |J-1|<=1/2 gives

    -81/rho_b^2<=E_W log J<=0,
    Var_W(log J)<=324/rho_b^2.

The target pin has |log J|<=18/rho_b. The full baseline log sum
therefore has absolute mean <=81n/rho_b^2+18/rho_b and variance
<=324n/rho_b^2. These vanish uniformly on the target-radius band,
since rho_b>=t/s and n/rho_b^2<=25x^2/n.

Fubini for the first profile gives exactly

    E_W J_1=E_fullsource phi=(2R/rho)H(a),
    a=1-g/(2R)>=a_s=(s+1)/(2s).

Since H decreases on [a_s,1] and R/q<=s,

    n E_W J_1 <=2s H(a_s)x=mu(s)x.                    (13)

The true pin contributes <=2/rho_b<=10x/n. Thus for large n the
sum has mean <=[mu(s)+d_mu/2]x. Its free-coordinate variance is
<=4n/rho_b^2<=100x^2/n. Chebyshev consequently gives

    P(sum_i J_1(w_i;qb)>[mu(s)+d_mu]x)
       <=400/(d_mu^2*n) ->0.                          (14)

Define A as A0 intersected with exactly these FIVE conditions:

1. |sum log J(w_i;qb)|<=e.
2. |Nsum(q_-)-E_-|<=e*n/4.
3. |Nsum(q_+)-E_+|<=e*n/4.
4. w belongs to the original D.
5. sum J_1(w_i;qb)<=[mu(s)+d_mu]x.

Here Nsum=sum_i[1+1_{|w_i|<a0(q)}] and
E_+/-=n+(n-1)a0(q_+/-)/r. Eventually a0<r, so the target pin
contributes exactly 1 and the other excess indicators are independent
Bernoulli(a0/r). Each count failure is <=4/(e^2*n). The baseline
failure tends to zero by its stated mean and variance. For the target
cap, negative pinned faces lie in D when t>rhi. On positive faces
the cap-complement free fraction is

    1+1/(2t)-1/(2r)
      <=1-1/(4rhi)<19/20

once t>=2rhi, since rhi<5. Its failure is at most (19/20)^(n-1).
All FIVE failures are therefore <=epsA/8 uniformly for sufficiently
large L. No independence between these five events is assumed. Thus

    lambda(A)>=(1-epsA)n log(rhi/rlo)>m*n.              (15)

Every probability above is a full-target ambient law used to choose A,
not a distribution assigned to actual increments.

## 5. Whole-matching degrees and genuine source-face retention

The matching q width is 2h/s. Since n/rho<=s*x, (11) transports
the baseline log product by at most 240xh<=240/L, which is <=e
eventually. Also |J_i-1|<=e eventually. Hence

    P0<=product_i J_i<=1/P0.

For ellpin=(1+a0/R)/2, ellpin'=-(1-epsilon)h/(2R^2)<0 and Nsum
is nondecreasing. At q_+ the mean ratio Nsum/(2n ellpin) is <=1;
at q_- it is >=1-2h-1/n, since a0<=1 and R,r>=1. Count
fluctuations alter this ratio by at most e/4. Once 2h+1/n<=3e/4,
the WHOLE matching interval satisfies

    1-e<=Nsum/(2n ellpin)<=1+e,
    U=[product_i J_i/(2ellpin)]sum_i N_i/J_i,
    Ulo*n<=U<=Uhi*n.                                  (16)

This sums all captured signed source faces and omitted-coordinate
terms. The first-profile transport from (12) is <=120xh; once
120h<=d_mu, it gives throughout matching

    sum_i J_1(w_i;q)<=[mu(s)+2d_mu]x.                 (17)

Fix ANY captured signed source face j before source-cap deletion and
soft clipping. Its pin loss is deterministic and <=log(5/3)<p=3/5.
The free-coordinate normalized inverse-capture joint density is the
product of 1_capture(z_i)/(2rho*J_i*ell_R(z_i)), i!=j. Thus the
coordinates are independent, with means J_1i/J_i and raw second
moments J_2i/J_i. Equations (5), (12), (17) yield

    E(V|j)<=p+m_face*x,
    Var(V|j)<=sum_(i!=j)J_2i/J_i
             <=s*x/(1-e)<10x,
    V=-log G.

The pin adds no variance. With a=chi-m_face>0,

    P_inv(V>chi*x|j)<=10x/(a*x-p)^2 ->0               (18)

once a*x>p. Therefore every face retains probability at least Ret
for all sufficiently large L, uniformly over retained targets and the
complete matching interval. Mixing faces only after this facewise
bound is valid requires no mixture-variance assertion. Since
bsoft/(theta/G)=min(1,G/theta) is one on V<=chi*x, the fullsource
shell row lies between Ret*theta*Ulo*n and theta*Uhi*n.

## 6. Weighted cap, complete column and the same global tau

From a0/R>=1/s and (16),

    #{i: |w_i|<a0}>=nu*n,
    nu=1/s-(1+1/s)e>0.

The source cap complement is exactly all z_i>=-rho+(1-q).
Negative source pins are excluded. Every counted FREE coordinate
captures the whole lower strip [-rho,-rho+(1-q)), with pre-cap
inverse probability >=(1-q)/[2rho(1+e)]. A positive source pin
omits at most one counted coordinate. Independence, then 1-y<=e^-y,
gives, on every captured source face,

    P_inv(D^c|j)<=exp[-((nu-1/n)*eta/(1+e))*x],
    eta=(1-qmax)/(2qmax)>0.                            (19)

Pointwise bsoft<=theta/G. Multiply (19) by the full face inverse
coefficients, sum them and integrate only over matched dq/q. The
radial relation

    (r/rhi)log((r+h)/(r-h))
       in [2(1-e)h/rhi,2(1+e)h/rhi]

holds uniformly for small h. Dividing by k cancels theta and bounds
the complete weighted source-cap loss by

    Lcap(w)/k<=M exp[-((nu-1/n)*eta/(1+e))*x] ->0.       (20)

For large n, nu-1/n>=nu/2, so this has a fixed positive exponential
rate in L and is eventually <=eS=e. It requires no chi<eta condition
and asserts independence only before cap and soft conditioning.

Before cap deletion raw rows lie between 2Ret(1-e)Ulo*k and M*k.
For a fixed matched source and target radius r>=R, its target capture
fractions ell_i(r)>1/2 are <=ell_i(R). If ell_i(r)<1, at most
one target pin sign is captured; otherwise at most two. The sum of ALL
target faces is <=n product_i ell_i(r)<=nG(z). Multiplying by arad
cancels dr/r, the ordinary target matching width is <=2h, and
bsoft*G<=theta. The complete column is <=2k<M*k; inactive columns
are zero.

The original scalar cap identity with Icap<=log t implies
B_n>=(n/2)log t>=nL/4>gamma*n for large L. The latter inequalities
follow from e^(L/2)>=x and fixed gamma, and do not require gamma<1.
Set the SAME global selection threshold

    tau=sqrt[A_n(B_n-gamma*n)]>0.

A_n<=exp x=n exp(-sigma*L) and B_n<=nL give

    tau/k<=rhi*L^(5/2)*exp(-cT*L) ->0,
    cT=sigma/2-chi(1-sigma)>0.                         (21)

For large L it is <=eT=e. The final row minus this SAME tau is
therefore >=b*k, both ambient degrees are <=M*k and target mass
is >m*n.

## 7. One common eventual threshold and original finite-history transfer

Every guard above has an explicit envelope tending to its required
limit: h and polynomial(L)/n vanish, t grows, target failures vanish,
the log product transport is <=240/L, profile transport divided by x
is <=120/L^2, (18) tends to zero and the cap/tau payments have fixed
positive exponential rates. The target-radius bounds are uniform on
[rlo,rhi]. These are finitely many strict inequalities with fixed
positive budgets, so their thresholds have a finite maximum. Positive
exponential-series terms prove all polynomial/exponential limits used
here. No uniform threshold over different c, gaps or primary tuples
is asserted or needed.

For EVERY original algorithmic a_alg>=8n^4 normalize by the last
minimum positive-contribution radius and reflect the last center into
the positive orthant. Let alpha=a_alg/(a_alg+1),
ccap=(a_alg*beta-1)/(a_alg+1). Since a_alg>2t, ccap>1. An actual
target item j whose increment meets A obeys

    ccap*r_j<=||u-w||<alpha+rhi<6,

so r_j<6. An actual source increment has r_i>=rho>20>6 and precedes
this SAME target. The displacement margin is

    alpha*beta-kappa-2/(a_alg+1)
       =zeta-(beta+2)/(a_alg+1)>0.

Indeed beta<2 and a_alg>=8n^4. Actual scaled center distance is
less than kappa*rho+(r_i+r_j)/(a_alg+1), at most
[kappa+2/(a_alg+1)]r_i<alpha*beta*r_i. Included sources really
capture the genuine target center in its own original eligibility test.
Both weights are <=1, so its actual weighted row is <=arad*tau<=tau.

On A and P the finite corner envelope is exactly lambda_a=c_a*lambda,
c_a=alpha^-n>=1, since

    alpha-ccap*m-(1-beta*m)=[(beta+1)m-1]/(a_alg+1)>0.

The true uncaptured-increment measure mu_a is dominated by lambda_a.
Energy, degree and target mass on the kernel support scale as
c_a^2,c_a,c_a. A and P are disjoint, so the two missing-mass terms
of the energy-difference identity are bounded by c_a*M*k times their
missing masses, whose sum is at most the TOTAL finite-cap deficit
Delta_a. Actual energy is <=tau*c_a*lambda(A). Writing
E=integral K d lambda d lambda gives

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/[c_a*M*k]
       >=[E-tau*lambda(A)]/[M*k]
       >=lambda(A)*b/M>=Gamma*n>gamma*n.               (22)

This holds for each original finite history uniformly in a_alg. No
actual-increment probability law, new selector or unpaid row factor
is introduced. gamma and Gamma may exceed one.

For general f in L1 use |f|. Points in a strict maximal level set have
actual positive witness radii; their open shrunk cubes finitely cover
each compact subset and selected dilates are still true witnesses.
Normalization applies to every positive terminal radius. Inner
regularity yields the full strict level set without compact-support,
scale-truncation, selector convergence or independent-history assumptions.

## 8. Entire original cover/root cost for the fixed tuple

Keep the captured budget, terminal self-weight and ENTIRE cover factor:

    C_n <=[1+(beta+1/a_alg)^n/tau]
       *[tau+a_alg^-n+alpha^-n J_a_alg(t)-gamma*n].

Only in these uniformly valid NUMERIC bounds let a_alg tend to infinity,
obtaining

    C_n<=(1+A_n/tau)(tau+B_n-gamma*n)
        =F(A_n,B_n-gamma*n).                           (23)

The exact scalar price is P=-log(1-sigma):

    B_n<=T_n+P*n+1,
    [A_n+2sqrt(A_n B_n)]/n
       <=exp(-sigma*L)+2sqrt(L)exp(-sigma*L/2).

The full positive-root remainder plus 1/n tends to zero. Eventually it
is <kappa_scalar=gamma-P-c, which was fixed positive in (6). This
proves F(A_n,B_n-gamma*n)<T_n-c*n with the root cross-term paid.
Thus the fixed-tuple assertion (4) is complete, subject to independent
review of the supplied proof. If gamma<=0, the weaker formula bound
follows directly from the original nonnegative-deficit bound
C_n<=F(A_n,B_n)<=F(A_n,B_n-gamma*n); the optimization below uses
only positive gamma.

## 9. Every c<C(s) uses a strict finite parameter choice

Fix s in [4,5], c<C(s), and d=C(s)-c>0. The following order makes
all choices concrete and preserves positive slack:

1. Choose chi>mu(s) so close to mu(s) that
   log(1+2chi)<log(1+2mu(s))+d/4.
2. Choose sigma>2chi/(1+2chi) so close to this value that
   P=-log(1-sigma)<log(1+2chi)+d/4.
3. Choose qmax<1 so close to 1 that qmax>1/s and log qmax>-d/8.
4. Choose positive gamma<log(s*qmax) so close to that upper endpoint
   that gamma>log(s*qmax)-d/8.

These choices exist by continuity and log(s*qmax)>0. They give

    gamma>log s-d/4,
    P<log(1+2mu(s))+d/2,
    gamma-P-c>d/4>0.                                  (24)

Now choose trim, mass m, Ret, d_mu and e in the order of Section 1,
then enlarge L to the common eventual threshold. The scalar root
remainder is paid only after the strict allowance (24) has been fixed.
No boundary parameter is inserted into a finite construction and no
parameter is allowed to drift with n. This proves (1) for every
strictly smaller c.

Conversely the sufficient primary constraints (4) have

    gamma-P<log s-log(1+2chi)<C(s),

so C(s) is the supremum of coefficients approached by this particular
strict-parameter recipe. This is a statement about the recipe and its
estimates, not an impossibility theorem for improved kernels, tail
estimates or other covering geometries.

## 10. Analytic maximization on [4,5]

Set l(s)=log((s+1)/(2s)). Direct differentiation gives

    mu'(s)=1+l(s)-1/s,
    C'(s)=[1+2mu(s)-2s mu'(s)]/[s(1+2mu(s))]
          =[1+2l(s)]/[s(1+2mu(s))].                  (25)

The denominator is positive. l(s) is strictly decreasing, since
l'(s)=-1/[s(s+1)]<0. Thus the derivative changes sign exactly where
l(s)=-1/2, namely

    s_*=1/(2a-1), a=exp(-1/2).

It lies strictly in [4,5]. One elementary proof uses

    exp(1/2)>13/8>8/5,
    exp(1/2)<=13/8+(1/48)/(1-1/8)=277/168<5/3.

The second line bounds the tail beginning with the cubic term by its
geometric ratio upper bound 1/8. Consequently 3/5<a<5/8, which is
equivalent to 4<s_*<5. The numerator in (25) is positive before
s_* and negative after it. Hence this is the unique maximum on [4,5].

At s_* one has l=-1/2 and mu=(s_*-3)/2. Therefore

    C(s_*)=log[s_*/(s_*-2)]
       =-log(1-2/s_*)=-log(3-4 exp(-1/2))=c_*.

Its log argument lies in (0,1), as also follows from 3/5<a<5/8.
This proves (2) analytically. The accompanying bounded exact-interval
check encloses s_* in (4.69348449,4.69348450) and c_* in
(0.55533956,0.55533957); the displayed extra digits are numerical
orientation, not a separate proof of optimality.

## 11. The quantifier step and exact scope of the limit statement

For any epsilon>0 fix s=s_* and c=c_*-epsilon. Section 9 supplies
ONE fixed strict primary tuple, followed by fixed trim/Ret/errors and
a finite threshold depending on that epsilon. Sections 1--8 then give
the first assertion of (3). By the definition of limsup this yields

    limsup_(n->infinity)(C_n-T_n)/n<=-c_*.

There is no exchange of the parameter limit with a finite-history or
operator limit. Each epsilon has its own proven eventual theorem.
No common threshold for all epsilon, or finite-domain inequality
C_n<T_n-c_*n, follows from this argument. The analytic maximum concerns
this fixed-ratio family with s in [4,5] and its specified sufficient
estimates; it is not a global optimum for C_n or all covering methods.

## 12. Actual local calculation and handoff status

The proof uses no parameter scan. A bounded local Python calculation
with Fraction arithmetic checks the newly derived uniform derivative
constants, baseline moments, count/variance bounds, and an alternating-
series enclosure of exp(-1/2). A 20-term positive atanh-log enclosure
with exact tails certifies the displayed s_*,c_* intervals. An initial
bounded Decimal evaluation supplied only orientation; SymPy was not
available in that Python, so the derivative identity (25) is proved by
the explicit algebra in the text rather than attributed to a tool.

REGISTRATION.json fixes this one bounded check before execution.
verify.py and verification.json preserve the actual exact-interval run.
The final script was actually run once and returned PASS on 18 exact
constant/interval checks; it did not evaluate a dimension threshold.
The script does not establish a dimension threshold or replace the
uniform local, probability, source-face, weighted-cap, same-history
and full-cost proofs supplied above. Independent whole-chain review
and main-agent acceptance remain pending.

Only upper/round085/parameter_family in this independent clone was
written. No earlier accepted result, paper, original source, shared
MATLAB session, n=7/lower-bound task, commit or push was changed.
