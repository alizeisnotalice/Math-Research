# R079 threshold transfer — independent ACCEPT

This is independent temporary staging, **not a Git clone**. The original
checkout was read only. Nothing was committed or pushed. All new files
are inside this `audit/round079/threshold` directory. No shared MATLAB,
production search, old 53-gate baseline, or producer program was run.

The accepted statement is, for every integer n with L=log n>=4096,

    C_n <= F(A_n(t), B_n(t)-(7166/25000)n)
        < n log n-n log log n-(17/100)n,
    t=n/[(89/100)L], A_n(t)=(1+1/t)^n,
    F(A,B)=(sqrt A+sqrt B)^2.

Here B_n is the original limiting corner-envelope mass from the R073
finite-cover framework. This is the original centered axis-parallel cube
maximal operator on Lebesgue measure, for arbitrary L1 inputs, all strict
levels and all positive witness radii. The argument below transfers the
accepted R073 local and finite-cover framework at the NEW parameters;
it is not an inference from a passed gate list. It uses only the old
UNWEIGHTED source-cap deletion. No weighted-cap theory is included.

## 1. Inputs, stale fields, and actual new calculation

Sources read: the complete R073 and R078 proofs, the R073 independent
full-chain acceptance, and all eight requested R079 certificate,
registration, amendment and stronger-guard files. The source hashes are
recorded in `independent_result.json`; the stronger receipt's three hash
links and the amendment's parent-registration hash agree with the
unchanged files. No original hash or receipt was rewritten.

The authoritative geometric tuple is `chosen_parameters`, namely

    s=3/2, qmin=2/3, qmax=89/100, sigma=11/100,
    chi=51/1000, h=1/(1536L), delta=2^-19,
    eB=eN=eT=2^-19, eJ=eRad=eS=eD=2^-20,
    epsA=2^-16, x=(1-sigma)L.

The amendment sets target mass lower coefficient mL=2889/10000 and
Ret=999999/1000000, and supplies the Bernstein replacement. Its gamma,
c and price are superseded for this audit by `stronger_registration`:

    gamma=7166/25000, price=5827/50000, c=17/100,
    gamma-price-c=1/10000.

The registration's `original_model_parameters`, `planned_gamma`,
`planned_final_c`, and Cantelli `retention` prose are copy remnants.
They do not define this theorem. In particular no old threshold,
old sigma, old qmax or Cantelli constant is imported. The amendment
and stronger choice are sequential refinements of the same geometry.

`independent_check.py` was actually run with local Python 3 standard
library Fraction arithmetic, PID 27558. It constructs its fixed inputs
independently, uses N32(L)=L^32/32!<=exp L, 24 positive atanh terms plus
their exact geometric logarithm remainder, and P48(y)=sum_0^48 y^j/j!.
It imports or executes no producer helper. Nine grouped arithmetic
checks pass. This is bounded new arithmetic, not a parameter search.

Useful displays (all decision inequalities were rational):

    x0=3645.44, h0=1/6291456, N32(L0)>1.49*10^80,
    (1-epsA)log(rhi/rlo)>0.2889235470880792>mL,
    H0>14.43057787285, P48(H0)>1849781>1000000,
    (1-eD)Gamma-gamma>0.0000031206583167,
    price+log(1-sigma)>0.0000061837440484,
    cap/k<2.114*10^-11<eS,
    tau/k<4.677*10^-9<eT,
    complete scalar remainder<1.470*10^-50<1/10000.

These displays support the fixed base. The following analytic argument
is what extends the conclusion to every L>=4096 and to the operator.

## 2. New active domain and unchanged local lemma interfaces

Define beta=1+1/t, zeta=n^-3, epsilon=zeta*t=1/(x n^2),
kappa=beta-zeta, rlo=1+delta and rhi=s*qmax-delta. The corner envelope
is d lambda=2^-n m(z)^-n 1_D dz, where m=max(1,||z||_infinity),

    D={||z||_infinity<t: some z_i<1-beta*m(z)}.

Let A0 be the full target annulus [rlo,rhi], and Pfull the complete
source band q=||z||/t in [qmin,qmax]. The source set is P=D intersect
Pfull. Matching means |r-sq|<h. For the complete interval set

    q_-=(r-h)/s, q_+=(r+h)/s, qb=q_+,
    rho=tq, R=sq-h, a0=(1-epsilon)q, g=R-a0.

At the new base, 2h<delta, h<=2^-20, epsilon*qmax<=h/s,
N32*qmin/x>4rhi, N32/x>rhi, x/N32<1 and x*N32^2>1.
Consequently [q_-,q_+] lies inside the source band: s*qmin=1,
rlo-h>1 and rhi+h<s*qmax. All active R satisfy

    r-2h<=R<=r, R>rlo-2h>1,
    g-R/3=epsilon*q-2h/3<=0, 0<g<=R/3,
    rho>4rhi, a0>0, beta<2, kappa>1, 2<=x<=L.

Both annuli lie inside ||z||<t. Targets have norm<3/2 and sources
norm>4, so A0 and P are disjoint. Inactive sources may have R<1;
their final columns are zero, and no R>=1 lemma is applied to them.
All the relative-gap R064/R066/R070 local interfaces accepted in R073
therefore apply on their genuine domains at this new h. No old finite
h-coefficient gate is substituted.

On the full source cube define

    ell_R(z_i)=min(1,(R+kappa*rho-|z_i|)/(2R)),
    phi_i=-log ell_R(z_i), G=product ell_R(z_i),
    theta=exp(-chi*x), bsoft=min(1,theta/G),
    arad(w)=||w||/rhi.

The fields are Borel; 0<ell_R<=1 and 0<bsoft<=1. The final kernel is

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
      *1_{||z-w||<kappa*rho}*arad(w)*bsoft(z).

Here A is specified below; bsoft is defined before source-cap removal
on ALL Pfull, including unmatched edge sources. The two bounded weights
are energy devices applied after the original selections.

## 3. Target laws, whole matching interval, and true source variance

For fixed target coordinate w_i define J,J1,J2 by the full-source
captured integrals of 1/ell, phi/ell, phi^2/ell, normalized by 2rho.
J2 is a second integrand profile, not J1 squared or E(V^2).
The R073 local identities and AC derivative bounds give

    |J-1|<1/(3rho), |d_q log J|<3/rho,
    0<=J1<3/(50rho), |d_q J1|<1/(2rho),
    0<=J2<1/(125rho), |d_q J2|<1/(10rho).

Positive-part kinks are harmless because the profiles are absolutely
continuous, and these derivative bounds hold almost everywhere.
At qb, R=r, and full-target Fubini yields the baseline mean-one J law,
the R073 log mean/variance bounds, and the full-source phi/phi^2 means.
Writing Q=log(6/5), the limiting profile constants are

    mu1=1/2-(5/2)Q<221/5000,
    mu2=1-5Q-(5/2)Q^2.

The finite-h structure, not merely g/R<=1/3, is used:

    0<=a-5/6<=h/(3R), a=1-g/(2R), R/q=s-h/q.

Since K2'(a)=-(log a)^2, this gives the baseline second-profile mean

    (1-1/n)[mu2-h(mu2+1/25)] <= E(S2)/x
        <= mu2+1/(250nqmin).

The n-1 free target coordinates and the deterministic target pin are
kept separate. The new rational enclosures place this mean strictly
inside (0.00528,0.005295)x at L0, and hence throughout the tail below.

On each full signed target face impose the same six conditions as R073:
the baseline |sum log J|<=eB, both endpoint count conditions with
tolerance eN*n/4, membership in D, sum J1<=(221/5000+1/100000)x,
and (527/100000)x<=sum J2<=(53/10000)x. All six probabilities are
calculated before conditioning on A or D. Free coordinates on each
such face are independent; the six events need not be independent.

For the baseline log event use mean bound
ML=(3/4)L^2/(qmin^2 N32)+(1/2)L/(qmin N32) and variance bound
VL=(3/2)L^2/(qmin^2 N32). New checks give ML<=eB/2 and
4VL/eB^2<=epsA/8. Each endpoint Bernoulli count fails with probability
at most 4/(eN^2 N32)<=epsA/8. The complement of D is absent on negative
pinned faces and on positive faces has free fraction at most
1+x/(2N32)-1/(2rhi)<=3/4. Thus its probability is at most
(3/4)^(n-1)<=(3/4)^74<=epsA/8.

For J1, the pin is <=half of the 10^-5*x mean slack at n>=9000.
Chebyshev on the remaining gap gives failure <=324000000/N32<=epsA/8.
For J2, the mean is at least x/200000 from both excluded endpoints.
Its free-coordinate variance is <=n/(125^2 rho_b^2), so a SINGLE
two-sided Chebyshev event costs <=5760000/N32<=epsA/8. No extra factor
two is missing. The union costs <=6epsA/8<epsA, and therefore

    lambda(A)>=(1-epsA)n log(rhi/rlo)>mL*n.

For the whole matching interval, ellpin=(1+a0/R)/2 has derivative
-(1-epsilon)h/(2R^2)<0, whereas the count Nsum is nondecreasing.
Endpoint calibration and 2h+1/N32<=3eN/4 give
1-eN<=Nsum/(2n ellpin)<=1+eN for every matching q. The log-product
shift is <=6Lh=6/1536; with the baseline event it gives

    P0<=product J_i<=1/P0, P0=1-eB-6/1536.

The pointwise bound x/(3qmin N32)<=eJ separately controls every
J_i in [1-eJ,1+eJ]. Summing ALL captured source signs and omitted
coordinates yields Ulo*n<=U<=Uhi*n, where

    Ulo=P0(1-eN)/(1+eJ), Uhi=(1+eN)/[P0(1-eJ)].

AC integration transports sum J1 by <=xh and sum J2 by <=xh/5.
The new h satisfies h<=10^-5 and h/5<10^-6, so uniformly on matching q

    sum J1<=2211*x/50000,
    5269*x/1000000<=sum J2<=5301*x/1000000.

Fix any captured signed SOURCE face before cap deletion or clipping.
The pin loss is deterministic and <=Q<183/1000. The inverse-capture
law is a PRODUCT on its free coordinates, because capture and inverse
weights factor and J_i are their normalizers. Their exact first and
raw second moments are J1_i/J_i and J2_i/J_i. Thus V=-log G obeys

    E V<=p+m*x, Var V<=v*x,
    p=183/1000, m=(2211/50000)/(1-eJ), v=531/100000.

For the variance, subtract squared means, then add the omitted
nonnegative J2 term and divide by 1-eJ; the new bound
(5301/1000000)/(1-eJ)<v pays it. No pin variance or mixture variance
is inserted. These are auxiliary ambient laws, never laws for mu_a.

## 4. New finite Bernstein retention and its whole-tail proof

Since kappa*rho=rho+a0 and |z_i|<=rho, every source coordinate has
ell_R>=1-g/(2R)>=5/6. Hence 0<=phi_i<=Q<p on every active full cube.
For each free coordinate Y_i=phi_i-E phi_i is centered and |Y_i|<=p.
Independence and k!>=2*3^(k-2), k>=2, give for 0<=lambda<3/p

    log E exp(lambda*sum Y_i)
        <=lambda^2*v*x/[2(1-p*lambda/3)].

This follows directly from the bounded moment series and log(1+y)<=y;
no unverified literature inequality is needed. For u=(chi-m)x-p>0,
the admissible choice lambda=u/(vx+pu/3)<3/p gives

    P_inv(V>chi*x | fixed face)<=exp[-H(x)],
    H(x)=u^2/[2(vx+(p/3)u)].

At x0=(89/100)*4096, the new exact P48(H0)>1000000, so the failure
is strictly less than 10^-6. Put a=chi-m>0 and d=p/3. Directly,

    H'(x)=u[a*v*x+p*v+d*a*u]/[2(v*x+d*u)^2]>0.

Thus u>0 and Ret=999999/1000000 hold for every x>=x0, uniformly in
target, matching q, omitted face and sign. Only after proving this
uniform estimate do we mix the source-face probabilities. The
bsoft/(theta/G) ratio is min(1,G/theta), equals one on V<=chi*x and
is nonnegative elsewhere. Therefore the full soft source-shell row
coefficient lies between Ret*theta*Ulo*n and theta*Uhi*n. This step
precedes source-cap removal. It requires no independence after cap
conditioning or soft clipping.

## 5. Complete degrees, UNWEIGHTED cap payment, tau, and all L>=4096

Set k=n*h*theta/rhi and M=2(1+eRad)Uhi. Exact radial integration gives
arad*log((r+h)/(r-h)); atanh(h/r)/(h/r) is between 1 and
1/[1-(h/r)^2]<=1+eRad. Raw rows are therefore at least
2Ret(1-eRad)Ulo*k and at most M*k.

For the complete column, on each matching target shell r>=R, all target
signs obey T_i/(2ell_i(r))<=1, and ell_i(r)>1/2. Hence the full
target-face coefficient is <=n*product ell_i(r)<=n*G(z). The arad
factor cancels dr/r to dr/rhi. The matched ordinary width is <=2h,
and bsoft*G<=theta, giving column <=2k. No target sign is omitted.
Inactive columns vanish. Target restriction only decreases degrees.

For source-cap removal define

    eta=(1-qmax)/(2qmax)=11/178,
    Clog=log(qmax/qmin),
    cS=(eta-chi)(1-sigma)=961/100000,
    cT=sigma/2-chi(1-sigma)=961/100000.

The equality of the NEW rates is calculated, not imported. The
unweighted cap-complement mass is <=Clog*n*exp[-(1-1/n)eta*x].
Both weights and the capture indicator are <=1, so this total
unweighted mass bounds its removal from each row. Since
eta*x/n<=eta*x/N32<=1/1025,
exp(eta*x/n)<=1/(1-eta*x/n)<=1025/1024. Therefore

    caploss/k <=(1025/1024)*Clog*rhi*1536*L*exp(-cS*L).

It is <eS by the new P48 base bound and monotonicity below.
This is exactly the old unweighted-cap chain, with its complete
logarithmic band mass, rhi, h denominator and positive finite-n factor.

The original cap identity and Icap<=log t imply B_n>=n log(t)/2.
At L>=4096, exp(L/2)>=(L/2)^2/2>=L>=x, so log t>=L/2 and
B_n>=nL/4>gamma*n. Also B_n<=nL and A_n<=exp x. Set the SAME
numeric threshold in the original finite algorithm to

    tau=sqrt[A_n*(B_n-gamma*n)]>0.

It obeys tau<=n*sqrt(L)*exp(-sigma*L/2), and hence

    tau/k<=rhi*1536*L^(3/2)*exp(-cT*L)<eT.

The new base bound uses sqrt(4096)=64 exactly; it does not bound
sqrt L by 64 throughout the tail. Put

    b=2Ret(1-eRad)Ulo-eS-eT, Gamma=mL*b/M.

Final row minus this same tau is >=b*k; final rows and columns are
<=M*k. New rational arithmetic gives b>0, M>2 and
0<gamma<(1-eD)Gamma<Gamma<1.

Here is the entire-tail extension, including the cases hidden by no
finite check. With N32=L^32/32!, x/N32 is proportional to L^-31,
epsilon*qmax/h is bounded by a multiple of L^-64, baseline log mean
and variance by sums of L^-30,L^-31, and every target count/profile
bad probability by L^-32. Source lower radii grow as L^31. Count,
individual-J, radial, pin and cap finite-n errors decrease. The J2 mean
lower envelope is an increasing product of positive increasing factors
(1-1/N32) and [mu2_lower-h(mu2_upper+1/25)]; its upper envelope
decreases. Constant product loss 6Lh does not vary. All fixed rational
geometry, target mass and degree inequalities thus hold on the tail.
Bernstein has its explicit monotonicity proof above.

For any positive prefix P of degree>=2,

    yP'(y)-P(y)>=-1+y^2/2>=0 for y>=2,
    yP'(y)-(3/2)P(y)>=-3/2-y/2+y^2/4>=0 for y>=4.

All omitted coefficients are nonnegative. Thus L/P48(cS L) and
L^(3/2)/P48(cT L) decrease once cS L,cT L>=4, which is true already
at L0 (both equal 39.36256). Using exp y>=P48(y) proves both complete
payments on every L>=4096. No sampled dimension-tail inference enters.

## 6. Every finite dilation and the same genuine selection history

Fix ANY original finite true-witness family and run the original
largest-eligible-radius algorithm with the numeric threshold tau above.
Take ANY algorithmic dilation a>=8n^4. Normalize a positive load-support
point by its last minimum positive-contribution radius and reflect its
last center into the positive orthant. Let alpha=a/(a+1),
c_cap=(a*beta-1)/(a+1), and u be the last normalized center.

The genuine uncaptured-increment measure mu_a is dominated by the full
finite envelope lambda_a, whose mass is alpha^-n J_a(t). On A and P
the envelope is exactly alpha^-n lambda: D lies inside the finite cap,
as for m>=1

    alpha-c_cap*m-(1-beta*m)=[(beta+1)m-1]/(a+1)>0.

Because x>1, t<n, so a>2t and c_cap>1. An actual target item meeting A
has c_cap*r_j<=||u-w||<alpha+rhi<3 and therefore r_j<3. An actual
source increment has r_i>=rho>4. It precedes THIS target item in the
SAME largest-radius history. No separately selected source history is
spliced in. Moreover beta<2 and a>=8n^4 imply

    alpha*beta-kappa-2/(a+1)
      =zeta-(beta+2)/(a+1)>0.

Indeed (beta+2)/(a+1)<4/(8n^4)<n^-3 for n>=1. An included pair's
true centers differ by less than

    kappa*rho+(r_i+r_j)/(a+1)
      <[kappa+2/(a+1)]r_i<alpha*beta*r_i.

Thus every included source increment genuinely captures that target
center in its own original eligibility test. Only included pieces of
the actual disjoint increments are integrated, costing at most the
original item weights. Since both new kernel weights are <=1, the
actual row is <=arad(w)*tau<=tau on mu_a-almost every target. No
auxiliary target/source probability law is imposed on mu_a.

Write c_a=alpha^-n>=1 and E=integral K d lambda d lambda. Ambient
energy scales by c_a^2, each degree by c_a and target mass by c_a.
Let nu_a=lambda_a-mu_a>=0 and Delta_a=nu_a(full space). The energy
difference identity, using row/column support A times P and the
ambient degree bound, is at most

    c_a*M*k*[nu_a(A)+nu_a(P)]<=c_a*M*k*Delta_a,

because A and P are disjoint. The actual energy is at most
tau*mu_a(A)<=tau*c_a*lambda(A). Hence

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
       >=[E-tau*lambda(A)]/(M*k)
       >=lambda(A)*b/M>=Gamma*n>gamma*n.

This proves the genuine full finite-cap deficit for EVERY original
finite history and EVERY a>=8n^4. The factor c_a is not dropped in
the wrong direction. No convergence of selectors or increment
measures is required.

## 7. General L1, every positive scale, full cover factor and root cost

For arbitrary f in L1 use |f| and any strict maximal level. Each point
has an actual positive witness radius. Its open a-shrunk witness cube
covers that point; compact subsets have finite subcovers. The selected
a-dilates remain those true witnesses. The normalization in Section 6
works for any positive terminal radius, so these finite families can
have arbitrary positive radii. Lebesgue inner regularity then gives
the whole strict level set. This requires neither compact support nor
any lower/upper scale cutoff. All levels are covered by homogeneity.

Insert the genuine deficit into the accepted original finite
cover/load/witness formula, retaining the ENTIRE cover factor and
terminal self-weight:

    C_n <=[1+(beta+1/a)^n/tau]
       *[tau+a^-n+alpha^-n J_a(t)-gamma*n].

For each fixed n, let a tend to infinity only in these numeric bounds
already uniform over every finite witness family. Its limit is

    (1+A_n/tau)*(tau+B_n-gamma*n)
       =F(A_n,B_n-gamma*n).

No common selected history across a, inputs or unrelated families is
needed. Let T=n log n-n log log n. The original scalar estimates give

    B_n<=T-n log(1-sigma)+1,
    A_n+2sqrt(A_n B_n)
       <=n exp(-sigma L)+2n sqrt L exp(-sigma L/2).

The new rational log enclosure gives -log(1-sigma)<5827/50000.
Consequently the full F value is less than T-c*n provided

    1/n+exp(-sigma L)+2sqrt L exp(-sigma L/2)<1/10000.

At L0 the true bracket is bounded by

    1/N32(L0)+1/P48(sigma L0)+128/P48(sigma L0/2)
       <1.470*10^-50<1/10000.

The N32 term decreases. The positive-prefix inequalities above also
prove decrease of L^0/P48(sigma L) and
L^(1/2)/P48(sigma L/2), since sigma L0/2>4. Thus this is the entire
tail remainder, not just a base or asymptotic check. Finally
gamma-price-c=1/10000 exactly, proving the strict stated bound on
every L>=4096 with all A, cross and self terms paid.

## 8. Audit conclusion and scope

ACCEPT for the displayed R079 threshold branch. The finite arithmetic
supports its margins; Sections 2--7 provide the new-parameter all-tail,
same-history, finite-a and general L1/all-positive-radius transfer.
This does not claim a global optimum or dominance over all existing
branches; any accepted envelope still takes the minimum of branches on
their respective valid domains. The R073/R078 historical source files
and receipts remain untouched. The original 53 checks were not rerun,
and their PASS alone was never treated as an operator proof.
