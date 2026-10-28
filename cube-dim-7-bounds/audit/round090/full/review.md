# R090 arbitrary fixed s>=4 and limsup — independent ACCEPT

Frozen source reviewed in full: `/Users/zhengzhihao/Desktop/T/audit/round090/general_s/PROOF.md`,
SHA256 `12605f79fc17c83516fd4423c69fd58b1340a2a780588490591795ee3ec59631`.

ACCEPT: for EVERY fixed finite real s>=4 and every real

    c<C_free(s)=log[s/(1+mu(s))],
    mu(s)=(s-1)+(s+1)log((s+1)/(2s)),

one fixed strict tuple and a finite L_*(s,c) give the original
Lebesgue centered axis-parallel cube maximal-operator bound

    C_n<=H_tau=(1+A_n(t)/tau)(tau+B_n(t)-gamma*n)
       <T_n-c*n,
    tau=n exp(-dL), t=n/[(1-sigma)L], L=log n,
    T_n=n log n-n log log n.

It follows that

    limsup_(n->infinity)(C_n-T_n)/n<=log(1-log2).

The coefficient -log(1-log2) is reached only as a limiting consequence
of actual fixed-finite-s theorems. No finite-domain theorem with that
boundary coefficient, uniform threshold over s or epsilon, explicit
dimension threshold, optimized F bound or global covering optimum is
accepted or claimed.

This directory is independent temporary staging, not a clone. All
original sources, old rounds, paper and Git state were read only.
No shared MATLAB, commit or push was used. The author's s=10 four
regimes, quadrature, finite differences and old grids were NOT run.

## 1. Fresh local proof without the old g/R<=.8 domain

Fix one FINITE real s>=4. Hold s,t,h,epsilon,w constant during q
differentiation. The actual eventual matching domain supplies

    q>=1/s, h<=.001, rho>4(s+1), 1<=R=sq-h<=s,
    a0=(1-epsilon)q>0, g=R-a0 in [0,R],
    epsilon*q<=h/s, |w|<=r<=min(s,R+2h).

The stronger geometric relation g/R<=1-1/s is retained for the
profile mean and count. It may exceed .8 when s>5, and is never
replaced by the old .8 bound. For the local value and derivative
estimates, the larger domain g/R<=1 suffices.

On source [-rho,rho] the two boundary layers have width g<=s<rho.
There ell=1-u/(2R), 0<=u<=g; in the central interval ell=1,
phi=0. Edge integration gives g+F, RQ^2 and (2R/3)Q^3 for
1/ell, phi/ell and phi^2/ell. By reflection take w>=0. If w>a0,
capture deletes the left interval of length m_w=w-a0<=s. Since
2rho-g>7s+8>s, deletion cannot reach the opposite edge. When
m_w>g it removes the whole near layer plus a central interval,
whose phi profiles vanish but whose inverse mass is its length.
Thus, with v=(g-m_w)_+,

    J=1+[F(R,g)+F(R,v)-m_w]/(2rho),
    J1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).

These are exact throughout matching, including |w|>R. The -m_w
term is retained. Null coordinate endpoints do not change integrals.

For 0<=u<=R, set z=u/(2R)<=1/2. Direct calculus gives

    Q=-log(1-z)<=log2<1,
    |Q_R|<=1/R, Q_u<=1/R,
    0<=F=2RQ-u<=2R,
    F_R=2Q-2u/(2R-u) in [-2,0],
    F_u=u/(2R-u) in [0,1].

F_R<=0 follows by integrating 1/(1-v)<=1/(1-z), and its lower
bound follows from Q>=0. No s<=5 estimate is used.

For A=F(R,g)+F(R,v)-m_w, -s<=A<=4s. Almost everywhere,
R'=s, 0<=g'<=s, |m_w'|<=1 and v'=0,g',or g'-m_w', hence
0<=v'<=s. Therefore

    |A'|<=6s+1,
    |J-1|<=2s/rho,
    |J'|<=[4s^2+6s+1]/(2rho).

The local radius guard gives 2s/rho<s/[2(s+1)]<1/2, so J>1/2.
The margin over 1/2 depends on s; no fixed endpoint margin such as
1/10 is asserted for arbitrary s. Consequently

    |d_q log J|<=[4s^2+6s+1]/rho=D_s/rho.

Also J1<=s/rho, J2<=2s/(3rho)<=s/rho. For u=g or v,
|d_q Q(R,u)|<=2s/R. Using the EXACT prefactor derivative
(R/rho)'=h/(q*rho) and differentiating both squared terms yields

    |J1'|<=h/(q*rho)+4s/rho<=5s/rho.

The positive parts m_w,v are Lipschitz, while Q,F are smooth in
a neighborhood of the closed relative domain with R>=1 and
1-u/(2R)>=1/2. Thus J,J1,J2 are AC; J>1/2 also makes log J
AC. All derivative estimates integrate across onsets. The local
constants depend only on the already fixed s, rather than on n or
an old bounded-s grid.

## 2. Five target events and true source-face variance

At baseline R=r, full-target Fubini still cancels ell exactly and
gives E_W J=1. With |J-1|<=2s/rho<=1/2, elementary log bounds
give free-coordinate log mean in [-4s^2/rho^2,0] and variance
<=16s^2/rho^2. The deterministic target pin contributes at most
4s/rho in absolute value. Hence the full-face log sum has envelopes

    absolute mean<=4s^2*n/rho_b^2+4s/rho_b,
    variance<=16s^2*n/rho_b^2.

Using rho_b>=n/(s*x), these are bounded by positive fixed-s
multiples of x^2/n and x/n, uniformly in target radius and sign.

The stronger ACTUAL geometry supplies a=1-g/(2R)>=a_s=(s+1)/(2s).
Full-target Fubini for J1 gives (2R/rho)H(a). Since H decreases
and R/q<=s, the free-sum mean is <=mu(s)x. Its pin contributes
<=s/rho_b<=s^2*x/n, eventually <=d_mu*x/2. Its free variance
is <=n*s^2/rho_b^2<=s^4*x^2/n. Chebyshev therefore costs at most
4s^4/(d_mu^2*n) for the threshold (mu+d_mu)x.

At both count endpoints, a0<qmax<1<rlo implies that the actual
target pin contributes exactly one to Nsum. Free excess indicators
are independent Bernoulli(a0/r), yielding failure <=4/(e^2*n)
for each count condition. Negative target pins lie in D when t>rhi.
On positive faces, t>=2rhi gives cap-complement fraction
<=1-1/(4rhi)<1-1/(4s)<1. For THIS fixed s its (n-1)st power
tends to zero. There is no false .95 bound independent of s.

These are exactly five events: one log, two counts, D and J1.
Each is eventually paid by e/8 uniformly over the fixed band; no
independence between events is needed. Their union gives
lambda(A)>=(1-e)n log(rhi/rlo)>m*n. All laws are full-target
ambient laws before conditioning, not laws for actual increments.

The q interval width is 2h/s and n/rho<=s*x. Hence product-log
transport is <=2D_s*x*h=O_s(1/L), and first-profile transport is
<=10s*x*h. The latter is <=d_mu*x once 10s*h<=d_mu.
Individual J errors vanish. The endpoint mean/count argument, with
ellpin'=-(1-epsilon)h/(2R^2)<0, gives the whole-interval count
ratio in [1-e,1+e] once 2h+1/n<=3e/4. Thus the complete signed
face sum has Ulo*n<=U<=Uhi*n and sum J1<=(mu+2d_mu)x.

Each pre-cap SOURCE face has the genuine product inverse law with
free means J1_i/J_i and raw second moments J2_i/J_i. Its pin is
deterministic and <=log2<1. The TRUE centered variance satisfies

    Var V=sum_(i!=j)[J2_i/J_i-(J1_i/J_i)^2]
       <=n*s/[rho(1-e)]<=s^2*x/(1-e)<2s^2*x.

The last inequality follows from e<1/12; no pin or mixture variance
is added. The mean is <1+m_face*x. For fixed a_face>0, Chebyshev
gives 2s^2*x/(a_face*x-1)^2, which tends to zero uniformly over
targets, matching q, faces and signs. Only after this control may
faces be mixed. min(1,G/theta)=1 on V<=chi*x then yields the
uniform raw Ret lower row. No old s<=5 variance constant, second-
profile event, normal approximation or actual-increment law is used.

## 3. The s-dependent count and complete weighted degrees

The chosen error is FIXED after s, with e<1/[2(s+1)]. Hence

    nu=1/s-(1+1/s)e>1/(2s)>0.

This is essential: a fixed error independent of arbitrarily large s
would not provide the same conclusion. The whole-matching count
therefore supplies at least nu*n counted coordinates.

Source D-complement is exactly all z_i>=-rho+(1-q). Each counted
FREE coordinate captures its full lower strip; the lower endpoint
uses w_i<a0 and the other endpoint uses
2rho+w_i+a0>2rho-s>7s+8>1. Negative source pins are excluded and
an arbitrary positive pin omits at most ONE counted coordinate.
Pre-cap inverse independence gives

    P_inv(D^c | face)<=exp[-(nu-1/n)*eta*x/(1+e)].

For fixed nu,eta>0, nu-1/n>=nu/2 eventually. Multiplying
bsoft<=theta/G, summing every face coefficient and integrating
ONLY over actual matched dq/q gives caploss/k<=M times that
exponential. theta cancels exactly. The complete cap payment tends
to zero without an unweighted mass or chi<eta restriction.

For columns, rho>4(s+1)>r ensures the one-boundary target capture
formula is valid even for large s. Its fraction is >1/2 and decreases
as r>=R grows, hence is <=ell_i(R). All target pin signs satisfy
T_i/(2ell_i(r))<=1; summing complete faces gives <=nG. arad
cancels dr/r, ordinary target width is <=2h, and bsoft*G<=theta.
The ENTIRE column is <=2k<M*k. Inactive columns vanish. Exact
radial log-width gives the usual raw row bounds. Neither argument
has a hidden target-radius cutoff at 5 or 6.

## 4. Same tau, actual finite history and general L1/all scales

The original R027 finite theorem and v3 cover audit permit EVERY
positive tau. For each fixed n, fix tau=n exp(-dL) before inputs,
witnesses, histories and finite a. The same value is used everywhere.
The new rate d-chi(1-sigma)>0 gives

    tau/k=rhi L^2 exp[-(d-chi(1-sigma))L] ->0.

After weighted cap and THIS tau are paid, row-tau>=b*k and both
ambient degrees are <=M*k. No optimizing square-root tau is used.

For every a_alg>=8n^4, t<n eventually ensures ccap>1. An actual
target item meeting A has radius<alpha+rhi<s+1. A source increment
has radius>=rho>4(s+1), so it precedes that SAME target in the
original largest-radius history. The finite-s dimensions replace all
old bounded-s constants explicitly. For actual n>=2 and beta<2,

    n^-3-(beta+2)/(a_alg+1)
      >n^-3-4/(8n^4)=(n-1/2)/n^4>0.

This and r_j<r_i make every included source capture the true target
center in its original eligibility test. Genuine increment pieces cost
at most their original item weights. Bounded kernel weights leave
the actual weighted row <=arad*tau<=tau.

The finite envelope equals c_a*lambda on A and P by the same
positive cap margin. Energy, degree and mass scale by c_a^2,c_a,c_a.
Disjoint A,P allow both missing-mass terms to be paid from the TOTAL
finite-cap deficit, giving

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
       >=[E-tau*lambda(A)]/(M*k)
       >=lambda(A)*b/M>Gamma*n>gamma*n.

The factors and c_a>=1 direction are correct; no gamma<1 or
deficit<=n restriction is needed. Bounds hold for every original
finite history and sufficiently large fixed n.

For arbitrary L1 input use |f|. True positive witnesses produce open
shrunk covers; compact subsets have finite subcovers, and selected
dilates remain true witnesses. Common normalization admits every
positive terminal radius and leaves weights and tau unchanged.
Inner regularity and level homogeneity give the entire strict level
set without compact support, radius truncation, scale lattice,
actual-increment independence or convergence of selected histories.

There are finitely many eventual envelopes for each FIXED tuple:
polynomial(L)/n, the fixed-s geometric target-cap base, O_s(1/L)
and O_s(1/L^2) transports, the face Chebyshev bound and positive-rate
cap/tau exponentials. All tend to the needed limits with fixed
positive budgets. Their thresholds have a finite maximum, although
these thresholds can grow with s and with inverse parameter gaps.

## 5. Complete H_tau and strict fixed-parameter quantifiers

The entire finite cover factor and terminal self-weight remain. For
fixed n,t,tau take only the scalar a limit in numeric bounds already
uniform over finite witness families. This yields C_n<=H_tau for
the same original operator; selected histories need not converge.

Eventually Q=B_n-gamma*n>0 because B_n>=n log(t)/2>=nL/4.
The exact expansion is Q+tau+A+A Q/tau. With
A<=n exp(-sigma L), B<=T_n+P*n+1 and B<=nL,

    H_tau<=T_n+(P-gamma)n
       +n[exp(-dL)+exp(-sigma L)
          +L exp(-(sigma-d)L)+1/n].

All terms tend to zero since 0<d<sigma. The fixed allowance
gamma-P-c>0 pays the FULL remainder. F(A,Q)<=H_tau does not
imply C_n<=F(A,Q), so no such inference or optimization is made.

For any fixed finite s>=4 and c<C_free(s), the strict primary
choices in the source exist by continuity: first chi>mu, then
sigma>chi/(1+chi), then d in (chi(1-sigma),sigma), then qmax<1
close to one and positive gamma<log(s*qmax). With the displayed
D/4,D/8 allowances they give gamma-P-c>D/4>0. Trim, m, Ret,
d_mu and the s-dependent e are chosen NEXT, all fixed before n
grows. The e=0 limit has strict mean and deficit slack, so the
additional smallness e<1/[2(s+1)] can be imposed simultaneously.
No old s<=5 constant, boundary tuple or n-dependent s enters this
finite-s theorem.

## 6. Limiting coefficient and the required order of quantifiers

Exact algebra gives

    (1+mu(s))/s
      =1-(1+1/s)log2+(1+1/s)log(1+1/s).

The ratio is positive because mu>0. Since 0<log2<1 and
0<=log(1+1/s)<=1/s, its difference from 1-log2 has absolute
value <=(2+1/s)/s and tends to zero. The limit 1-log2 is strictly
positive, so continuity of log proves

    C_free(s)->c_infty=-log(1-log2).

Now fix epsilon>0, choose ONE finite s with
C_free(s)>c_infty-epsilon/2, and set c=c_infty-epsilon.
This is strictly below C_free(s). The fixed-s theorem supplies
ONE fixed tuple and its own finite threshold. Hence eventually
(C_n-T_n)/n<-c_infty+epsilon. The definition of limsup yields

    limsup (C_n-T_n)/n<=-c_infty=log(1-log2).

This is not a substitution s=infinity into a kernel or a uniform
threshold argument. The order is epsilon, finite s, all fixed
parameters, THEN n. It supports every strict coefficient below the
limit eventually; it does not produce the boundary coefficient on
a finite domain or a threshold common to epsilon.

## 7. Actual finite auxiliary, identities and decision

This audit saved REGISTRATION.json before one independent local
Python Fraction auxiliary, PID 33534. It uses exact polynomial
operations for the s-dependent derivative/moment coefficients and
independent log16 rigorous intervals, rather than the author's log24,
s=10 quadrature regimes or finite differences. The auxiliary verifies
D_s=4s^2+6s+1, its positive shift at s=4, the new s^4 moment and
2s^2 variance bookkeeping, and the sign/interval of c_infty. Its
finite floor arithmetic is only a consistency check; the actual
all-s J>1/2 proof uses the s-dependent margin in Section 1, not
a numerical margin uniform over unbounded s.

Independent exact intervals give approximately
1.18138706178865<c_infty<1.18138706185749. The exact expression,
not these decimal digits, is the accepted coefficient. The finite
calculation does not prove local uniformity or quantify a dimension
threshold; Sections 1--6 supply the analytic proof.

All original source identities and author receipt links were checked,
and source hashes are rechecked at handoff. No mathematical gap was
found. ACCEPT the arbitrary fixed finite-s theorem, its complete
same-tau original-operator transfer and the stated limsup consequence,
with the explicit scope limits at the start of this report.
