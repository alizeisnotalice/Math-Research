# R083 complete s=9/2 eventual theorem — independent ACCEPT

Acceptance object:

    There exists finite L_* such that for all integers n with log n>=L_*,
    C_n <= Phi83=F(A_n(t),B_n(t)-(1501/1000)n) < T-(11/20)n,
    t=n/[(387/1000)log n], T=n log n-n log log n,
    F(A,B)=(sqrt A+sqrt B)^2.

This is the original Lebesgue centered axis-parallel common-radius cube
maximal operator, for arbitrary L1 inputs, positive strict levels and
all positive witness radii. No explicit threshold is certified, including
L=4096. The statement gives <=Phi83; strictness is proved for the
comparison Phi83<T-.55n. It does not claim a globally optimal branch.

The frozen author PROOF.md has SHA256
6d2dff6fedf7cde05d4823830bca674235bd84be4608fc35ea129aa44d649c0e.
All four R083 files and the three authorized R081 background files
were read. The final explicit inverse-product-density and actual
19-check receipt paragraphs were included in this identity and reviewed.
REGISTRATION.json records all seven source hashes. This directory is
independent temporary staging, not a Git clone. No author file, original
repository file, paper, accepted ledger, shared MATLAB or Git state was
modified; no commit/push occurred. Only files in this directory were
written. Neither the author's 19-check program nor any old baseline was
executed. The accepted R081 finite-cover framework is used as supplied
background; every interface changed by s=9/2 is checked below.

## 1. New local formulas and explicit sufficient domain

local_s45.md independently rederives the s=9/2 coordinate lemma from
the one-dimensional integrals. Its sufficient local guards include
h<=.001, rho>=16, 1<=R<=4.5, |w|<=r<=min(4.5,R+2h), and
0<=g/R<=7/9. They are eventual guards of the actual geometry.

The disjoint source boundary strips have ell>=11/18. Integration
gives exactly the author's J,J1,J2 formulas. When |w|>R, deletion
continues beyond the near boundary layer into the central ell=1,
phi=0 region. The -m term remains in J; setting v=0 correctly
removes the near-edge J1,J2 terms. Deletion cannot reach the opposite
edge. Thus the formula is valid throughout the complete matching
interval, rather than only at its baseline endpoint.

The independent bounds are stronger than needed:

    |J-1|<=9/(2rho)<7/rho,
    |d_q log J|<=3873/(88rho)<90/rho,
    J1<=9/(8rho), |d_q J1|<6.547/rho<40/rho,
    J2<=3/(8rho).

All differentiations hold t,h,epsilon,w fixed. Positive parts are
Lipschitz and denominators are bounded away from zero. J>=2/3 on
the sufficient domain, so log J is AC as well. Integrating a.e.
derivative bounds across onset kinks is legitimate. No s=4 profile
constant is borrowed, and the cubic positive prefix 79/48 is sufficient
for exp(.5)>18/11, whereas the quadratic prefix is insufficient.

## 2. Active geometry, full-target laws and five events

The complete source interval is q_±=(r±h)/(9/2). Once 2h<delta it
lies inside [2/9,.9999]. On it R=(9/2)q-h lies in [r-2h,r] and
exceeds one. The exact relation

    g-(7/9)R=epsilon*q-(2/9)h<=0

is eventual because epsilon*qmax/h is a fixed multiple of L/n^2.
Also rho=tq grows uniformly, and a0>0, g>0, beta<2, kappa>1,
t>rhi and rho>4rhi hold eventually. The full source band has
R>0 and G>0, but the R>=1 local bounds are needed only on matched
columns; unmatched edge columns are zero. Both bands lie below t
and are disjoint.

At qb, R=r. Averaging the full-target coordinate capture indicator
gives exactly ell_R(z). Fubini cancels 1/ell_R, proving E_W J=1.
With |J-1|<=7/rho<=.5, log inequalities give

    -49/rho^2<=E_W log J<=0,
    Var_W(log J)<=196/rho^2.

The true target pin is deterministic with |log J|<=14/rho. The
sum has absolute mean <=49n/rho_b^2+14/rho_b and variance
<=196n/rho_b^2. Both tend uniformly to zero as n/rho_b^2<=s^2*x^2/n.

The same cancellation for J1 gives

    E_W J1=(2R/rho)H(a), H(a)=1-a+a log a, a>=11/18.

H decreases in a and R/q<=s, yielding free-sum mean <=mu*x,
mu=7/2-(11/2)log(18/11)<.7914. The actual pin adds at most
9/(8rho_b)<=x/100000 for n>=506250. The mean is therefore
<=.79141x, the free variance <=(6561/256)x^2/n, and Chebyshev
across the remaining x/50000 gap gives 64072265625/n.

Each endpoint count consists of a true pin of count one and n-1
independent Bernoulli(a0/r) excess indicators. Its failure is at most
4/(e^2 n). For target D, a negative pinned face belongs to D; on a
positive pinned face the cap-complement free fraction is
1+1/(2t)-1/(2r)<=8/9+.01<.9 once t>=50. Its failure is
<=.9^(n-1). These are exactly five target events: log product,
two endpoint counts, D, and first profile. Each failure is eventually
<=epsA/8, uniformly in radius and signed face. Their mutual
independence is unnecessary. Union bound gives

    lambda(A)>=(1-epsA)n log(rhi/rlo)>1.503n.

These ambient target laws are not assigned to actual increments.

## 3. Whole matching, true source variance and narrow retention gap

The interval length is 2h/s. With n/rho<=s*x, the log-product
transport is at most 180xh=69.66/L. Together with the baseline
event, it yields product J_i in [P0,1/P0] eventually, where
P0=1-2e. Individual J_i in [1-e,1+e] follow separately from the
local value bound, not from the product estimate.

ellpin=(1+a0/R)/2 has derivative -(1-epsilon)h/(2R^2)<0,
and Nsum is nondecreasing. At q_+ the mean ratio
Nsum/(2n ellpin) is <=1. At q_- its shortfall from one is
at most 2h+1/n, using a0<=1 and R,r>=1. Endpoint deviations
change the ratio by at most e/4. Consequently the ratio lies in
[1-e,1+e] on the whole interval once 2h+1/n<=3e/4.

Summing every captured signed source face gives exactly

    U=(product J_i)/(2ellpin) * sum_i N_i/J_i,
    Ulo*n<=U<=Uhi*n.

The omitted coordinate and both signs are accounted for. J1 transport
costs at most 80xh; once 80h<=1e-5, the whole-matching profile
sum is <=.79144x.

Fix a captured signed source face before cap or soft conditioning.
Its pin has deterministic loss <=log(18/11)<.5<.6. The free
joint density is the product over i!=j of
1_capture/(2rho J_i ell_R). Thus these coordinates really are
independent, with mean J1_i/J_i and raw second moment J2_i/J_i.
The true centered variance, with the deterministic pin omitted, is

    Var V=sum_(i!=j)[J2_i/J_i-(J1_i/J_i)^2]
         <=n*(3/(8rho))/(1-e)
         <=27*x/[16(1-e)]<2x.

No mixture variance or whole-sum raw second moment has been used.
The mean is <=.6+[.79144/(1-e)]x<.6+.79147x. Therefore

    chi*x-EV >=(3/100000)x-3/5=:u(x),
    P_inv(V>chi*x | face)<=2x/u(x)^2

when u>0. This positive gap survives the inverse normalization.
For x>.6/(3e-5), the last bound decreases, as its derivative is
-2(ax+p)/(ax-p)^3<0. It tends to zero and is eventually <=1e-5,
giving Ret=.99999 on every face and sign. Only after this uniform
facewise estimate may faces be mixed. The soft ratio
bsoft/(theta/G)=min(1,G/theta) equals one on V<=chi*x.
No Gaussian, compound-Poisson, finite-dimension model, lower-variance
assumption or extra second-profile target event enters.

## 4. General-s weighted cap and complete degrees

From a0/R>=2/9 and the whole-matching count,

    mcount>=nu*n, nu=2/9-(11/9)e=.22221>0.

The source cap complement is exactly all
z_i>=-rho+(1-q). Negative source pins cannot lie there. For each
counted free coordinate the entire lower strip of length 1-q is
captured, using |w_i|<a0 and rho>4rhi. Under its pre-cap inverse
law that strip has probability at least (1-q)/[2rho(1+e)].
Any positive source pin omits at most one counted coordinate.
Product-law independence and 1-y<=exp(-y) give

    P_inv(D^c | face)<=exp[-(nu-1/n)*eta*x/(1+e)],
    eta=(1-qmax)/(2qmax)=1/19998>0.

This is applied before conditioning. Multiplying bsoft<=theta/G,
summing all nonnegative signed-face coefficients, and integrating
only the actual dq/q interval gives, after the exact radial
normalization and division by k,

    caploss/k<=M exp[-(nu-1/n)*eta*x/(1+e)] ->0.

The theta factors cancel; there is no chi subtraction in this rate
and no condition chi<eta. Once nu-1/n>=nu/2, the L exponent has a
fixed positive rate, independently checked as >2.15e-6. The loss
is eventually <=e. Raw rows are between
2Ret(1-e)Ulo*k and M*k.

For a fixed source and every matched target radius r>=R, the
coordinate fraction ell_i(r)>1/2 and decreases with r. Hence
product ell_i(r)<=G. If ell_i(r)<1, at most one target pin sign
is captured; if ell_i(r)=1, at most two are. Therefore each
signed-face contribution satisfies T_i/(2ell_i(r))<=1 and the
complete target-face coefficient is <=nG. arad changes dr/r to
dr/rhi; the ordinary matching width is <=2h. Since bsoft*G<=theta,
the entire column is <=2k<M*k. Restrictions to A and P decrease
upper degrees, and inactive columns are zero.

## 5. One eventual threshold, positive roots and same tau

All geometry, local and good-target bounds above are uniform in the
full compact target band and complete matching domain. The bad
target probabilities are bounded by polynomial(L)/n, O(1/n), or
.9^(n-1); product/profile transports are O(1/L), O(1/L^2).
Source-face failure tends to zero as 1/x, and cap payment decays
exponentially. These are finitely many inequalities with fixed
positive budgets, so their finite thresholds have a finite maximum.
No old threshold is transferred and no selector/history convergence
is asserted. Positive exponential-series terms prove the stated
polynomial-versus-exponential limits.

The original scalar cap identity gives B_n>=n log t/2. Eventually
log t>=L/2, whence B_n>=nL/4>1.501n. Thus gamma>1 is harmless:
the same global tau=sqrt[A_n(B_n-gamma*n)] is positive.
Using A_n<=exp x and B_n<=nL gives

    tau/k<=rhi*L^(5/2)*exp[-(379/2000000)L] ->0.

The exponent is strictly positive despite its narrow value. This
very same tau is eventually <=e*k. Paying cap loss and tau gives
row minus tau >=b*k. Both ambient degrees remain <=M*k. The
independent rational enclosure verifies (1-e)Gamma>1.501, with
Gamma=mL*b/M. No gamma<1 or Gamma<1 restriction is needed.

## 6. Original finite history and arbitrary L1/all-scale transfer

For every a_alg>=8n^4, normalize by the last minimum positive
contribution radius and reflect the last center. Eventually
a_alg>2t, ccap>1. A true target item meeting A obeys

    ccap*r_j<=||u-w||<alpha+rhi<6,

so r_j<6. A true source item satisfies r_i>=rho>4rhi>6 and
precedes that same target item in the original largest-radius
history. The displacement guard is

    alpha*beta-kappa-2/(a_alg+1)
      =n^-3-(beta+2)/(a_alg+1)>0.

Because r_j<r_i, included actual centers have distance
<kappa*rho+(r_i+r_j)/(a_alg+1)
<alpha*beta*r_i and are eligible in that target's original test.
Both kernel weights are <=1, so actual weighted row<=arad*tau<=tau.
The target separation constant is 6 and the source bound exceeds
6 explicitly; neither the old s=4 constant 5 nor an altered history
is substituted.

The finite envelope equals c_a*lambda on A and P, c_a=alpha^-n>=1,
because the cap margin is exactly
[(beta+1)m-1]/(a_alg+1)>0. Energy, degree and target mass scale by
c_a^2,c_a,c_a. The two missing-mass terms in the energy difference
cost at most c_a*M*k times their missing masses on A and P.
Those sets are disjoint, so their sum is bounded by the total
finite-cap deficit Delta_a. Actual energy<=tau*c_a*lambda(A), hence

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
            >=[E-tau*lambda(A)]/(M*k)
            >=lambda(A)*b/M>Gamma*n>1.501n.

The c_a>=1 direction is correct. There is no deficit capped by n,
division by 1-gamma, or missing finite-a factor. The estimate is
uniform over all the original finite histories.

For any L1 input use |f|. Every point of a strict maximal level
set has a true positive-radius witness; its open shrunk cube
covers the point. Compact subsets admit finite such covers and
selected dilates remain true witnesses. Every positive terminal
radius is allowed by normalization. Inner regularity and level
homogeneity then give the whole strict level set. No compact
support, radius truncation, scale lattice or increment iid law is
introduced, and no limit of histories is taken.

## 7. Complete root cost and accurate old-function comparison

At every finite a_alg the entire original cover factor, captured
budget and terminal self-weight are retained. Only in the numerical
inequality, already uniform over finite witness families, take
a_alg to infinity. The exact result is

    C_n<=(1+A_n/tau)(tau+B_n-gamma*n)=Phi83.

The original scalar bounds are
B_n<=T+n log(1000/387)+1 and
[A_n+2sqrt(A_n B_n)]/n<=exp(-.613L)+2sqrt L exp(-.613L/2).
The complete positive-root remainder, including 1/n, is eventually
<.001. Since log(1000/387)<.95 and 1.501-.95-.55=.001,
this gives Phi83<T-.55n. The root cross-term is paid explicitly.

For the old complete functions, changing variables
v=1/rho-1/t in I_n gives the author's exact integral inequality
I_n<=2t/(n-1). The nonnegative core and n/(n-1)<=2 imply
B_n>=n log t-2n/x. At R081's x=.4L this is a 5/L loss,
not 4/L. The accurate price log(5/2)>.916 yields, on its true
positive-root domain,

    Phi81>=B81-1.3n>T-n[.384+5/L].

Using merely the weaker .6 price would not separate the branches.
For R080, -log(.89)>=.11 and 2/(.89L)<=4/L give

    Phi80>=T-n[.28029+4/L].

At L>=64, .55>.384+5/64 and .55>.28029+4/64. Thus after
enlarging only the unknown new threshold to >=64, Phi83 is strictly
below both old complete upper-bound formulas on their actual common
validity domains. These are lower bounds for old formulas, never
lower bounds for C_n; they do not activate old branches outside their
accepted domains or compute the new threshold.

## 8. Actual independent evidence and final decision

One registered bounded Python/SciPy run executed check.py with
PID 29940. It uses independent exact Fraction log20 enclosures,
fresh derivative arithmetic and sensitive full-cost margins; it
imports no author code. Twelve direct quadrature cases spanning
full capture, edge deletion, central deletion and signs match all
three local formulas to 2.22e-16 maximum absolute error. Finite
differences are diagnostic only. receipt.json records exact rational
values, executable/script identities, display-only values and all
case results. All seven frozen source hashes were rechecked at the
end of the actual run and at final audit packaging.

The computation certifies only bounded arithmetic consistency.
The analytic local proof and the whole-chain review supply the
uniform/eventual reasoning. Neither these checks nor the author's
receipt are a threshold certificate. No mathematical gap was found
in the changed local, five-event, source-face, weighted-cap,
same-history, all-L1/all-scale or complete-cost interfaces.

Decision: ACCEPT the stated eventual theorem and common-domain
complete-function comparisons, with finite unquantified L_*.
