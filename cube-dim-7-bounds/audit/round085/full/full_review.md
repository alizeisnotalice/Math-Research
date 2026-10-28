# R085 fixed parameter family and limiting coefficient — independent ACCEPT

The frozen source reviewed is PROOF.md with SHA256
ded160a902cdd543a7be9ebbb9c1f4d076c33f866d8aa96e9481f92d69bb41fa.
This is independent staging, not a clone. No author file, old result,
original repository, paper or Git state was modified. Only this owned
audit directory was written. No author baseline/check program, parameter
scan, shared MATLAB, commit or push was used. The supplied accepted
finite-cover framework is used as background; every new local, probability,
weighted-cap, history and parameter-limit interface is checked below.

Accepted statement: for each fixed s in [4,5] and each real c<C(s),
one fixed strict tuple gives a finite L_*(s,c) such that for every integer
n with log n>=L_*(s,c), the original centered axis-parallel Lebesgue
common-radius cube maximal operator, for arbitrary L1 inputs and all
positive witness radii, satisfies

    C_n<=F(A_n(t),B_n(t)-gamma*n)<T_n-c*n,
    t=n/[(1-sigma)log n], T_n=n log n-n log log n.

The unique recipe maximum is

    s_*=1/(2 exp(-1/2)-1),
    c_*=-log(3-4 exp(-1/2))=0.5553395620244... .

For every epsilon>0 the eventual bound with c_*-epsilon follows,
and limsup (C_n-T_n)/n<=-c_*. No finite-domain inequality with
coefficient exactly c_* or common threshold over epsilon is asserted.
This is optimal only for the stated sufficient-parameter recipe on
s in [4,5], not globally over cube covering methods.

## 1. Strict choices and the uniform local geometry

For fixed primary s,chi,sigma,qmax,gamma satisfying (4), chi>mu(s)>0
and cT=sigma/2-chi(1-sigma)>0. Positive gamma<log(s*qmax)
allows a fixed trim with log(rhi/rlo)>gamma. The subsequent choices
gamma<m<log(rhi/rlo) and gamma/m<Ret<1 have strict slack.
At e=0, the three conditions (5) are strictly satisfied:
the mass exceeds m, mu+2d_mu=(mu+chi)/2<chi, and Gamma=m Ret>gamma.
Continuity therefore supplies a sufficiently small positive FIXED
e<1/12 satisfying all of them simultaneously. Then m_face can be
chosen strictly between (mu+2d_mu)/(1-e) and chi. Every constant
is fixed before n grows. Uniformly in s in [4,5],
nu=1/s-(1+1/s)e>1/10.

The complete matching interval q_±=(r±h)/s lies inside [1/s,qmax]
once 2h<delta. Its R lies in [r-2h,r] and exceeds one. The exact
identity g-(1-1/s)R=epsilon*q-h/s gives the new relative domain
g/R<=1-1/s<=4/5 when epsilon*qmax<=h/s. This guard is eventual
because epsilon*qmax/h is a fixed multiple of L/n^2. Uniformly on
matching, rho grows beyond 20, r<rhi<5, and beta<2, kappa>1,
t>2rhi, a0>0 hold. Bands lie below t and are disjoint. Unmatched
edge sources can have R<=1 but their final columns vanish; local
R>=1 bounds are not applied there. R>0 and G>0 on the full source
band eventually.

local_uniform.md independently derives the exact J,J1,J2 formulas
and AC bounds on this enlarged domain. Central deletion |w|>R is
explicitly included: v=0, the central portion has phi=0, and the
necessary -m_w remains in J. The deletion cannot reach the opposite
edge. With h<=.001 and rho>=20, stronger independent estimates are

    |J-1|<=6/rho<9/rho,
    |d_q log J|<=243/(4rho)<120/rho,
    J1<=9/(5rho)<2/rho,
    |d_q J1|<9.002/rho<60/rho,
    J2<=18/(25rho)<1/rho.

The profiles and log J are AC through all positive-part kinks.
Differentiation holds s,t,h,epsilon,w fixed. The author's coarser
constants 17/(2rho),178/(3rho),583/(10rho) are also arithmetically
safe; none is imported from a smaller s domain.

## 2. Exactly five full-target events

At baseline R=r, averaging a full uniform target coordinate capture
indicator gives ell_R(z). Fubini cancels the inverse factor, hence
E_W J=1. From the author's safe |J-1|<=9/rho and log inequalities,

    -81/rho^2<=E_W log J<=0,
    Var_W(log J)<=324/rho^2.

The target pin has deterministic absolute log contribution <=18/rho.
Thus the baseline sum's absolute mean is bounded by
81n/rho_b^2+18/rho_b and its variance by 324n/rho_b^2. Both vanish
uniformly in target radius and signed face because n/rho_b^2<=25x^2/n.

For J1, exact full-target Fubini yields (2R/rho)H(a), where
a>=a_s=(s+1)/(2s). H decreases on [a_s,1] and R/q<=s, so
the free-sum mean is <=mu(s)x, mu(s)=2s H(a_s)>0. The actual pin
adds at most 2/rho_b<=10x/n. For large n its addition is
<=d_mu*x/2. Free variance is <=4n/rho_b^2<=100x^2/n. Chebyshev
across the remaining d_mu*x/2 gap proves exactly 400/(d_mu^2 n).

The true target pin contributes one to each endpoint Nsum because
a0<r. The n-1 excess indicators are independent Bernoulli(a0/r),
so each endpoint count failure is <=4/(e^2 n). For D, negative
pinned faces are inside D; positive pinned faces have complement
free fraction 1+1/(2t)-1/(2r)<=1-1/(4rhi)<19/20 once t>=2rhi.
Their failure is <=(19/20)^(n-1).

These are precisely five events: baseline log, two counts, target D
and first profile. Each failure is eventually <=epsA/8, uniformly
in the complete target band and signed faces. No independence between
the five events is needed. Union bound gives target mass
lambda(A)>=(1-epsA)n log(rhi/rlo)>m*n. These laws apply to full
ambient target faces, not to actual algorithmic increments.

## 3. Complete matching and genuine source-face retention

The q interval has width 2h/s and n/rho<=s*x. The uniform local
log derivative gives shift <=240xh<=240/L. It is eventually <=e,
so the baseline event yields product J_i in [P0,1/P0]. Separately,
the pointwise local value bound yields every J_i in [1-e,1+e].

ellpin' =-(1-epsilon)h/(2R^2)<0 and Nsum is nondecreasing.
Its normalized mean ratio is <=1 at q_+ and >=1-2h-1/n at q_-.
The endpoint fluctuations change it by <=e/4. Thus when
2h+1/n<=3e/4 the ratio stays in [1-e,1+e] on the WHOLE interval.
Summing all captured signed source faces and their omitted-coordinate
terms gives exactly U=(product J_i)/(2ellpin) sum_i N_i/J_i,
with Ulo*n<=U<=Uhi*n. No source sign or omitted coordinate is lost.

The J1 derivative gives shift <=120xh. Once 120h<=d_mu, the
complete-matching sum is <=[mu+2d_mu]x. Fix any captured signed
source face before cap deletion and soft clipping. Its free joint
density really factors as product 1_capture/(2rho J_i ell_R), i!=j.
The pin is deterministic with loss <p=.6. Hence

    EV<=p+[(mu+2d_mu)/(1-e)]x<p+m_face*x,
    Var V=sum_(i!=j)[J2_i/J_i-(J1_i/J_i)^2]
         <=n/[rho(1-e)]<=s*x/(1-e)<10x.

This is the TRUE centered variance of the face loss, not mixture
variance or the raw second moment of the full sum. The pin adds no
variance. With fixed a=chi-m_face>0,

    P_inv(V>chi*x | face)<=10x/(a*x-p)^2 ->0.

Thus every face and sign retains at least the fixed Ret eventually,
uniformly over retained targets and complete matching. The fixed
positive slack is essential but need not be bounded below across
different parameter tuples. Only after facewise control may faces
be mixed. min(1,G/theta)=1 on V<=chi*x, giving the claimed raw
soft-row bounds. No model, Gaussian approximation, lower variance
or additional second-profile target event is needed.

## 4. Weighted cap and entire degrees

The whole-matching count and a0/R>=1/s give at least nu*n counted
coordinates. Source D-complement is exactly all
z_i>=-rho+(1-q). Negative source pins are excluded. Each counted
free coordinate captures its whole lower strip of length 1-q; the
endpoint conditions follow from |w_i|<a0, rho>20 and r<5. Its
pre-cap inverse probability is at least (1-q)/[2rho(1+e)]. A
positive pin omits at most one counted coordinate. Independence
before conditioning gives

    P_inv(D^c | face)<=exp[-(nu-1/n)*eta*x/(1+e)],
    eta=(1-qmax)/(2qmax)>0.

Multiplying bsoft<=theta/G, summing the full face coefficients,
and integrating only over the actual matched dq/q interval gives
caploss/k<=M exp[-(nu-1/n)*eta*x/(1+e)]. The radial log-width
is uniformly within 1±e of 2h/rhi. theta cancels against k.
Since nu-1/n>=nu/2 eventually, this loss decays at a fixed positive
exponential rate for each tuple and is eventually <=e. No condition
chi<eta or post-cap independence is imposed.

For a fixed source and any matched target radius r>=R, the
target-coordinate capture fractions exceed 1/2 and decrease in r.
At most one target pin sign is captured when its fraction is <1,
and at most two otherwise. Summing all target faces therefore gives
coefficient <=n product ell_i(r)<=nG. arad cancels dr/r; ordinary
matching width is <=2h and bsoft*G<=theta. The full column is
<=2k<M*k. Target/source restrictions only reduce upper degrees;
inactive columns are zero. Raw rows have the stated lower bound
2Ret(1-e)Ulo*k and upper bound M*k.

## 5. Eventual threshold, same tau and original finite history

All relevant envelopes are uniform on the fixed target band and full
matching domain: local and geometry guards, polynomial(L)/n target
errors, count errors, geometric target-cap errors, O(1/L) and
O(1/L^2) transports, source-face Chebyshev error, and exponentially
decaying cap/tau payments. Fixed positive budgets give finitely many
finite thresholds, whose maximum is finite. This does not assert a
uniform threshold over s,c or parameter gaps, and none is needed.

B_n>=n log t/2>=nL/4>gamma*n eventually, so the same global
tau=sqrt[A_n(B_n-gamma*n)] is positive even if gamma>1. With
A_n<=exp x and B_n<=nL, tau/k<=rhi L^(5/2)exp(-cT L), where
cT=sigma/2-chi(1-sigma)>0. This same selection tau is eventually
<=e*k. Paying cap and tau leaves row minus tau >=b*k, while both
degrees remain <=M*k and lambda(A)>m*n.

For every original finite a_alg>=8n^4, ccap>1 eventually. An actual
target meeting A has ccap*r_j<alpha+rhi<6, hence r_j<6. An
actual source has r_i>=rho>20>6 and precedes that SAME target.
The margin n^-3-(beta+2)/(a_alg+1)>0 and r_j<r_i ensure that
every included true source center lies within alpha*beta*r_i of
that target in its original eligibility test. Bounded kernel weights
give actual weighted row<=arad*tau<=tau. No changed selector or
history is introduced.

On A and P the finite envelope equals c_a*lambda exactly, with
c_a=alpha^-n>=1, by the displayed positive cap margin. Energy,
degree and mass scale by c_a^2,c_a,c_a. The two missing-mass terms
cost at most c_a*M*k times their missing masses on disjoint A,P,
whose sum is <=Delta_a. Actual energy<=tau*c_a*lambda(A). Thus

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
            >=[E-tau*lambda(A)]/(M*k)
            >=lambda(A)*b/M>Gamma*n>gamma*n.

The c_a>=1 direction is correct; no gamma<1 or deficit<=n
restriction is required. This holds uniformly for original histories.

For any L1 input use |f|. True positive-radius witnesses provide open
shrunk covers of compact subsets of each strict maximal level set;
finite subcovers are available and selected dilates remain genuine
witnesses. Normalization permits every positive terminal radius.
Inner regularity recovers the whole level set. No compact support,
truncated radius range, lattice, independent-increment law or history
limit is introduced.

## 6. Complete cost and fixed-parameter approximation

The entire finite cover factor, captured budget and terminal self-weight
are retained. Taking only the numeric a_alg limit after uniformity
gives exactly C_n<=F(A_n,B_n-gamma*n). The full positive-root
remainder, including the scalar +1/n term, is bounded by
exp(-sigma L)+2sqrt L exp(-sigma L/2)+1/n and tends to zero.
The exact price is P=-log(1-sigma). Since the fixed allowance
gamma-P-c is positive, the complete root cost is eventually <T_n-c*n.
The cross-term is paid, not discarded.

For c<C(s), let d=C(s)-c>0. The proposed order of primary choices
is valid by continuity. It gives strict chi>mu, strict sigma above
2chi/(1+2chi), qmax<1 but close enough to one, and positive gamma
strictly below log(s*qmax). The inequalities imply

    gamma>log s-d/4,
    P<log(1+2mu)+d/2,
    gamma-P-c>d/4>0.

Trim, mass, Ret, d_mu, e and m_face are then fixed with their
strict slacks before n grows. No boundary value or n-dependent
parameter is inserted into a finite theorem. Conversely the primary
constraints force gamma-P<log s-log(1+2chi)<C(s), so C(s) is
the supremum approached by this specific sufficient recipe.

## 7. Analytic maximum and limsup quantifiers

With l=log((s+1)/(2s)), direct differentiation yields
mu'=1+l-1/s and C'=[1+2l]/[s(1+2mu)]. The denominator is
positive, while l'=-1/[s(s+1)]<0. The sole stationary point is
l=-1/2, namely s_*. The elementary exponential enclosures in the
source imply 3/5<exp(-1/2)<5/8, so 4<s_*<5. The derivative
changes from positive to negative, proving a unique maximum.
At that point mu=(s_*-3)/2, whence
C(s_*)=log[s_*/(s_*-2)]=-log(3-4exp(-1/2)).

For any epsilon>0 fix s=s_* and c=c_*-epsilon. The already proved
strict-parameter theorem supplies one fixed tuple and its own finite
threshold. It gives eventually (C_n-T_n)/n< -c_*+epsilon.
By the definition of limsup, allowing arbitrary positive epsilon
yields limsup<=-c_*. This is not an exchange of parameter and
operator limits. It gives neither a common threshold for all epsilon
nor C_n<T_n-c_*n at a finite domain.

## 8. Complete-function comparison with R083 on the common domain

The main agent supplied a separate single-run exact Fraction receipt
for four fixed comparisons, not the author baseline. Its relevant
enclosures are c_*>5553/10000 and

    c_old=1501/1000+log(387/1000)<5517/10000.

The first is also implied by this audit's independent c_* interval.
The supplied old-price enclosure is not replayed here. The original
I<=2t/(n-1) bound and nonnegative core give, with x=.387L,

    Phi83>=B83-(1501/1000)n
         >T_n-[c_old+2000/(387L)]n.

The direction is correct: it is a lower bound for the old complete
upper-bound formula, not a lower bound for C_n. For L>=65536,
2000/(387L)<1/10000, verified directly by
20000000<387*65536=25362432. Therefore Phi83>T_n-.5518n.

Choose c=5553/10000<c_* in the R085 family theorem, fix its strict
tuple, and use its own eventual threshold. Its COMPLETE Phi85 then
satisfies Phi85<T_n-.5553n. Subtracting yields

    Phi83-Phi85>35n/10000

on the actual common validity domain, enlarged to L>=65536. Every
inequality uses an actual fixed R085 tuple; the supremum coefficient
is not inserted as a finite-construction coefficient. No old branch
is activated outside its true domain and no dimension threshold for
the new tuple is computed. The strict margin follows from
.5553-.5518=.0035; no root remainder has been omitted.

## 9. Actual evidence, provenance and decision

One independent registered Python/SciPy run, PID 31041, executed
check.py. Independent Fraction log24 bounds and alternating-series
enclosures for exp(-1/2) certify the displayed s_* and c_* intervals,
endpoint derivative signs, fresh local constants, positive uniform
count floor and true-variance bound. Twenty-four direct integral
cases have maximum error 2.22e-16; finite differences are diagnostics
only. receipt.json records exact rational values and actual run/script
identities. The frozen source hash was rechecked after execution and
at final audit packaging. No dimension threshold was evaluated.

These computations support bounded consistency; the analytic local
proof and whole-chain review supply uniformity and quantifier control.
No mathematical gap was found. ACCEPT the fixed-parameter eventual
family, its unique analytic recipe maximum, and the stated limsup
consequence, with the limitations above. Main-agent acceptance remains
the final proof responsibility.
