# R088 free algorithmic threshold — independent full-chain ACCEPT

Frozen author source: `/Users/zhengzhihao/Desktop/T/audit/round088/new_tau/PROOF.md`,
SHA256 `75c56ce9cbd59efedc1a6863d9a125b14f30088ab85106cb98760d0507c35151`.

ACCEPT the following EVENTUAL theorem. For every fixed s in [4,5]
and every real

    c<C_free(s)=log[s/(1+mu(s))],
    mu(s)=(s-1)+(s+1)log((s+1)/(2s)),

one fixed strict tuple and a finite L_* give, for every integer n with
L=log n>=L_*,

    C_n<=H_tau=(1+A_n(t)/tau)(tau+B_n(t)-gamma*n)
       <T_n-c*n,
    tau=n exp(-dL), t=n/[(1-sigma)L],
    T_n=n log n-n log log n.

Also ACCEPT the fully fixed s=9/2 tuple in the source with c=9/10.
No explicit dimension threshold is certified here. **The conclusion
is the complete H_tau bound. A square-root F conclusion at these new
parameters is not proved and is not asserted.** R089 is outside this
review and was not read or mixed into the acceptance.

This is independent temporary staging, not a Git clone. Original files
were read only. New files are confined to `audit/round088/full`; no old
round, paper, shared MATLAB or Git state was changed, and no commit or
push occurred. No author13 or historical baseline was executed.

## 1. Direct source check: the original algorithm permits every tau>0

The complete original R027 cube-cap proof and revision_v3/cover/audit
were read directly. Both formulate the finite theorem for FIXED n,t,
every a>2t and every tau>0, BEFORE any scalar optimization. The
selection is a largest-radius remaining eligible cube, with frozen
nonnegative weights w_i=|I_i|/(2a r_i)^n and actual test

    sum_(previous selected i) w_i 1_(a beta Q_i)(center)<=tau.

Increasing the prefix adds nonnegative frozen terms, so eligibility
can only be lost and the selected radii are nonincreasing. A finite
family terminates for any positive tau, including tau>1. The first
rejection prefix supplies the full covering factor

    1+(beta+1/a)^n/tau.

Every selected star remains a true input witness. At a positive load
point, choose the smallest contributing radius and its last selected
index among ties. The previous terms capturing that center are a
nonnegative subsum of its OWN actual eligibility test and cost <=tau.
The self-term is paid separately by a^-n. None of these arguments
uses a square root, tau<=1, an a-dependent tau or a relation to A,B.

R027 takes a scalar a limit with n,t,tau held fixed, THEN optionally
minimizes a numerical expression over tau. Thus choosing
tau=n exp(-dL)>0 is legal: for each fixed dimension n, fix this same
number before ANY witness family, input, level, selection history or
finite a. Dependence on dimension is permitted because C_n is a
dimension-specific constant. tau is independent of f, witnesses,
history and a. Weight invariance under the common terminal-radius
normalization ensures it is also the same tau in normalized histories.

This verifies the new mechanism's main authorization in the ORIGINAL
finite theorem, not by extrapolating an already optimized F formula.

## 2. Fixed-parameter assumptions and valid reuse of R085

R085's local, full-target, inverse-face, weighted-cap and finite-history
proofs were read previously and are independently accepted. Only those
portions independent of its old tau payment are reused. R088 does NOT
invoke the entire old fixed-tuple theorem with an invalid old primary
condition: for the new fixed instance, sigma/2>chi(1-sigma) fails.
That old square-root payment is replaced explicitly below.

The new strict primary conditions are

    chi>mu(s)>0, 0<sigma<1,
    chi(1-sigma)<d<sigma,
    1/s<qmax<1, 0<gamma<log(s*qmax),
    c<gamma+log(1-sigma).

The interval for d is nonempty exactly when
sigma>chi/(1+chi). It implies d>0 and d<1; in fact tau may grow
with n although tau/n and tau/k decay. No claim that absolute tau
tends to zero is required.

The trim, mass, retained fraction and error are then FIXED with
positive slack: log(rhi/rlo)>m>gamma, m*Ret>gamma,
d_mu=(chi-mu)/4>0 and e in (0,1/12) satisfying the three displayed
mass/mean/deficit inequalities. Existence follows by continuity at
e=0, where the face mean is (mu+chi)/2<chi and Gamma=m*Ret>gamma.
A fixed m_face strictly between the divided transported mean and
chi therefore exists. The count floor nu>1/10 and all degree/row
constants remain positive. No constant drifts with n.

## 3. Local profiles, five events and the genuine source law

For x=(1-sigma)L, h=L^-2 and epsilon=1/(x n^2), all active
geometry is eventually uniform: full matching q_±=(r±h)/s lies
within [1/s,qmax], R lies in [r-2h,r], R>1, r<5, rho>20,
h<=.001, a0>0, beta<2, kappa>1 and t>2rhi. The exact identity
g-(1-1/s)R=epsilon*q-h/s<=0 gives g/R<=4/5 and ell>=3/5.
Inactive edge columns vanish and receive no R>=1 local assertion.
Both bands lie below norm t and are disjoint.

The accepted exact J,J1,J2 formulas include the central-deletion
case |w|>R: v=0, extra deletion has phi=0 and the -m_w contribution
remains in J. Positive-part functions are Lipschitz and the profiles
are AC, including log J; the uniform derivative bounds integrate
across onsets. Their constants have no dependence on sigma,d,tau
except through the actual geometric x,rho bounds.

At baseline R=r, full-target Fubini gives mean-one J and the exact
first-profile loss mean. The log sum's uniform envelopes are
81n/rho^2+18/rho and 324n/rho^2. The first-profile free mean is
<=mu*x, its deterministic pin is <=10x/n and its free variance is
<=100x^2/n. These imply the stated baseline event and the first-
profile failure <=400/(d_mu^2 n) eventually. The two endpoint count
events cost <=4/(e^2 n) each. Negative pinned target faces lie in D;
positive faces cost <=(.95)^(n-1) once t>2rhi. Thus exactly five
failures are each eventually <=e/8, without mutual independence,
giving lambda(A)>m*n.

The full matching product-log drift is <=240xh<=240/L. The
first-profile drift is <=120xh. Together with individual-J and
endpoint count calibration, these give the stated P0,Ulo,Uhi and
whole-matching sum J1<=(mu+2d_mu)x. ALL captured source signs and
omitted-coordinate terms are retained in U.

On each captured SOURCE face before cap/soft conditioning, the free
joint density really factors into 1_capture/(2rho J_i ell_R).
Its pin is deterministic and <p=.6. Exact free means and raw second
moments are J1_i/J_i and J2_i/J_i, so the TRUE centered variance
is their raw-second sum minus squared means, bounded by
s*x/(1-e)<10x. No pin variance, face-mixture variance, Gaussian
law or actual-increment product law is used. With a_face>0,

    P_inv(V>chi*x | face)<=10x/(a_face*x-p)^2 ->0.

This is uniform over targets, full matching, faces and signs. Only
after that control are faces mixed. The soft ratio equals one on
V<=chi*x, giving the fullsource Ret lower row. Every step precedes
cap deletion and is independent of the new algorithmic tau.

## 4. Weighted cap and full ambient degrees remain valid

The whole-matching count supplies at least nu*n counted coordinates.
The exact source-cap complement requires all z_i>=-rho+(1-q).
Each counted FREE coordinate captures its full lower strip of length
1-q and has inverse probability at least (1-q)/[2rho(1+e)].
Negative source pins cannot meet this complement; an arbitrary positive
pin removes at most ONE count. Pre-cap product independence yields
the uniform exponent (nu-1/n)*eta*x/(1+e)>0 eventually.

Pointwise bsoft<=theta/G, the full signed face sum and actual matched
dq/q integration give

    caploss/k<=M exp[-(nu-1/n)*eta*x/(1+e)] ->0.

theta cancels exactly against k. The weighted captured deletion is
paid directly; no unweighted cap mass, eta>chi condition or post-cap
independence appears. It is eventually <=e*k.

The exact radial logarithmic width with arad gives raw row bounds
2Ret(1-e)Ulo*k and M*k. For each source and matched r>=R, target
capture fractions exceed 1/2 and decrease in r. Summing EVERY target
pin sign gives coefficient <=nG. arad cancels dr/r, ordinary matched
width is <=2h and bsoft*G<=theta. Thus the complete column is
<=2k<M*k; restrictions decrease upper degrees and inactive columns
are zero. These degree arguments have no hidden tau dependence.

## 5. New relative tau payment and every original finite history

For the SAME threshold fixed before every selection,

    tau/k=rhi L^2 exp[-(d-chi(1-sigma))L] ->0.

Its rate is strictly positive. Hence it is eventually <=e*k and,
together with the weighted cap, leaves row minus tau>=b*k. The
actual tau is not the square-root choice and no old square-root-
generated history is reused.

For every a_alg>=8n^4, t<n eventually gives ccap>1. An actual
target item meeting A has radius<alpha+rhi<6 and an actual source
increment has radius>=rho>20>6. It precedes that SAME target in
the selection built with the NEW tau. For actual n>=2, beta<2,

    n^-3-(beta+2)/(a_alg+1)
      >n^-3-4/(8n^4)=(n-1/2)/n^4>0.

This and r_j<r_i force each included true source center to capture
the true target in its own original eligibility test. Integrating only
included genuine increment pieces costs at most original weights;
arad,bsoft<=1 give actual weighted row<=arad*tau<=tau. The bound
is uniform over arbitrary finite families and all their actual histories.

The genuine uncaptured measure is dominated by the full finite corner
envelope, independently of tau. On A and P this envelope equals
c_a*lambda exactly by the positive cap-containment margin. Energy,
degree and target mass scale by c_a^2,c_a,c_a. The disjoint A,P
supports allow both missing-mass terms to be paid from the TOTAL
finite-cap deficit Delta_a. Thus

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
       >=[E-tau*lambda(A)]/(M*k)
       >=lambda(A)*b/M>Gamma*n>gamma*n.

This uses c_a>=1 in the right direction. gamma>1 causes no restriction;
there is no assumed deficit<=n. The new tau is used consistently in
actual eligibility, ambient-row subtraction and finite load payment.

All local, probability, transport, retention and new cap/tau bounds
have explicit uniform envelopes tending to their required limits.
They are finitely many inequalities with FIXED positive budgets, so
one finite maximum threshold exists. The proof uses no convergence
of histories and no dimension grid. It does not compute L_*.

## 6. General L1, every positive radius and scalar a limit

For arbitrary real or complex f in L1 use |f|. Each point of a
strict maximal level set has an actual positive witness radius. Its
open a-shrunk cube covers that point; compact subsets have finite
subcovers. Every selected dilate remains a true witness. The weight
normalization permits every positive terminal radius and leaves tau
unchanged. Inner regularity and level homogeneity recover the whole
level set, not only compactly supported or truncated-scale inputs.
No radius lattice, input product law, selected-history convergence or
independent-history splice is introduced.

The original finite cover factor and terminal self term are retained:

    C_n<=[1+(beta+1/a)^n/tau]
        [tau+a^-n+alpha^-n J_a(t)-gamma*n].

For each fixed n, t and this tau, the inequality holds for every
sufficiently large finite a and every finite true-witness family.
Only this family of NUMERIC constants is passed to its a limit.
J_a->B_n(t) and the remaining factors converge, giving precisely

    C_n<=H_tau=(1+A/tau)(tau+B-gamma*n).

No supremum over inputs is exchanged with a selection limit, because
each finite-a bound already holds for the same original operator.

## 7. COMPLETE H_tau expansion and the two different new rate tests

Eventually B>=n log t/2>=nL/4>gamma*n, so Q=B-gamma*n>0.
The exact expansion is

    H_tau=Q+tau+A+A Q/tau.

There are TWO distinct exponential requirements:

* Relative row payment uses d-chi(1-sigma)>0.
* The complete cross cost uses sigma-d>0.

Using A<=n exp(-sigma L), Q<=B<=nL and tau=n exp(-dL),

    tau/n=exp(-dL), A/n<=exp(-sigma L),
    A Q/(tau*n)<=L exp[-(sigma-d)L].

The fixed positive d,sigma,sigma-d imply that the complete normalized
cost remainder

    exp(-dL)+exp(-sigma L)+L exp[-(sigma-d)L]+1/n

tends to zero. Its fixed allowance is gamma-P-c>0, with the exact
price P=-log(1-sigma). Thus H_tau<T_n-c*n with every factor and
cross term paid. The different row and cost inequalities have not
been conflated.

Although F(A,Q)<=H_tau algebraically, C_n<=H_tau does NOT imply
C_n<=F(A,Q). The optimizing square-root threshold may fail the row
payment. This audit explicitly accepts H_tau alone at the new tuple
and performs no hidden minimization back to F.

## 8. Every strict coefficient below C_free(s)

Fix s in [4,5] and c<C_free(s), and let D=C_free(s)-c>0. The
source's order of choices is valid: first chi>mu close enough that
log(1+chi)<log(1+mu)+D/4; next sigma>chi/(1+chi) close enough
that P<log(1+chi)+D/4; then d strictly between chi(1-sigma) and
sigma; then qmax<1 close enough that log qmax>-D/8 and qmax>1/s;
finally positive gamma<log(s*qmax) with
gamma>log(s*qmax)-D/8. Continuity and log(s*qmax)>0 ensure
these choices exist. They give

    gamma>log s-D/4,
    P<log(1+mu)+D/2,
    gamma-P-c>D/4>0.

Trim, mass, Ret, d_mu, e and m_face can now be chosen with strict
fixed slack before n grows. All preceding proofs apply, giving one
finite threshold for each FIXED tuple. No boundary coefficient,
n-dependent parameter drift or common threshold over different gaps
is asserted. The theorem is restricted to s in [4,5]; neither s>5
nor a global maximum is part of this acceptance.

## 9. Fully fixed s=9/2 instance and actual new evidence

The exact fixed tuple s=9/2, chi=.8, sigma=.45, d=.445,
qmax=.9999, trim=10^-6, m=1.503, Ret=.99999, e=10^-5,
gamma=1.5, m_face=.796 and c=.9 meets all strict requirements.
Its inverse source mean upper is <.795697623<.796, its source
face gap is .004, and its true-variance Chebyshev envelope is
10x/(.004x-.6)^2, eventually <=1-Ret. The weighted cap rate is
positive and every uniform geometry/target event is eventual as above.

The NEW two rate margins are EXACTLY

    d-chi(1-sigma)=1/200,
    sigma-d=1/200.

Consequently tau/k=rhi L^2 exp(-L/200) tends to zero and the
full scalar remainder is bounded by

    exp(-89L/200)+exp(-9L/20)+L exp(-L/200)+1/n.

The fixed price P<299/500 and 1.5-299/500-.9=1/500 pay it
eventually. Therefore C_n<=H_{n exp(-89L/200)}<T_n-.9n at one
finite unquantified threshold. No claim at 4096 or any explicit L0
is made.

This audit saved REGISTRATION.json before its independent bounded
Python Fraction run, PID 32499, using log20 enclosures rather than
the author's log24/13-check program. Five grouped checks passed;
the author program and old baselines were not replayed. They confirm
fixed mean/mass/deficit/price gaps, both new rates, and the full-H
allowance, not a finite dimension threshold. The original R027
universal-tau proof and v3 audit were read directly, and source
identities are locked and rechecked. Computation is supporting exact
arithmetic; Sections 1--8 supply the proof.

## 10. Additional complete-function comparisons supplied by the main agent

The separate primary_comparison registration, script and actual PID
32481 receipt were read, not rerun. Its four new exact comparisons
give c87=619/400+log(37073/100000)<.556 and
200000/(37073*65536)<.0001. The original scalar integral
I<=2t/(n-1) and nonnegative core imply

    Phi87>=B87-(619/400)n
       >=T_n-[c87+200000/(37073L)]n.

Thus Phi87>T_n-.5561n for L>=65536 on Phi87's TRUE domain.
Combining with H88<T_n-.9n gives Phi87-H88>.3439n on the
actual common validity domain. In particular the known explicit
R087 branch requires L>=2^80, and the new H88 requires its own
uncomputed L_*; the comparison does not activate either at 65536.
These are lower bounds for an old upper-bound FORMULA, never lower
bounds for C_n. The direction is correct even though the new function
is H and the old one is F; neither is silently optimized or replaced.

The same receipt gives C_free(5)>.9493 by LOWER-bounding the log
after an UPPER bound on mu(5). This supports only the family
consequence that the fixed strict coefficient .9493 is eventually
attainable for some fixed s=5 tuple. It supplies no explicit tuple,
dimension threshold or result for the separate R089 task, which
was not reviewed here.

## 11. Final decision and boundaries

No gap was found. ACCEPT the new arbitrary-tau mechanism, complete
H_tau eventual theorem, fixed s=9/2 coefficient .9 and strict fixed-s
family c<log[s/(1+mu(s))]. Accept the stated full-function comparison
on its actual common domain. Keep all original run identities and
eventual scopes separate. No F bound, explicit threshold, s>5 result,
global optimum, paper edit or old-branch domain extension follows.
