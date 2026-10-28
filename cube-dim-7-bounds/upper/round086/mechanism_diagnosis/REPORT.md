# R086: three bounded mechanism diagnoses, with the threshold scope explicit

Status: diagnostic results, not a new C_n bound. Exactly three mechanisms
were considered and this batch stops here. The .55533956 obstruction below
belongs to the specified pointwise row/column deficit recipe with the
SQUARE-ROOT threshold. It does not apply to arbitrary choices of the same
global threshold, actual C_n, all energy certificates, or all covering
methods. The new non-root threshold mechanism in R088 already bypasses
that particular obstruction at the scalar-cost level; its operator
candidate remains outside this batch and under separate review.

## Common diagnostic class and complete cost

Let L=log n, x=(1-sigma)L, t=n/x, 0<sigma<1 fixed. Targets lie in a
fixed compact radial set 1<r<Rmax and sources have q=rho/t in fixed
compact bands [q0,q1] with q0>0 and q1<1. There are finitely many
branches, independent of n. Each branch f_j is C1, strictly increasing,
with 0<m_j<=f_j'<=M_j<infinity. Target radial bands are trimmed so
the complete inverse matching intervals exist:

    I_j(r)={q: |r-f_j(q)|<h}, h=L^-2,
    R_j(q)=f_j(q)-h<=r,
    a0=(1-epsilon)q, kappa*rho=rho+a0.

No turning-point map, n-dependent geometry, growing number of branches
or exponential matching multiplicity is excluded by this note. Such a
case needs a new analysis. For fixed compact domains the coordinate
formulas and AC estimates remain bounded, but the old sign of ellpin'
is NOT automatic: its numerator is f_j-h-q*f_j'. Whole-matching count
calibration would have to use its O(h) variation and endpoint counts,
not the old monotonic-ratio proof. This is a proof obligation, not a
new theorem being silently activated.

Let a_j(r)>=0 be radial target weights, with sum_j a_j(r)<=1. Let
theta_j(q)=exp(-chi_j(q)*x), with fixed Lipschitz chi_j, and

    G_j=product_i ell_(R_j)(z_i),
    b_j=min(1,theta_j/G_j).

The retained target conditions are assumed to give the genuine fullsource
inverse coefficient U_j between Ulo*n and Uhi*n, with constants of order
one, and facewise retention at least Ret. The COMPLETE precap row is
bounded by

    Ret*Ulo*n*a_j(r) integral_(I_j(r)) theta_j(q)dq/q
       <= Row0_j(w)
       <= Uhi*n*a_j(r) integral_(I_j(r)) theta_j(q)dq/q. (4)

Final rows subtract the nonnegative source-cap loss. For a source at q,
the true target-face sum is at most nG_j, since target radius r>=R_j,
ell_i(r)<=ell_i(R_j) and T_i/(2ell_i(r))<=1. The COMPLETE column of
the SUM of kernels therefore has the valid upper bound

    Col_sum(z)<=n sum_j theta_j(q)
       *integral_(targetband_j intersect [f_j(q)-h,f_j(q)+h])
             a_j(r)dr/r.                              (5)

It is a SUM at a shared source, not the maximum of separate columns.
Terms are zero outside their branch source bands.
Source caps or good-target restrictions can reduce it; any claimed
reduction needs its own proof. Separate degree budgets cannot be reused.

On an actual target the original eligibility argument gives row<=a_j(r)*tau
for each branch, hence total row<=A(w)*tau with A(w)=sum of the active
target weights <=1. Keeping this stronger bound, a valid common degree
upper bound D_n yields the pointwise-deficit certificate

    certified_deficit >= [E-tau integral A(w)d lambda(w)]/D_n,
    D_n>=max(sup Row_sum, sup Col_sum).                 (6)

All finite-a, source-cap and actual-history interfaces would still need
verification before a new construction became an operator theorem.

The entire scalar cover cost at a selected threshold is

    H_tau=(1+A_n/tau)(tau+B_n-gamma*n)
         =B_n-gamma*n+tau+A_n+A_n(B_n-gamma*n)/tau.     (7)

Its algebraic minimum over tau is F(A_n,B_n-gamma*n), but a deficit
proved using one tau cannot simply be moved to that minimizer. For fixed
sigma,gamma, the original scalar identities give

    A_n=n exp(-sigma*L)(1+o(1)), (B_n-gamma*n)/(nL)->1.

Thus the square-root threshold is
tau_sqrt=n sqrt(L) exp(-sigma*L/2)(1+o(1)). This exponent, rather than
the mere use of one global tau, is the source of the old .55534 barrier.

## Mechanism 1: nonlinear matching with one global source ceiling

For one branch and theta independent of q, define the exact log width
W_f(r)=log(q_+(r)/q_-(r)), where q_+/-=f^-1(r+/-h). The row
cost in (4) is proportional to a(r)*W_f(r); its column cost is
theta*n integral a(r)dr/r over an ordinary target width at most 2h.
As h tends to zero,

    W_f(r)=(2h/[q*f'(q)])(1+o(1)),
    row/column local ratio =r/[q*f'(q)].              (8)

This Jacobian must be paid. For f(q)=Kq^alpha and a(r)=r/Rmax,
row is asymptotic to 2nh*theta/(alpha*Rmax), while the complete
column budget is 2nh*theta/Rmax. The corresponding degree-normalized
row factor is min(1,1/alpha). At alpha=1/2 the apparent row factor 2
also doubles the maximum row degree and gives no gain; at alpha=2 the
column dominates and the normalized factor is 1/2. These are two exact
analytic examples, not a parameter scan.

There is also a strict relaxation obstruction. Let A be the UNION of
all target sets. The measure mu_relax=lambda*1_(A^c) satisfies every
actual-row constraint in this relaxation, because it has no targets in
A. Its total deficit is lambda(A). It need not be realizable by an L1
history; it proves only that this relaxed certification cannot force
more than lambda(A). Equivalently E<=D_n*lambda(A) in (6). Thus

    gamma_cert<=lambda(A)/n<=log Rmax.                 (9)

If smax=sup_active f(q)/q, then Rmax<=q1*smax<smax. The true inverse
mean rate is mu(s), with

    mu(s)=(s-1)+(s+1)log((s+1)/(2s)).

The near-full-retention recipe requires chi>mu(smax). With tau_sqrt,
the precap row divided by a(r)*tau is bounded by a constant times
h/sqrt(L)*exp([sigma/2-chi(1-sigma)]L). A positive certificate needs
sigma/2>chi(1-sigma). Therefore its mass-minus-price coefficient is

    gamma_cert+log(1-sigma)
       <log[smax/(1+2mu(smax))]<=c_*.                  (10)

This obstruction is for the stated square-root, near-full-inverse-
retention, pointwise-deficit recipe. Nonlinear matching alone gives no
certified improvement there. Degenerate maps, altered column estimates
or a different threshold are outside this conclusion.

## Mechanism 2: finite bands, local ceilings and radial reweighting

The exact costs are (4)--(6). With fixed Lipschitz chi_j, each matched
interval has width O(h), and its exponent changes by O(xh)=o(1).
Radial weights multiply BOTH the ambient row and the actual-row bound
a_j*tau; they cannot change the sign of their difference.

Assume chi_j(q)>=mu(f_j(q)/q) on contributing cells, as in the
near-full-retention recipe. Since q<1 and mu is increasing on s>1,
every target of radius r has precap row, divided by its total positive
radial weight, at most

    C*n*h*exp(-mu(r)*x+o(1)),                          (11)

with a fixed C for the finite collection. Source-cap deletion only
reduces this. If tau=n exp(-delta*L+o(L)), any fixed radius with
mu(r)>delta/(1-sigma) has negative row-minus-actual-threshold for
large L. Let S_delta>1 be the unique solution

    mu(S_delta)=delta/(1-sigma).

For every fixed eta>0, positive contributions to (6) are eventually
confined to radii below S_delta+eta. Indeed (11)/tau tends uniformly
to zero on r>=S_delta+eta. For the positive-contribution set Omega,

    E-tau integral A d lambda
       <=integral_Omega Row_sum d lambda
       <=D_n*lambda(Omega).

Hence the limiting certificate has gamma_cert<=log S_delta. At the
SQUARE-ROOT threshold delta=sigma/2 this yields

    gamma_cert+log(1-sigma)
       <=log[S_delta/(1+2mu(S_delta))]<=c_*.           (12)

This remains true if different bands use different chi_j, and it pays
their true joint columns. It is not an obstruction for every delta.

A finite counterexample shows why independent band gains cannot be
added at a shared source. Take two target vertices A1,A2 and one source
vertex P, each of ambient mass 1; each target has kernel 1 to P. Use
the SAME positive tau=1/4. Removing only P leaves both targets, satisfies
both actual-row constraints, and has deficit 1. The individual formulas
give (1-tau)=3/4 each, whose naive sum 3/2 is false for that feasible
measure. The complete joint column is 2, and (6) correctly gives only
(2-2tau)/2=3/4. The new script enumerates exactly the eight vertices
of this three-node relaxation to check this example. It is not an
operator counterexample.

The new threshold mechanism bypasses (10)--(12) while retaining ONE
global tau. If tau=n exp(-delta*L), the complete H_tau correction is
bounded by

    [H_tau-(B_n-gamma*n)]/n
       <=exp(-delta*L)+exp(-sigma*L)
                         +L exp(-(sigma-delta)L),     (13)

and the row payment is

    tau/k<=rhi*L^2*exp(-[delta-chi(1-sigma)]L).         (14)

The strict sufficient range is chi(1-sigma)<delta<sigma. It replaces
the old sigma/2 exponent; (12) must not be applied to it. As a CONDITIONAL
scalar example, s=9/2, chi=4/5, sigma=9/20, delta=89/200 has both
gaps in (13)--(14) equal to 1/200. Also log(20/11)<299/500 leaves
gamma=3/2,c=9/10 a scalar allowance 1/500. These facts alone do not
certify a C_n theorem; R088 supplies the separate whole-chain review.
The extra 1/n from B_n<=T_n+price*n+1 also tends to zero.
They already show why .55534 cannot be called a barrier for all
same-global-tau constructions.

## Mechanism 3: different clipping or thinning atypical targets

Fix a local ratio s and retain the same product column certificate
b(z)*G(z)<=Theta, with 0<=b<=1. Pointwise its largest feasible weight is

    b_opt(z)=min(1,Theta/G(z)).                        (15)

The current soft clip is therefore row-optimal under THAT pointwise
column constraint. A hard source cutoff, a different scalar soft profile
or source deletion is dominated by (15) at the same Theta. This does
not prove optimality when the true column can be certified by a different
source-dependent function instead of G.

For the full-target typical set, the baseline first-profile sum divided
by x converges in probability to mu(s), with target variance O(x^2/n).
J_i->1 uniformly, pin/x->0, and the matching transport is o(x).
Thus on a target set of ambient fraction tending to one, the genuine
source-face loss has mean mu(s)x+o(x) and variance O(x). The lower
mean assertion follows from the SAME Fubini identity and two-sided
Chebyshev for the target profile, not from assuming a law for actual
increments. All these assertions are uniform on fixed compact ratio
and radius bands; finitely many bands can be intersected.

If Theta=exp(-chi*x) with chi<mu(s) by a fixed gap, set
d=(mu(s)-chi)/2>0. Under each genuine inverse source-face law,

    E_inv min(1,exp(chi*x-V))
       <=P_inv(V<=(chi+d)x)+exp(-d*x)->0.              (16)

The probability tends to zero by the stated mean/variance bounds.
By (15), no alternative b satisfying the same product column budget
has a larger row. Therefore the current degree certificate, whose
column scale is Theta*n times the matched radial width, has normalized
row and O(n)-deficit coefficient tending to zero on those targets.

Trying to keep only targets whose profile mean is below mu(s) by a
fixed amount also fails in this recipe: their ambient fraction is
O(1/n) by target-profile concentration, so their target mass is o(n).
The relaxation cap (9) then forces gamma_cert=o(1). Fixed finite
band multiplicity does not restore O(n) mass.

The product column bound is not an arbitrary missing factor at central
sources. On a negative source face take its free coordinates in a central
box so that all corresponding target capture fractions equal one. If
ellpin is the fraction for the pinned source coordinate, the FULL target-
face coefficient is exactly (n-1)ellpin+1/2, while nG=n ellpin; their
ratio is 1-O(1/n). The negative pin belongs to the original source cap.
If each target face keeps ambient fraction at least 1-epsA, restricting
to the good set changes this by at most n epsA in coefficient, giving
a constant-factor lower bound when epsA is small. Thus a uniform
pointwise column reduction needs a real additional argument. Deleting
these sources is allowed, but remains dominated by (15) if the same
product column certificate is used.

The unresolved route is consequently precise: replace the product
column relaxation by a proved source-dependent true-column bound,
or exploit further actual-history constraints beyond (6). This batch
does not establish either improvement. The non-root threshold escape
in Mechanism 2 bypasses the old scalar restriction; its full operator
candidate is being reviewed separately.

## Analytic facts used beyond the old [4,5] interval

The smax or S_delta in (10)--(12) need not lie in [4,5]. Their use
requires a new algebraic check on ALL s>1, supplied here:

    mu(1)=0, mu'(s)=1+log((s+1)/(2s))-1/s,
    mu''(s)=1/[s^2(s+1)]>0.

Hence mu is strictly increasing on s>1. Also mu(s)/s tends to
1-log 2>0, so its range is (0,infinity) and S_delta exists uniquely.
For C(s)=log[s/(1+2mu(s))],

    C'(s)=[1+2log((s+1)/(2s))]/[s(1+2mu(s))].

The numerator decreases from 1 to 1-2log 2<0 and has its unique
zero at s_*=1/(2exp(-1/2)-1). Thus C(s) has the SAME unique
global maximum on s>1, namely -log(3-4exp(-1/2)). This is solely
the analytic function used in the square-root-recipe obstruction;
it is not a global maximum over all threshold choices or operators.

## Actual evidence, limitations and stop

verify.py was actually run once for a bounded collection of diagnostic
checks: exact power-map factors, the eight-state shared-column example,
an alternating-exponential/atanh enclosure of c_*, and the conditional
non-root scalar gaps. It executes no old producer or threshold program,
no dimension or decimal parameter grid, and establishes no new C_n bound.
The exact receipt is verification.json.
The final run returned PASS on 13 diagnostic checks.

All three diagnoses are within the explicitly stated pointwise deficit
relaxation. The measure and graph counterexamples need not arise from
actual L1 histories and therefore prove only limits of that relaxation.
Cases with a new true-column theorem, new actual-history coupling,
degenerate/n-dependent matching or other growing structures remain
unresolved here. This finite batch stops after these three mechanisms.

Only this owned R086 directory in the agent's independent clone was
written; no old proof, paper, original source, shared MATLAB session,
commit or push was changed.
