# 125: a genuine joint nested-count moment, with a complete positive tail

Prepared unrun. One new deterministic optimization and at most one
rationalized candidate's complete interval certificate require a primary
gate. No735471-state data, optimization or integral has run in preparation.
The sealed124 files and all eight single-h certificates are read only.

## 1. Same original moment, different legitimate scalar model

Keep the original body32 unit-weight common-time Cauchy model. As in
123/124, Sigma=sum S_i v_i v_i^T, d=det(Sigma)^(-1/2),
Y_i=1/(4S_i) are iid Gamma(1/2,1), and M=E d^(1/6).
The h thresholds concern these same auxiliary Gaussian mixture inputs;
they are not different semigroup times or a changed generator/source.

For the SAME ordered eight h values from124, let N_j=# {Y_i<=h_j}.
Read124's complete rational beta_jk and CDF endpoints p_j,L,p_j,U,
and its exact Phi_k(n), k0..6,n0..32. For k>=1,

 E[d Phi_k(N_j)]<=beta_jk,  E d<=beta0.

The beta0 values are checked equal across all eight saved files. In
particular B0 is the accepted exact002 full-fan input, not a new peak
or determinant sample. All rank-peeling/independent-subset proofs and
zero branches from124 remain in force; no single-h optimization is rerun.

For49 nonnegative weights with x0>0 and sum x=1 define

 W(N)=x0/beta0 + sum_(j=1)^8 sum_(k=1)^6
                        x_jk Phi_k(N_j)/beta_jk.

Linearity alone gives E[dW]<=1. No independence of the N_j is assumed.
Holder yields M^6 <= H^5, H=E W(N)^(-1/5). This model can use information
from several h simultaneously, so124's fixed single-h obstruction does
not apply to it.

The exact target is the complete123 rational
Theta=36893488147419103232/629200415790081024169921875.
Only H_U^5<Theta can pass this conditional moment/I gate. MC values
are not inputs. The global cube g0=1/65000 remains separately unproved.

## 2. The true joint law and the correct monotone endpoint order

For any increasing CDF vector p, put q1=p1,qj=pj-p_(j-1),q9=1-p8.
The raw category counts c1..c9 have the ONE multinomial law

 Pr(c)=32!/(product c_i!) product q_i^c_i, sum c_i=32,
 N_j=c1+...+cj.

The eight binomial counts are never multiplied as independent variables.
All lower-endpoint cumulative p values form one increasing vector,
and all upper endpoints form another; q differences are taken inside
each complete vector. Category probabilities are not independently
rounded in an order that would describe a different law.

Couple the32 Bernoulli counts using the same independent uniforms:
N_j(p)=sum_i1{U_i<=p_j}. Since every Phi is nondecreasing and x>=0,
W increases coordinatewise in N and W^(-1/5) decreases. Thus the exact
full expectation obeys

                  H(p_U)<=H(p_true)<=H(p_L).          (1)

The tail-relaxed upper function below need not itself be monotone at
its K16/17 seam. Equation(1) is applied to the EXACT H first; only then
is its fixed-endpoint expectation relaxed. No monotonicity of a truncated
or tail-relaxed integrand is required or asserted.

## 3. Complete prefix and the explicit tail

Set K=N_8 and freeze Kmax=16. Enumerate precisely
0<=N1<=...<=N8<=16, equivalently nonnegative c1..c8 with sum<=16.
The exact number is sum_(K=0)^16 binom(K+7,7)=binom(24,8)=735471.
For each state set c9=32-K and use its actual multinomial coefficient.
The positive prefix has its original subprobability; it is NOT normalized.

For K>16 positivity gives the pointwise lower bound

 W(N)>=Wtail(K):=x0/beta0+sum_(k=1)^6 x_8k Phi_k(K)/beta_8k>0.

Since K has Binomial(32,p8), all16 tail values K17..32 are retained.
With Pprefix(p) denoting the full positive prefix sum, a legal complete
upper and lower are

 H_L=Pprefix(p_U) <=H(p_true),
 H(p_true)<=H_U=Pprefix(p_L)
             +sum_(K=17)^32 BinomProb(K;32,p8,L) Wtail(K)^(-1/5). (2)

There is no discarded tail, missing label pattern, conditional
renormalization or independence assumption. The lower bound is on
the specified Holder objective H, not on actual M or I.

## 4. One fixed finite optimization, followed by one exact candidate

Build735471 prefix rows and16 tail rows, each with49 features, using
the midpoint CDF vector and its joint multinomial law. The tail features
retain only the common intercept and the six last-h features. This is
the complete upper expression of(2) at that midpoint, not a Monte Carlo
or a re-evaluation of an old33-state optimization. Its objective and
analytic gradient are new and convex in positive weights.

Use one SLSQP start: embed the exact124 best h=1/256 simplex weights
at its corresponding h block, with the same common x0. There is one
solver invocation, at most1000 objective/analytic-gradient calls and
500 iterations, with x0>=2^-24. No restart or second candidate search.
Nonconvergence/cap exhaustion retains the best feasible candidate only;
it does not prove an optimum or an obstruction. Float mass is checked
but never renormalized.

Rationalize the candidate to the2^-44 simplex grid: floor all nonnegative
non-intercept weights and assign the remaining exact mass to x0. Verify
all49 weights nonnegative, x0>0 and exact sum1. Release all large float
arrays before interval certification. The float matrix is only an
in-memory derived array (735487*49 float64); it is not an archive input.

## 5. Compact exact arithmetic for the complete saved certificate

Use fixed precision192 and scale S=2^192. The saved p endpoints already
lie on this grid. Convert each whole CDF endpoint vector to positive
integer category numerators q_i, with sum S. This describes two honest
multinomial laws. Precompute q_i^c for c0..32, not floating probabilities.

For the rationalized x, compute the264 group sums
T_j(n)=sum_k x_jk Phi_k(n)/beta_jk and the common intercept as exact
Fractions once. Enclose them outward on the absolute1/S grid, saving
their integer lower/upper tables. Each state then uses eight integer
additions and the common intercept to get wL/S<=W<=wU/S.

Arb192 encloses one reciprocal fifth root of this W interval. Round
outward to integers vL,vU on the1/S grid and verify EXACTLY

              vL^5 wU <= S^6 <= vU^5 wL.             (3)

Therefore vL/S<=W^(-1/5)<=vU/S. Identical integer W intervals may reuse
the same already-checked root through a bounded65536-entry cache; no
approximate root equality or stale candidate cache is used.

Every32-input multinomial probability has the COMMON denominator S^32:
its numerator is the integer multinomial coefficient times product
q_i^c_i. Multiplying by vL or vU gives denominator S^33=2^6336.
Thus all prefix lower/upper numerators can be accumulated as exact
positive integers. Lower uses the complete p_U law; upper uses p_L.
The16 Binomial tail upper numerators have this same denominator.
No floating sum or cumulative rounding error enters the final certificate.

Save every state's pair(vL,vU), ordered lexicographically by the complete
nondecreasing count tuple, in deterministic gzip chunks of8192 states.
Chunk headers record coverage indices, first/last tuples, state hash,
per-K counts and EXACT lower/upper weighted sums. No single file may
exceed the project100MB limit. Checkpoints are atomic at completed
chunks; the raw float feature matrix is not kept.

The final positive rational H_L,H_U give the exact tests and certificate

 M^6<=H_U^5, and I^5<=(3/16)^6 H_U^5/384^5 via123.

The single saved candidate is marked strict pass only if H_U^5<Theta.
Otherwise report its saved interval/margins as unseparated or
noncompetitive. No family-wide, multi-h-optimality or actualM exclusion
is inferred, even if this one candidate's H_L exceeds the target.

The full saved124 h_4 certificate upper and all eight exact dual lower
bounds are separate read-only comparison inputs. If the new H_U^5 is
below their minimum dual lower, this proves a strict improvement over
those eight fixed single-h cap-budget families. It is still no lower
bound on actualM and no C6 conclusion. No comparison is built from the
displayed0.068475 value, and no pending stronger moment table is adopted.

## 6. Saved-only verification, resources and stop

The independent checker uses no Arb, SciPy, objective or determinant
evaluation. It rechecks source weights/caps, exact group tables, every
inequality(3), lexicographic full coverage, all true multinomial/binomial
integer numerators and the exact positive sums. It therefore verifies
the COMPLETE candidate interval rather than trusting chunk decimals.
The 124 CDF/beta validity remains pinned to its passed pure series/root
certificates; no old peak or single-h optimization is rerun.

One new candidate, one finite735471-state prefix plus16-value tail,
fixedK16, fixedeight h, fixed192bit bounds; no precision/K/h escalation,
MC, blind enumeration of76904685 total states, native, new agent,
MATLAB, lower operator work, stage or push. At most one science CPU,
BLAS threads1, nice10, no wall cap;121's15 native jobs have priority.
First genuine error stops, preserving committed chunks and candidate;
no automatic restart. After this fixed candidate and its saved check,
seal and stop for the primary's next decision.

The unchanged original whole-space positive-source mechanism still
handles generalL1 and all positive common cube radii. No conditional
Gaussian generator replaces it. Passing this M/wholeI slot would leave
globalg0 and final same-model S/C6 acceptance for separate certification.
