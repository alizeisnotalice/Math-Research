# 124: rank peeling and one-threshold scalar certificates for the123 moment

Prepared unrun. The only proposed science is eight fixed one-h scalar
optimizations and their33-term exact/Arb certificates, after one primary
proof/source gate. No determinant Monte Carlo, native, old peak/fan/tuple
program, spatial grid, or new agent is run.

## 1. Fixed original moment and complete exact target

Keep the original six-dimensional32 body directions v=(1,+/-1,...,+/-1),
unit weights, common Cauchy time, Lebesgue measure. Use independent
S_i>0 with E exp(-lambda S_i)=exp(-sqrt(lambda)),

 Sigma=sum_i S_i v_i v_i^T, d=det(Sigma)^(-1/2),
 M=E[d^(1/6)]=E[det(Sigma)^(-1/12)].

123 proves a complete all-positive-time/full-face budget
I<= (3/16)^(6/5) M^(6/5)/384, using boundary Minkowski AFTER the exact
positive Gaussian time integral. Thus with the globally UNPROVED target
g0=1/65000, the sufficient strict test is

 M_U^6 < Theta =36893488147419103232/629200415790081024169921875
              =(8/9)^6*24576^5/65000^6.                (1)

The complete rational value is copied from123's actual results, and its
symbolic form is copied from123 inputs/theorem. They are checked equal
using exact arithmetic, never derived from a displayed decimal. The123
MC numbers are not certificate inputs. No conditional Gaussian generator
is used in the original whole-space source transfer.

## 2. Determinant peeling, independence and the exact cut moment

Put Y_i=1/(4S_i), iid Gamma(1/2,1). Fix an independent subset J of k
body columns, U=span(v_j:j in J), W=Uperp, r=6-k. Write C as the
restriction to W of Sigma. All selected columns project to zero, so C
depends ONLY on the other Y_i, and is the actual original-law projected
Gaussian covariance mixture. In orthogonal U/W blocks, let A,B,C be
the covariance blocks of the unselected columns. C is positive definite
on W because the complete projected frame is tight32I_W. Its Schur
complement A-BC^-1 B^T is positive semidefinite. Adding the selected
positive rank-k block therefore gives

 det Sigma >= D_J product_(j in J) S_j * det C,
 D_J=det Gram(v_j:j in J)>0.                           (2)

This also follows by positive Cauchy-Binet terms. For k6 the empty
quotient has det C=1. There is no conditioning on a zero-probability
fiber, no independence assumed between determinants, and no un-restored
omitted-column law: independence holds because the selected columns are
EXACTLY zero in C, leaving disjoint independent input sets.

For h>0 set z=1-exp(-h). The one-dimensional cut moment is exactly

 E[S^(-1/2) 1{Y<=h}]
  = (2/sqrt(pi)) integral_0^h exp(-y)dy=2z/sqrt(pi).

The actual rank-r projected heat density satisfies
Q_r(0)=(4pi)^(-r/2) E det C^(-1/2). Combining this with(2), positive
Tonelli and the independent inputs gives

 E[d product_(j in J)1{Y_j<=h}]
 <= D_J^(-1/2) (2z/sqrt(pi))^k (4pi)^(r/2) Q_r(0).    (3)

For k1..3 the exact033 full-law peak is
Q_r(0)=sqrt(D_J) R_J/(2pi)^r. Hence(3) becomes

       E[d 1_J] <= (2^k R_J/pi^3) z^k.                (4)

All times, constants and measures are those of the original law.
Finiteness follows from the accepted projected peaks; no new peak is
computed. Arbitrarily small/large Y_i are included through the expectation.

## 3. Exact B0..3 and legal coarse B4..6

Let N_h=# {i:Y_i<=h}. Its UNCONDITIONAL law is
Binomial(32,p), p=erf(sqrt(h)). It is not a count under normalized
Dirichlet coordinates. Every1,2,3 distinct body columns are independent.
Summing(4) yields E[d binom(N_h,k)]<=B_k z^k, with

 pi^3 B1=64R1,
 pi^3 B2=4(96R21+240R22+160R23),
 pi^3 B3=8 sum_j n_j R3,j,
 (n_j)=(480,1920,480,640,1440), sum n_j=4960.          (5)

These are the COMPLETE body32 unordered sums, not the no-ones26970
ordered-triple remainder. All R values are read from exact033/036.

There is an already independently checked exact full fan input002:
R0=integral_R6 exp(-psi)=520066141/48769479475200. Its saved producer
and independent fan audit agree. Since p1(0)=R0/(2pi)^6,

 B0=E d=(4pi)^3 p1(0)=R0/pi^3.                        (6)

The science checks R0<27/2^21, the weaker108 full6D coefficient. The
saved fan is read only; no cone or peak computation is rerun.

For k4..6, do NOT replace the number of independent selected subsets
by binom(N_h,k). Any rank-j subspace contains at most2^(j-1) unoriented
body lines: some j-coordinate projection is injective on that subspace,
so it contains at most2^j oriented +/-1 vertices, paired by negation.
Greedy choices consequently give at least

 N(N-1)(N-2)(N-4)...(N-2^(k-2))/k!

independent unordered k-subsets when N exceeds the last threshold.
Define Phi_k(N) as its CEILING in that range and0 when N<=2^(k-2).
The zero branch is essential: multiplying two negative factors below
the threshold would give a false positive count. Since the actual
independent subset count is integral, it is >=Phi_k. For k0 use1 and
for k1..3 use binom(N,k). Every Phi_k is nonnegative and nondecreasing.

For any independent k columns a nonzero k-by-k +/-1 minor is divisible
by2^(k-1) (subtract its first row from all other rows). Cauchy-Binet thus
gives sqrt(D_J)>=2^(k-1). Accepted108 on the TRUE original rank-r
projection in(3) gives, with r=6-k, M_r=2+2r/3,

 E[d Phi_k(N_h)]<= B_k z^k,
 pi^3 B_k <= binom(32,k) 2^(1-4r) M_r^(r/2).          (7)

Use the right side of(7) as the fixed coarse cap, without enumerating
any old k-tuples. Explicitly its pi^3 coefficients are

 k4: binom(32,4)*5/192,
 k5: binom(32,5)*sqrt(8/3)/8,
 k6: 2binom(32,6).                                    (8)

For r0 the zero-dimensional peak is1, not an108 formula at a
nonexistent positive dimension. Dependent subsets are never peeled;
binom(32,k) is used ONLY as an upper bound for the total number of
independent subsets when summing valid bounds.

## 4. Holder scalar model and rigorous rational normalization

For a_k>=0 with a0>0 put W(N)=sum_(k=0)^6 a_k Phi_k(N)>0.
Holder(6,6/5), applied to d^(1/6)=(dW)^(1/6)W^(-1/6), proves

 M^6 <= E[dW] * (E[W(N_h)^(-1/5)])^5.                (9)

Let b0=B0, bk=B_k z^k for k>=1. Then E[dW]<=sum a_k bk.
The exact objective H(a)=sum_(n=0)^32 BinomProb(n;p) W(n)^(-1/5)
is convex, being a positive sum of convex functions of a positive
linear form. The normalization sum a_k bk=1 fixes scaling. At a0=0,
the n0 term is infinite because p<1; no boundary0 workaround is valid.

For computable rational certificates the code encloses bk with Arb512,
rounds DOWN/UP to the fixed2^-240 grid as l_k, beta_k, and uses

 W_x(n)=sum_k x_k Phi_k(n)/beta_k,
 x_k>=0, sum x_k=1, x0>0.                            (10)

This gives E[dW_x]<=sum x_k=1 EXACTLY. Eight separate finite float
optimizations find candidate x; their objective values are not evidence.
Each final x is floored to the2^-44 simplex grid, negative floating
components are discarded ONLY while constructing these nonnegative
weights, and all leftover rational mass is assigned to x0. The result
is rechecked exactly for nonnegativity, strict x0 and sum1. No random
sample, event, state or probability is dropped or renormalized.

The full33 terms of(10) are then certified with positive rational W(n).
Arb log/exp gives a positive fifth-root reciprocal interval v_n, rounded
outward to2^-192; exact comparisons v_L^5<=1/W(n)<=v_U^5 certify each.
p=erf(sqrt h) is enclosed and rounded outward to2^-192. Since W is
nondecreasing, E W(N)^(-1/5) is nonincreasing in p (couple the32
Bernoulli variables with common uniforms). Therefore the exact rational
sums at p_U,v_L and p_L,v_U give H_L,H_U. The formal success test is

                   H_U^5 < Theta.                   (11)

No displayed threshold, float objective, MC estimate or statistical
confidence bound enters this comparison. A positive exact residual
certifies a true upper bound for M, and through123 a COMPLETE I bound.
Global g0 remains a separate unproved requirement for C6.

## 5. Optional rigorous diagnosis of the fixed scalar family

A float non-pass cannot exclude its exact optimum. To distinguish this,
the same33 saved root intervals give lower endpoints g_j,L for the
gradient of H(x). The binomial probability at each n has its maximum
on[p_L,p_U] at an endpoint or at p=n/32; these are exact rational checks.
Convexity and degree-1/5 homogeneity give the support-plane lower bound

 min_(simplex) H >= (6/5)H(x)+min_j partial_j H(x)
                >= Ldual:=(6/5)H_L+min_j g_j,L.       (12)

This is a finite exact lower certificate on the EXPLICIT scalar budget
optimum, never a lower bound on actual M. Define gamma=min_k l_k/beta_k.
For normalization with the exact bk instead of beta, a·beta<=1/gamma;
scaling shows its optimal H^5 is at least

                  gamma max(0,Ldual)^5.               (13)

If(13)>Theta, this fixed h and fixed B-table Holder relaxation is
rigorously obstructed. Otherwise report unseparated; do not infer an
obstruction from optimization status. Tighter actual mixed moments,
different h, multiple h, or the original M are not excluded.

## 6. Finite batch, independent saved checking and stop

The ONLY h values, in order, are
1/4096,1/2048,1/1024,1/512,1/256,1/128,1/64,1/32.
Each has one fixed SLSQP start x=(1,0,...,0), at most500 iterations
and2500 objective/analytic-gradient evaluations; at most20000 objective
evaluations total. No restarts, new h, precision refresh or MC. Solver
nonconvergence/budget exhaustion merely returns its best candidate for
exact certification. A genuine nonfinite/math/runtime error stops the
batch, preserving completed per-h checkpoints. At most one CPU, no
wall/CPU cap, no child science workers. An observer records complete exit.

All eight fixed candidates and certificates are saved, then stop. The
separate saved-only checker uses Fraction,128-term alternating exp/erf
series and384-bit rational square-root bounds with a pinned768-bit pi
interval to validate the beta/p enclosures. It checks every fifth-root,
33-term sum, gradient support and rational comparison without rerunning
optimization, Arb, old fan/peaks or determinant samples.

No multiple-h count enumeration is run or generated. Its potential
9-bin multinomial reduction remains future work requiring new input
and release. GeneralL1, all positive common cube radii and the same
whole-space source/generator interface from123/001/027 remain unchanged.
No complete g, unconditional positive S or new C6 is claimed here.
