# A069: full Green C2 off the true S3 and qualified derivative finite parts

Let a_1,...,a_m be nonzero vectors spanning R6, with positive spectral weights
already absorbed in the columns. Let c_t(z)=t/[pi(t^2+z^2)], let p_t be the
density of sum_j a_j Z_(j,t) for independent scale-t Cauchy variables, and set
G(x)=integral_0^infinity p_t(x)dt. All statements apply to the original unit-weight
body32 model. No symmetry is needed for the general theorem.

Define the **true** exceptional set

    S3=union_(J subset {1,...,m}, |J|<=3) span(a_j:j in J).

This finite union is closed. It includes every span of any rank-at-most-three
subset, however many dependent columns that subset contains. Equivalently one
can take the union over all unordered three-element subsets, interpreting each
span at its actual rank. No generic-independence assumption is permitted.

## 1. Four independent large parameters give a t^2 density bound

Choose any invertible six-column matrix N. Set

    C=||N^-T||_F^2,        c=(1+C)^-1>0,
    S=sum_j ||a_j||_1,     R_delta=delta/(2S).

The symbol psi(xi)=sum_j|a_j.xi| satisfies psi(xi)>=c||xi||2. Indeed,
||xi||<=sqrt(C)||N^T xi||2<=sqrt(C)psi(xi), and 1+C>=sqrt(C).

Fix x with dist(x,S3)>=delta>0. In any representation x=sum_j a_j z_j, call
an index large when |z_j|>R_delta. If its columns had rank at most three, their
sum would lie in S3, while the remaining sum has Euclidean norm at most
R_delta sum_j||a_j||2<=delta/2. This contradicts the distance hypothesis.
Thus the large columns have rank at least four. A rank-greedy selection gives
four independent large columns. It does not assert that every four large
columns are independent.

For each independent four-element set J let U_J=span(A_J) and
D_J=det(A_J^T A_J)>0. Project the remaining independent-variable sum to U_J^perp.
The quotient is two-dimensional. Its characteristic exponent is precisely the
restriction of psi, since selected columns vanish there; it is still at least
c||eta||. Fourier inversion gives its density q_(t,J) the uniform peak bound

    ||q_(t,J)||infinity
      <= (2pi)^-2 2pi integral_0^infinity r exp(-ct r)dr
       = 1/(2pi c^2 t^2).                                  (1)

On the event that all selected four parameters exceed R_delta in modulus,
their joint density is at most t^4/(pi^4 R_delta^8). The map from those four
parameters to U_J has Jacobian sqrt(D_J). Positive coarea/disintegration of
the remaining variables by their orthogonal projection therefore bounds the
restricted pushforward density by

    t^4 q_(t,J)(proj_(U_J^perp)x)/(pi^4 R_delta^8 sqrt(D_J)).

This is a pushforward of positive measures, not conditioning on a zero-probability
fiber. The union over independent J covers every parameter fiber above the
qualified x. Summing the positive restrictions, then using (1), proves

    p_t(x) <= M_delta t^2,      t>0, dist(x,S3)>=delta,
    M_delta= [sum_(|J|=4,D_J>0) D_J^-1/2]
                  /(2pi^5 c^2 R_delta^8).                    (2)

The coarea bound initially holds almost everywhere. For a given x at distance
at least delta, apply it on a neighborhood using delta'<delta, use continuity
of the fixed-t density, and let delta' increase to delta. This proves the stated
pointwise bound, including the boundary of the distance condition. Fixed-time
smoothness follows from Fourier inversion and coercivity of psi.

For rational columns let d_A be a common positive denominator. The Gram
determinant of the integer matrix d_A A_J is a positive integer, so
D_J>=d_A^-8. Using pi>2 gives a wholly rational, unoptimized majorant

    M_delta <= Mbar_delta
       = binomial(m,4) d_A^4 /(64 c^2 R_delta^8).              (3)

No independent-four-set enumeration is required for this coarse majorant.
For body32, S=192 and d_A=1, so R_delta=delta/384. This is a qualitative
uniformity certificate; no useful covering radius or affordability is claimed.

## 2. Derivative domination uses the complete positive mixture

Reorder the frame as [N,M]. With nu the product law of m-6 standard Cauchy
variables w, the full positive density is

    p_t(x)=integral k_(t,w)(x)nu(dw),
    k_(t,w)=|det N|^-1 product_i=1..6
                 c_t((N^-1(x-tMw))_i).

For directions u,v write b=N^-1u, d=N^-1v and define

    C0=1,
    C1(u)=||b||1,
    C2(u,v)=||b||1||d||1+2 sum_i |b_i d_i|.

The scalar score l'_t(z)=-2z/(t^2+z^2) has modulus at most1/t, and
|l''_t(z)|<=2/t^2. Direct product differentiation therefore gives pointwise

    |D_u k| <= C1(u) k/t,
    |D_uv k| <= C2(u,v) k/t^2.

The global six-factor peak |det N|^-1 pi^-6 t^-6 supplies a fixed-t dominator
for these derivatives under nu. Integration yields the crucial **relative**
bounds for the original complete density:

    |D_u p_t(x)|<=C1(u)p_t(x)/t,
    |D_uv p_t(x)|<=C2(u,v)p_t(x)/t^2.                       (4)

No large-parameter event, rank selector, or indicator function is differentiated.
The event decomposition proves (2) only; the complete mixture independently
proves (4), after which the two inequalities may be combined.

## 3. Uniform cutoff tails and G in C2

Let H_epsilon(x)=integral_epsilon^infinity p_t(x)dt. These functions are smooth
for epsilon>0. At large t, the complete six-factor peak and (4) bound derivatives
of orders0,1,2 by constants times t^-6,t^-7,t^-8, all integrable. On a set at
distance at least delta from S3, (2) and (4) give for k=0,1,2 the exact tail bounds

    0<=G-H_epsilon <= M_delta epsilon^3/3,
    |D_u G-D_u H_epsilon| <= C1(u) M_delta epsilon^2/2,
    |D_uv G-D_uv H_epsilon| <= C2(u,v) M_delta epsilon.       (5)

To justify the derivatives on the left, first integrate the bounded derivatives
on the right and take their uniform limits on a small closed neighborhood whose
distance to S3 is positive. Uniform convergence of the functions and their first
two derivatives, or differentiation under these integrable local dominators,
shows that the limit is C2 and has those derivatives. Every off-S3 point has
such a neighborhood. Therefore

    G is finite and C2 on R6 minus S3.                     (6)

Equation (5) then follows by integration of t^(2-k), and Mbar_delta may replace
M_delta throughout. Mixed derivatives use the displayed C2(u,v), not an
unexplained directional-to-coordinate conversion.

No conclusion is drawn about actual differentiability or Hessian divergence
on S3. A040's positive-moment obstruction there is compatible with (6), and
does not decide signed derivative cancellation on S3. The older scalar theorem
already supplies finiteness/continuity off S2; this theorem uses the stricter
off-S3 region for two derivatives.

## 4. Positive-cutoff copied differences come before constant extraction

Assume a complete finite simplicial fan has been proved for the same psi. For
each cone use its exact rays r_i, q_i=psi(r_i)>0,
alpha=|det(r_1,...,r_6)|/product q_i, b_i=r_i.x/q_i and
J_(i,k)=r_(i,k)/q_i. Any actual symmetry multiplicity is included exactly once.
The full positive-cutoff formula is

    H_epsilon=(2pi)^-6 sum_fullfan alpha DD5(f_epsilon;b_1,...,b_6),
    f_epsilon(b)=atan(b/epsilon).

At epsilon>0, f is analytic in a strip around the real line. Its Cauchy contour
divided-difference formula justifies labelled differentiation, including
zero and repeated node values:

    partial_(b_i) DD5=DD6(base plus b_i),
    partial_(b_i b_j) DD5=DD7(base plus b_i,b_j),       i!=j,
    partial_(b_i b_i) DD5=2 DD7(base plus b_i,b_i).          (7)

Distinct labels do not gain the diagonal factor2 merely because their values
coincide. The node map is linear, so the physical derivatives are

    D_k H_epsilon=(2pi)^-6 sum alpha sum_i J_(i,k) DD6(base+b_i),

    D_kl H_epsilon=(2pi)^-6 sum alpha [
       sum_i 2 J_(i,k)J_(i,l) DD7(base+b_i+b_i)
       +sum_(i<j)(J_(i,k)J_(j,l)+J_(j,k)J_(i,l))
                                      DD7(base+b_i+b_j)].  (8)

These are derivatives of actual cutoff functions, not derivatives of an
artificial sign/zero-jet prescription.

## 5. Complete Laurent cancellation gives the derivative constant formulas

For a fixed tuple of n+1 real nodes, DDn is a finite linear functional of
normalized jets through order n. At nonzero b, the value of f_epsilon has
constant (pi/2)sign(b), and every positive-order b-jet is O(epsilon). At zero,

    f_epsilon^(j)(0)/j! = (-1)^((j-1)/2)/(j epsilon^j), j odd,
                         0,                             j even,

where the even case includes j=0. Thus for n=5,6,7,

    DDn(f_epsilon;tuple)
       =(pi/2) C_n(tuple)
        +sum_(odd 1<=j<=n) L_(n,j)(tuple) epsilon^-j
        +O(epsilon).                                      (9)

C_n is the formal Hermite functional with value sign(b) at nonzero b and all
other data zero (including all data at zero). It is a constant coefficient,
not a claim that sign is differentiable at a collision. For rational nodes,
C_n and every L_(n,j) are rational. Residue-series coefficients supply an exact
implementation without division by a zero difference.

At a fixed x off S3, Section3 proves finite limits for H_epsilon and all
physical first and second derivatives. Insert (9) into the **complete finite
fan** sums (8). Every negative Laurent channel in each complete physical field
must vanish exactly, since otherwise its highest nonzero negative power would
prevent a finite limit. DD6 has channels epsilon^-1,-3,-5; DD7 additionally
has epsilon^-7. These cancellations are full-fan identities, not separate-cone
regularizations, and no uniform cone-by-cone expansion is required.

Consequently G, D_kG and D_klG are obtained from the corresponding DD5/DD6/DD7
expressions by replacing each divided difference with C_5/C_6/C_7 and replacing
(2pi)^-6 with

    (2pi)^-6(pi/2)=1/(128pi^5).                            (10)

The order of operations is essential: differentiate the positive-cutoff formula
first, prove uniform derivative convergence, then extract the complete constant.
One may not first use B052's scalar artificial jets and differentiate their
point-dependent assignment rule. Neither a finite formal constant nor tiny
numerical principal residues replace the off-S3 and complete-fan hypotheses.

## 6. Exact qualification and future derivative contract

For rational inputs, inspect every unordered triple of the actual frame. Use
exact rank-greedy independent columns within that triple, with matrix U of its
actual rank r=1,2,3. Its exact squared residual is

    rho_J=x.x-(U^T x)^T (U^T U)^-1(U^T x).

The minimum over all binomial(m,3) triples is dist(x,S3)^2. A zero residual
rejects this derivative contract; it does not assert derivatives do not exist.
Negative residuals, incomplete enumeration, zero columns or failed rank-six
qualification are invalid inputs. With rho>0, delta0=min(rho,1)/2 satisfies
delta0^2<=rho. A rational ball radius r<delta0 gives the uniform neighborhood
margin delta=delta0-r. A rectangle may use any certified Euclidean radius bound,
for example the sum of its coordinate halfwidths, before using (3)--(5).

`derivative_contract.py` implements this exact qualification and returns explicit
rational remainder budgets through order2 from one fixed six-column basis. Its
low-level Laurent tuple function is explicitly formal. A future genuine oracle
must additionally bind the complete fan proof, actual frame normalization and
multiplicities, exact point and requested fields. It must retain all zero/repeated
nodes and exact signs. Unsupported qualification cannot be silently weakened
to off-S2, and a request crossing S3 cannot use these full-G C2 budgets.

This establishes a future local interface on qualified regions. It does not
prove an affordable whole-face continuous-I certificate, delete S3 regions,
or replace the continuing positive-band g-cover method.
