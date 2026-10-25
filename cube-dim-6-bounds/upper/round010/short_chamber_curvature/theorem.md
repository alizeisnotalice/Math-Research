# The body32 short potential is not convex on the strict ordered face

## Statement

Let V={(1,s2,...,s6): sj in {-1,1}}, let the 32 variables Z_v be
independent centered Cauchy variables with density

    c_t(u)=t/[pi(t^2+u^2)],

and let p_t be the density of sum_(v in V) v Z_v. For any fixed T>0 put

    K_T(x)=integral_0^T p_t(x) dt,
    D={x: x1=1>x2>x3>x4>x5>x6>0}.

Then K_T is C2 on D and **is not convex on D**. More precisely, for
some y in D and some i in {2,3,4,5},

    D^2_(e_i-e_(i+1)) K_T(y)<0.                       (A)

Consequently an entire sufficiently short symmetric segment about y is
in D and has a strictly negative midpoint Jensen gap. This includes
T=4/125. The proof gives existence of an ordinary interior negative
curvature point; it does not numerically locate that point. No new
G/H/p point evaluations are used.

Write S_r for the finite union of the spans of at most r columns of V.
We prove all regularity needed below; in particular, we do not assume
that G or K_T is differentiable everywhere off S2.

## 1. Exact density representation

For any six independent columns B and the 26 remaining columns A_N,
independence and a nonsingular change of variables give, at every x,

    p_t(x)=1/|det B| integral_(R^26)
      product_(j=1)^6 c_t((B^-1(x-A_N u))_j)
      product_(j in N) c_t(u_j) du.                  (1)

The right side is a continuous density: the six factors are bounded by
(pi t)^-6 and the other 26 factors integrate to one. Hence this is a
pointwise formula, with its full determinant factor. It also supplies
fixed-time derivatives by dominated differentiation. The direction set
is preserved up to signs by every signed coordinate permutation; Cauchy
symmetry and independence therefore give the same invariance for p_t
and K_T.

## 2. Local uniform bounds away from S3

Let E be a compact subset of R6 minus S3, let
eta=dist(E,S3)>0, and take delta=eta/(64 sqrt(6)). For a parameter vector
a with A a=x in E, the columns with |a_j|>=delta must span a space of
dimension at least four. Indeed, if they span a space L of dimension at
most three, then L is one of the subspaces included in S3 and

    dist(x,L) <= sum_(|a_j|<delta) |a_j| ||v_j||
              <32 delta sqrt(6)=eta/2,

which contradicts the definition of eta. Thus there are four independent
columns I with all four parameters at least delta in absolute value.

There are finitely many such I. For each I extend it to a basis B_I from
V. On the corresponding event in (1), the four indicated Cauchy factors
are at most t/(pi delta^2), and the two other basis factors are at most
1/(pi t). Discarding the event after this bound and integrating the 26
remaining probability densities gives at most

    t^2/[|det B_I| pi^6 delta^8].

The events cover every fiber over E. Summing their positive contributions
is legitimate either directly by the fiber change of variables or first
as a positive density inequality and then by continuity. Since the
nonzero integer determinants have absolute value at least one,

    p_t(x) <= binom(32,4) t^2/(pi^6 delta^8)  (x in E).       (2)

This argument does not substitute a conditional probability on a
zero-probability fiber.

For a fixed basis B0 and direction h put b=B0^-1 h. The elementary
Cauchy bounds

    |c_t'(u)| <= t^-1 c_t(u),
    |c_t''(u)| <= 2 t^-2 c_t(u)

applied to (1) show

    |D_h p_t(x)| <= ||b||_1 t^-1 p_t(x),
    |D_h^2 p_t(x)| <= 2 ||b||_1^2 t^-2 p_t(x).       (3)

The same product calculation bounds mixed derivatives. The constants
are independent of x and t. Combining (2) and (3) supplies locally
uniform integrable bounds O(t^2), O(t), and O(1) as t decreases to zero.
For t in a compact subinterval of (0,T], ordinary dominated
differentiation in (1) suffices. Therefore

    K_T is finite and C2 on R6 minus S3.              (4)

Only this region is used for derivatives in the rest of the proof.

## 3. Coordinate-value count removes the wall regularity issue

If x is in a span of r<=3 sign vectors, each coordinate of x is one of
the at most 2^r values sum_(j=1)^r epsilon_j a_j. Opposite sign patterns
give opposite values, so there are at most 2^(r-1)<=4 distinct coordinate
absolute values (the zero vector causes no exception).

Every x in D has six distinct positive coordinate values and hence is
off S3. If M_i x replaces x_i,x_(i+1) by their average, for any
i=2,3,4,5, then M_i x has exactly five distinct positive coordinate
values: the average is strictly between its two original values and
strictly separated from all other coordinates. Thus M_i x is also off
S3. Every point on the segment from M_i x to x, except the wall endpoint,
is in D. All of this segment has the C2 regularity in (4).

No genericity claim, exceptional-ray argument, or extension across S2
is needed.

## 4. A rank-two point has infinite short potential

Let

    v=(1,1,1,1,1,1), w=(1,-1,-1,-1,-1,-1),
    c=(1,1/2,1/2,1/2,1/2,1/2)=(3/4)v+(1/4)w.

Complete v,w by the four columns obtained from v by flipping respectively
coordinates 2,3,4,5. They form a basis B. The exact verification file
checks its determinant and B^-1 c=(3/4,1/4,0,0,0,0).
Let C=max_j sum_k |(B^-1 A_N)_(j,k)|, a finite positive number, and
t0=min(1,1/C). Restrict the 26 variables u in (1) to |u_k|<=t.
Then |(B^-1 A_N u)_j|<=Ct. For 0<t<=t0, the first two solved parameters
have absolute value at most 2, and the remaining four at most Ct.
Consequently their six density factors have product at least

    [t/(5 pi)]^2 [1/(pi(1+C^2)t)]^4
      = t^-2/[25 pi^6 (1+C^2)^4].

The integral of each remaining Cauchy density over [-t,t] is exactly
1/2. Formula (1) therefore proves the explicit pointwise estimate

    p_t(c) >= 2^-26 t^-2 /
                 [25 |det B| pi^6 (1+C^2)^4]  (0<t<=t0).     (5)

It follows that K_T(c)=infinity for every T>0. For any x_n tending to c,
fixed-time continuity and Fatou's lemma give

    liminf_n K_T(x_n)>=integral_0^T p_t(c)dt=infinity.         (6)

Thus approaching this point inside D forces divergence. This is an
independent positive proof, consistent with the stronger quantitative
A009 S2 barrier but not dependent on that barrier's numerical constants.

## 5. A strictly negative adjacent-coordinate curvature must occur

Use the exact rational vector

    z=(0,51,22,2,-24,-51)/200,
    x(s)=c+s z,   0<s<=1,
    F(s)=K_T(x(s)).

These points all belong to D, sum_(j=2)^6 z_j=0, and the four initial
partial sums of z2,...,z6 are (51,73,75,51)/200, all positive.
By (4), F is C2 on (0,1]. By (6), F(s) tends to infinity as s decreases
to zero, whereas F(1) is finite. The mean value theorem therefore gives
some s_* in (0,1) with F'(s_*)<0.

Set x=x(s_*) and g_j=partial_j K_T(x). Summation by parts gives the
exact identity

    F'(s_*)=[51(g2-g3)+73(g3-g4)+75(g4-g5)+51(g5-g6)]/200.

At least one g_i-g_(i+1) is consequently strictly negative. For that
i put m=M_i x, a=(x_i-x_(i+1))/2>0, d=e_i-e_(i+1), and

    q(r)=K_T(m+r d),   0<=r<=a.

Section 3 proves C2 regularity along this entire segment and in a
neighborhood of m. Swapping coordinates i and i+1 fixes m and changes
r to -r, so permutation symmetry makes q even near zero. Hence
q'(0)=0, while q'(a)=g_i-g_(i+1)<0. Applying the mean value theorem to
q' gives r_* in (0,a) with q''(r_*)<0. The point y=m+r_*d lies strictly
inside D and proves (A).

By continuity of q'', a sufficiently short interval about r_* has
q''<0. The midpoint Jensen gap there is strictly negative, as follows
also from the integral second-difference formula. Thus this is a strict
analytic disproof of chamber convexity, not merely a failure to prove it.

Equivalently, had K_T been convex on D, symmetry and differentiability
at each single-equality wall would force g2>=g3>=...>=g6. The same
summation-by-parts identity would force F to be nondecreasing, in direct
conflict with its left-end divergence and finite interior values.

## Scope

The result excludes the proposed single stationary-point/tangent
certificate based on convexity of the whole strict ordered chamber.
It leaves open convex subregions, different lower-bound constructions,
and the required global essential lower bound m_T. It proves no new C6
bound. In particular, the positive B009 Jensen gap at its three points
is entirely compatible with this theorem: curvature need not have one
sign throughout the chamber.
