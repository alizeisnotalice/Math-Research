# Independent second review: cube-cap operator bound

Reviewer: fresh native Sol; date 2026-10-27. Frozen source commit:
`79f8acb20c418e60c5196a3c830103253e34c900`.

**Decision: ACCEPT the analytic operator theorem below.** No unresolved core gap
was found. The R028 decimal is a corollary conditional on a separately verified
enclosure of the scalar expression at its registered rational input. This review
does not assert global optimality, sharpness, novelty, or a lower bound.

This review uses the definitions and the proof itself, not ledger entries,
historical PASS labels, or author/auditor conclusions as mathematical premises.
The full derivation below is the review's mathematical evidence. The source
receipts were read only to identify scope and provenance. No research extension,
optimization, scan, new candidate, or lower-bound computation was performed.

## 1. Definitions and theorem accepted

Write B(y,r)=y+(-r,r)^n, with n a positive integer and r>0. For s>0,
sB(y,r)=B(y,sr). All cubes have sides parallel to the fixed coordinate axes and
use one common radius in every coordinate. Define

    M_n f(x) = sup_{r>0} (2r)^(-n) integral_{B(x,r)} |f(y)| dy,
    C_n = sup_{0<||f||_1<infinity, lambda>0}
            lambda |{M_n f>lambda}| / ||f||_1.

Lebesgue measure and general real or complex L1 functions are intended. Thus C_n
is the full centered cube weak-(1,1) constant, with every positive common radius.

For every t>1, define

    B_n(t)=1-((1+1/t)/2)^n+n log t
            -(n/2) integral_1^t
              [1+1/(2t)-1/(2rho)]^(n-1) d rho/rho.

The accepted result is

    B_n(t)>0,
    C_n <= [(1+1/t)^(n/2)+sqrt(B_n(t))]^2.                 (T)

The proof first establishes, for each finite a>2t and tau>0,

    b=a(1+1/t), alpha=a/(a+1), c=(b-1)/(a+1)>1,
    P=((a+1)/a)^n, A_a=((b+1)/a)^n,
    ell_a(rho)=(1+c-alpha/rho)/2,
    J_a=1-ell_a(1)^n+n log t
          -(n/2) integral_1^t ell_a(rho)^(n-1) d rho/rho,
    E_a=a^(-n)+P J_a,
    C_n <= (1+A_a/tau)(tau+E_a).                          (F)

All quantities in (F) are finite. Notice b is finite and the covering factor
really is b+1; replacing b+1 by b at finite a would be unjustified.

## 2. Finite selection: source increments and rejection prefix

Take any finite supplied family Q_l=B(y_l,r_l). Initially every member is
eligible. At each step choose an eligible remaining cube of greatest radius;
arbitrary tie breaking is allowed. If none is eligible, stop. For the jth chosen
cube B_j, freeze

    I_j=B_j minus union_{i<j} B_i,
    v_j=|I_j|, w_j=v_j/|aB_j|=v_j/(2ar_j)^n.

Eligibility of an unchosen cube with center y means

    sum_{i already chosen} w_i 1_{bB_i}(y) <= tau.         (S)

The frozen w_i are nonnegative. Thus a remaining cube can only lose eligibility,
never regain it. Each I_j is measurable, the I_j are pairwise disjoint, and
their union equals the union of the chosen small cubes. The chosen radii are
nonincreasing: each subsequent eligible set is a subset of the previous eligible
set after removing the chosen member, so its maximum cannot grow. Zero-volume
increments cause no problem and can be retained with weight zero.

Here is the precise prefix argument needed for rejected cubes. If a supplied
cube Q=B(y,R) is unchosen, it must become ineligible at some step; otherwise the
finite procedure could not terminate without choosing it. Let m be the first
step after which it is ineligible. Immediately before every step j<=m, Q was
still eligible and unchosen. The greatest-eligible-radius rule therefore gives
r_j>=R for every j<=m. Furthermore

    sum_{j<=m} w_j 1_{bB_j}(y)>tau.

For each j in this sum and every x in Q, the triangle inequality gives

    ||x-y_j||_infinity < R+b r_j <= (b+1)r_j.

Thus Q is contained in (b+1)B_j for every such j, and the final enlarged load

    L(x)=sum_{j chosen} w_j 1_{(b+1)B_j}(x)

exceeds tau on Q. Consequently, if U is the union of all supplied cubes and V
is the union of the selected small cubes,

    U subset V union {L>tau}.

Integration and |V|=sum_j v_j now give

    |U| <= sum_j v_j + tau^(-1) integral L
         = (1+A_a/tau) sum_j v_j,                         (C)

since |(b+1)B_j|/|aB_j|=((b+1)/a)^n=A_a.
No assumption about disjoint selected cubes is used. Disjointness belongs to
the frozen small increments, whose exact definition is essential.

## 3. The actual selected-star load and its earlier terms

Define W(x)=sum_j w_j 1_{aB_j}(x). If W(x)=0 there is nothing to prove. Otherwise
among its positive contributing terms choose one of least radius, and among
those choose the last selected index k. All other positive contributors have
index j<k: a larger radius cannot have been selected later, while an equal
radius is handled by the last-index rule. This establishes eligibility payment
for every possible other term, not merely the ones with strictly larger radius.

Translate x to 0 and divide every coordinate and radius by r_k. The weights
are unchanged because both numerator and denominator scale by the same r_k^n.
Thus r_k=1, ||y_k||_infinity<a, and all other contributing radii are >=1.

A previous contributing term is called captured if y_k belongs to bB_j.
Their total weight is at most tau by the actual condition (S) at the time B_k
was selected. They form a sub-sum of that condition; other previous terms can
only increase the full nonnegative selection load. The self term is separate:

    w_k <= |B_k|/|aB_k|=a^(-n).                           (self)

For an uncaptured previous contributor,

    ||y_k-y_j||_infinity>=b r_j,
    ||y_j||_infinity<a r_j, ||y_k||_infinity<a.

If r_j>=t then

    ||y_k-y_j||_infinity<a+a r_j
       <= a r_j+(a/t)r_j=b r_j,

which is impossible. Every unpaid term therefore satisfies

    1<=r_j<t.                                            (R)

This step pays all large-radius terms through selection. It does not replace
the eligibility condition by any geometric model or omit equal-radius terms.

## 4. A cap domain that contains every unpaid increment

For each y in an unpaid I_j subset B_j, set

    z=y/(a+1), u=y_k/(a+1), m(z)=max(1,||z||_infinity).

Then ||u||_infinity<alpha. Also

    ||y||_infinity <= ||y_j||_infinity+||y-y_j||_infinity
                    <(a+1)r_j,

so ||z||_infinity<r_j<t and m(z)<=r_j. Uncapturedness gives

    ||y_k-y||_infinity
      >=||y_k-y_j||_infinity-||y-y_j||_infinity
      >(b-1)r_j.

Hence ||u-z||_infinity>c r_j>=c m(z). Each such increment maps inside

    D(u)={z: ||z||_infinity<t, ||u-z||_infinity>=c m(z)}.

Enlarging a strict inequality to a non-strict one gives a valid upper domain.
For the weight, r_j^(-n)<=m(z)^(-n). Thus disjointness of I_j and the common
change of variables dy=(a+1)^n dz prove

    sum_{unpaid j} w_j
      = (2a)^(-n) sum_j integral_{I_j} r_j^(-n) dy
      <= ((a+1)/a)^n 2^(-n) integral_{D(u)} m(z)^(-n) dz
      = P H(u).                                          (H)

The Jacobian and 2^(-n) normalization are both explicit. This is an upper
bound for actual source increments, regardless of overlaps among the cubes or
their stars. Combining it with selection payment and (self) gives

    W(x)<=tau+a^(-n)+P H(u).                              (W)

## 5. Coordinate monotonicity of the cap integral

Simultaneous reflection of z and u in any coordinate preserves the integral,
so H depends on |u_i|. It suffices to compare 0<=u_i<=alpha<1.

On the inner cube ||z||_infinity<1, m=1. The captured set is its intersection
with u+(-c,c)^n. The overlap in coordinate i has length

    L_i(1)=min(2,1+c-u_i).

Indeed u_i+c>1, so the upper endpoint is 1; the lower endpoint is max(-1,u_i-c).
These lengths are positive and nonincreasing in u_i. The normalized uncaptured
inner contribution 1-2^(-n) product_i L_i(1) is therefore coordinatewise
nondecreasing.

For 1<=rho<t partition the shell by the unique coordinate attaining maximal
absolute value. Tied-maximal sets have n-dimensional measure zero. On each
face z_i=+rho or -rho the remaining coordinates range in [-rho,rho]; the map
(rho,z_other) to z has absolute determinant 1. The captured coordinate overlap
length in any other coordinate is

    L_l(rho)=min(2rho,(1+c)rho-u_l).

The positive face always satisfies |rho-u_i|<c rho because 0<=u_i<1<=rho and
c>1. The negative face satisfies |-rho-u_i|<c rho exactly when
u_i<(c-1)rho. Thus the captured shell density is

    S(rho,u)=sum_i [1+1_{u_i<(c-1)rho}]
                         product_{l!=i} L_l(rho).

For any fixed coordinate u_j, every summand is nonincreasing: the jth summand
contains the nonincreasing indicator count, while every other summand contains
the nonincreasing overlap length L_j. All other factors are nonnegative.
Therefore S is coordinatewise nonincreasing, and

    H(u)=1-2^(-n) product_i L_i(1)
       +2^(-n) integral_1^t rho^(-n)
            [2n(2rho)^(n-1)-S(rho,u)] d rho

is coordinatewise nondecreasing. This is a face parameterization with a fixed
radial weight, so it does not differentiate a moving intersection volume.
Threshold equalities affect at most finitely many radial values and hence no
integral. The same reasoning is valid when n=1, with empty products equal to 1.
We conclude

    H(u)<=H(alpha,...,alpha).                             (corner)

The closed comparison corner need not be realized by any selected center.

## 6. Exact finite-a corner integral

With b=a(1+1/t), a>2t,

    c-1=(a/t-2)/(a+1)>0,
    alpha/(c-1)=at/(a-2t)>t.

It follows that alpha>(c-1)rho for all 1<=rho<=t. At the comparison corner
exactly one face in each pair is captured and

    L(rho)=(1+c)rho-alpha,
    L(rho)/(2rho)=ell_a(rho), 0<ell_a(rho)<1.

The inner cube term is 1-ell_a(1)^n. The normalized full shell is

    2^(-n)rho^(-n) 2n(2rho)^(n-1) d rho=n d rho/rho.

The n captured faces give

    2^(-n)rho^(-n) n L(rho)^(n-1) d rho
       =(n/2)ell_a(rho)^(n-1) d rho/rho.

Therefore H at the corner equals precisely J_a, including the factor n/2.
Equations (W) and (corner) yield the pointwise inequality

    W(x)<=tau+E_a.                                       (L)

Since ell_a(1)<1, the inner term is positive. The shell integrand is
n[1-ell_a^(n-1)/2]/rho>=n/(2rho), also for n=1. Thus J_a>0.

## 7. From finite geometry to the full maximal operator

Let f be any L1 function and lambda>0. For each fixed R>0 the cube average is
continuous in its center: the change after translating by h is bounded by
(2R)^(-n)||f(.+h)-f||_1, which tends to zero. Consequently
E={M_n f>lambda}, the union of all fixed-radius strict superlevel sets, is open.

For each x in a compact K subset E choose one actual positive radius R_x whose
centered cube average exceeds lambda. Supply the small cube B(x,R_x/a), whose
a-star is exactly that bad cube. These small open cubes cover K; choose a finite
subcover. Apply the finite selection construction to that family. Every chosen
star retains its original average greater than lambda, so, term by term,

    lambda v_j <= w_j integral_{aB_j}|f|.

For a zero increment both sides vanish. Summing the finite family and using
(L) proves

    lambda sum_j v_j <= integral |f|W <=(tau+E_a)||f||_1.

Since K is contained in the union of all supplied cubes, (C) now gives

    lambda |K| <= (1+A_a/tau)(tau+E_a)||f||_1.

Inner regularity takes the supremum over compact K subset E and gives (F).
This argument permits arbitrarily small or large witness radii, finitely many
different radii in each cover, and completely general L1 input. There is no
product, finite-source, compact-support, smoothness, or minimum-radius assumption
on f. Open and closed cube averages agree, since each boundary is null.

## 8. Scalar limit and the square-root bound

Keep n,t,tau fixed. At every finite a>2t, (F) is already valid for the same full
operator and for all f,lambda. As a tends to infinity,

    alpha ->1, c ->1+1/t, P ->1, a^(-n)->0,
    A_a ->A=(1+1/t)^n,
    ell_a(rho) ->ell(rho)=1+1/(2t)-1/(2rho).

On the fixed compact interval [1,t], these affine functions of 1/rho converge
uniformly. Their nonnegative integer powers do so as well, proving J_a->B_n(t).
Therefore a scalar limit of (F) yields

    C_n <= (1+A/tau)(tau+B_n(t)).                         (SL)

No selected family, increment, f, or geometric domain is claimed to converge.
It is solely a limit of constants already valid for the entire operator.

For t>1, ell(1)=(1+1/t)/2<1, while 0<ell(rho)<=1 on [1,t]. The inner term
1-ell(1)^n is positive and the shell term is at least (n/2)log t>0. In particular
B_n(t)>0. Expanding the right side of (SL) gives

    A+B_n(t)+tau+A B_n(t)/tau.

The arithmetic-geometric mean inequality minimizes the last two terms at
tau=sqrt(A B_n(t)), proving (T). This is a scalar choice, independent of f.

## 9. n=7 expression and connection to the recorded decimal

Put q=1+1/(2t) and h=(1+1/t)/2. Expanding a sixth power gives the exact identity

    integral_1^t (q-1/(2rho))^6 d rho/rho
      =q^6 log t+sum_{j=1}^6 binom(6,j)q^(6-j)(-1/2)^j
                                      (1-t^(-j))/j.

The six coefficients multiplying q^(6-j)(1-t^(-j)) are

    -3, 15/8, -5/6, 15/64, -3/80, 1/384.

This follows by integrating rho^(-j-1) exactly; it is not quadrature. Hence
B_7(t)=1-h^7+7 log t-(7/2) times this integral. This recovers the exact scalar
formula implemented by R028, including the possibly signed log coefficient
7-(7/2)q^6.

R028 specifies t=3391749/524288. A bounded local standard-library Fraction run
in this review independently checked all six coefficients, the recorded exact
rational coefficient and constant for this t, and the exact comparison between
the stored dyadic upper endpoint and 24.554960525965. It made zero calls to
log, sqrt, quadrature, optimizers, scans, or candidate generators. An attempted
SymPy import failed because it is absent from this python3 environment; the
check then used Python's existing exact Fraction arithmetic, with no install.

That finite check does not validate the transcendental enclosure stored in R028.
If an independent scalar review establishes that its interval contains
F(t)=[(1+1/t)^(7/2)+sqrt(B_7(t))]^2, then the accepted theorem and the strict
recorded endpoint comparison establish

    C_7 <= F(3391749/524288) <24.554960525965.

The operator theorem itself has no such numerical dependency.

## 10. Issues and paper-ready scope

No issue requiring mathematical repair was found. Three proof details should
be explicit in a paper: use the first rejection prefix in Section 2; use the last
selected index among least-radius positive contributors in Section 3; and keep
the self a^(-n) term and finite b+1 covering factor until the scalar limit.
These are supplied above rather than assumed from a previous audit.

The exact acceptable statement is (T) for all integer n>=1 and all real t>1,
with general L1 input on Lebesgue R^n and all positive common cube radii. Its
n=7 rational-plus-log specialization and the implication from any correctly
enclosed scalar value are accepted. A paper may present this as an internally
derived theorem after the parent integrates independent reviews. This review
does not establish a publication priority claim, a best known global record,
a global minimizer in t, a sharp constant, any rectangular maximal theorem
with independently varying side radii, or any lower-bound result.

Only the assigned independent clone's `second_review/operator/` was authored.
Source files, main repository, other projects, original three/four-dimensional
materials, and lower-bound originals were not modified. No commit or push was
made. Work stops at this review; no paper is generated here.
