# A085: exact local data for the actual body32 rank3 layer

This package concerns only the original unit-weight body32 frame
V={(1,s2,...,s6):si=+/-1} and the **new** point

    x*=(1,1,1,1,1/2,0).

It computes the actual maximal in-plane four-column law and the actual
28-column quotient data needed by accepted076. It supplies no complete
body32 Green semiconcavity constant Gamma, scalar vertex certificate,
continuous-I certificate, global g cover, or new Cn bound. The parent retains
final proof acceptance. All finite claims below are separated from the
analytic inequalities used to turn the computed rational data into bounds.

## 1. Exact coordinates, selected span, and local scales

Use the orthogonal coordinates (y,z), with

    y0=(x1+x2+x3+x4)/2,  y1=x5,  y2=x6,
    z1=(x1+x2-x3-x4)/2,
    z2=(x1-x2+x3-x4)/2,
    z3=(x1-x2-x3+x4)/2.

The selected span is U={z=0}; y*=(2,1/2,0). The complete in-plane index set
is exactly the four body columns whose first four entries are1. Their
tangential columns, in the frozen (+1,-1) enumeration, are

    v1=(2,1,1), v2=(2,1,-1), v3=(2,-1,1), v4=(2,-1,-1).       (1)

No other original column lies in U. The orthogonal map and all32 transformed
columns were checked exactly; results retain original zero-based body indices
and full original column tuples, so these labels do not depend on an implicit
kernel ordering.

The parent's independently sealed `primary_unique_U/proof.md` proves that U
is the **only** incident rank3 frame span through x*. Its manifest SHA256 is
137668810b3687776f22f16763b1316f5802565bb340594e90ee04b363b6b5cc; the proof
SHA256 is551616307a7d47c8bd73499e24f15ddfcb253d44f496b5bdea2e97a098d0bbd1.
That proof was read, not rerun here. Its status at this package's preparation
is pending independent review. The present seven-configuration batch does not
claim a separate all-triple classification or duplicate its sign verifier.
Once that appendix is accepted, no second incident amplitude remains for
this particular point. Until then, the single-U claim is explicitly dependent
on that appendix.

This batch independently verifies the stronger exact distance

    dist(x*, true S2)^2=1/8.                                  (2)

For any pair of body columns, a vector in its span has each coordinate in
{+s,-s,+d,-d}. Since x* has nonnegative coordinates, its squared distance to
such a vector is at least its squared error when assigned to the two free
nonnegative levels |s|,|d|. All two-level partitions of the six values
(1,1,1,1,1/2,0), modulo swapping the labels, number32. Their exact minimum
within-cluster squared error is1/8, attained by grouping the four ones and
the last two values. This calculation is a relaxation lower bound for every
body pair. Equality is attained by the actual columns (1,1,1,1,1,1) and
(1,1,1,1,-1,-1) with coefficients5/8 and3/8. Thus(2) is an exact certificate,
not a numerical nearest-plane search or a496-pair run.

For every nonincident independent body triple J, the positive integer Gram
determinant for [A_J,2x*] and Hadamard's bound det(A_J^T A_J)<=6³ give

    dist(x*,span A_J)^2>=1/(4*6³)=1/864>1/32².                (3)

This requires no triple enumeration. Dependence is handled by choosing an
actual basis. With sum_v|v|=32sqrt6<80, the fixed scales are

    rho=1/512,          R0=rho/81=1/41472,
    Q_y: |y-y*|_infinity<=h=1/32,       T=1/64.               (4)

They satisfy the076 geometric construction: 16rho<=1/32, 2R0 sum|v|<2rho,
and the in-plane ball of radius4rho fits inside Q_y. The radius is a local
proof parameter; it does not certify the values of G on that ball by itself.

## 2. Exact actual in-plane density, including normalization

Write c_t(s)=t/[pi(t²+s²)]. The first three columns of(1) have determinant
of absolute value8. Solving their parameters with the fourth parameter w
gives

    w1=w+(y1+y2)/2,
    w2=y0/4-y2/2-w,
    w3=y0/4-y1/2-w,
    w4=w.

The true three-dimensional density is consequently

    p_t^U(y)=1/8 integral_R product_(j=0)^3 c_t(w-zeta_j(y)) dw,
    zeta=(0, -(y1+y2)/2, y0/4-y2/2, y0/4-y1/2).             (5)

The change-of-variable Jacobian1/8 is essential. This is the law of all four
original in-plane variables, not a selected three-column subframe. Put
d_(j,l)=zeta_j-zeta_l. For separated centers, all four upper-half-plane poles
are simple, and contour integration gives

    p_t^U(y)=t³/(8pi³) sum_(j=0)^3 product_(l!=j)
                                      [d_(j,l)(d_(j,l)+2it)]^-1. (6)

The sum is real, and is even in t after dividing by t³. Taking its constant
term gives the actual nonnegative amplitude

    a_U(y)=1/(8pi³) sum_j product_(l!=j) d_(j,l)^-2>0,
    p_t^U(y)=t³ a_U(y)+t⁵ r_U(t,y).                          (7)

Thus a_U is homogeneous of degree-6 in y. In the original6D Green expansion,
a_U(y)N(z) has degree-5 because N has degree1. Neither a_U nor p_t^U is the
full6D Green function.

At y*, the scalar centers are (0,-1/4,1/2,1/4). The complete four-residue
calculation yields

    pi³ p_t^U(y*)/t³
       =2560/[(1+16t²)(1+64t²)(9+64t²)],
    a_U(y*)=2560/(9pi³),                                    (8)
    pi³ r_U(t,y*)
       =-40960(64t²+7)²/[9(1+16t²)(1+64t²)(9+64t²)],
    pi³ r_U(0,y*)=-2007040/81.

These are exact local symbolic results. Their executed work is one complete
3D density path and four residues; no6D/body32 full kernel was called.

## 3. Uniform C4 amplitude package on Q_y

The six absolute center separations have these exact uniform bounds on Q_y:

| Pair j,l | Minimum | Maximum |
|---|---:|---:|
| 0,1 | 7/32 | 9/32 |
| 0,2 | 61/128 | 67/128 |
| 0,3 | 29/128 | 35/128 |
| 1,2 | 93/128 | 99/128 |
| 1,3 | 61/128 | 67/128 |
| 2,3 | 7/32 | 9/32 |

All are positive. For a linear form d(y)=d(y*)+ell.h and any positive integer
p, its exact Taylor coefficient at multiindex alpha, |alpha|=k, is

    [h^alpha] d(y)^-p
      =(-1)^k binomial(p+k-1,k) (k!/alpha!)
                                 product_i ell_i^alpha_i / d(y*)^(p+k). (9)

For a uniform derivative bound replace d(y*) by the positive lower bound
delta for |d| and replace ell by -|ell|, so all coefficients are positive.
Applying the product formula to(7) gives35 exact partial derivatives at y*
and35 uniform absolute partial-derivative bounds through total order4. These
are stored with their multiindices in results Q003. Every entry is scaled by
pi³; no floating-point input or approximate denominator is used.

If B_alpha is such a scaled bound, then

    A_k=sum_(|alpha|=k) (k!/alpha!) B_alpha                  (10)

is a valid majorant for pi³||D^k a_U|| as a multilinear Euclidean operator:
each component of a unit direction has absolute value<=1. This deliberately
simple conversion is rigorous but not optimized. The exact A0,...,A4 are in
the same result record. In particular,

    757700950491136/4365571466025
       <= pi³ a_U(y) <=652147834224640/1326226321161,       y in Q_y. (11)

The lower bound uses all upper separation bounds in the positive formula(7).

For the global bounded extension used in076, set a_ext(y)=chi(y)a_U(y), with
chi=1 on |y-y*|<=3rho and0 on |y-y*|>=4rho. On the transition let

    v=(|y-y*|²-9rho²)/(7rho²),
    chi=h4(v)=1-126v⁵+420v⁶-540v⁷+315v⁸-70v⁹.

This is nonnegative and C4, and all first four derivatives match at both
ends. With H_r bounding |h4^(r)|, use H0=1, H1=315/128, and for r=2,3,4 the
finite coefficient absolute sums recorded in Q003. Since

    ||Dv||<=8/(7rho)=d1,        ||D²v||<=2/(7rho²)=d2,

the cutoff operator bounds are

    E0=1, E1=H1 d1, E2=H2 d1²+H1 d2,
    E3=H3 d1³+3H2 d1d2,
    E4=H4 d1⁴+6H3 d1²d2+3H2 d2².

Q003 records the exact scaled extension bounds

    pi³||D^k a_ext|| <=sum_(r=0)^k binomial(k,r) E_r A_(k-r). (12)

These global-extension bounds are much larger than the natural-box bounds.
They are provided honestly as finite interface data, not asserted competitive
full-G curvature budgets.

## 4. Uniform C2 package for the t5 correction

The residue formula allows the time remainder to be bounded without a time
grid, numerical integration, or a full3D fan. For each residue j let delta_jl
be its uniform separation lower bound and g_jl=|gradient d_jl| componentwise.
Define the positive multivariate generating function

    M_j(h,tau)=product_(l!=j)
       [(delta_jl-g_jl.h)(delta_jl-g_jl.h-2tau)]^-1.

Its coefficients majorize the absolute spatial/time coefficients of the
corresponding complex residue product. Here 2T=1/32<min delta_jl=7/32, so
the following evaluated series and all needed derivatives converge. Odd
time powers do not enter the real density. Hence, for |alpha|<=2,

    pi³ |D_y^alpha r_U(t,y)| <=
      (alpha!)/(8T²) [h^alpha] sum_j
           { [M_j(h,T)+M_j(h,-T)]/2-M_j(h,0) },             (13)
             y in Q_y, 0<=t<=T.

The right side retains precisely the nonnegative even time powers of order
at least2, divided by T². For t<=T it dominates the corresponding divided
series at t. Derivatives can be taken termwise because of the strict margin.
Thus this is a uniform analytic bound, not validation only at t=0 or t=T.

All10 bounds in(13) and the operator sums R0,R1,R2 obtained by(10) are exact
rationals in Q004. The center value from(8) lies within the zeroth-order
bound. Applying the same cutoff chi gives the extension bounds

    pi³||D^k r_ext(t,.)|| <=sum_(r=0)^k binomial(k,r) E_r R_(k-r),
                                                 k=0,1,2.       (14)

Their normalization is p_t^U=t³a_U+t⁵r_U throughout; the factors t³,t⁵ and
pi^-3 are not included a second time in(12)-(14).

## 5. Actual quotient law and sharp Euclidean coercivity

The28 outside projections form the following seven nonzero directions,
each with four original labels:

    2e1, 2e2, 2e3,
    (1,1,1), (1,1,-1), (1,-1,1), (-1,1,1).                 (15)

For the **marginal quotient density only**, independence and Cauchy stability
allow the multiplicity4 to be summed in the exponent. For L1=|xi|_1 and
Linf=|xi|_infinity,

    psi(xi)=8L1+4 sum_(four tetrahedral signs) |sign.xi|
           =8L1+max(16Linf,8L1).                           (16)

To check the last identity, use sign/permutation invariance and write
a>=b>=c>=0. The four absolute values sum to4a when a>=b+c and2(a+b+c)
otherwise. It follows that the exact Euclidean constants are

    16sqrt2 |xi| <=psi(xi)<=16sqrt3 |xi|.                  (17)

For the lower bound in the first region, set r=b+c<=a. Then

    (3a+r)²-8(a²+r²)=(a-r)(a+7r)>=0,

and b²+c²<=r². In the second region max(a,b,c)<=L1/2, so
a²+b²+c²<=L1²/2. For the upper bound, the first region has linear
coefficient norm8sqrt11<16sqrt3, and the second uses L1<=sqrt3|xi|.
Equality in the lower bound occurs at(1,1,0)/sqrt2, and in the upper bound
at(1,1,1)/sqrt3. The checker verifies the finite polynomial identities and
the exact transformed direction list; this prose supplies the uniform
ordered-cone argument.

By076/072 the deficit norm and its Lipschitz constant therefore satisfy

    |z|/(786432pi) <= N_psi(z) <=3|z|/(1048576pi),
    Lip(N_psi)<=3/(1048576pi).                             (18)

The upper bound for a single renormalized signed first-moment layer's
distributional Hessian is, from076,

    6/[pi(16sqrt2)^5]=3/(2097152sqrt2 pi).                 (19)

Equation(19) concerns the all-time **linearly subtracted** moment. The actual
finite-time layer still requires its linear term and its upper-tail budget.
Neither is silently omitted when constructing a future complete remainder.

## 6. All28 weighted labels remain separate

The tangential outside matrix C has

    CC^T=diag(16,28,28),          ||C||²=28.

Fix as quotient differentiation basis the actual outside body columns

    (1,1,-1,-1, 1, 1),
    (1,-1,1,-1, 1,-1),
    (1,-1,-1,1,-1, 1).

Their quotient matrix is2I3. Their tangential transport C_J(2I)^-1 has
operator norm1, checked from its exact Gram matrix. Its use in derivatives
is fixed; it is not chosen as a function of z or w.

For each weighted original label j, a separate quotient basis containing
that column is registered. Axis directions have absolute determinant8 and
squared norm ||B_J^-1||_(2->1)^2=3/4. Tetrahedral directions use that label
and the second/third axis representatives; their absolute determinant is4
and the squared norm is9/2. There are seven distinct inverse matrices but
28 explicit label mappings. The exact norms are verified by the eight
dual signs: ||L||_(2->1)=max_(epsilon in {+/-1}³)|L^T epsilon|.

For the full28-variable positive coarea integrals, these bases give

    integral_(Bw=z) product c_t(w_l) d_Bw <=1/(8pi³)t^-3,
    integral_(Bw=z) |w| product c_t(w_l) d_Bw <=11/(4pi³)t^-2,
    integral_(Bw=z) |w|² product c_t(w_l) d_Bw <=11/(2pi³)t^-1. (20)

The latter two constants sum the original per-label W1,W2 from076. Sixteen
tetrahedral labels have |det|^-1=1/4 and twelve axis labels have1/8, giving
sum_j |det B_(J(j))|^-1=11/2. Every original tangent c_j and weighted index
is retained. The marginal-law grouping in(16) is **not** reused to merge
these signed or tangentially weighted fields without a separate theorem.

## 7. A genuine single-profile budget, not the full Gamma

At the center, equations(8) and(18) imply

    [5/(13824pi⁴)]|z| <= a_U(y*)N_psi(z)
                            <=[5/(6144pi⁴)]|z|.           (21)

Thus the actual selected rank3 layer has an explicitly quantified negative
normal profile. Its norm is anisotropic; (21) gives lower/upper Euclidean
comparisons, not an isotropic equality.

On |x-x*|<=rho, the singular term f_U(y,z)=-a_U(y)N_psi(z) alone has a
valid distributional Hessian upper constant

    C_cusp=3/(1048576pi⁴) [rho A2+2A1]
       =3465893065395297320960
          /[65143564997786663569923 pi⁴],                  (22)

using the unextended natural-box A1,A2 of(10). Its Lipschitz constant is at
most3(A0+rho A1)/(1048576pi⁴). The stored orientation-only decimal values
are respectively0.0005461905056499521 and1.4964826464204757e-5. The exact
rational-times-pi formulas are the certificate; a purely rational upper
bound for(22) also follows by replacing pi⁴ by81.

These constants pay **only f_U**. They do not pay the positive curvature of
the C1,1 remainder in076. A large negative singular Hessian does not cancel
an unpaid positive Hessian at nearby points in a certified upper-bound
argument. In particular C_cusp must never be supplied as the Gamma for G.

## 8. Remaining interface and feasibility decision

Conditional on acceptance of the parent's unique-U appendix, the missing
items for a full local body32 Gamma are now concrete:

* Integrate and sum the selected-U cutoff-removal, weighted quadratic, and
  signed-linear remainders, including the cutoff derivatives, the finite-time
  linear subtraction and upper tail. This package supplies the required local
  in-plane and quotient constants but has not assembled that total.
* Import the independently owned086 rank-at-least-four aggregate bound. This
  package neither proves it nor enumerates2^32 cutoff patterns. Its release
  and source hash must be matched before use.
* Pay the original full-density large-time contribution t>=T, with the same
  model and normalization, and match the resulting Gamma to actual scalar
  vertex enclosures. No native or old32-vertex run was authorized here.

The coarse global h4 extension in(12) already makes the fourth derivative
majorant very large; Q003 records its exact size. Formula(14) likewise has
a large extension penalty. These are existence-level constants and are
not claimed to yield a competitive full-G box penalty. The present round
stops at the fixed package rather than expanding parameters or rerunning
kernels to disguise that limitation. No usable complete Gamma was obtained.
This is a limitation of the current bounds and unpaid ledger, not a claim
that the actual Green curvature is large or that076 fails.

## 9. Execution and source boundary

Seven new configurations ran once under the pre-execution source/input hashes.
All passed. The exact work ledger is32 column transforms,32 scalar two-level
partitions, one3D density path/four residues,35 point amplitude partials,
35 amplitude box bounds,10 remainder box bounds, seven quotient inversions,
56 dual-sign values, and28 label mappings. The cutoff/operator bounds and
coercivity identities are included in those registered configurations.

Guard PID91773 and child PID91774 completed with no timeout or live child.
Checker CPU0.233814 seconds, checker wall0.338497917 seconds, full guard
wall0.402595291 seconds. The whole-program limit was180 seconds including
cleanup, with165 seconds of child wait and10 seconds reserved for killing
and waiting only its own process group. CPU limits were guard soft2/hard29,
child soft/hard27, within the30-second task budget.

No6D/body32 full kernel, native routine, old061 point, old baseline, MATLAB,
recursive agent, automatic retry, staging, commit, or push occurred. The
source and runtime pins are in registration.json; exact output is in
results.json and execution_receipt.json. The primary unique-U appendix is
bound separately after its read-only arrival, leaving the seven original
configuration sources and registration unchanged.
