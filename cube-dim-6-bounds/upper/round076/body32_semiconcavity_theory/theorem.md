# A076: a fixed-cutoff route to local semiconcavity for finite Cauchy frames

This is a new analytic proof submission, pending the parent's independent
review. Its mathematical scope is the general finite-frame statement below.
The finite checks accompanying it do not constitute a body32 calculation, a
body32 curvature constant, a whole-face integral, a global g certificate, or
a Cn improvement. No claim about differentiability on the exceptional sets
is made. A072's three-dimensional deficit norm and finite-overlap product
lemma are accepted inputs, cited explicitly rather than rerun.

## 1. Statement and conventions

Let a_1,...,a_m be nonzero vectors spanning R6, with all positive spectral
weights absorbed in these columns. Put

    c_t(s)=t/[pi(t²+s²)],
    p_t=density of sum_i a_i W_(i,t),     W_(i,t) independent with density c_t,
    G(x)=integral_0^infinity p_t(x)dt,
    S2=union_(|J|<=2) span(a_j:j in J).

The set S2 is the true frame set, not the union of all abstract two-dimensional
subspaces, and not the set of pairwise intersections of rank3 frame spans.
Any subset of actual columns of rank at most two has its span in S2, by
selecting a basis from that subset. Repeated or parallel columns are allowed.

**Proposed theorem.** Every x0 outside S2 has a ball B on which G is finite,
locally Lipschitz, and satisfies D²G<=C_B I as distributions for a finite,
in-principle computable constant C_B. More precisely, on B there is an exact
representation

    G(x)=R(x)-sum_(U incident) a_U(Y_U x) N_(psi_U)(Z_U x),          (1)

where the sum is over finitely many distinct rank3 frame spans containing
x0; Y_U and Z_U are orthonormal tangential/normal coordinates; a_U>=0 is C4
on the required local set; R is C1,1 locally; and N_(psi_U) is the convex
three-dimensional deficit norm in A072. Empty sums are permitted. In
particular this proposed theorem specializes to the original body32 model,
without symmetry or general-position assumptions. No numerical value of
C_B is supplied for that specialization.

The proof below constructs (1), including its variable-amplitude remainder.
It does not infer semiconcavity from a leading asymptotic alone, and does not
use the divergent positive J2 bound from040 as a signed-Hessian estimate.

## 2. Fixed smooth parameter cutoffs and the exact rank partition

Write S=sum_i |a_i|. Let delta=dist(x0,S2)>0. Let F3 be the finite set of
distinct rank3 spaces spanned by actual independent triples. Let gamma be
the minimum distance from x0 to members of F3 not containing x0; take gamma=1
if this list is empty. Choose

    rho=min(delta,gamma,1)/16,       R0=rho/(1+S),
    B={x:|x-x0|<rho}.

Then 2 R0 S<2rho. Use an even C2 cutoff chi=chi_(R0), equal to1 for |s|<=R0,
zero for |s|>=2R0, and between0 and1. On the transition write

    chi(s)=H((|s|-R0)/R0),
    H(v)=1-10v³+15v⁴-6v⁵,       0<=v<=1.

Its first two derivatives match zero at both endpoints. We use
|chi'|<=2/R0 and |chi''|<=6/R0². Set eta=1-chi. Expand the scalar product law
using c_t=(c_t chi)+(c_t eta). A pattern L specifies the indices with eta;
its factors are fixed C2 functions of their scalar arguments. There are
finitely many patterns and their pushforward densities sum exactly to p_t.
The phrase "large columns" below means the columns with indices in L. A
nonzero integrand in this pattern has |w_i|>R0 for i in L and |w_i|<2R0
for i outside L. No indicator depending on x is differentiated.

If rank(A_L)<=2, its large sum is in true S2 and its remaining sum has norm
at most 2R0 S<2rho. This cannot reach B, whose distance from S2 is at least
delta-rho>=15rho. The entire pattern density is identically zero there.
If rank(A_L)=3 but U=span(A_L) does not contain x0, the same argument uses
dist(B,U)>=gamma-rho>=15rho. Such a pattern is also zero throughout B.

For an incident U define the **maximal in-plane index set**

    I_U={i:a_i in U},          O_U={1,...,m}\I_U.

In orthonormal U and U-perpendicular coordinates write a_j=(c_j,b_j) for
j in O_U and define Cw=sum c_j w_j, Bw=sum b_j w_j. Every b_j is nonzero,
and these b_j span R3: the full original frame spans R6. These two facts
will be needed separately later.

For x=(y,z) in B and all outside w in the chi support, u=y-Cw lies within
3rho of x0 in U. At this u, a representation by the in-plane frame cannot
have its in-plane large columns of rank<=2: its small in-plane sum has norm
less than2rho, whereas dist(u,S2)>=delta-3rho>=13rho. Thus summing every
in-plane cutoff pattern of rank3 is exactly the *full* in-plane density
p_t^U(u). Consequently the aggregate of all patterns whose large-column
span is U has the following exact local density:

    K_U(t,y,z)= integral_(Bw=z) p_t^U(y-Cw)
                       product_(j in O_U)[c_t(w_j)chi(w_j)] d_B w. (2)

Here d_Bw is positive coarea measure dH^(n-3)(w)/sqrt(det(BB^T)), n=|O_U|.
Equivalently one may solve any three independent quotient columns J and
integrate the remaining variables, with Jacobian |det B_J|^-1. This is an
ordinary pushforward-density formula, not conditioning on a probability-zero
event.

The patterns have a unique rank and, at rank3, a unique span U. Distinct U
may intersect geometrically away from true S2, so their resulting densities
may both be nonzero at the same spatial point. This is allowed. If a *single
actual large-column subset* were contained in two distinct rank3 spaces,
its rank could be at most two; it is that actual subset's span which lies
in true S2. We do not assert that the entire geometric intersection lies
in true S2.

## 3. Derivatives of rank-at-least-four patterns keep four tails

For 0<t<=1 the large factor l_t=c_t eta and its derivatives obey

    |l_t|<=t/(pi R0²),
    |l_t'|<=4t/(pi R0³),
    |l_t''|<=20t/(pi R0⁴).                                      (3)

Each derivative in (3) is supported in |s|>=R0. For example use
|c_t'|<=2t/(pi|s|³), |c_t''|<=6t/(pi|s|⁴) and the two cutoff bounds. Thus
derivatives of the multiplier do not lose the large-tail support. Either
type of factor f_t=c_t chi or c_t eta also satisfies the global bounds

    |f_t^(k)|<=B_k(R0) t^(-1-k),       k=0,1,2,                  (4)

where one may take B0=1/pi, B1=(1+2/R0)/pi and
B2=(2+4/R0+6/R0²)/pi. Both factors have total mass at most1.

For a fixed pattern with rank(A_L)>=4, choose four independent columns from
L, and complete them to an invertible six-column basis N using two more
actual columns. In the fixed-basis integral, x derivatives act only on the
six selected scalar factors through N^-1. The other factors are integrated
against subprobability measures. Every product-rule term still contains
four large factors or their derivatives. Equations(3)-(4) give, uniformly
in x,

    |D^k p_(t,L)(x)|<=K_(L,k) t^(2-k),       k=0,1,2.           (5)

The constants are finite explicit sums of products of the bounds in(3)-(4),
entries of N^-1, and |det N|^-1. These fixed-t densities and their first two
derivatives are continuous, by dominated differentiation of that same fixed
integral. Integrating(5) from0 to any T<=1 gives a C2 function with bounded
Hessian. Finitely many patterns can be summed. No derivative is applied to
a rank selector, an event indicator, or a moving integration domain.

## 4. Full in-plane density: t³ amplitude and t⁵ correction

Consider the three-dimensional frame consisting of all columns in I_U.
On the closed U-ball |u-x0|<=4rho it is separated from its own true S2, since
that S2 is a subset of the original one. Its density has a finite full-fan
Fourier formula

    p_t^U(u)=(2pi)^-3 sum_cones alpha Re product_(j=1)^3
                                      (t-i beta_j(u))^-1,         (6)
    alpha=|det(r1,r2,r3)|/product_j psi_in(r_j),
    beta_j(u)=r_j.u/psi_in(r_j).

One constructs this fan by cutting R3 with the hyperplanes a_i.xi=0 and
triangulating each pointed polyhedral cone using its existing extreme rays.
In dimension3 this only triangulates a polygon and requires no new rays.
Each extreme ray is perpendicular to two independent in-plane columns, so
r_j.u cannot vanish on the stated off-S2 compact set. Thus there is a
positive uniform lower bound beta_* for all |beta_j| there. The finite
formula(6) is real analytic in t for |t|<beta_* and smooth in u with uniform
derivative bounds of every fixed order.

The real part of a product of three factors in(6) is odd in t. Moreover,
three independent large scalar parameters are necessary to reach the
off-S2 compact set. Selecting them and bounding their densities by tail
bounds as in the rank argument gives p_t^U(u)=O(t³), uniformly in u. This
last estimate uses a positive event decomposition only for a value bound;
it does not differentiate the events. The analytic t coefficient therefore
vanishes identically on the open neighborhood. Its u derivatives vanish
there as well. It follows that, for some fixed 0<T<=min(1,beta_*/2),

    p_t^U(u)=t³ a_U(u)+t⁵ r_U(t,u),       0<t<=T,                (7)

where a_U is smooth and nonnegative (take the positive density's t^-3
limit), and r_U and its required u derivatives are uniformly bounded.
Shrink T to work for all incident U; there are only finitely many.

For later full-fiber integrals extend a_U to a globally bounded nonnegative
C4 function on U, agreeing with it for |u-x0|<=3rho and vanishing for
|u-x0|>=4rho. An explicit cutoff can be formed by applying

    h4(v)=1-126v⁵+420v⁶-540v⁷+315v⁸-70v⁹

to v=(|u-x0|²-9rho²)/(7rho²) between0 and1 and clamping to1 and0 outside.
Its first four derivatives match and it is decreasing from1 to0. Extend
r_U with this same cutoff. The extensions retain finite C4 bounds, uniform
in t for r_U. Equation(7) is used only where the extension agrees with the
original density expansion.

In (2), the t⁵ r_U term contributes a C2 time-integrated remainder. Indeed,
choose any fixed quotient basis of three columns. The selected factors
and their derivatives satisfy(4), the residual scalar laws have masses<=1,
and the amplitude and its derivatives are bounded. Its derivatives of
total order k<=2 are O(t^(5-3-k))=O(t^(2-k)), uniformly on B. These bounds
also justify differentiating the fiber integral in fixed coordinates.

The integral from T to infinity of the original p_t is smooth with bounded
local derivatives: a fixed six-column product and its score bounds give
O(t^(-6-k)) for derivatives of order k. This tail is another C2 remainder.

## 5. Remove outside cutoffs before expanding the amplitude

Only the leading part t³ a_U(y-Cw) of(2) remains. Replace the product of
outside chi factors by1. The error uses the exact telescoping identity

    1-product_j chi(w_j)=sum_j eta(w_j) product_(l<j) chi(w_l).    (8)

For the j-th term, b_j is nonzero, so choose a quotient basis containing
this column. In its fixed-coordinate formula the j-th selected factor is
large, and its derivatives retain the support and O(t) bounds in(3).
The other two selected factors obey(4). The bounded amplitude a_U(y-Cw)
and all of its derivatives through order2 are included in the product rule.
Residual factors have masses<=1. Consequently this error's derivative of
total order k<=2, including its leading t³, is O(t^(3+1-2-k))=O(t^(2-k)).
Its time integral is C2.

This step is deliberately performed on the **unexpanded bounded amplitude**.
There are no scalar w_l weights in(8). We never require an arbitrary weighted
column and the large column to belong to a common independent basis; they
may be parallel. Only after this C2 cutoff error is paid do we consider the
unrestricted full outside law and write

    a_U(y-Cw)=a_U(y)-Da_U(y).Cw+E_U(y,w).                       (9)

The quadratic remainder and the signed linear term need separate arguments.

## 6. Weighted quotient moments and the quadratic remainder

Let B:R^n->R3 be any full-rank quotient with every column b_j nonzero, and
write P_t(w)=product_j c_t(w_j). For each j choose a three-column basis
J(j) containing j. Set D_j=|det B_(J(j))|^-1. Positive coarea and the elementary
scalar suprema give, for all z,

    M_(j,k)(t,z):=integral_(Bw=z)|w_j|^k P_t(w)d_Bw
                         <= W_(j,k) t^(k-3),       k=0,1,2,      (10)
    W_(j,0)=D_j/pi³,   W_(j,1)=D_j/(2pi³),   W_(j,2)=D_j/pi³.

Indeed sup c_t=1/(pi t), sup |s|c_t=1/(2pi), and sup s²c_t=t/pi.
Use the weighted bound on the j-th solved coordinate, the ordinary peaks
on the other two, and integrate all remaining scalar probability densities.
The selected basis may depend on j for this **bound**, and no two arbitrary
columns are asserted to be simultaneously independent. The unconditional
Cauchy moments need not exist; (10) is a pointwise density estimate for a
weighted positive pushforward. In particular, with finite constants W0,W1,W2,

    integral_(Bw=z) P_t <= W0 t^-3,
    integral_(Bw=z) |w| P_t <= W1 t^-2,
    integral_(Bw=z) |w|² P_t <= W2 t^-1.                         (11)

The last two follow by |w|<=sum|w_j| and |w|²=sum w_j².

Now a is the bounded C4 extension in(9). Taylor's theorem and bounded
derivatives give the following uniform estimates, including y derivatives
through the orders that make total order at most2:

    |D_y^r E(y,w)| <= K |w|²,                 r=0,1,2,
    |D_w D_y^r E(y,w)| <= K |w|,              r=0,1,
    |D_w² E(y,w)| <= K.                                          (12)

For example E=a(y-Cw)-a(y)+Da(y).Cw. Its w derivative is
C^T[Da(y)-Da(y-Cw)], its w Hessian is C^T D²a(y-Cw) C, and two y
derivatives in the Taylor error require a through order4. This is why
C4, rather than merely C2, was retained for the extension.

Fix **one** quotient basis J0 for all differentiations, put L=B_(J0)^-1,
and solve w_(J0)=L(z-B_rest w_rest). The moving fiber has become a fixed
Lebesgue integral. Each z derivative acts as a fixed directional w derivative
on E and on the three selected Cauchy factors. Their relative derivative
bounds through order2 are O(t^-1) and O(t^-2). The second z derivative of
the integrand is therefore bounded in integral by

    K[ t^-2 integral|w|²P_t + t^-1 integral|w|P_t
                                      + integral P_t ]
        <= K' t^-3.                                             (13)

Mixed y,z and y,y derivatives are bounded by K(t^-1 M2+M1) and K M2,
respectively, with M1,M2 as in(11); they are even smaller as t decreases.
Different bases were used only to estimate the three positive weighted
integrals after the fixed-J0 differentiation. This is legitimate equality
of coarea integrals, not a derivative of a basis selector.

For fixed t>0, these differentiations are justified locally by dominated
convergence. On a bounded z neighborhood, shifts of a selected c_t factor
are bounded by a finite constant times its unshifted value; the polynomial
weights in(12) are controlled by(11) at a fixed reference z. The derivatives
are continuous. Multiplying(13) by t³ gives a constant small-time dominator
for the second derivatives; lower orders and mixed derivatives are also
integrable. Hence

    integral_0^T t³ integral_(Bw=z) E(y,w)P_t(w)d_Bw dt            (14)

is C2 jointly in y,z, with a locally bounded, explicitly majorizable Hessian.
No logarithmic divergence is inferred from the rough unweighted score bound.

## 7. The signed first-moment layer is locally C1,1

Define

    m_j(t,z)=integral_(Bw=z) w_j P_t(w)d_Bw,
    psi(xi)=sum_l |b_l.xi|,       psi(xi)>=c|xi|, c>0.            (15)

The integral is absolutely finite at each z by(10). Choosing a quotient
basis containing j proves fixed-t smoothness and derivative bounds

    ||D_z^k m_j(t,.)||infinity <= K_(j,k) t^(-2-k),
                                      k=0,1,2,3.                (16)

Derivatives can fall on the affine solved w_j, in which case a zeroth
weighted moment from(10) replaces a first one. The selected scalar score
bounds and the product rule still give exactly the powers in(16). Constants
are explicit in |det B_J|^-1 and B_J^-1. For instance, with
Lambda=||B_J^-1||_(2->1), one safe third-derivative bound is
K_(j,3)=48 Lambda³/(pi³|det B_J|), using
|c_t^(r)|/c_t <=(r+1)!t^-r for r<=3. The function m_j is odd in z and obeys
m_j(t,z)=t^-2 m_j(1,z/t).

With q_t(z)=(2pi)^-3 integral exp(-t psi(xi))cos(xi.z)dxi, its signed moment
has the exact Fourier representation

    m_j(t,z)=t/(2pi)^3 integral sign(b_j.xi) sin(xi.z)
                                                  exp(-t psi(xi))dxi. (17)

One rigorous justification differentiates
exp[-t|b_j.xi+lambda|-t sum_(l!=j)|b_l.xi|] at lambda=0 in tempered
distributions. Pairing the physical derivative with a Schwartz test is
dominated by the weighted bound(10); the Fourier difference quotient is
bounded by t exp(t|lambda|)exp(-t psi). The hyperplane b_j.xi=0 has measure
zero because b_j!=0. Thus differentiation gives i t sign(b_j.xi)exp(-t psi)
for the Fourier transform with positive phase, and Fourier inversion gives
the positive sine sign in(17). This avoids assuming a finite unconditional
first moment.

The required all-time expression is **renormalized**:

    F_j(z)=integral_0^infinity t³[m_j(t,z)-D_zm_j(t,0).z]dt.       (18)

At t=0 the absolute bounds in(16) make t³|m_j|=O(t) and
t³|D_zm_j(t,0).z|=O(|z|), both integrable. At infinity oddness eliminates
the quadratic Taylor term, and (16) gives a cubic remainder
O(|z|³t^-5); after multiplying by t³ this is O(t^-2). Thus(18) converges.
Neither integral of its two terms separately over all positive time is
asserted convergent.

Subtract the linear term inside(17), and integrate t using
integral_0^infinity t^4 exp(-t psi)dt=24/psi^5. Absolute Fourier Fubini is
valid: sin(xi.z)-xi.z is O(|xi|³) at zero and O(1+|xi|) at infinity, so
its product with psi^-5 is integrable in R3. Polar coordinates then give

    F_j(z)= -3/(4pi²) integral_S²
        sign(b_j.theta) (theta.z)|theta.z| / psi(theta)^5 dtheta. (19)

Here integral_0^infinity[sin(ar)-ar]/r³ dr=-pi a|a|/4 follows by an ordinary
integration by parts from A072's accepted integral of (1-cos(ar))/r²; both
boundary terms vanish. This use does not repeat its earlier numerical fixture.

The real function s|s| is C1 with Lipschitz first derivative and distribution
second derivative 2 sign(s), with no atom at zero. Consequently(19) is C1,1
and has the global two-sided Hessian bound

    ||D² F_j|| <= [3/(2pi²)] integral_S² psi(theta)^-5 dtheta
                 <=6/(pi c^5).                                 (20)

The finite-time layer L_j(z)=integral_0^T t³m_j(t,z)dt is related to it by

    L_j(z)=F_j(z)+T D_zm_j(1,0).z - Q_j(T,z),
    Q_j(T,z)=integral_T^infinity
                       t³[m_j(t,z)-D_zm_j(t,0).z]dt.             (21)

Q_j is C2 on every bounded z set. For second derivatives, oddness gives
D²m_j(t,0)=0 and(16) bounds |D²m_j(t,z)| by K_(j,3)|z|t^-5, leaving
an integrable t^-2 tail. The analogous cubic and quadratic Taylor bounds
handle orders0 and1. In particular ||D²Q_j||<=K_(j,3) Z/T on |z|<=Z.
Equations(20)-(21) prove the desired local C1,1 statement with a quantitative
Hessian bound. The linear term in(21) has zero Hessian but must be retained
for the exact identity. Multiplication by the smooth coefficients
Da_U(y).c_j preserves C1,1 on bounded local sets, including mixed derivatives.

## 8. Normal deficit, actual remainder, and all overlaps

After Sections5-7, the only rank3 contribution not already in a C1,1
remainder is

    a_U(y) integral_0^T t³ q_t(z)dt,

where q_t is the full quotient density for all columns outside U. A072 proves

    N_psi(z)=integral_0^infinity t³(q_t(0)-q_t(z))dt
            =3/(8pi²) integral_S² |theta.z|/psi(theta)^4 dtheta.  (22)

Hence the truncated layer equals

    integral_0^T t³q_t(z)dt
       =T q_1(0)-N_psi(z)
                       +integral_T^infinity t³(q_t(0)-q_t(z))dt. (23)

The last term is C2 on bounded z sets: evenness, scaling, and the global
second-derivative bound O(t^-5) give O(t^-2) integrands for its derivatives
through order2 after the indicated Taylor cancellations. Its multiplication
by a_U is C2. The coefficient Tq_1(0) uses q_t(0)=t^-3q_1(0).

Combining (2), (5), (7), (8), (14), (21), and (23) yields the **exact** local
representation(1). Its remainder comprises finitely many C2 terms and the
explicit C1,1 signed-moment terms; it is therefore C1,1 locally. Each a_U is
nonnegative, as required. The construction includes every incident rank3
span once and remains valid when different normal variables tend to zero
at different rates. All remainders were bounded on the same full ball B,
not only along an isolated normal line or one stratum.

Apply A072's finite-overlap product lemma. With L_U a Lipschitz constant for
N_(psi_U), M_U=sup_B N_(psi_U)(Z_Ux), and C_R a local Hessian norm bound
for the constructed R, a permissible semiconcavity constant is

    C_B=C_R+sum_U [M_U ||D²(a_U o Y_U)||infinity
                         +2L_U ||D(a_U o Y_U)||infinity].         (24)

The dropped term -a_U D²N is nonpositive as a distribution because a_U>=0
and N is convex. Positive mollification gives the same conclusion when N
is nonsmooth. No disjointness hypothesis is used. A072 also gives
L_U<=3/(4pi c_U^4) when psi_U>=c_U|xi|.

This proof therefore supplies the previously missing variable-amplitude and
overlap argument at the analytic level, subject to independent review of
the proof submission. It does not claim that a computed body32 C_B is useful
or small, or that G is C2 on rank3 strata.

## 9. Quantitative contract and vertex certificate

All constants above have a finite constructive majorant. A future actual
box certificate must supply the following data for the *same* frame:

1. Certified positive off-S2 distance and the finite incident/nonincident
   rank3 partition, with a spatial ball and fixed R0 satisfying Section2.
2. Four-large/six-basis choices for the finite pattern bounds, or a rigorously
   proved aggregate majorant which dominates their entire finite sum; each
   relevant determinant and inverse bound must be certified.
3. For each incident U, the complete in-plane frame, a full three-dimensional
   fan without lost cones or multiplicities, positive ray-node margins,
   a common T, and certified C4 amplitude / C2 correction bounds for(7).
4. The complete quotient frame, certified coercivity c_U>0, a quotient basis
   for differentiation, and for every outside j a basis containing j for(10).
   Parallel outside columns require separate choices, not a false joint basis.
5. Certified bounds for every part of R, including the cutoff derivatives,
   the unexpanded cutoff-removal error, the quadratic weighted remainder,
   and the renormalized signed layer with its linear term and upper tail.
6. A final certified upper bound for(24), together with true G scalar vertex
   enclosures and exact input/window/model identity. No copied finite-part
   Hessian on S3 or old positive J2 bound substitutes for this constant.

For a box Q=center+product_i[-h_i,h_i] contained in such a neighborhood,
let l_v<=G(v) be verified lower bounds at every vertex. Concavity of
G-C_B|x|²/2, with the ordinary multilinear vertex weights lambda_v(x), gives

    G(x)>=sum_v lambda_v l_v
          -(C_B/2)sum_i [h_i²-(x_i-center_i)²]
        >=min_v l_v-(C_B/2)sum_i h_i²,           x in Q.            (25)

On an orthonormal five-dimensional face use its five side coordinates and
32 vertices; a nonorthonormal face parameterization must pull back the
quadratic form rather than silently reuse Euclidean side lengths. A finite
cover of a convex box by qualified neighborhoods may use the maximum of
their distributional constants if the entire box is covered. Formula(25)
is a conditional whole-box lower interface, not an already executed cover.
The G^(6/5) and continuous-I follow-up additionally needs the matching
positive bounds and integration contract. No such whole-face result is
claimed in this round.

## 10. Evidence scope

The executable preparation registers only three new small symbolic
configurations: fixed cutoff algebra, one non-axis four-column three-dimensional
density at u=(1,3,7), and the weighted/signed-layer formal identities. The
tiny density uses four scalar upper-half-plane residues. It is not a six-
dimensional full kernel, not a body32 kernel, and not an old axis baseline.
Its numerical status is reported separately from the theorem's proof status.
No computation can substitute for the uniform estimates and exact partition
above; no successful tiny fixture is presented as a real body32 certificate.
