# A072: negative normal cusps, independent six-column frames, and the body32 gap

The target is local semiconcavity of the original body32 Green function on
R6 minus its true S2. This note does not prove that target and does not provide
a body32 counterexample. It proves a normal-profile identity, a local sufficient
condition allowing finitely many overlapping rank3 strata, and semiconcavity
for the full six-independent-column Cauchy family. It also gives an abstract
counterexample to inferring semiconcavity merely from a negative cusp leading
term. The finite fixtures support only their explicitly stated scopes.

Here semiconcavity means: every point has a ball on which D²f<=C I in the
distributional sense, for a finite C>=0. Equivalently f-C|x|²/2 is concave on
that ball. A finite local cover gives a uniform local constant on each compact
subset of the open domain. No segment crossing the removed S2 is asserted.

## 1. Exact three-dimensional normal deficit

Let q_t be a symmetric three-dimensional Cauchy stable density with exponent
psi, a positive continuous even degree-one homogeneous function satisfying
c|xi|<=psi(xi)<=C|xi|. This includes every true three-dimensional orthogonal
quotient of a finite spanning Cauchy frame. With Fourier normalization
q_t(z)=(2pi)^(-3) integral exp(-t psi(xi)) cos(xi.z) dxi, define

    N_psi(z)=integral_0^infinity t³(q_t(0)-q_t(z))dt.

The integrand is nonnegative because Fourier inversion expresses the
difference using 1-cos>=0. Tonelli, integral t³exp(-t psi)dt=6/psi^4, and polar
coordinates in R3 give the exact identity

    N_psi(z)=[3/(8pi²)] integral_S2 |theta.z|/psi(theta)^4 dtheta.       (1)

For the radial integral, integral_0^infinity (1-cos(ar))/r² dr=pi|a|/2.
One nonoscillatory proof uses r^-2=integral_0^infinity u exp(-ur)du,
then Tonelli and the elementary Laplace transform: the inner result is
a²/(u²+a²), whose u integral is pi|a|/2. These steps also prove finiteness;
they do not interchange a conditionally convergent signed integral.

Equation(1) is a positive integral of absolute linear functionals, hence an
even convex degree-one norm. It has finite Lipschitz constant; for instance

    [3/(4pi C^4)]|z| <=N_psi(z)<=[3/(4pi c^4)]|z|.                    (2)

Thus -N_psi has negative-semidefinite distributional Hessian. This is a
renormalized normal-profile identity. It does not by itself identify the
leading coefficient, remainder, or summed strata of the full body32 G.

## 2. A sufficient condition that permits every finite overlap

On a ball B suppose an actual function has an exact representation

    f(x)=R(x)-sum_(j=1)^J A_j(x) N_j(P_j(x-x_j)),                     (3)

where J is finite, P_j are orthogonal projections, N_j are norms with Lipschitz
constants L_j, A_j>=0 are C² on a neighborhood of the closed ball, and R has a
finite distributional Hessian upper bound C_R I there. A C1,1 remainder is a
sufficient special case. Put M_j=sup_B N_j(P_j(x-x_j)). Then

    D²f <=[C_R+sum_j(M_j||D²A_j||infinity
                         +2L_j||grad A_j||infinity)] I.             (4)

To prove this without assuming smooth norms, mollify each N_j with a positive
mollifier. Convexity, nonnegativity and the same Lipschitz bound are retained.
For a unit direction h, the Hessian product identity is

    D²(-A N)[h,h]= -A D²N[Ph,Ph]
                  -2 DA[h] DN[Ph]-N D²A[h,h].

The first term is nonpositive. Bounding the other two terms gives(4), with
M_j+o(1) before passing to the locally uniform mollification limit. This limit
preserves the distributional upper inequality.

No disjointness of the strata, projections, or normal neighborhoods is needed.
Coincident spans can be combined by adding their nonnegative amplitudes.
Nevertheless an actual representation(3), its signs, and its uniform remainder
bound must first be proved; a formal asymptotic expansion is not enough.

If additionally 0<m<=f<=M and f is locally Lipschitz with constant L, then for
1<q<2 the chain rule under smooth approximation gives a semiconcavity constant

    q M^(q-1) C_f + q(q-1)m^(q-2)L²                                  (5)

for f^q, where C_f>=0 is a valid constant for f. Thus the proposed route toward
G^(6/5) integration is mathematically coherent if such quantitative bounds are
available. This note supplies no body32 values of m,M,L,C_f and no area integral.

## 3. Complete proof for axis6 and any six independent columns

For axis6,

    p_t(y,z)=product_(i=1)^3 c_t(y_i) product_(j=1)^3 c_t(z_j),
    c_t(s)=t/[pi(t²+s²)].

Any point outside its true S2 has at least three nonzero coordinates. Relabel
those coordinates as y and choose a closed neighborhood with |y_i|>=delta>0;
the remaining three coordinates are z in a bounded set. Define

    A(y)=pi^-3 product_i y_i^-2>0,
    H_y(t)=product_i(1+t²/y_i²)^-1,
    q_t(z)=product_j c_t(z_j).

Then the exact factorization is p_t=A(y) H_y(t) t³q_t(z). The function
F(y)=G(y,0)=A(y)pi^-3 integral_0^infinity H_y(t)dt is C² (in fact smooth)
on this neighborhood: H_y and its first two y derivatives are bounded at0
and O(t^-6) at infinity, uniformly on compact y sets away from zero.

Subtracting the normal deficit(1) gives the exact representation

    G(y,z)=F(y)-A(y)N_axis3(z)+R(y,z),
    R=A(y) integral_0^infinity t³(q_t(0)-q_t(z))(1-H_y(t))dt.          (6)

Here R is C² jointly in y,z with uniformly bounded Hessian on the chosen
closed neighborhood. The following bounds justify all required derivatives.

* At small t, 1-H_y and every y derivative of order1 or2 are O(t²), uniformly.
  For |alpha|<=2, |D_z^alpha q_t|<=C_alpha t^(-3-|alpha|).
  Thus the differentiated integrands in(6) are bounded by O(t^(2-|alpha|)),
  integrable even at |alpha|=2. No signed cancellation is needed here.
* At large t, evenness gives grad q_t(0)=0. On bounded z sets,
  q_t(0)-q_t(z)=O(t^-5), grad q_t(z)=O(t^-5), and D²q_t(z)=O(t^-5).
  These follow from the global O(t^-5) Hessian bound and Taylor's theorem.
  After multiplying by t³ the bounds are O(t^-2). The y derivatives of
  1-H_y are bounded, as are A and its local derivatives, so mixed derivatives
  are also uniformly integrable.

Dominated differentiation now proves the asserted C² remainder. Applying(4)
to(6) proves:

    Axis6 Green is locally semiconcave on R6 minus true S2.            (7)

It need not be differentiable on S3: the negative norm cusp is real. For any
invertible six-column matrix N, G_N(x)=|det N|^-1 G_axis6(N^-1x), so the same
statement follows under a linear change of coordinates, with the appropriate
finite Jacobian and operator-norm factors in its semiconcavity constant.
This theorem concerns exactly six independent columns, not an arbitrary
frame with additional coupled columns. No extension to body32 is inferred.

## 4. New true path and a genuine overlap geometry

The fixed new axis6 point is x0=(1,2,4,0,0,0), distinct from069's old fixture.
It is off S2 and on S3. The one symbolic path x(s)=(1,2,4,s,0,0), 0<s<1/2,
has the true density integral

    G(x(s))=pi^-6 integral_0^infinity
       t²/[(t²+1)(t²+4)(t²+16)(t²+s²)]dt.                            (8)

Exact partial fractions evaluate(8); no formal DD derivative is used at S3.
The one-sided cusp coefficient predicted by(6) is

    pi^5 G(x(s))=7/1440-|s|/128+O(s²).                              (9)

The derivative jump is therefore -1/(64pi^5), a negative atom in the second
distributional derivative, compatible with semiconcavity. The finite program
also checks the finite one-sided second-derivative limit directly from(8).

The fixed8-column frame consists of I6 plus v7=e1+e4 and v8=e2-e4. At
x*=(1,1,1,0,0,0), the selected distinct rank3 spaces

    W1=span(e1,e2,e3), W2=span(v7,v8,e3)

both contain x*. Their intersection is span(e1+e2,e3), of dimension2.
The complete28-pair qualification proves x* is nevertheless off this frame's
true S2. Hence rank3 strata can intersect away from S2, even for a tiny finite
frame. This is an exact geometry stress fixture only: no Green or Hessian of
this8-column frame is computed, and no semiconcavity conclusion for it follows.
Lemma(4) says finite overlaps would be harmless if the full signed decomposition
and uniform remainder were established there; it does not establish them.

## 5. A precise counterexample to the leading-cusp shortcut

On |s|<1/2 set f(0)=0 and

    f(s)=-|s|+s² log(1/|s|), s!=0.

It is continuous, C² away from0, and f(s)=-|s|+o(|s|). Its first derivative
jump at0 is -2. Yet for s>0,

    f''(s)=2log(1/s)-3 -> +infinity.

No finite Hessian upper bound exists on any neighborhood of0. The negative
Dirac atom cannot compensate for arbitrarily large positive classical curvature
at nearby nonzero points. Thus even a negative cusp, a negative gradient jump,
and C² off the singular set do not imply local semiconcavity without a controlled
remainder. This is an abstract analytic counterexample, not a Cauchy-frame or
body32 counterexample. In particular it does not refute the target conjecture.

## 6. Exact missing steps for body32

Accepted052 gives finite continuous G off S2; accepted069/B069 gives C² only
off true S3. Accepted040's positive J2 divergence on S3 proves neither signed
Hessian divergence nor a positive gradient jump. None of these facts fills
the following gaps:

1. Construct an exact local decomposition such as(3) for the actual coupled
   body32 law near every off-S2 point meeting S3, with all directions sharing
   a rank3 span accounted for and with proved nonnegative amplitudes.
2. Prove that every amplitude has the needed uniform local derivative bounds.
   It is not enough to identify a pointwise positive leading coefficient.
3. Prove an upper-Hessian or C1,1 remainder uniformly through all rank3 overlaps
   and changing relative normal scales. The separated axis identity(6) uses
   independence of the tangential and normal blocks and does not survive the
   additional body32 columns without further work.
4. Supply computable semiconcavity/Lipschitz constants and compatible scalar
   enclosures before using(5) for a continuous full-G^(6/5) area certificate.

No body32 positive-gradient-jump or positive-Hessian-divergence counterexample
has been found by this bounded task. No body32 semiconcavity theorem, global g,
complete I or Cn improvement is claimed. The parent retains final proof review.
