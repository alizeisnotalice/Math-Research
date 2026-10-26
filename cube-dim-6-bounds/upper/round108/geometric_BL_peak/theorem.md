# 108 — universal original-law quotient peaks via geometric BL

Proof submission; parent review remains separate from arithmetic checks.
The symbol and all quotient measures are exactly those fixed in scope.md.
The Brascamp–Lieb inequality is classical, not a new theorem claimed here.

## 1. The special geometric BL inequality, with analytic details

Let finitely many unit vectors u_i in R^r and c_i>0 satisfy

    sum_i c_i u_i u_i^T=I_r.                                  (1)

Then sum_i c_i=r by taking trace. For the nonnegative integrable functions
needed below, f_i(s)=exp(-alpha_i |s|), alpha_i>0, we prove

    integral_R^r product_i f_i(u_i.x)^c_i dx
                   <= product_i (integral_R f_i)^c_i.         (2)

First replace f_i by f_i,R=f_i 1_[-R,R], R>0. These are bounded,
nonzero, compactly supported functions. Set g_i(t,s)=H_t f_i,R(s), where
H_t has kernel (4pi t)^(-1/2)exp(-s^2/(4t)), and let
h(t,x)=product_i g_i(t,u_i.x)^c_i and F(t)=integral h(t,x)dx.
For t>0 all g_i are positive and smooth.

Here are sufficient estimates for differentiating and integrating by
parts, rather than a formal IBP argument. On a compact time interval
0<a<=t<=b, write m_i=integral f_i,R. From the heat kernel,

 g_i(t,s)<=C_i exp(-(|s|-R)_+^2/(4b))
          <=C_i' exp(-s^2/(8b)),
 |(log g_i)'(t,s)|<=(|s|+R)/(2a),
 |g_i''(t,s)/g_i(t,s)|<=(|s|+R)^2/(4a^2)+1/(2a).          (3)

The first constant may depend on a,b,R,m_i. The last two estimates
follow by representing the derivative ratios as weighted averages over
source points y in[-R,R]. Equation(1) now bounds h by
C exp(-|x|^2/(8b)); its required time and space derivatives are bounded
by that Gaussian times a fixed polynomial. Dominated differentiation
is therefore valid. With ell_i=log g_i evaluated at(t,u_i.x), a smooth
radial cutoff equal to1 on the ball of radius L and0 outside radius2L
has all boundary terms tend to0 as L→infinity by those same estimates.
Consequently the divergence integral is exactly zero and

 F'(t)=integral h [sum_i c_i (ell_i')^2
                    -|sum_i c_i u_i ell_i'|^2] dx >=0.     (4)

Indeed the matrix A whose columns are sqrt(c_i)u_i has AA^T=I_r.
Thus ||A||=1, and applying it to the vector(sqrt(c_i)ell_i') proves
the bracket nonnegative. The identity preceding(4) uses
partial_t log g_i=ell_i''+(ell_i')^2 and
div(h sum c_i u_i ell_i')=h[sum c_i ell_i''+|sum c_i u_i ell_i'|^2].

For the large-time limit, substitute x=sqrt(t)y. Exactly,

 sqrt(t) g_i(t,sqrt(t)s)
 =(4pi)^(-1/2) integral f_i,R(z) exp(-(s-z/sqrt(t))^2/4) dz.

For t>=1 the inequality (s-z/sqrt(t))^2>=s^2/2-z^2/t
bounds the scaled product by C exp(-|y|^2/8). Pointwise its limit is
(4pi)^(-r/2) product m_i^c_i exp(-|y|^2/4); the powers of t cancel
because sum c_i=r. Dominated convergence gives

    lim_(t→infinity) F(t)=product m_i^c_i.                 (5)

At t↓0 heat approximation converges to f_i,R outside the two endpoint
hyperplanes u_i.x=±R; their finite union is Lebesgue null. Fatou gives
integral product f_i,R(u_i.x)^c_i <=liminf_(t↓0)F(t).
Monotonicity(4) and(5) bound this by product m_i^c_i. Finally let R↑infinity;
both sides converge monotonically to the two sides of(2). This proves
the exact BL case used here, including all truncation and limit steps.

The author-hosted Barthe–Cordero-Erausquin–Ledoux–Maurey paper,
*Correlation and Brascamp-Lieb inequalities for Markov semigroups*, p.2,
states the geometric constant1 under the two projection/frame identities
and discusses the heat-flow method. Official arXiv math/0505065v5,
Bennett–Carbery–Christ–Tao, is another primary background reference.
The proof above is self-contained and does not rely on unverified
equality cases or general BL finiteness theorems.

## 2. Project the complete body frame, keeping zero columns correctly

Let E be any r-dimensional subspace of R^6, 1<=r<=6, with orthogonal
projection P. Set a_v=P v and z_v=|a_v|^2 for v in V. Exact sign
orthogonality gives sum_v vv^T=32I_6, hence

    sum_v a_v a_v^T=32I_E, sum_v z_v=32r.                 (6)

For z_v>0 put u_v=a_v/sqrt(z_v), c_v=z_v/32, and
f_v(s)=exp(-32|s|/sqrt(z_v)). These obey(1), and their product to the
c_v powers is exactly exp(-psi(xi)) on E, since
c_v 32/sqrt(z_v)=sqrt(z_v). Each one-dimensional integral is

    integral f_v(s)ds=sqrt(z_v)/16.                      (7)

When z_v=0 its original symbol term is identically zero on E. Omit that
index entirely; it contributes nothing to(6) and no u_v is defined.
There is no 0^0 convention inside the BL application. Later z log z is
defined to be0 at z=0 by continuity. Thus(2) gives

 integral_E exp(-psi(xi))d xi
 <= product_(z_v>0)(sqrt(z_v)/16)^(z_v/32)
 =16^(-r) exp((1/2) E[Z log Z]),                        (8)

where Z=epsilon^T P epsilon and epsilon has all six independent uniform
signs. The average over32 representatives equals the average over all64
signs because Z is even. The same is true for Z log Z, including zeros.

## 3. Exact quadratic moment and a dimension-uniform peak bound

Write d_i=P_ii. For every symmetric P the finite sign expansion gives

 E Z=tr P,
 E Z^2=(tr P)^2+4 sum_(i<j)P_ij^2.

For an orthogonal rank-r projection tr(P^2)=r. Therefore

 E Z=r,
 E Z^2=r(r+2)-2 sum_i d_i^2
       <=r(r+2-r/3)=r M_r,   M_r=2+2r/3,                (9)

using sum d_i=r and sum d_i^2>=r^2/6. This is a finite exact Rademacher
identity, not a Gaussian approximation. Define a probability on the
positive Z atoms by weight Z/r relative to the uniform sign measure.
Jensen for log gives

 E[Z log Z]/r <=log(E Z^2/r)<=log M_r.                  (10)

Equivalently use log z<=log M_r+z/M_r-1 and(9); multiplying by z
is harmless at z=0 by continuity. Combining(8)-(10),

 Q_E(0):=(2pi)^(-r) integral_E exp(-psi(xi))d xi
            <= M_r^(r/2)/(32pi)^r.                    (11)

This is valid for every E, without requiring E to be a body-triple
orthogonal complement. The original projected columns span E by(6).
Also psi(xi)>=sqrt(sum_v(v.xi)^2)=sqrt(32)|xi| on E, so its Fourier
integral is absolutely convergent and defines a continuous density.
It is exactly the time-one density of the projection of the original
body32 law: its characteristic function is exp(-psi|_E). No column
has been deleted from this law except columns whose actual projection
is zero. An arbitrary omitted-column law need not satisfy(6), (9), or(11).

For time t>0 and every x in E, positive density and Fourier inversion give

    q_(E,t)(x)<=q_(E,t)(0)=t^(-r)Q_E(0).                (12)

The first inequality is cos(x.xi)<=1 against the positive Fourier
multiplier; the second is degree-one scaling. Equation(12) is a peak
bound with no x-decay. It is not a pointwise isotropic-Cauchy shape bound.

## 4. The 3D normal-potential coefficient and exact comparisons

For r=3, M_3=4. Thus

    Q_E(0)<=1/(4096pi^3).                               (13)

The accepted098 identity concerns the actual original4D quotient q_t
after projecting away a body pair and its actual3D line marginal:
H(z)=integral_0^infinity t^2 q_t(z)dt=Q_E(0)/(2|z|), z!=0.
Thus every one of those original-law normal potentials satisfies

    H(z)<=1/(8192pi^3 |z|).                             (14)

The normal-potential identity is an098 dependency; it is not proved by
Fourier cos<=1 or by the BL peak alone. In104 the retained scalar Cauchy
law is restored before the time-first identity, so(13) applies there.
It would be incorrect to insert(13) into a deleted-law peak before that
restoration.

The old098 coefficient is Q_old=1/[pi^2(32/sqrt3)^3]. Their ratio is

    lambda=Q_new/Q_old=8/(3sqrt3 pi).                   (15)

The universal sharp-Khintchine coercivity32/sqrt2, if independently
accepted in105, gives ratio sqrt(8/27) relative to098. The new ratio
divided by that ratio is2sqrt2/pi<1 (pi>3>sqrt8). This compares the
universal coefficients only; it neither evaluates nor orders every
105 directional-cap estimate, which can be smaller in its qualified cone.

Every saved104 P norm and B norm is linear in the same old3D peak.
Replacing it by(13) multiplies their sum by lambda and each q=6/5
spatial integral allocation by lambda^(6/5). All104 source-event,
pair-distance, face-Jacobian, disjoint-core, and completeness arguments
remain explicit dependencies. The three-configuration script rescales
only complete saved104 dyadics, preserves the paidCEH allocation, and
does not regenerate any old geometry, function, radial integral, or log.
The unproved candidate g0=1/65000 remains conditional. Any resulting
positive core allowance omits unpaid regions; any failed complete-middle
allocation rejects that expression's competitiveness, not the true model.

## 5. The full6 special case and remaining gaps

For r=6, P=I, so every z_v=6 and (10) has no moment loss. Directly
u_v=v/sqrt6 and c_v=3/16 satisfy(1), giving

    p_t(x)<=p_t(0)<=216/(2^30 pi^6)t^(-6).               (16)

The saved001 determinant160 peak is1/(160pi^6)t^-6; the exact ratio
of(16)'s coefficient to it is135/4194304. This comparison is between
two peak upper bounds. For r=4, (11) gives M_4=14/3 and the coefficient
(14/3)^2/(32pi)^4, potentially useful only where a genuine original4D
peak is the required input. No new tail or spatial budget is supplied.

The actual upper C6<=64 is unchanged. A global g_L, a complete competitive
I_U including far, and positive global S=64(g_L/6)^(6/5)-I_U are still
required. This is a new application/constant and a lawful budget transfer,
not a new classical BL theorem, a lower-bound result, or a C6 breakthrough.
