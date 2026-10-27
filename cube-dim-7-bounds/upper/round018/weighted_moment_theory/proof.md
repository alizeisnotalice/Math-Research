# R018: determinant-weighted directional moments and a two-node Gaussian lower bound

2026-10-27. Pure Sol theory/specification. No new scalar, moment, phase, sample, MC or MATLAB evaluation by this owner. Main has authorized the separate numerical owner to run only the finite7constant/3gate contract after source locking. Original7D body64, unit weights and common Cauchy time remain fixed.

## 1. Weighted probability law and exact normalization

For P~Dirichlet(1/2)^64 put A=sum_v vv^T/P_v, B=A^-1, d=(detA)^(-1/2), and Q_x=x^T Bx. Define

    a=71/2,
    b=(65/2)(67/2)(69/2)=(a-3)(a-2)(a-1),
    K0=Gamma(a)/(pi^(7/2)Gamma32),
    KD=Gamma(a-3)/(pi^(7/2)Gamma32).

The accepted radial formulas are

    p1(x)=K0 E_P[d(1+Q_x)^(-a)],
    G(x)=KD E_P[d Q_x^(-3)]

with the latter extended positive. Write p0=p1(0)>0. The finite positive normalizer is Ed=p0/K0, so

    dnu(P)=d(P)dPr(P)/(Ed),
    G(x)=(p0/b) E_nu Q_x^(-3).                              (1)

The weight is d, not the raw radial moment or the R017 trace proxy. The scalar ratio KD/K0=1/b is exact. The R016 p0 is a finite time-one density; G(0) remains infinite.

Define the directional density Fourier moments

    J_(2k)(x)=(-1)^k D_x^(2k) p1(0)
       =(2pi)^(-7) integral_R7 (x·xi)^(2k)exp(-psi(xi))dxi.   (2)

These are NOT the R011 H_epsilon moments, which have an additional1/psi factor. The Euclidean psi floor gives all moments in(2) finite. Also A>=2048I, so0<=Q_x<=|x|²/2048 uniformly over P. Expansion or dominated differentiation of p1(t x) gives

    mu_k(x):=E_nu Q_x^k
      =k! J_(2k)(x)/[(a)_k p0(2k)!], k0,1,2,3.              (3)

Thus the proposed normalization is correct precisely with the density moments(2). Explicitly,

    mu1=(J2/p0)/(2a),
    mu2=(J4/p0)/[12a(a+1)],
    mu3=(J6/p0)/[120a(a+1)(a+2)].

No new Gamma integration or shape sample is needed to obtain these three moments once the even Fourier coefficients are known.

## 2. Two-node Gauss quadrature is a strict lower bound

For x!=0, Q_x>0 almost surely. Put

    m=mu1, V=mu2-m², W=mu3-2m mu2+m³,
    H=m mu3-mu2².

In the actual model V>0, W>0 and H>0. Here is a proof of nondegeneracy, rather than an assumption from interval numerics. The Dirichlet law has positive density on the simplex interior, and the weighted law nu has the same full support because d is positive continuous there. At uniformP, A=4096I and Q_x=|x|²/4096>0. Let one P_i approach1 and all other63 approach0 at the same rate. Their direction Gram matrix is64I-v_i v_i^T, which is SPD (its smallest eigenvalue is57). Hence A tends to infinity in every direction and Q_x tends to0. Continuity along an interior path shows that Q_x has infinitely many support values. Thus V=Var(Q)>0; W=E[Q(Q-m)²]>0; and H>0 by strict Cauchy-Schwarz for sqrtQ and Q^(3/2). All moments are finite by the uniform upper bound for Q.

Define c=W/V, S=m+c, T=H/V=mc-V. The monic polynomial

    P2(z)=z²-Sz+T

is orthogonal to1,z under the Q law. It has two distinct positive roots r1<r2. One elementary argument is that an orthogonal quadratic must have two sign changes within the convex support: multiplying by the product of any fewer interior sign-change factors would produce a polynomial of degree<2 with nonzero fixed-sign inner product, contradicting orthogonality. Thus the roots lie strictly inside the positive support interval. Their exact weights

    lambda1=(r2-m)/(r2-r1), lambda2=(m-r1)/(r2-r1)

are positive and sum1. Indeed P2(m)=-V<0 places m strictly between the roots. Dividing any cubic by P2 proves this rule is exact for degrees0,...,3. Stable optional root evaluation uses r2=(S+sqrt(S²-4T))/2 and r1=T/r2, not subtraction of close positive numbers.

For f(z)=z^-3, f''''(z)=360z^-7>0. The cubic Hermite polynomial matching f,f' at both r1,r2 satisfies

    f(z)-Hermite(z)=f''''(theta)/4!*(z-r1)²(z-r2)²>=0,

for each positive z, with a suitable theta in the positive interval containing z and the nodes. Integrating the polynomial is legal and finite. Even if E f(Q) is infinite, the inequality remains valid. Therefore

    E_nu Q^-3>=lambda1/r1³+lambda2/r2³.                     (4)

For this actual non-two-point law the inequality is strict, although no positive numerical remainder is claimed. This method retains x through its directional moments and is outside the two named R017 proxy chains.

## 3. Root elimination and a rational expression

The rule is a probability measure of meanm, so it is at least m^-3 by Jensen. Its inverse-third moment can be evaluated without roots or weights. The polynomial relation z²-Sz+T=0 at both nodes gives

    R2(x)=lambda1/r1³+lambda2/r2³
         =[c³+(2c+m)V]/T³
         =[W³+2WV³+mV^4]/H³.                               (5)

For example the first expression follows by the recurrence for negative moments, or by cubing the inverse of the2by2 Jacobi matrix with diagonal(m,c) and off-diagonal sqrtV. Substitution c=W/V and T=H/V yields the second expression with no lost powers. All terms in its numerator are positive, and H>0. The Gaussian lower function is

    L2(x)=(p0/b)R2(x)<=G(x).                                (6)

For exact rational moment coefficients and rational x, R2 is exact rational, so only the common pi^-7 scale needs ball arithmetic. Positive denominators are certified exactly. A naive interval evaluation through nearly coincident roots is unnecessary.

For a general probability law not known to be this model, V=0 means Q=m a.s.; if m>0, the exact answer is m^-3, while m=0 implies Q=0 and the inverse moment is infinite. A two-atom positive law makes(4) exact. Under the current full-support model, degeneracy occurs only at x=0; there define G=L2=infinity rather than use0/0. Failure of an interval to separate V or H is not a proof of degeneracy and does not permit substituting an equality case.

## 4. Exactly six B7-invariant Fourier coefficients

Let s_i=x_i², p_j=sum_i s_i^j, e2=sum_(i<j)s_i s_j and e3=sum_(i<j<k)s_i s_j s_k. For partitions lambda of k<=3 define

    T_lambda=(2pi)^(-7) integral xi1^(2lambda1)...
                           xi_l^(2lambda_l) exp(-psi(xi))dxi,
    rho_lambda=T_lambda/p0.

The six partitions are(1),(2),(11),(3),(21),(111). Sign symmetry removes odd coordinate powers and permutation symmetry identifies their positions. Consequently

    J2/p0=rho1 p1,
    J4/p0=rho2 p2+6rho11 e2,
    J6/p0=rho3 p3+15rho21 sum_(i!=j)s_i²s_j+90rho111 e3.   (7)

Equivalently the last two use only p1,p2,p3:

    J4/p0=3rho11 p1²+(rho2-3rho11)p2,
    J6/p0=15rho111 p1³+(15rho21-45rho111)p1p2
                             +(rho3-15rho21+30rho111)p3.

Thus the whole directional moment function is determined by six constants, not a63Dquadrature or a Green-point grid. Together with the p0 coefficient the finite extraction has seven outputs.

## 5. Extract the seven constants without a full group orbit

For an accepted ordered simplex S, let u_r=xi_r/q_r and alpha_S=|detW|/product q. After rescaling its cone coordinates, the seven coordinate parameters are independent rate1exponentials E_r in the elementary integration formula. Define linear forms ell_j=sum_(r in S)(u_r)_j E_r and P2=sum_j ell_j², P4=sum_j ell_j^4, P6=sum_j ell_j^6. The notation E below means only integrating a polynomial by E product E_r^n_r=product n_r!, not numerical sampling.

Put R_lambda=pi^7 T_lambda. Exact sign/permutation averaging gives

    R0=pi^7p0=7!sum_S alpha_S,
    R1=6!sum_S alpha_S E P2,
    R2=6!sum_S alpha_S E P4,
    R11=5!sum_S alpha_S E(P2²-P4),
    R3=6!sum_S alpha_S E P6,
    R21=5!sum_S alpha_S E(P2P4-P6),
    R111=4!sum_S alpha_S E(P2³-3P2P4+2P6).                  (8)

Here R2 with a partition subscript means the fourth-order coefficient; it is distinct from the Gauss inverse-moment R2(x) in(5). A numerical implementation must use distinct variable names, for example R_part2 and gauss_inverse3. Finally rho_lambda=R_lambda/R0. All seven quantities in(8) are exact positive rationals.

To verify the group factors, a monomial involving l fixed even-powered coordinates, summed over B7, equals2^7(7-l)! times the sum over ordered distinct coordinate assignments. Division by(2pi)^7 leaves(7-l)!/pi^7. The ordered distinct sums for partitions11,21,111 are respectively P2²-P4, P2P4-P6, and P2³-3P2P4+2P6. This proves(8) and avoids explicit645120group images.

A common positive integer denominator L=lcm_r q_r permits ell_j=L^-1 sum_r kappa_(r,j)E_r, with integer kappa=(L/q_r)xi_r. Compute integer polynomials first and divide degree2kresults by L^(2k) only after integration. At degree6 in seven variables there are at most924homogeneous monomials. The concrete sparse schedule in registration gives at most12,621,420integer multiplication/accumulation pairs over212simplices, plus additions, factorial contractions and a small rational final sum. This is an a priori algorithmic count, not observed timing or a hardware saturation claim.

## 6. General x, continuity and what a global certificate would require

For every x!=0, L2(x) is a positive smooth rational function of its polynomial moments times p0. It is B7 invariant and homogeneous of degree-6. In particular it is continuous on every closed cube face. A legal floor obtained solely from this bound must satisfy

    g_cert<=inf_(x in boundary cube)L2(x).

A single upper enclosure of L2(x*) below a proposed target excludes reaching that target with this particular two-node lower function. It is not an upper bound on G(x*) and cannot refute the true Green floor. This remains valid at a face/corner boundary, since continuity gives a relative neighborhood of positive face measure where the lower function also fails the target. If all three initial points pass, that is only a necessary test and gives no global floor.

For a subsequent global analysis no6Dgrid is necessary at the level of formal variables. On the positive face x1=1 set u_j=sum_(i=2)^7 s_i^j, j1,2,3; then p_j=1+u_j and(7) determines L2 from only these three power sums. The actual domain of six equal-weight atoms s_i in[0,1] lies in the compact Hausdorff moment outer set

    [[u1,u2],[u2,u3]]>=0,
    [[6-u1,u1-u2],[u1-u2,u2-u3]]>=0.

The constraints follow by integrating t(a+bt)² and(1-t)(a+bt)² over the six atoms. Certifying a positive inequality over this larger set would suffice, but failure on the larger set need not imply failure on the actual domain. For a target g the exact inequality is p0*(W³+2WV³+mV^4)-b*g*H³>=0, with H>0. Any future polynomial/interval/SOS analysis must separately prove coverage and signs; it is not part of this initial batch.

The simpler weighted Jensen function is(p0/b)mu1(x)^-3, giving a uniform value at |x|²7. It is no larger than L2, but no claim that it or L2 reaches the target can be made without the new moments. No general analytic inevitability of failure has been established here. The finite coefficient-plus-corner gate is therefore meaningful, while an automatic full-domain optimization would be premature.
