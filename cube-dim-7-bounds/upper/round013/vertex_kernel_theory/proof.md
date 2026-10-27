# R013: a true-vertex formula for the original truncated Green kernel

2026-10-27. Sol theory only. Fixed dimension7, body64 directions with first coordinate1, unit weights, common Cauchy time. No mathematical numerical run, MC, MATLAB, commit or push by this owner. The accepted complete R010 ordered fan and the accepted complete R003 spectral-vertex catalog are inputs; their bytes are locked in sources.json. Claims about their completeness remain dependent on their independent certificates.

## 1. Exponential integration over a simplex

For an affinely independent simplex Delta=conv(v0,...,v7), let d_Delta=|det(v1-v0,...,v7-v0)|. If lambda_j=z·v_j are distinct, then

    integral_Delta exp(z·y)dy
      =d_Delta sum_(j=0)^7 exp(z·v_j)
                       /product_(k!=j) z·(v_j-v_k).              (1)

There is no extra factor7! in(1). To prove it, integrate the exponential power series on the standard barycentric simplex. The degree-m contribution is h_m(lambda)/(m+7)!, where h_m is complete homogeneous. This follows from the elementary monomial integral product alpha_j!/(m+7)!. On the other hand the order7 divided difference of t^n is0 for n<7 and h_(n-7)(lambda) otherwise, by the generating identity product_j(1-lambda_j u)^(-1). Thus the same series is DD7 exp, whose distinct-knot formula is the sum in(1). Coincident lambda values are removable singularities of the full simplex expression.

Put P={psi<=1}. A top R010 cone with primitive rays xi_i=Mw_i and q_i=psi(xi_i) yields the simplex conv(0,xi_1/q_1,...,xi_7/q_7), with d_Delta=alpha=|detW|/product q_i because detM=1. Its B7 images positively triangulate P, up to zero-measure boundaries. This is an actual complete signed-permutation expansion; no individual simplex is multiplied by the group order without transforming its phase.

## 2. Aggregation and why artificial vertices disappear

For every node v of this triangulation define the homogeneous rational function

    T_v(z)=sum_(Delta incident to v) d_Delta
                     /product_(u other vertex of Delta) z·(v-u).

Equation(1), summed over simplices, gives

    integral_P exp(z·y)dy=sum_v exp(z·v)T_v(z).                  (2)

Take the local tangent cone C_vDelta of every closed simplex containing v, even when v is not one of that simplex's vertices. These local cones have disjoint interiors and cover the geometric tangent cone C_vP. Indeed any sufficiently short segment from v inside P belongs to one of these simplices, since the finitely many closed simplices not containing v stay a positive distance from v. Conversely each local cone lies in the tangent cone. Overlapping cone interiors would give overlapping simplex interiors near v. A local cone for which v is a simplex vertex is cone(u-v:u other vertex of Delta). If v is not a simplex vertex, it has lineality along the smallest positive-dimensional simplex face containing v.

This distinction avoids an unnecessary assumption that all signed-permutation images of the ordered pulling triangulation agree face-to-face on sorting walls. Only the accepted positive full-dimensional cover with disjoint interiors is used.

Here is a self-contained rational-cone cancellation argument, including cones that have lines. For a simplicial cone K=E R_+^7 with generator columns e_j, distributional differentiation gives

    product_j D_(e_j) 1_K=|detE| delta_0.                       (3)

This follows by the linear coordinate change y=Et and differentiation of the product of seven one-dimensional Heaviside functions. For a finite simplicial partition of a cone, take the product P(D) of all generator directional derivatives, retaining enough multiplicities for every simplicial term. Applying it to the partition identity expresses P(D)1_C as a polynomial in derivatives of delta_0, with polynomial symbol

    P(z) sum_K |detE_K|/product_j(z·e_(K,j)).

If C has a nonzero lineality vector e, D_e1_C=0. Apply D_e as well. Injectivity of polynomial differential operators on delta_0 gives

    (z·e)P(z) sum_K |detE_K|/product_j(z·e_(K,j))=0

as a polynomial identity; therefore the rational sum itself is identically0. Signs from v-u=-(u-v) contribute the same common factor(-1)^7 and do not affect cancellation. If C is pointed instead, the rational sum is its ordinary exponential integral on the nonempty convergence region Re z·h<0 for h in C\{0}, and extends rationally. This also proves independence of any positive simplicial subdivision.

A general rational polyhedral cone, including one with lines, admits a finite positive pointed simplicial partition: subdivide by the coordinate orthants first and triangulate the resulting pointed cones. This permits the same differential-operator argument on the local cones above. Local cones at a point that is not a simplex vertex contribute0 by lineality. A nonvertex node of P, including0 and artificial ordered-wall/fan nodes on a positive-dimensional face, also has lineality in its full tangent cone. The partition identity therefore makes the sum of the remaining simplex-vertex coefficients, namely T_v, identically0. Only the true vertices remain in(2). The accepted R003 catalog has21B7 orbits and504868vertices. This reasoning proves cancellation; it does not claim that floating-point summation of unsimplified artificial terms would be stable.

For a true vertex v, T_v is precisely the rational continuation of integral_(C_vP) exp(z·h)dh. It is homogeneous of degree-7. A triangulation of this pointed cone using only its true extreme rays shows that its possible poles lie on x·e=0 for actual edge directions e of P incident to v. Artificial triangulation diagonals may introduce additional displayed denominators, but their poles must cancel exactly. Generic evaluations must either use reduced templates or test every denominator of the representation actually used.

## 3. Layer integration, regulator and the sign of the logarithm

For epsilon>0,

    exp(-epsilon psi)/psi
       =integral_(psi)^infinity exp(-epsilon r)(epsilon/r+1/r²)dr.

The absolute integral in frequency and r is finite, since psi>=c|xi| for c>0 and dimension7. Thus Fubini is justified before splitting rational vertex terms. Set b_v=x·v. Homogeneity and i^(-7)=i give, for generic real x,

    integral_(rP) exp(i x·xi)dxi=i sum_v T_v(x)exp(i r b_v).

The left side is real. Its real part is -sum_v T_v(x)sin(r b_v). It is O(r^7) as r tends to0. Expanding the finite sum therefore proves the exact cancellation identities

    sum_v T_v(x)b_v^m=0 (m0,...,6),
    sum_v T_v(x)b_v^7=7!vol(P).                                (4)

In particular sum_v T_v(x)b_v=0. Do not integrate individual vertex terms down to0 without this regulator. For delta>0 integration by parts gives

    integral_delta^infinity exp(-epsilon r)(epsilon/r+1/r²)sin(br)dr
      =exp(-epsilon delta)sin(b delta)/delta
          +b integral_delta^infinity exp(-epsilon r)cos(br)/r dr.

After summing with -T_v, the first term tends to0 by(4). The common cos0 term in the second cancels by the same identity. Finally

    integral_0^infinity exp(-epsilon r)(cos(br)-1)/r dr
      =-0.5log(1+b²/epsilon²),

proved by differentiating in b, using the elementary Laplace integral of sine, and integrating from b0. Define C_v(x)=b_vT_v(x). The result is

    H_epsilon(x)=(2pi)^(-7) sum_(v true vertex) C_v(x)
                                  *0.5log(epsilon²+b_v²).       (5)

The common logepsilon term cancels because sum_v C_v=0. The sign in(5) is positive. Each C_v is homogeneous of degree-6 in x and need not be positive for real generic x. No clipping is valid.

There is an independent algebraic check directly from the accepted R011 kernel. For distinct projected boundary nodes the cone DD6 expands as alpha sum_j h_epsilon(b_j)/product_(k!=j)(b_j-b_k). The coefficient at its boundary vertex v is exactly b_v times the simplex coefficient T_v, since the eighth simplex node is0 and contributes the extra denominator b_v. Thus grouping the R011 DD coefficients gives precisely C_v and(5), including its normalization and positive sign. The rational identities extend beyond these generic distinct-node tests. This is a second derivation, not a second numerical run.

Central symmetry gives T_(-v)(x)=T_v(-x)=-T_v(x), so C_(-v)=C_v and the logarithms agree. Hence(5) may be evaluated over252434opposite pairs, with one factor2. This is an outer-term count, not a proof that every coefficient is cheap.

As an independent normalization identity, the large-epsilon expansion has its first nonzero term at b^6:

    sum C_v=0, sum C_v b_v²=0, sum C_v b_v^4=0,
    sum C_v b_v^6=7!vol(P).

It recovers H_epsilon(0)=7!vol(P)/[6(2pi)^7epsilon^6] after scaling x to0, consistent with R011. This limit involves cancellation of degree-6 coefficients and should not be evaluated by substituting x0 into the singular rational terms.

## 4. Exact tangent cone with six zero forms

At a true unit-ball vertex v put Z(v)={z in V:z·v=0} and

    g=sum_(z not in Z(v)) sign(z·v)z.

Then g·v=1 and, for each fixed h and sufficiently small positive t,

    psi(v+t h)=1+t[g·h+sum_(z in Z(v))|z·h|].

Thus the tangent cone is exactly {h:g·h+sum_Z|z·h|<=0}. If Z consists of exactly six independent rows z1,...,z6, take L with rows(-g,z1,...,z6), and define (a,b1,...,b6)=L^(-T)x. Since Lv=(-1,0,...,0), a=-x·v. The cone in (t,s)=Lh coordinates is t>=sum_j|s_j|. On the open convergence domain a<-max_j|b_j|,

    T_v(x)=1/|detL| integral_R6 integral_(sum|s|)^infinity
                            exp(a t+b·s)dt ds
          =-64 a^5/[|detL| product_j(a²-b_j²)],                 (6)

because each one-dimensional integral of exp(a|s|+bs) is -2a/(a²-b²), and the t integral contributes -1/a. Rational continuation gives(6) everywhere its denominator is nonzero. Therefore

    C_v(x)=64 a^6/[|detL| product_j(a²-b_j²)].                   (7)

All constants and signs are thus fixed. Rescaling v or choosing a different six-row basis without simultaneously changing g,L is forbidden. Here v is the true psi-unit vertex; g remains the unnormalised sign gradient of the original64directions.

At a0 with all b_j nonzero, formula(7) is0, not a pole. Conversely x·edge=0 can make displayed denominators vanish even though H_epsilon is smooth. Neither kind of individual-node behavior classifies the physical Green singular set.

## 5. More than six zero forms: exact remaining problem

Choose any six independent zero rows Z0. Each zero row is uniquely c_z Z0. With L=(-g;Z0), define F(s)=sum_(z in Z(v))|c_z·s|. Then F is a norm on R6 and the tangent cone is t>=F(s). On its convergence domain,

    T_v(x)=-1/[a|detL|] integral_R6 exp(aF(s)+b·s)ds,
    C_v(x)=1/|detL| integral_R6 exp(aF(s)+b·s)ds.                (8)

If the sign fan of F is positively simplicially triangulated in6D, with a simplicial cone generated by columns u1,...,u6 and F linear on it, then its contribution to C is exactly

    |detU|/[|detL| product_j(-aF(u_j)-b·u_j)].                  (9)

There is no positivity assertion for the rational continuation at an arbitrary physical x. Formula(9) is a finite, exact alternative to the full7D DD table; an efficient compact rational form still has to be compiled and certified.

A second exact compiler uses the already accepted R010 star. Let v be the ordered representative, r its unique ordered catalog ray, G_v its exact signed-permutation stabilizer, and S range over ordered top simplices containing r. Then

    C_v(x)=sum_(Q in G_v) sum_(S containing r) alpha_S
                         /product_(u in S, u!=v) x·(v-Qu).      (10)

The factor x·v in the simplex T coefficient has canceled. Every global incident simplex appears exactly once; ordered-cone interiors have trivial signed-permutation stabilizer. Formula(10) is a certified source identity, but evaluating it unsimplified can defeat the intended speedup. It may also contain artificial diagonal poles.

These formulas reduce compilation to21prototypes, followed by exact orbit unfolding and opposite pairing. They do not prove that all prototypes have a small numerator/denominator. An advertised <=500000outer vertex terms is possible because the paired count is252434; a corresponding fast, stable operation bound requires measured/certified template sizes and actual ball widths. No such computation is performed in this theory batch.

## 6. Exceptional directions, continuation and epsilon0

For epsilon>0 the true H is globally smooth and holomorphic on |Im x|<c epsilon, c=32sqrt2 from the primary R011 sharp-floor receipt. Equality(5) on the complement of finitely many real pole hyperplanes therefore determines its unique analytic continuation. This is not permission to divide a ball interval containing0, replace a pole by0, or claim an interval proof via a sampled limiting direction. An exceptional point needs a certified confluent/symbolic extension or the already accepted full-fan DD algorithm.

If x avoids all true-edge pole hyperplanes, each T_v(x) is finite. Taking epsilon down to0 in(5), monotone convergence of the actual H to G and finiteness of the vertex sum give

    G(x)=(2pi)^(-7) sum_v C_v(x)log|x·v|.                      (11)

The sufficient simpler implementation contract may require every b_v nonzero. More generally if b_v=0 but T_v is finite, C_v=b_vT_v=0 and that term is defined by its zero limit, so it causes no divergence in(11). This extension must be established in the chosen reduced representation, not inferred from an uncanceled artificial denominator.

On true-edge pole directions the regulated formula may need cancellation across vertices. An individual zero projection, tangent coefficient pole, or artificial triangulation pole does not prove G infinite. R001's exact S3 classification is not replaced. No point H/G formula here proves a global face lower bound g, and the existing all-scale/general-L1 operator comparison still requires its original global hypotheses.
