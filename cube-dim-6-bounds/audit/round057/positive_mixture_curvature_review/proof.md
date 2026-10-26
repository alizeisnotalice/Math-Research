# A057 independent proof audit

Verdict: **ACCEPT_THEORY_AND_CONDITIONAL_BOX_INTERFACE** for the frozen files
identified in `input_manifest.json`. No mathematical correction to the reviewed
theorem or its conditional checker is required. This is an independent audit,
not a claim that the supplied synthetic inputs describe body32.

## 1. Fixed integration variables and finite weighted mass

Let A have exactly the 32 columns (1,s2,...,s6), with all sign vectors occurring
once, and let every scalar Cauchy factor have scale t and unit weight. Choose any
invertible six-column submatrix N and write A=[N,M]. For every t>0, the change of
variables in the six selected scalar coordinates gives the density

    p_t(x) = |det N|^-1 integral_R26
      product_i=1..6 f_t((N^-1(x-My))_i) product_j=1..26 f_t(y_j) dy.

This is a positive convolution integral. It is not a conditional probability on
the zero-measure event AZ=x. The product of the last 26 factors is a probability
density, while each of the first six factors is at most (pi t)^-1. Therefore

    0 < p_t(x) <= |det N|^-1 pi^-6 t^-6.

For I=[a,T] or [a,infinity), with a>0 and T>a when finite, integration proves

    0 < F(x)=integral_I p_t(x)dt
      <= |det N|^-1 pi^-6 /(5a^5) < infinity.

The same argument proves the density is continuous, using bounded domination
locally in x for each t. The strict positivity is pointwise, including at all
projection-zero/collision strata.

Fix a right inverse R once and for all, AR=I. Partition its rows as R_N,R_M,
and define K=(-N^-1 M; I26). Since NR_N+MR_M=I, the translation

    y = R_M x+w

has Jacobian one and transforms the first six coordinates into

    N^-1(x-M(R_M x+w)) = R_N x-N^-1 Mw.

Thus, with P_x(t,w)=product_k=1..32 f_t((Rx+Kw)_k),

    F(x)=|det N|^-1 integral_(I x R26) P_x(t,w) dt dw.

All subsequent derivatives hold w fixed in this final representation. Its
domain and base measure dt dw are independent of x. The x-dependence of the
preceding coordinate substitution creates no omitted derivative of a domain,
Jacobian, or probability law. Infinite base-measure mass is harmless: Tonelli
and the previous bound give integral P_x=|det N|F(x), a finite positive number.

The independent matrix check uses a different legal basis from the author:
the all-positive column and the five columns with one sign flip. It verifies
the affine coefficient identity and AK=0 exactly. Its determinant is -32; no
assertion about a replacement kernel normalization is made.

## 2. Differentiation under the infinite mixture

For l_t(z)=log f_t(z), direct differentiation gives

    l'_t(z)=-2z/(t^2+z^2),
    l''_t(z)=2(z^2-t^2)/(t^2+z^2)^2,
    |l'_t(z)|<=1/t,       |l''_t(z)|<=2/t^2,
    l''_t(z)>=-2/t^2.

The first bound follows from 2|zt|<=z^2+t^2. For the lower second-derivative
bound, multiplication by the positive denominator leaves
2z^2(z^2+3t^2)>=0. The absolute bound follows also from
|z^2-t^2|<=z^2+t^2. The author's sharper upper bound 1/(4t^2) is valid.

Fix x0. For h in any bounded neighborhood, integration of the score along the
coordinate shifts yields

    P_(x0+h)(t,w) <= exp(||Rh||_1/t) P_x0(t,w)
                  <= exp(||Rh||_1/a) P_x0(t,w).

The neighborhood supremum of the exponential factor is finite. In a fixed
direction u, set

    sigma_x = sum_k (Ru)_k l'_t((Rx+Kw)_k),
    tau_x = sum_k (Ru)_k^2 l''_t((Rx+Kw)_k).

Then

    D_u P_x=P_x sigma_x,
    D_u^2 P_x=P_x(sigma_x^2+tau_x),
    |sigma_x|<=||Ru||_1/a,
    |tau_x|<=2||Ru||_2^2/a^2.

Consequently P_x and both of these derivatives have a common local integrable
dominator, a fixed constant times P_x0. For mixed partials one may either use
the corresponding product formulas or polarization in the finitely many
coordinate directions. Dominated convergence establishes F in C^2 locally at
every point and permits both differentiations. This argument includes the
entire t-to-infinity tail; a numerical tail cutoff is neither used nor needed.

Normalize on the fixed space by

    dmu_x(t,w)=P_x(t,w)dt dw/(|det N|F(x)).

This is a probability measure. The score and log-second-derivative bounds
make their first and second moments finite under mu_x. Quotient differentiation
therefore gives the genuine identity

    D_u^2 log F = E_mu_x tau_x + Var_mu_x(sigma_x)
                >= -2||Ru||_2^2/a^2.

The weights depend on x; their derivative is precisely the variance term, not
something that may be discarded before differentiation. Our symbolic check of
the three-mass identity validates this algebra only. The local dominator above
is the proof for the actual infinite mixture.

## 3. The full-frame constant and the six-basis fallback

The diagonal entries of AA^T equal 32. An off-diagonal entry is the sum of a
nonconstant sign character over all 2^5 sign vectors, hence zero. Thus the
single matrix R0=A^T/32 satisfies AR0=I and R0^T R0=I/32. Substitution gives

    Hess log F(x) >= -I6/(16a^2),                  for all x in R6.

For x=c+s u, integrate the directional second derivative with weight 1-s on
[0,1]. Since integral_0^1(1-s)ds=1/2,

    log F(c+u) >= log F(c)+grad log F(c).u-||u||^2/(32a^2).

This factor is consistent with the source code's Q32=I5/32 on x1=1. The matrix
was fixed before choosing u; there is no illicit direction-dependent mixture.

If R=R0+D is another right inverse, then AD=0 and R0^T D=AD/32=0, so

    R^T R=I/32+D^T D >= I/32.

This proves exactly the claimed optimality within the fixed-right-inverse
argument. It does not prove sharpness of the actual Hessian bound. Our distinct
matrix configuration checks this decomposition for a nonzero D=KE with two
rational nonzero entries in E; the general identity is the proof.

The six-factor bound likewise follows by keeping the residual scaled Cauchy
law independent of x and taking b=N^-1u. Its upper bound is
(1/4+6)||b||^2/a^2=25||b||^2/(4a^2), because Var sigma<=E sigma^2 and
||b||_1^2<=6||b||^2. For the source's specified
N=(J6-2I6)diag(-1,1,1,1,1,1), NN^T=4I+2J and its inverse is I/4-J/32.
Restricting to u1=0 gives exactly Q6=I5/4-J5/32. Q6-Q32=(7I5-J5)/32 is
positive definite: its eigenvalues are 1/16 and 7/32. The fallback is valid but
uniformly less favorable for this tangent payment.

The source's centered-box Q6 maximum follows because a positive semidefinite
quadratic attains its box maximum at a corner. At a corner the sum of squares
is fixed; minimizing the squared signed sum gives the stated expression. For
five equal halfwidths the least absolute signed sum is h, giving 39h^2/32.
The 39/5 comparison with Q32 is a curvature-cost ratio only.

## 4. Entire boxes and the shared value denominator

Let d range over a rectangle D in R5. For a fixed positive semidefinite Q,
the established tangent inequality is

    F(c+(0,d)) >= f exp(g.d/f-d^TQd/a^2),

where f=F(c)>0 and g is the matching five-component gradient of F at c. The
exponent is concave in d. Every rectangle point is a convex combination of
vertices; applying concavity and then the minimum proves that every point is
bounded by the lowest vertex exponent. Degenerate sides cause only duplicate
vertices and do not invalidate the argument. Neither chamber convexity nor
projection separation is used.

For interval data f in [L,U], L>0, and independent g_i in [g_i^-,g_i^+], at a
fixed vertex d the smallest g.d is

    b_d=sum_i min(g_i^- d_i,g_i^+ d_i).

The logarithm of the remaining factor to minimize is log f+b_d/f. Its derivative
is (f-b_d)/f^2. If b_d<=L it is nondecreasing on the interval; if L<b_d<U it
decreases then increases; if b_d>=U it is nonincreasing. These statements include
b_d<=0 and both endpoint equalities. Therefore the minimizer is exactly

    f_d=min(U,max(L,b_d)),
    C_d=f_d exp(b_d/f_d-d^TQd/a^2).

The extrema over the compact interval jet and the finite vertex set commute.
Thus min_d C_d is the exact robust minimum of this tangent expression over the
independent interval jet. It is not asserted to be the exact minimum of F itself.
Genuine correlations between value and gradient can improve, never invalidate,
this conservative lower bound.

The shared f in the prefactor and denominator matters. Our new synthetic band
configuration has F interval [2,5]. At its four distinct corners, the frame32
triples (b_d,f_d,exponent) are

| d1,d2 | b_d | f_d | exponent |
| --- | --- | --- | --- |
| -1,0 | -1/3 | 2 | -13/72 |
| -1,7 | 20/3 | 5 | 23/36 |
| 12,0 | 3 | 3 | -1 |
| 12,7 | 10 | 5 | -49/72 |

The global robust minimum is 3/e at the interior choice f_d=3. The independent
exact exponential enclosure establishes 1.103<3/e<1.104. This is an interface
example, not a claimed value of any body32 F.

A pointwise minimum of several Q matrices need not be convex, so its negative
need not satisfy the same vertex rule. The source correctly prohibits that
construction. Taking the maximum of already valid whole-box lower certificates
from separately fixed Q choices remains valid.

## 5. Exact arithmetic, metadata, and proof boundaries

The checker computes d, b_d, f_d, q_d, and the exponent using `Fraction`. Its
top-level rational parser rejects floats and booleans. Arb rational division
and exponential evaluation provide ball enclosures; conversion of the outward
endpoints to rational numbers preserves the endpoints, and acceptance requires
the lower endpoint to exceed or equal the exact target. The minimum of the
corner lower endpoints is a rigorous lower bound. The minimum of their upper
endpoints is an upper bound on this *robust tangent minimum*, not on inf_box F.

The independent exponential check uses S_n(x)=sum_k=0..n x^k/k!, x>=0 and n=200.
For x<n+2, the omitted-term ratios are bounded above by x/(n+2), giving

    S_n(x) <= exp(x)
      <= S_n(x) + [x^(n+1)/(n+1)!]/[1-x/(n+2)].

For a negative exponent, take reciprocals in reverse order. All operations in
this check are exact rational arithmetic. For each of the eight distinct tested
corner exponent/factor pairs, the resulting rational bracket lies inside the
source Arb enclosure at 64 bits. The point-box equality is exactly 2 and correctly
passes target 2 while returning UNRESOLVED at an exact target only 10^-19 larger.

The root model literal and its full normalization string are checked, as are
positive a, finite T>a for B, the absence of T for H, face center/box shape,
interval order, and matching jet model, window and center. Window metadata use
Python dictionary equality, rather than a semantic rational canonicalization;
equivalent rational strings can conservatively mismatch. This is not an
authentication mechanism. The nonempty `source_proof` string identifies a claim
to be checked externally; it does not establish the jet's correctness.

There is no lower-curvature certificate for full G obtained by setting a=0.
The domination and the curvature cost require a>0. Zero/repeated Fourier
projections, S3, S2, ties and zeros do not affect this real-space proof, but a
future numerical jet provider must still evaluate its correct Hermite limits.

Finally, Cauchy scaling gives p_t(r x)=r^-6 p_(t/r)(x), so for r>0,

    H_a(r x)=r^-5 H_(a/r)(x),
    B_(a,T)(r x)=r^-5 B_(a/r,T/r)(x).

Differentiating with respect to x also gives the corresponding gradient scaling
with factor r^-6. A fixed cutoff is not spatially homogeneous. Changing cutoff,
window, factor weights, or model therefore requires matching certified jets or
a proved exact transformation. A band value with an H gradient, full-G jets,
and sparse/augmented-frame jets are not substitutions. Positivity gives G>=F,
but a global g floor still needs a complete paid cover and valid full-G symmetry/
homogeneity transfer. No such coverage or production kernel calculation appears
in this audit.
