# R013 independent mathematical audit: true-vertex kernel identity

2026-10-27. **The main identities are accepted with one explicit local sign correction below.** They establish an exact true-vertex representation and sufficient regular-direction epsilon0 formula, conditional on the accepted complete fan/catalog. They do not establish a fast certified implementation, all-exceptional-direction evaluator, global floor, new integral budget or operator bound. This owner performed pure mathematical review and source hashing only; no mathematical numerical/symbolic program was run. Parent retains final proof responsibility.

## 1. Local correction in the frozen text

At proof.md:44 the auxiliary rational sum with denominator product_j(z dot e_j) is called the ordinary exponential cone integral. Literally this omits(-1)^7. For K=E R_+^7 on its convergence domain,

    integral_K exp(z dot h)dh
       =|detE|/product_j(-z dot e_j)
       =(-1)^7 |detE|/product_j(z dot e_j).

The displayed differential-operator sum uses the latter denominator, while the actual simplex coefficient T_v uses v-u=-e. Thus T_v has precisely the correct sign for the ordinary tangent-cone integral. The homogeneous rational cancellation is unchanged by this common nonzero sign. The final T_v,C_v, logarithm and tangent formulas all have the correct signs independently. If the sentence's 'rational sum' was intended to mean the already sign-adjusted T_v sum, it should say so explicitly. This review does not silently accept the literal unsigned identification and does not edit the frozen source.

## 2. Simplex coefficient and global cover

For an affine7simplex the determinant d_Delta, not its volume d_Delta/7!, multiplies DD7exp. Expanding the exponential on the standard barycentric simplex gives degree-m coefficient h_m(lambda)/(m+7)!, exactly DD7exp. Hence the distinct-node expression

    integral_Delta exp(z dot y)dy
       =d_Delta sum_v exp(z dot v)/product_(u!=v) z dot(v-u)

is correct. Coincident projected nodes are removable in the summed simplex expression; individual fractions are not separately defined there. The normalized fan simplex has vertices0,xi_i/q_i and determinant alpha=|detW|/product(q_i), because the primitive coordinate map has determinant1. No further q or factorial belongs here.

The required input is the full positive signed-permutation cover of P={psi<=1}, with disjoint full-dimensional interiors; the ordered simplex is transformed along with its phase. A single ordered cone may not be multiplied by the group order while retaining its phase. These completeness properties are accepted R010 dependencies, not established by the vertex catalog alone. The catalog's21orbits,504868vertices and252434opposite pairs are source-attributed counts; this audit did not recount them.

## 3. Artificial-node cancellation without face-to-face assumptions

Let v be any node of the full cover. Every simplex containing v contributes its geometric tangent cone at v, including simplices for which v is not a vertex. Those finitely many cones cover the tangent cone of P at v up to boundaries and have disjoint interiors. A simplex not containing v stays a positive distance away. Directions in the interior of two local cones would give overlapping simplex interiors sufficiently near v, contradicting the input cover. This argument works even when signed-permutation subdivisions disagree on a sorting wall.

If v lies in a positive-dimensional face of a simplex, that local cone has lineality. Any polyhedral cone can be split into finitely many pointed simplicial cones by coordinate orthants and triangulation. For a pointed simplicial cone K with generator columns e_j,

    product_j D_(e_j) 1_K=|detE| delta_0.

Choose a common product of all generator directional derivatives, with multiplicities, for a finite cone partition. Applying it to the indicator partition gives a polynomial in derivatives of delta_0 whose symbol is the common denominator times the sum of simplicial rational terms. For a cone with a nonzero lineality vector e, D_e1_C=0. Applying D_e and using injectivity of polynomial differential operators on delta_0 forces that rational sum to vanish identically. All identities of indicators are distributional; lower-dimensional boundaries have zero ambient measure. No common convergent Laplace domain for a cone with lines is assumed.

For a pointed cone, Laplace integration on the strict dual convergence domain proves independence of positive subdivisions, with the sign convention in section1. Thus local cones of simplices in which v is not a vertex contribute zero rational continuation. At any nonvertex of P the full tangent cone has lineality, so the remaining simplex-vertex coefficient T_v is zero. This rigorously eliminates artificial ordered-wall nodes and0. Non-face-to-face intersections cause no missing contribution: omitting the non-vertex local cones would have been invalid, but the source includes them and proves their zero contributions.

For a true vertex v, T_v is the rational continuation of its ordinary pointed tangent-cone Laplace integral and has degree-7. Triangulating using only true extreme rays shows that possible poles are among hyperplanes normal to actual incident edge directions. Denominators introduced by other triangulations can be artificial and must cancel exactly. This statement concerns rational functions, not floating residuals of large signed sums.

## 4. Regulated layer integral, sign and normalization

For epsilon>0,

    exp(-epsilon psi)/psi
       =integral_psi^infinity exp(-epsilon r)(epsilon/r+1/r^2)dr.

The corresponding absolute frequency/r integral is finite in dimension7: near0 the frequency density behaves like1/|xi|, and at infinity it is exponentially damped. Fubini is therefore legitimate before vertex splitting. Homogeneity and i^(-7)=i give

    integral_(rP) exp(i x dot xi)dxi=i sum_v T_v(x)exp(i r b_v),
    b_v=x dot v.

The left side is real by central symmetry; its real part is-sum_v T_v sin(r b_v). Its expansion at r0 yields sum_v T_v b_v^m=0 for m<7 and sum_v T_v b_v^7=7!vol(P). In particular sum_v C_v=0 for C_v=b_vT_v.

One must keep the common lower regulator delta until after summing. Integration by parts leaves a boundary sin(b delta)/delta term, whose sum vanishes in the limit, and an integral of cos(br)/r. Subtracting its common cos0 term is justified by sum C_v=0. The elementary regulated identity

    integral_0^infinity exp(-epsilon r)[cos(br)-1]/r dr
        =-1/2 log(1+b^2/epsilon^2)

then yields the **positive** logarithm formula

    H_epsilon(x)=(2pi)^(-7) sum_(v true) C_v(x)
                                  [1/2 log(epsilon^2+b_v^2)].

The logepsilon term cancels. Individual vertex terms cannot be integrated separately to r0 without this cancellation. Coefficients need not be positive; clipping changes the formula. Independently, each order6cone-DD coefficient at v equals b_v times its order7simplex coefficient, since the extra simplex node0 supplies denominator b_v. This gives the same sign and constants.

Central symmetry gives T_(-v)=-T_v and C_(-v)=C_v, hence exactly one factor2 for opposite pairs. The large-epsilon first nonzero coefficient uses sum C_v b_v^6=7!vol(P) and gives H_epsilon(0)=7!vol(P)/[6(2pi)^7epsilon^6], consistent with the accepted fan-origin normalization. Evaluating at x0 requires the limit of the whole expression, not substitution into singular rational coefficients.

## 5. Six-zero and nonsimple tangent cones

At a unit vertex, g=sum_nonzero sign(z dot v)z satisfies g dot v=1. The true tangent cone is

    {h:g dot h+sum_(z dot v=0)|z dot h|<=0}.

With six independent zero forms Z0 and L=(-g;Z0), invertibility follows because g dot v=1 and every zero row annihilates v. In coordinates Lh=(t,s), L^(-T)x=(a,b), one has a=-x dot v. If these are all the zero forms, the tangent cone is t>=sum|s_j|. It is generally **not a simplicial7cone**: its cross-section is the six-dimensional cross-polytope, with generators(1,+/-e_j). The explicit product formula exploits its separable Laplace integral; it is not an assumption of a simple polytope vertex.

On a<-max|b_j|, direct integration gives

    T_v=-64a^5/[|detL| product_j(a^2-b_j^2)],
    C_v=64a^6/[|detL| product_j(a^2-b_j^2)].

The sign, determinant and homogeneity degrees are correct. The original64direction sign gradient and true unit vertex must be kept when changing the six-row basis. At a0 with every b_j nonzero, C_v is zero, not a pole.

For more zero forms, each is c_z Z0 and F(s)=sum_z|c_z dot s| is a norm because it contains the independent basis forms. The exact cone is t>=F(s), giving T_v=-(a|detL|)^(-1)integral exp[aF(s)+b dot s]ds and C_v=|detL|^(-1) times this integral. A positive simplicial fan of linearity for F gives contribution |detU|/[|detL|product_j(-aF(u_j)-b dot u_j)]. These are valid finite source identities and rational continuations; compact templates or efficient certified evaluation have not been produced by the theory alone.

The alternative R010star identity is also correct. For the ordered representative v there is only one ordered representative in its B7orbit, so every global incident simplex comes from an ordered simplex containing v and a group element in the stabilizer of v. Full-dimensional simplex interiors have trivial chamber stabilizer, preventing duplication. The node0 denominator is canceled by b_v, leaving six differences x dot(v-Q u). This does not require the global cover to be face-to-face. Unsimplified star denominators may nevertheless contain artificial poles and many terms.

## 6. Exceptional directions and exact acceptance boundary

For epsilon>0 the actual H is smooth, with holomorphic continuation in the accepted tube. The vertex formula initially holds off a finite union of representation denominator hyperplanes. It determines a unique extension, but this is not a numerical instruction to divide by a zero-containing interval, assign a pole zero or perturb the point. A certified confluent/symbolic limit or the accepted complete DD evaluator is required at exceptions.

If every true-edge projection is nonzero, every T_v is finite. The finite epsilon0 limit exists and MCT of actual H_epsilon↑G identifies

    G(x)=(2pi)^(-7) sum_v C_v(x)log|x dot v|.

When x dot v=0 but T_v is finite, C_v=0 and the regulated term is identically zero for epsilon>0, so its limit is zero. This is different from an uncanceled pole in a chosen template. The simpler contract requiring all vertex projections nonzero is sufficient only together with denominator regularity; neither condition may silently replace the other. If true-edge or artificial denominators vanish, cancellation can occur across terms; a pole does not diagnose G=infinity. The physical singular classification remains the accepted S3 theorem.

Strictly accepted scope: exact regulated vertex identity, artificial-node cancellation, explicit six-zero formula, general norm-fan/star compilation interfaces, and the stated regular-direction trueG limit. Not accepted as completed: all21compact templates, tentative six-zero classification counts, runtime/width improvement from the outer pair count, a new H/G numeric value, a global floor or a C7bound. finite_contract.md correctly leaves those computations unreleased and requires actual exact identity/coverage certificates rather than sample agreement. Pure mathematical review is now complete and stops.
