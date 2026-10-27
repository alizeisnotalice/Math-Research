# R020 independent audit: psi7-004 tangent coefficient

2026-10-27. **The one-prototype coverage, ten-pair coefficient formula and generic star identity are accepted.** No mathematical correction to the formula is needed. One terminology caveat is recorded below. This is a pure theoretical audit; it does not execute the prototype compiler or any of the three gates. Parent retains responsibility for integrating the independent numerical/source witnesses.

## 1. Frozen source, dependencies and terminology

This review binds to proof SHA256724a8a92c850c709db9101e77e2bff440a068312e644e6ffef09cbe09bdca343, prototype SHA256974ad51bc0eb2a6ba7de569f6db27cda73db9145b7bc9f5cedeb88f70e1645ff and registration SHA2560145563bbe7d2b083b413d94c152f951add61b488fbbbe87c51edc112a6b7f43. The entire proof, literal prototype, contract, sources and receipt were read. Complete R003catalog, R009rays and R010ordered cover are accepted prior dependencies; their declared source bytes are hash-verified, not regenerated. The prior R013 auxiliary odd-dimensional sign lapse is explicitly corrected in the new source and in this review.

The selected literal prototype is p=(4,3,3,2,2,1,1), psi(p)=344, v=p/344, with zero rows7,11,12,17,18,33,34 and basis the first six. Its seven-zero classification and psi rate are catalog/finite-compiler inputs. No64row computation, gradient or detL was executed by this auditor.

Terminology: source section1 calls this the smallest 'nonsimple' type. Precisely, it is the first selected class with **more than six zero forms**, beyond the explicit separable six-zero class. A six-zero tangent cone already need not be simplicial: the epigraph of the6D l1norm has twelve extreme generators. Thus the label must not be read as a claim that every six-zero polytope vertex is geometrically simple. This wording caveat does not affect the present template.

## 2. Rank, coordinates and the actual cone

The literal relation z34=z18+z33-z17 is correct by row addition. Rank6 of the displayed basis can also be seen without a numerical determinant. In its nullspace, row differences force h4=h5, h5=h6+h7, h6=h7 and h2=h3. The rows then force h1=4h6 and h2=3h6, so the common nullspace is exactly span(p). This independently licenses the six-coordinate system once g dot v=1 is supplied by the true unit-vertex normalization.

With g the original unnormalized sign gradient and L=(-g;Z0), one has Lv=(-1,0), so a=-x dot v in dual coordinates(a,b)=L^(-T)x. The tangent cone is t>=F(s) with

    F=|s1|+|s2|+|s3|+|s4|+|s5|+|s6|+|-s4+s5+s6|.

Integrating t from F(s) to infinity supplies-1/a in T_v; multiplying by x dot v=-a yields

    C_v=|detL|^-1 integral exp[aF(s)+b dot s]ds.

No further vertex scale or q factor belongs here. The transformation y=(-s4,s5,s6), y4=-y1-y2-y3 has absolute coordinate Jacobian1 and gives F3=sum_(i=1)^4|y_i|, c=(-b4,b5,b6,0). Integration uses dy1dy2dy3 inherited from s, not the induced Euclidean surface measure of the zero-sum hyperplane in R4. Inserting that hyperplane's area factor would give a wrong normalization.

## 3. Positive14to20cone cover and unimodularity

On the zero-sum hyperplane, a nonzero point has a nonempty proper positive-index set I and negative complement J, up to zero-coordinate walls. There are the stated fourteen sign cones. Nonnegative transportation with these marginals writes every such point as a positive combination of roots r_ij=q_i-q_j, i in I,j in J. Conversely every such combination has the required signs. Thus no root cone is missing.

The singleton or co-singleton cases have three independent root generators. For a two-by-two split, label the four roots A=r_ij, B=r_il, C=r_kj, D=r_kl. Their only linear relation is A+D=B+C. Given a nonnegative representation, subtract min(coefficient B,coefficient C) from both off-diagonals and add it to both diagonals. The result lies in cone(A,B,D) or cone(A,C,D). Their common part is cone(A,D); the unique relation rules out an interior overlap. Hence the split is positive and complete, not justified merely by matching a volume.

Every resulting triple is a spanning tree on four vertices. Its reduced incidence determinant is plus or minus1 by deleting a nonroot leaf repeatedly. These are determinants in the actual three coordinates, so each cone has absolute Jacobian1. Eight three-ray sign cones plus two triangles for each of six four-ray sign cones give the complete twenty-cone decomposition. The literal ten representatives are exactly the four singleton cones and two triangles for each of the three two-by-two splits containing index1; their negatives give the remaining ten. The prototype's fixed diagonal and negation rule agree with this proof. No general6Dfan or full group unfolding is needed.

## 4. Independent pairing and every dimensional sign

On a<-max_i|b_i| the positive cone integrals converge. Since F3(r_ij)=2, a representative root triple with d_j=c dot r_j contributes

    1/product_j(-2a-d_j),

and its opposite contributes1/product_j(-2a+d_j). This uses the **negative** of the linear exponent in each denominator, as required for an exponential cone integral. In dimension3 the conversion from positive exponent denominators has sign(-1)^3; it cannot be dropped.

Set u=-2a. Adding the opposite terms algebraically gives

    2u[u^2+d1d2+d1d3+d2d3]/product_j(u^2-d_j^2).

The three independent absolute-value coordinates integrate to

    product_(i=1)^3[-2a/(a^2-b_i^2)]
        =-8a^3/product_(i=1)^3(a^2-b_i^2).

Multiplying these expressions yields exactly

    C_v=32a^4/[|detL| product_(i=1)^3(a^2-b_i^2)]
          *sum_(ten triples)
            [4a^2+d1d2+d1d3+d2d3]
                  /product_(j=1)^3(4a^2-d_j^2).

Thus the factor32, sign, a^4, detL, quadratic scales and single opposite-pair factor are independently verified. The two odd-dimensional minus signs cancel. This is the parent's proposed expression derived as a conclusion, not used as a premise. It is homogeneous of degree-6, as a vertex C coefficient must be.

The genuine tangent extreme generators are the independent-coordinate ones(1,+/-e_i) and the root ones(2,r_ij), with the other block zero. Transportation shows the F3unit ball is the convex hull of r_ij/2; a linear functional with values1at i,-1at j and0elsewhere uniquely exposes that root. The unit ball of the sum norm is the convex hull of its two block unit balls. These facts verify the asserted extreme generators. The displayed linear pole factors are therefore a+/-b_i and2a+/-(c_i-c_j). A general rational evaluation is allowed only when every denominator used by its representation is nonzero. Rational continuation may give signed contributions outside convergence; positivity is not asserted there. At a0 with all displayed factors nonzero the simplified C is zero. Dividing separately by a is neither needed nor licensed by this C-only batch.

## 5. Star identity, stabilizer and absence of missing factors

All coordinates of v are positive. Its sole physical stabilizer choices independently swap the equal pairs(2,3),(4,5),(6,7), hence the stated eight elements; no sign flips or swaps of unequal magnitudes fix v. Every B7orbit has one ordered-chamber representative. Therefore every global simplex incident to this v arises from an ordered simplex containing its unique ordered ray and a member of this stabilizer. Different full-dimensional images have disjoint interiors, so there is no extra divisor from a simplex stabilizer.

At a true vertex, any simplex containing v must have v as one of its own vertices, by extremality in P. Its tangent cones cover the actual pointed tangent cone even if the global signed-permutation cover is not face-to-face. For one simplex the seven denominator directions are-v and Q u_r-v for its other six boundary vertices. Its determinant is the accepted alpha. The ordinary7D exponential integral has(-1)^7 in positive-exponent denominator notation; multiplying by x dot v cancels the node0factor. Six remaining difference signs give

    alpha/product_(r!=r_v) x dot(v-Q u_r).

Summing the incident ordered simplices and eight stabilizer images is exactly the source C_star. No2, factorial, group order, q or normalization factor is missing. Both template and star are integrals of the same true tangent cone on the strict exposing domain. Their equality on this nonempty open set proves rational identity. Agreement at three gates alone would not establish that identity; the gates test the compiler implementation only.

The claimed bound8*212on source terms is a bound from the accepted input size, not a measured incidence count. Actual selected incidence, detL, g and gate values remain the separate finite compiler's responsibility. This auditor did not generate or check those numerical witnesses.

## 6. Generic denominator and physical-point contract

For each literal beta, every absolute component is strictly less than1. The physical point x=L^T(-1,beta)=g+sum beta_i z_i therefore has a=-1,b=beta and x dot v=1 exactly, without floating conversion. Since F(s)>=sum|s_i|,

    aF(s)+b dot s<=-(1-max|beta|)sum|s_i|<0 for s!=0.

The remaining positive tdirection is also strictly negative in the exponent, so x strictly exposes v on its full tangent cone. Thus every star difference x dot(v-Q u_r) for r!=r_v is positive, and every template quadratic denominator is positive: |b_i|<1 and |c_i-c_j|<2. No adaptive pole search is required. The precise genericity is for the representations used, not for every unrelated full-frequency ray projection.

The fixed three gates must use exact physical Fraction coordinates from this L, not the literal dual beta as if it were physical x. They should record actual term counts, all signs and exact equality C_template=C_star>0. The theorem predicts these outcomes but is not execution evidence. Their separate compiler must still verify the catalog rate/zero set, fixed64row gradient, matrix data, cone witnesses and actual gate totals as authorized.

Strict accepted range: this one prototype's rational coefficient formula, coverage/sign/scale proof, generic denominator contract and rational identity with the accepted R010star. Not obtained: all21templates, full H/G, epsilon0 Green evaluation, derivatives, globalg, an I budget or C7bound. The three coefficient gates do not test any such claim. No new mathematical numerical/symbolic software, sample, gate, phase/DD, MATLAB, source edit, commit or push occurred in this audit. Source hashes and metadata receipt accompany the report; the audit stops here.
