# R009: exact ordered-cone arrangement rays and conditional Green interface

2026-10-27. Sol theory only: this owner executes no mathematical script, cofactor enumeration, numerical quadrature, MC or MATLAB. The main owner may release only the finite ray-catalog batch below. No chamber enumeration, triangulation or Green/g/I certificate is claimed from a ray list.

## 1. Fixed process and integer coordinates

The original characteristic exponent is

    psi(xi)=sum_(v in V)|v·xi|,
    V={(1,s2,...,s7): sj in{-1,+1}}.

Work on the fundamental ordered cone xi1>=...>=xi7>=0. Put

    w_i=xi_i-xi_(i+1)>=0, xi8=0,
    xi_i=sum_(k=i)^7 w_k.

The map xi=Mw is integer upper triangular with diagonal1, so det M=1 and its inverse is the integer difference map. Define lambda(w)=sum_(k=1)^7 k w_k=sum_i xi_i. A nonzero nonnegative w has lambda>0. The section lambda=1 is a compact6-simplex. Its standard-simplex coordinates are u_k=k w_k; sumu=1.

For v define its prefix walk s_v(k)=sum_(i=1)^k v_i. Then

    v·xi=sum_k s_v(k)w_k=sum_k[s_v(k)/k]u_k.                       (1)

Integer w normals are therefore preferable for exact kernels. If one integrated in u coordinates, the frequency Jacobian would be1/7!, not1; converting rays to xi or using w with det M=1 avoids that factor.

## 2. Exactly35 noncut forms and29 distinct internal cuts

All prefix walks start at0 with first step+1. A form in(1) cuts the relative interior of the simplex iff some prefix coefficient is negative. Its first coefficient is+1; hence a negative coefficient gives both signs and an interior zero by a segment between nearby interior points of opposite sign. If every prefix is nonnegative, the form is strictly positive on the interior and does not cut it.

Reflection gives an independent count: nonnegative length7 walks ending at7-2j, j0,...,3, number binom(7,j)-binom(7,j-1). Summing telescopes to binom(7,3)=35. They all automatically start with+1 and belong to V. Thus64-35=29 internal cuts.

A constructive first-hit classification independently lists all29, and fixes their order in `fan_normals.json`:

- First hit−1 at step3: prefix+--, arbitrary4-step tail,16 cases.
- First hit at step5: prefix++--- or+-+--, arbitrary2-step tail,8 cases.
- First hit at step7: the five walks+++----,++-+---,++--+--,+-++---,+-+-+--.

Tails are lexicographic with−1 before+1; the two step5 templates and five step7 templates are taken in the displayed order. The resulting prefix normals are the29 explicit integer rows in the JSON. Facet normals e1,...,e7 come first, rows0..6, followed by these29, rows7..35. This is the frozen36-row ordering.

Every internal normal begins with1. If two linear hyperplanes in w have the same nonempty section on lambda=1, their restricted affine linear functions are nonzero scalar multiples. Their difference then vanishes on the entire affine section; since it is linear on R7, it must vanish identically (a multiple of lambda that also vanishes at lambda=1 has coefficient0). Thus their normals are proportional. The first coordinate1 fixes the scalar to1, and consecutive prefix differences recover the same v. No distinct internal forms have the same whole section. None is a facet normal because an internal form has both positive and negative entries.

This does not forbid coincident restrictions on a **boundary face**. Those create degeneracies handled by exact rank and primitive-ray deduplication; one must not delete globally distinct normals because their traces agree on a face.

A discarded nonnegative-prefix form causes no missing boundary subdivision: on w>=0 its zero set is exactly the coordinate face where every w_k with s_v(k)>0 vanishes. Its active normal is then a linear combination of those active coordinate-wall normals. Thus the7walls plus29internal forms describe the complete restricted arrangement, including its boundary.

## 3. Complete vertex/ray enumeration from six active normals

Let N be the frozen36x7 integer row matrix. Consider every6-subset I of its rows, exactly

    binom(36,6)=1,947,792.

For its6x7 matrix N_I, the signed cofactors

    z_j=(-1)^j det(N_I with column j deleted), j=0,...,6,

satisfy N_I z=0. If all are0, rank<6 and this row choice does not define a unique ray. Otherwise rank=6. Accept iff all z are nonnegative or all nonpositive, with at least one nonzero; orient to nonnegative and divide by the gcd of nonzero absolute entries. Mixed-sign kernels are rejected. Zero components are retained.

The accepted primitive w uniquely represents a ray. Store xi=Mw, lambda=sum k w_k and the exact normalized section point xi/lambda (or u_k=k w_k/lambda). Since M and M^(-1) are integer, gcd(w)=gcd(xi); a primitive w gives primitive xi. No extra coordinate rounding or arbitrary scale is needed.

**No accepted extraneous ray.** For an accepted w>=0, take generic positive points approaching w that avoid all29hyperplanes. A finite sign subsequence belongs to one full-dimensional arrangement chamber and its closure contains w. The chosen6 independent rows are active valid boundary inequalities of that chamber. On lambda=1 they isolate w/lambda, so it is a vertex of that chamber; equivalently w is an extreme ray of its cone. A row need not be an irredundant facet by itself for active-rank sufficiency.

**No missing arrangement vertex.** Every full-dimensional chamber in the orthant has only the36 possible supporting hyperplanes. At a section vertex, the active homogeneous normals have rank6: together with lambda they must have rank7, and they annihilate its nonzero ray so cannot have larger rank. Choose any6 independent active rows. Their subset occurs in the enumeration and its cofactor kernel is that ray. This also covers boundary/degenerate vertices with more than6 active normals. A rank<6 row subset cannot itself certify a unique vertex, but some rank6 subset exists at every vertex.

Deduplication is by the oriented primitive integer w, not by floating normalized coordinates or an arbitrary incidence subset. Many6-subsets can produce the same degenerate ray. Keep all36 exact dot products/zero incidence (or a mask), the signs of all64 original body forms, and at least one rank6witness with a nonzero6-minor.

For an independent finite coverage check, stratify by c coordinate-wall rows among the chosen6:

|c|number|
|---:|---:|
|0|475020|
|1|831285|
|2|498771|
|3|127890|
|4|14210|
|5|609|
|6|7|

The formula is binom(7,c)binom(29,6-c), and the sum is1,947,792. Wall rows allow reduction to a(6-c)x(7-c) kernel after setting their coordinates0; a correct implementation can exploit this without changing the finite set.

## 4. Additional rays are not the old spectral-vertex catalog

The arrangement includes the sorting walls. For example w=e1, hence xi=(1,0,...,0), is certainly one of its rays, obtained from six coordinate walls. But no body dot product vanishes there, so psi is locally one linear function. Its point xi/psi(xi) lies in a relative facet interior of{psi<=1}, not a vertex. More generally odd-length prefix axes give no zero body dot product. Thus the arrangement must contain rays not classified as true spectral-unit-ball vertices; it is not the earlier21-orbit catalog under a new name.

A useful optional exact label is r_body=rank{v: v·xi=0}, using **all64** body normals, including discarded noncut forms that can vanish on a wall. The true spectral-body vertex criterion is r_body=6. Indeed active supporting gradients at xi/psi(xi) span one nonzero gradient with nonzero dot product against xi, plus the zero-body-normal span; their rank is1+r_body. A polytope vertex requires rank7. This label must not be inferred from the rank6of wall+cut normals.

The signed-permutation group has order2^7·7!=645120. Sorted nonnegative xi is a unique ray-orbit representative, but stabilisers matter for lower-dimensional faces. No catalog orbit size or true-vertex count is asserted before the new enumeration.

## 5. What a future complete chamber certificate must add

The ray catalog alone is **not** a list of full-dimensional fan cones and is not a Fourier decomposition. For each full-dimensional sign chamber one still needs a rational exact H/V description, an interior witness, boundary incidences and a positive triangulation of its6-dimensional lambda-section using catalog vertices. A common pulling order can produce compatible triangulations, but the actual simplices/coverage must be certified.

An incremental construction starting from the simplex and splitting by the29cuts can prove full coverage inductively. Volume sums are useful diagnostics but do not alone prove coverage or disjointness. All chamber inequalities must have constant weak signs on every simplex ray, and the simplices' interiors must partition the chamber. No automatic2^29sign-vector search or triangulation is included in the released ray batch.

On a sign chamber, psi=sum_(noncut)L_v+sum_(cut)sign_v L_v is exactly linear. After producing a complete positive simplicial fan, generate its true signed-permutation images to cover R7 up to measure-zero boundaries. A single cone cannot be assigned full group symmetry or merely multiplied by645120 for a general x. Only the full process after the sum is invariant.

## 6. Exact finite Fourier kernel once those inputs exist

Suppose a certified full simplicial fan has cones C=R[0,infinity)^7, with7 independent frequency rays r_i in one sign chamber. Let q_i=psi(r_i)>0. Positive linearity on that chamber gives psi(Ra)=sum_i a_iq_i. For t>0, absolute Fourier integrability yields

    p_t(x)=(2pi)^(-7)sum_C |det R|
                         product_(i=1)^7(tq_i-i x·r_i)^(-1).       (3)

Keep every term; individual real parts may be negative. Nonunimodular determinants are not dropped. Scaling a ray by a positive factor changes its determinant and denominator by cancelling factors. If rays were generated in w coordinates, use xi=Mw in the phase and detR_xi=detR_w. If using u rays, the1/7!Jacobian must be included.

For epsilon>0, integrate all positive times from epsilon. The double absolute frequency/time integral is

    integral e^(-epsilon psi(xi))/psi(xi) dxi<infinity,

because psi>=8|xi|, giving radial power r^5 near0 and exponential decay at infinity. Thus cone/time/Fourier exchanges are justified despite their complex phases. With

    alpha_C=|det R|/product_iq_i,
    b_i=x·r_i/q_i,

one gets

    H_epsilon(x)=(2pi)^(-7)sum_C alpha_C
               DD6[log(epsilon-i·)](b1,...,b7).                   (4)

For distinct knots, partial fractions have coefficients i^(-6)/product_(k!=j)(b_j-b_k)=−1/product(b_j-b_k). The upper logarithmic boundary at infinity cancels, so the integral has the **positive** DD6 sign. At repeated knots use Hermite divided differences, with

    f_epsilon^(k)(b)/k!=−i^k/[k(epsilon-i b)^k], k>=1.             (5)

The logarithm uses the right half-plane epsilon>0, so zero/repeated real knots are safe. For all b=0, DD6=f^(6)(0)/6!=1/(6epsilon^6), agreeing with integral_epsilon^infinity t^(-7)dt. This verifies the dimension7 sign and time power analytically.

Pairing opposite cones gives2Re of the paired expression. If all cones are retained already, no further factor2 is inserted. On rational x and rays, equal knots are decided exactly before confluent evaluation; overlapping interval balls do not prove exact equality. Signed cancellation is retained, and arithmetic interval error is separate from any short-time error.

Formula(4) is a finite deterministic expression for the **true** truncated kernel if, and only if, the complete fan/triangulation inputs have been certified. It removes a shape-expectation approximation at that stage; it does not make those missing combinatorial inputs free.

## 7. Analytic and original-G transfer interfaces

The original Fourier expression for H_epsilon is globally C_infinity and extends to |Im x|<8epsilon. Its kth directional derivative has the deterministic bound

    Gamma(6+k)/[120pi^4·8^(7+k)epsilon^(6+k)].

Within that tube every cone denominator has positive real part for t>=epsilon, so the log branch in(4) remains compatible. Certified AD/interval evaluation plus these global bounds can support a whole closed face cover or polynomial certificate, including all original S3 locations.

A continuous H_epsilon face lower certificate gives a valid lower g for the original G, since H_epsilon<=G and G is homogeneous of degree−6. The original full-G finite source and all-common-radius/general-L1 comparison remain the legal transfer; one must not substitute H_epsilon's spread-out source LH_epsilon=−p_epsilon for the original delta source.

If a future continuous computation gives B>=integral_[0,1]^6H_epsilon(1,y)^(7/6)dy, R001's whole-face short-time bound gives

    I_G<=[B^(6/7)+E_epsilon]^(7/6),
    E_epsilon=(7/8)2^(9/7)pi^(-13/7)epsilon^(8/7).

Alternatively a separately certified full-G I bound may be reused with the H-derived g. No Minkowski shape budget is intrinsically required by the fan method, but a full g and I certificate still needs substantial continuous arithmetic. No fan count, spectral volume, truncated node sample or algebraic cone formula alone is a C7 bound.

For epsilon decreasing to0, only after global cone cancellation may one extend the formula to every offS3 point. Some individual cone time integrals have divergent zero-knot terms although true G is finite there. R001's local p_t=O(t) offS3 supplies the limit/continuity interface. On S3 the true G is infinite. The safe first computational target is epsilon>0, not an unproved assignment to individual zero-knot limits.

## 8. Status

The exact35/29classification, frozen36normal list, ray completeness/deduplication proof and conditional Fourier/DD6 interface are established here analytically. This owner has not enumerated a single candidate combination or produced full chambers/triangulations. The accompanying finite registration authorises only the separately released1,947,792-combination ray catalog, then stops for main review. No spectral catalog count, H value, g, I or C7 result is claimed.
