# R010: exact chamber splitting and pulling triangulation from the complete ray catalog

2026-10-27. Sol theory only. This owner runs no rank, cell, face, determinant, quadrature or MC computation and no MATLAB. The separately accepted R009 arrangement catalog is the input; its nearly two-million6-subset enumeration is not rerun here. No numerical release follows from this file alone. No H, g, I or C7 certificate is claimed from structural output.

## 1. Input and coordinates

Use R009's complete catalog of primitive nonnegative w rays and its frozen36 integer normals: seven coordinate walls e_i, followed by29internal prefix forms c_j. The64 body directions and their prefix normals are unchanged. Write

    xi=Mw, xi_i=sum_(k=i)^7 w_k, det M=1,
    N(w)=lambda(w)=sum_(k=1)^7 k w_k=sum_i xi_i>0,
    u=diag(1,...,7)w/N(w), sumu=1.

All catalog rows, their exact witnesses, incidences and source hashes must be verified before use. The actual accepted R009 input has46rays; the proofs below apply to any complete finite catalog. Global vertex order is its canonical primitive-w lexicographic order. Store all cells/faces as exact catalog-index bitsets.

## 2. The catalog contains every intermediate vertex

For a partial sign vector sigma=(sigma1,...,sigmak), k<=29, let

    P_sigma={w>=0, lambda(w)=1, sigma_j c_j·w>=0 for j<=k}.

Every vertex of this bounded polytope has six independent active homogeneous normals chosen from its seven walls and firstkcuts. The R009 complete6-subset catalog therefore contains that vertex direction, even though the other29-k cuts have not yet been imposed. Conversely every compatible catalog point belongs to P_sigma. Thus

    P_sigma=conv{w_i/N_i: catalog i compatible with sigma}.        (1)

Some catalog points can be redundant interior points at intermediate levels; this causes no problem for(1). The equality follows by including all true vertices, not by assuming every compatible point is already a vertex.

Let I_sigma denote that entire compatible index set. Normalising positive columns by N_i does not change their linear rank. For a nonempty set on lambda=1, its affine dimension is linear rank minus1: differences have lambda0, while any base point has lambda1. Consequently

    P_sigma has dimension6 iff rank{w_i:i in I_sigma}=7.           (2)

No floating LP, edge intersection or convex-hull routine is required for this feasibility/full-dimension test.

## 3. Finite incremental sign splitting proves complete chambers

Start with the full catalog bitset I_empty, representing the standard lambda simplex. At levelk+1 form

    I_plus={i in I_sigma:c_(k+1)·w_i>=0},
    I_minus={i in I_sigma:c_(k+1)·w_i<=0}.

Keep precisely the children of rank7. Formula(1) also applies to each child because its true vertices belong to the global catalog. Therefore this filtering is an exact halfspace split of P_sigma.

At least one child has full dimension. If the cut crosses the interior, both do; their relative interiors are disjoint. If it does not cross, the lower-dimensional child is a boundary face of the retained child and need not be separately retained. Thus discarding rank<7 children loses no point from the union of the closed full-dimensional children. Induction through29levels proves full coverage of the closed simplex and pairwise interior disjointness of the final cells. Distinct final cells disagree on at least one strict interior cut sign and meet only on lower-dimensional faces.

At every level, a full-rank bitset uniquely determines its partial polytope. Duplicate bitsets are safely identified; if two records with the same full-rank bitset carry conflicting signs, this is an inconsistency, since no nonzero cut can vanish on a rank7span. Do not silently merge a contradictory certificate.

The final compatible catalog points are exactly the extreme vertices of their cell. Each catalog point has a6-independent-normal witness; all those normals are active valid inequalities for any final sign cell containing it. Hence it is a vertex of that cell. This is stronger than the intermediate statement, and removes any concern about redundant vertices in the final pulling step.

An integer interior witness is simply sum_(i in I)w_i. Every wall or signed cut that is nontrivial on the full-dimensional cell is positive at at least one vertex and nonnegative at all. Its value on the sum is strictly positive. All cuts and walls are nontrivial on a rank7cell. Thus a rational LP interior witness is unnecessary.

## 4. Exact facet identification for every encountered face

Fix a final cell C. For a face represented by a vertex bitset F, let r=rank{w_i:i in F}; its cone dimension isr and its lambda-section dimension isr-1. For each of the original36 normals n form

    F_n={i in F:n·w_i=0}.

The sign of n on C is fixed; zero sets are independent of that sign. Hence cone(F_n) is an exposed face of cone(F). If F_n is proper and has rankr-1, it is a facet. If n vanishes identically on F, its set equals F and is not a proper facet. Duplicate facet bitsets from different normals are deduplicated.

These are **all** facets. At a relative-interior point of a proper facet G of F, some original inequality that is not identically zero on F must become equality; otherwise that point would be relative interior to F. Its zero face is a proper face of F containing G, so by maximality of a facet it equals G. This proves completeness without a numerical hull or a volume argument.

Every face encountered recursively remains an intersection of C with original zero hyperplanes. It has exactly the subset of C's extreme rays on it, so its vertex bitset and rank fully specify its geometry. A global memo table keyed by bitset is valid across cells: the same bitset has the same cone, and the same zero intersections/ranks. Conflicting cached face data are an error.

## 5. Positive pulling recursion, coverage and disjointness

Define Pull(F) recursively using the same global vertex order:

- If rank(F)=1, F contains one primitive ray (otherwise the catalog has duplicate collinear rays). Return that singleton simplex.
- Otherwise choose a=min(F). Find every facet G of F as in Section4 that does **not** contain a. For each simplex S in Pull(G), return S union{a}.

Every returned r-ray simplex is independent. By induction S spans its rankr-1facet, while a is outside that facet's supporting span. No Steiner ray is introduced. Canonicalise simplex tuples by sorted catalog indices; duplicates must agree in provenance and are not double-counted.

For coverage, take the compact lambda-section of F. Rays from its vertex a through generic interior points exit on exactly one facet not containing a; facets/ridges where the exit is nonunique have lower dimension. Coning that facet from a therefore partitions the section up to its boundary, and the recursive facet triangulation partitions it further. Closedness gives coverage of boundary points as well. Intersections of distinct full-dimensional returned simplices have empty relative interior, because their generic exit facet or recursive facet simplex differs. This proves coverage and interior disjointness inductively, not merely a matching volume.

A common vertex order is the usual pulling rule. The recursion is defined on the geometric face bitset, so the same face is triangulated identically wherever encountered. The construction gives compatible face restrictions by the same induction (faces containing a inherit cones over their opposite facets; faces not containing a use their own recursion). For the Fourier application, the proved positive full-dimensional coverage/interior disjointness would already suffice even if only that property were needed.

No two distinct final cells contain the same independent7-ray simplex: such a simplex has nonempty7-dimensional cone interior and cannot lie in the closures of two cells with opposite cut signs. Thus global top-simplex tuples can be uniquely assigned to their cell; a duplicate across distinct cells is an error, not something to hide by a volume sum.

## 6. Exact volume normalisation is a strong independent check

For a top simplex let W be the7x7 matrix of primitive w columns, N_i=lambda(w_i), and

    U=diag(1,...,7)W diag(1/N_i).

All columns of U have coordinate sum1. Its absolute determinant is exactly the volume ratio of their6-simplex to the standard u simplex: subtract one column from the other six and use the sum-of-rows functional to obtain the ordinary6-dimensional determinant ratio. Hence

    volume_ratio=|det U|=7!|det W|/product_iN_i.                    (3)

The positive triangulation and chamber coverage prove

    sum_all_top_simplices 7!|det W|/product_iN_i=1.                 (4)

Every contribution is positive rational and the sum is checked exactly. The full equality is powerful, but(4) alone would not prove no overlaps/missing cells. Sections2–5 supply that proof. Record per-cell sums as well as the global sum, determinants and denominators.

## 7. Exact psi coefficients and primitive-ray kernel weights

For a final cell with cut signs sigma_j, the w-coordinate exponent is the integer linear form

    ell_C=sum_(35 noncut v)s_v+sum_(j=1)^29 sigma_j c_j.             (5)

Noncut prefixes are nonnegative, so their absolute values equal their forms on the entire ordered cone. Every catalog ray in C must satisfy

    q_i:=psi(Mw_i)=ell_C·w_i>0.                                   (6)

The left side is checked directly as the64-body sum of integer absolute dot products. Because all simplex rays satisfy the same weak signs, equality extends linearly over their positive cone. Store ell_C, all qi and the sign/interior witness; a sample sign guessed in floating arithmetic is not a substitute.

Let xi_i=Mw_i and R_xi=MW. Since detM=1,

    |det R_xi|=|det W|,
    alpha=|det W|/product_iq_i,
    b_i=(x·xi_i)/q_i.                                             (7)

These primitive-ray quantities are the recommended kernel representation. If section-normalised xi_i/N_i are used, both determinant and q_i acquire matching productN_i factors which cancel in alpha; b_i is unchanged. If one uses u columns U instead, the frequency transformation has determinant1/7!, which cancels the7! in(3). Mixing detU with primitive qi would be wrong.

## 8. Signed-permutation expansion and true truncated Green formula

The signed permutation group G=B7 has2^7·7!=645120 elements. Its images of the ordered cone coverR7, with disjoint interiors and measure-zero sorting walls. Positive full-dimensional simplicial cones from Sections2–5 have interiors in the strict ordered chamber, so no nonidentity group element fixes such an interior. Sum the actual transformed phase vectors g xi_i; do not merely multiply a general-x single-cone contribution by the group order or divide by a ray stabiliser.

For t>0, the complete decomposition gives the exact finite identity

    p_t(x)=(2pi)^(-7)sum_(C,S)sum_(g in G)|det W_S|
                         product_(i=1)^7[tq_i-i x·g xi_i]^(-1).    (8)

All terms are retained. The orthogonal group maps preserve psi and absolute determinants. Opposite g images may be paired into2Re if a genuine one-per-opposite-pair index set is used; otherwise no extra2 is inserted. Only the full process is symmetric, not an arbitrary single cone.

For epsilon>0, with alpha_S,b_i(g) from(7),

    H_epsilon(x)=(2pi)^(-7)sum_(C,S),g alpha_S
                 DD6[log(epsilon-i·)](b1(g),...,b7(g)).            (9)

This is R009's verified positive DD6 sign. The logarithm lies in the right half-plane; repeated or zero knots use Hermite values f^(k)(b)/k!=−i^k/[k(epsilon-i b)^k]. All-zero knots yield1/(6epsilon^6). Rational knot equality is exact; interval overlap is not an equality test.

The cone/time/Fourier interchange is justified by integral e^(-epsilon psi)/psi<infinity, since psi>=8|xi|. Thus the finite deterministic H expression is for the original common-time Cauchy process and avoids the shape-expectation error. Its evaluation still requires controlled signed cancellation and complete group phases; no H value is computed in this structural batch.

R001/R003 provide the remaining legal interfaces: H_epsilon<=G, a full-face continuous H lower certificate gives g for original G; a certified H^(7/6) integral plus the whole-face short-time bound pays the original I. Original G retains its whole-space delta finite source/generalL1/all-common-radius comparison. A chamber/triangulation count or exact spectral volume is not a C7 improvement.

## 9. Finite implementation and scope boundary

The proposed registration freezes29bitset-split levels, rank7child tests, exact face/facet rank recursion, top-simplex determinants/qi and exact volume1. Natural finite bounds are2^29sign cells,2^Ncatalog face subsets and binom(N,7)independent top tuples, with N=46 for the accepted input. A sharper general hyperplane-region bound is sum_(j=0)^6binom(29,j)=621616, but no actual cell/face/simplex count is asserted before the run. These bounds are not instructions to allocate all subsets: generate only encountered cells/faces and cache them sparsely.

All29levels are sequential in dependence, but cells within a level and final-cell triangulations are independent and may use15single-thread processes. Shared-face memoisation is an optimisation; if workers produce the same face key their exact facet and pulling data must agree on merge. No old cofactor scan, new catalog ray, automatic Green evaluation or I/g mesh is in this batch.

This owner ends at the proof/finite contract. Numerical release and execution belong to the main owner/Terra, and no such release is implied by the mathematical feasibility argument.
