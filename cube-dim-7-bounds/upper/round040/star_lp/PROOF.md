# R040 · A laminar obstruction to the common-star index-cover LP, and containment compression

## Status and scope

This note proves one exact obstruction and one exact replacement interface. It does **not** prove a new upper or lower bound for the centered cube maximal-operator constant C_n. All cubes are open, centered, axis-parallel, and have one common radius in every coordinate; measure is Lebesgue measure. Boundaries have no effect on the statements. The input class remains arbitrary L1 and all positive scales.

The R039 common-star index-cover LP has optimum `1+d/2` on a depth-d binary laminar family in every dimension. The family maps to an actual original selection history, with every increment strictly positive, whenever `d a^(-n) <= tau`. In particular its optimum can grow as `tau a^n`, despite a simple root-witness payment of the actual selected union. Thus a uniform small K theorem for that exact LP is false in this parameter regime. This is a limitation of its cluster restriction and payment rule, not of the operator.

A finite containment-compression LP is then given. It uses actual increment volumes, has at most 2N columns, contains the old pointwise load certificate, and compresses the whole tree into one root column. On this example its optimum equals the exact selected-union/common-mass ratio.

## 1. The precise R039 LP being tested

For a finite actual history put

\[
 Q_i=Q(y_i,r_i),\quad S_i=aQ_i,\quad V_i=|S_i|,
 \quad I_i=Q_i\setminus\bigcup_{h<i}Q_h,\quad v_i=|I_i|.
\]

A column is a nonempty set J with nonempty common open intersection of its stars S_i, together with a chosen maximum-radius member j(J). For nonnegative coefficients theta_J the constraints are

\[
 \sum_{J\ni i}\theta_J\ge1\quad(i=1,\ldots,N),\qquad
 \sum_J\theta_J\mathbf1_{S_{j(J)}}\le K\quad\text{a.e.}       \tag{1}
\]

Let K_star be the minimum K. All possible common-star columns are allowed; no pivot grid, center-only pivot restriction, or chosen clustering algorithm is imposed. The R039 arrangement-cell interface has this same optimum because it preserves the representative when enlarging a column. Its payment bound is `lambda sum_i v_i <= ((a+1)/a)^n K_star ||f||_1` for every original common witness f.

## 2. Exact abstract laminar calculation

**Theorem 1.** Suppose the stars form a complete binary rooted tree of depth d, with root depth zero, and satisfy:

1. each child star is strictly contained in its parent star;
2. stars from different branches are disjoint;
3. radii strictly decrease with depth;
4. every leaf star has positive volume.

Then the full LP (1) has

\[
                         K_{\rm star}=1+\frac d2.              \tag{2}
\]

**Lower bound, with an explicit nonnegative dual certificate.** Every common-star cluster is an ancestor chain. If its representative has depth k, that node is the highest node of the chain. Put a probability measure nu on the leaf stars by giving each of the 2^d leaves mass 2^(-d), uniformly inside that leaf. A depth-k representative contains exactly its descendant leaves, so

\[
                 \nu(S_j)=2^{-k}.                              \tag{3}
\]

For a node i of depth l define

\[
 p_i=\begin{cases}2^{-l-1},&0\le l<d,\\2^{-d},&l=d.\end{cases} \tag{4}
\]

(When d=0 the second case applies.) For every admissible cluster J whose representative has depth k,

\[
 \sum_{i\in J}p_i
 \le\sum_{l=k}^{d-1}2^{-l-1}+2^{-d}
 =2^{-k}=\nu(S_{j(J)}).                                       \tag{5}
\]

Consequently every feasible theta satisfies

\[
 1+d/2=\sum_i p_i
 \le\sum_J\theta_J\sum_{i\in J}p_i
 \le\sum_J\theta_J\nu(S_{j(J)})\le K.                          \tag{6}
\]

The first inequality uses precisely the item-cover constraints; the last integrates the paid representative load. These are exact rational inequalities and prove the lower bound for **all** fractional clusters. They do not assume symmetry of an optimizer.

**Upper bound.** Write q_root=1 and q_v=1/2 for every other node. For each node v of depth k and each descendant leaf w put one column on the full path from v to w, with representative v and coefficient

\[
                       q_v\,2^{-(d-k)}.                       \tag{7}
\]

The total representative coefficient at v is q_v. A node at depth l receives coverage

\[
 2^{-l}+\frac12\sum_{k=1}^{l}2^{-(l-k)}=1.                     \tag{8}
\]

At any point the representative stars containing it lie on one ancestor chain. Their total coefficient is at most `1+d/2`, with equality throughout each leaf star. This proves (2). The case d=0 is the one-column certificate of value one. No appeal to strong LP duality is needed.

The same proof for an m-ary tree gives `K_star=1+d(1-1/m)`. Only the binary case is needed below.

## 3. Real centered-cube history with strictly positive increments

**Theorem 2.** Fix n>=1, `a>8/3`, `1<=beta<=2`, and tau>0. For every integer d>=0 satisfying

\[
                            d a^{-n}\le\tau,                  \tag{9}
\]

there is an actual original history of `N=2^(d+1)-1` cubes, with all small cubes pairwise disjoint, whose full common-star LP has value (2).

**Construction.** Index nodes by words sigma in {-,+} of length at most d. Let

\[
 R_\sigma=8^{-|\sigma|},\quad y_\varnothing=0,\quad
 y_{\sigma\pm}=y_\sigma\pm\frac{R_\sigma}{2}e_1,
 \quad S_\sigma=Q(y_\sigma,R_\sigma),
 \quad Q_\sigma=Q(y_\sigma,R_\sigma/a).                         \tag{10}
\]

Here S_sigma is exactly aQ_sigma; it has not been changed into a rectangle or uncentered witness. The same construction works in every n, with all other center coordinates zero.

The child-star containment follows from `1/2+1/8<1`; sibling stars are disjoint because their center separation is R_parent and the sum of their radii is R_parent/4. Thus Theorem 1 applies.

More strongly, the inflated stars beta S of siblings are disjoint because their radius sum is `beta R_parent/4 <= R_parent/2 < R_parent`. Also

\[
 \tfrac12R_\sigma+\tfrac\beta8R_\sigma<\beta R_\sigma
 \qquad(1\le\beta\le2),                                      \tag{11}
\]

so every inflated child is contained in its inflated parent. Therefore inflated stars of unrelated nodes are disjoint.

The small parent cube is disjoint from each child star: in the first coordinate their half-width sum is `(1/a+1/8)R_parent < R_parent/2`, by `a>8/3`. All later descendants lie in a child star. Cubes in different branches are disjoint as well. Hence **all** original Q_sigma are pairwise disjoint. In the actual history,

\[
 I_\sigma=Q_\sigma,\quad
 v_\sigma=(2R_\sigma/a)^n>0,\quad
 w_\sigma=v_\sigma/V_\sigma=a^{-n}.                           \tag{12}
\]

**All eligibility tests.** Order the input by decreasing radius, breaking same-level ties arbitrarily, and use the original algorithm with `b=a beta`. At the center of a depth-k node, precisely its k ancestors among earlier nodes have their bQ=beta S containing that center. Unrelated inflated stars are disjoint from its own inflated star, hence cannot contain its center. Every ancestor has already been selected and there are no earlier descendants. Thus its actual eligibility sum is exactly

\[
                                k a^{-n}\le\tau.              \tag{13}
\]

Induction over the order proves every cube stays eligible and is selected. Equal-radius tie choices do not matter. This is a finite instance of the original largest-eligible-radius algorithm, not an abstract admissible overlap graph.

**One original L1 witness.** All S_sigma lie in the root star. For any lambda>0 and epsilon>0 take

\[
 f_\epsilon=(1+\epsilon)\lambda\mathbf1_{S_\varnothing}.
                                                                    \tag{14}
\]

Then every S_sigma is a strict original witness, with integral exactly `(1+epsilon)lambda V_sigma`. The original input is nonnegative, integrable, and has finite support; it is an allowed member of the general input class. No claim that the whole selected union is a strict level set at level lambda is needed.

For the original parameterization `beta=1+1/t`, `t>1`, `a>2t`, impose also `a>8/3`; then all conditions above hold unchanged.

## 4. What is obstructed, and what is actually cheap

At the legal depth `d=floor(tau a^n)` the LP obeys

\[
                 K_{\rm star}=1+\tfrac12\lfloor\tau a^n\rfloor.
                                                                    \tag{15}
\]

Thus no bound independent of a holds for this LP over these legal histories when n,tau are fixed and a is allowed to increase. At fixed admissible a it gives an explicit lower obstruction of this size; it does **not** assert unbounded depth at fixed a while retaining (12). In a regime `tau a^n/(n log n) -> infinity`, this rules out an `o(n log n)` or fixed-proportion `n log n` upper estimate for this **exact K_star certificate**, uniformly over all selected histories.

The obstruction is not caused by zero increments: every v_sigma in (12) is strictly positive. Nevertheless fine-scale indices have very small actual volumes, whereas (1) requires a full unit of index coverage and pays each cluster by the same common-star volume constant regardless of its actual volume.

Let

\[
 A_{n,d}=\sum_{l=0}^d(2\,8^{-n})^l.
\]

The actual selected volume and the exact common-witness mass cost are

\[
 U=\sum_\sigma v_\sigma=(2/a)^n A_{n,d},\qquad
 \mathfrak m=V_\varnothing=2^n.                               \tag{16}
\]

For the second equality the root constraint gives the lower bound, and `g=1_(S_root)` attains it because every other star lies inside the root. Thus

\[
                  U/\mathfrak m=a^{-n}A_{n,d}.                \tag{17}
\]

This stays below `(4/3)a^(-n)`, for all n>=1 and d, while K_star grows linearly in d. There is also a direct cheap payment: the actual selected union lies in S_root, so `lambda U <= lambda V_root <= integral_(S_root) f`. Neither a common pivot for the full family nor a common-star cluster is necessary for that payment. Statements (15)--(17) are selected-union/certificate comparisons, not operator lower bounds.

## 5. A restricted finite repair: containment compression with actual volume costs

This interface keeps explicit geometric columns rather than returning to the unrestricted exact witness-mass LP.

For every index j define the **containment column**

\[
 H_j=\{i:r_i\le r_j,\ Q_i\subset S_j\}.                       \tag{18}
\]

It includes j. Containment can be checked directly by

\[
               \|y_i-y_j\|_\infty+r_i\le a r_j.               \tag{19}
\]

Include these N columns and the N singleton columns `{j}`, each anchored at j. Duplicates can be removed. A column h=(J,j) has the exact actual-volume cost

\[
                  c_h=\frac{\sum_{i\in J}v_i}{V_j}.           \tag{20}
\]

For containment columns, the real increments are disjoint and lie in S_j, so `0<=c_h<=1`. For singleton columns, `c_h=v_j/V_j<=a^(-n)`. Neither statement replaces the actual increments by independently chosen pieces.

Consider the finite LP

\[
 \min K,\quad\theta_h\ge0,\quad
 \sum_{h:i\in J_h}\theta_h\ge1\ (\forall i),\qquad
 \sum_h c_h\theta_h\mathbf1_{S_{j_h}}\le K\quad\text{a.e.}     \tag{21}
\]

The last inequalities are imposed on the finitely many positive-volume arrangement cells of the original stars. Thus (21) has at most 2N columns and no unspecified pivot selection or measurability requirement. Zero-cost columns are allowed.

**Theorem 3 (original-instance payment).** Every feasible certificate in (21), for every original common witness f, satisfies

\[
                          \lambda U\le K\|f\|_1.             \tag{22}
\]

**Proof.** The index coverage, the exact costs (20), and the original witness inequalities give

\[
 \lambda U
 \le\sum_h\theta_h\lambda\sum_{i\in J_h}v_i
 =\sum_h c_h\theta_h\lambda V_{j_h}
 \le\sum_h c_h\theta_h\int_{S_{j_h}} f
 \le K\|f\|_1.                                               \tag{23}
\]

All summands are nonnegative and finite. This proves the direction needed for an upper certificate. The coefficients depend on the true geometric history, not on independently optimized witness masses. ∎

Taking every singleton with coefficient one recovers exactly the old load `W=sum_i(v_i/V_i)1_(S_i)`. Hence (21) is never worse than that certificate on a fixed history. It is generally a restricted LP; no claim is made that the 2N columns represent every useful volume-payable cluster. In particular enlarging a column changes its cost (20), so arbitrary enlargement is **not** automatically an equivalent transformation.

For the tree, H_root is the whole family. Its cost is exactly `a^(-n)A_(n,d)`. The one-column choice theta_root=1 proves

\[
                  K_{\rm contain}\le a^{-n}A_{n,d}.           \tag{24}
\]

Conversely apply (22) to `lambda=1`, `f=1_(S_root)` with weak witness inequalities (the proof allows them), or to (14) and let epsilon decrease to zero. It gives `K_contain >= U/V_root`. Consequently (24) is equality. Replacing the exact containment cost by the coarser geometric bound one yields a still-valid K=1 certificate on the tree.

Common-star columns can be appended with their R039 cost `((a+1)/a)^n`, or their exact actual-volume cost when computed. The proof (23) remains valid whenever a column's charged cost is at least `sum_(i in J) v_i/V_j`. This observation does not itself give a uniform global bound for (21).

## 6. Limits and remaining proof responsibility

1. The all-n obstruction is rigorous for the **unmodified** R039 common-star/unit-index-coverage LP. It survives arbitrary pivots and fractional columns. It is not a C_n statement and does not refute mass-incidence, containment, or other geometric payment schemes.
2. The fixed-a legal depth limit (9) is explicit. Unbounded d at fixed n and arbitrary a concerns general geometric families; the positive-increment legal-history construction supplies all d only when a is allowed to increase accordingly.
3. The repair (21) is an exact finite payment interface, with a genuine strict improvement on the exhibited tree. A dimension-dependent uniform upper bound over all legal histories is **not proved**. No improvement to the n log n operator upper bound is claimed.
4. A finite numerical diagnostic cannot establish such a uniform theorem. Any useful next argument must bound (21), or a clearly specified extension, using genuine history geometry and actual increment volumes.
5. No production optimization or MATLAB call was used. The supplied small check uses Python standard-library rational arithmetic and checks all columns needed for the exact tree certificates. The proof is analytic; the check is supporting evidence only.

## 7. Sources and actual verification

This note is an internal derivation. No unpublished internal result is described as a publication, and no external literature theorem is needed. Actual read sources and SHA-256:

| Source | SHA-256 |
|---|---|
| `/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round039/joint_mass_dual/PROOF.md` | `3fa6117e186343466f3a1cbeb18ed4084a12f4b303e966f14480687a1edd24f6` |
| `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/upper/round032/joint_cover_theory/proof.md` | `6bfd95576ef401f27c23be1eaf3e7b9b19293883f9b2b524a329246ee7c3fdd2` |
| `/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round037/general_cover/PROOF.md` | `6261341afbdafe868616e3e365f3484e0dab7e9e8e6ca09f83c56c93751df90e` |
| Main repository `AGENTS.md` | `18025c63fdc4225e0c326bc535b810c9a9e54ab683b990dce818758f3ca732d5` |

Independent sparse clone HEAD: `326124055562b85b1df9c9b922328bcce9d6f644`. Its checkout reported an unavailable object in an unrelated excluded project, while `AGENTS.md` and `cube-dim-7-bounds` were materialized and the initial status was clean. Only this owned directory was written. No fetch, rebase, commit, or push was performed.

Reproduce the small exact check with `python3 verify.py` from this directory. The output file records tested parameters, exact rational values, source hashes, interpreter version, and the actual check counts.
