# r077: convex source-allocation baseline and the grouping obstruction

## Scope and model

Fix a bounded cap-two source \(0\le F\le2\), finitely many fixed-corner
rectangles \(R_j\), and use the closed success constraints

\[
 {1\over |R_j|}\int_{R_j}G_j\ge1,\qquad 0\le G_j\le F.    \tag{1}
\]

The closed problem is the natural convex relaxation of strict success:
strict inequalities have the same infimum whenever a feasible solution has
slack, by multiplying slightly below the cap or by a local repair. Define

\[
 {\cal P}_{\rm sep}(F;\{R_j\})
 =\inf\left\{\sum_j\int\Phi_N(G_j):\ (1)\hbox{ holds}\right\}. \tag{2}
\]

This is the cost of assigning one source function to each listed rectangle.
It is a useful baseline, but is not automatically the cost of a globally
optimal allocation: one assigned function can in principle satisfy several
rectangle constraints at once.

For one rectangle R, values outside R are useless. The KKT condition is the
water-filling rule

\[
 G(y)=\min\{F(y),h_\lambda\}\mathbf1_R(y),                 \tag{3}
\]

up to the subgradient convention at the kink of \(\Phi_N\), where
\(h_\lambda\) is a constant selected to make the mean equal to one. This
gives a simple analytic benchmark. If \(F\ge1\) on R, then
\(G=\mathbf1_R\) is the Jensen minimizer and the cost is
\(|R|\Phi_N(1)=|R|\).

If instead one permits a single G to satisfy a group \({\cal I}\) of
constraints, its KKT profile is generally zero outside
\(\bigcup_{j\in{\cal I}}R_j\), and on that union has the form

\[
 G(y)=\min\left\{F(y),
  (\partial\Phi_N)^{-1}\!
  \left(\sum_{j\in{\cal I}}{\lambda_j\over|R_j|}
                      \mathbf1_{R_j}(y)\right)\right\},   \tag{4}
\]

with the usual zero branch when the multiplier is below
\(\partial\Phi_N(0+)\). The overlap levels in (4) show why grouping is a nontrivial combinatorial
allocation problem rather than a consequence of the separately minimized
program (2).

## A two-dimensional cap-two KKT example

Let

\[
 R_1=[0,1]^2,\qquad R_2=[0,3/4]\times[0,1],\qquad
 F=2\mathbf1_{R_2}.                                      \tag{5}
\]

These are the rectangles from the two receivers \((1,1)\) and
\((3/4,1)\), respectively, with the same fixed corner \(0\).
For the separate program, the R_1 constraint can receive mass only from
R_2. Jensen on its area \(3/4\) gives the exact minimum

\[
 \inf_{G_1}\int\Phi_2(G_1)
 ={3\over4}\Phi_2(4/3),                                   \tag{6}
\]

attained by \(G_1=(4/3)\mathbf1_{R_2}\). The R_2 constraint separately has
minimum

\[
 \inf_{G_2}\int\Phi_2(G_2)={3\over4}\Phi_2(1)={3\over4},   \tag{7}
\]

attained by \(G_2=\mathbf1_{R_2}\). Hence

\[
 {\cal P}_{\rm sep}
 ={3\over4}\Phi_2(4/3)+{3\over4}.                         \tag{8}
\]

By contrast, the single grouped source
\(G=(4/3)\mathbf1_{R_2}\) has average one on R_1 and average \(4/3\) on
R_2. Its cost is only

\[
 {3\over4}\Phi_2(4/3)<{\cal P}_{\rm sep}.                 \tag{9}
\]

Replacing \(4/3\) by \(4/3+\varepsilon\), with small positive epsilon,
makes both successes strict while preserving the strict gap in cost.
Thus even two nested rectangles in the same fixed-anchor family can defeat
an unconditional claim that the allocation optimum is additive.
Here the gain is legitimate within one fixed-anchor Hardy family. Across
different anchors, identical values of G do not make their Hardy charges
free: the distinct-anchor count or allocation rule must still be paid for.
Crossing same-anchor rectangles give the same grouping phenomenon through
their common overlap.

## A sufficient separated-copy additivity mechanism

Let \(F=\sum_{k=1}^K F^{(k)}\) have pairwise disjoint support packets
\(Q_k\), and partition the rectangle constraints into families
\({\cal I}_k\). Assume the geometric isolation condition

\[
 R_j\cap\operatorname{supp}F\subset Q_k
 \qquad(j\in{\cal I}_k).                                  \tag{10}
\]

Then a constraint in family k cannot receive useful mass from another
packet. Splitting any feasible assigned source into its restrictions to the
disjoint Q_k neither raises its modular nor weakens any constraint. The
convex allocation problem, even when arbitrary grouping is allowed inside
each family, decomposes as

\[
 {\cal P}(F;\{{\cal I}_k\}_{k=1}^K)
 =\sum_{k=1}^K {\cal P}(F^{(k)};{\cal I}_k).               \tag{11}
\]

For translated copies of one packet and a translation-invariant rectangle
family satisfying (10), this is exactly \(K{\cal P}_0\), while
\(J(F)=KJ(F^{(0)})\). The normalized allocation cost is therefore unchanged.
This is a sufficient condition, not a general property of separated
physical packets: if a constraint rectangle can see two packets, a grouped
source can lower the optimum as in (9).

## Practical implication

The r076 separated-copy argument legitimately preserves a chosen localized
allocation because its translated successful rectangles are isolated from
the other packets. It does not prove that this allocation is globally
optimal among all cross-rectangle groupings. Any future finite-grid source
allocation must either establish an isolation/Carleson condition like (10),
or optimize the grouped KKT system (4); simply summing one source cost per
corner is an upper construction, not a lower characterization.

## Whole receiver-cell and rational-certificate interface

Suppose G is constant on a finite axis-parallel source grid, a receiver cell
P does not cross a source-grid hyperplane in any coordinate, and a fixed
anchor a is outside the open coordinate interval of P in every coordinate.
For a fixed orientation, the signed surplus

\[
 D(x)=\int_{R_a(x)}G-\sigma|R_a(x)|                       \tag{12}
\]

is multi-affine on the closed cell P: each source-grid box contributes a
product of coordinate lengths, each affine in its one receiver coordinate.
The orientation is fixed because no coordinate of x crosses a_i.

The multilinear interpolation formula expresses D at every point of P as a
convex combination of its vertex values. Therefore vertex inequalities
\(D(v)\ge0\) imply D(x)\(\ge0\) on all of P. If \(\sigma>1\), every
nondegenerate rectangle in P has average at least sigma and is hence
strictly successful at the original level one. Vertices where a side length
vanishes may have D=0 and do not cause a problem; they are boundary
degeneracies. This is a valid whole-cell replacement for a point-cloud
check. One must still subdivide whenever a receiver cell crosses a source
grid face or an anchor coordinate.

For a rational LP certificate, the safe directions are as follows.

1. A rational primal source must satisfy all vertex inequalities exactly.
   Mixing a rounded floating primal with a known strict source repairs
   feasibility because the constraints are linear; the mixing amount must
   be bounded from an explicit strict margin. Convex mixing preserves
   \(0\le G\le F\le2\).
2. An upper chord \(U\ge\Phi_2\), evaluated on that rational feasible primal,
   gives a certified upper bound on the true modular.
3. There are two valid lower-bound routes. A standard LP dual requires its
   usual stationarity equations. The present finite-grid certificate can
   instead use a box-residual lower bound, which deliberately does not force
   stationarity.  Write the success constraints as
   \(-Ag\le-b\), the epigraph tangent constraints as
   \(m_k g_l-z_l\le-b_k\), with \(0\le g_l\le F_l\) and
   \(z_l\ge0\).  For rational multipliers
   \(\lambda_r,\theta_{lk}\ge0\) satisfying
   \(\sum_k\theta_{lk}\le w_l\), put
   \[
   d_l=\sum_k\theta_{lk}m_k-\sum_r\lambda_rA_{rl}.
   \]
   The Lagrangian infimum over the explicit box gives the certified bound
   \[
   {\rm LB}=\lambda\cdot b+\sum_{l,k}\theta_{lk}b_k+
             \sum_l F_l\min\{0,d_l\}.                     \tag{13}
   \]
   Indeed the z coefficient is
   \(w_l-\sum_k\theta_{lk}\ge0\), so its infimum over \(z_l\ge0\) is zero,
   while \(\inf_{0\le g_l\le F_l}d_lg_l=F_l\min(0,d_l)\).
   Thus exact post-quantization checks of nonnegativity and
   \(\sum_k\theta_{lk}\le w_l\) suffice for this residual certificate;
   \(d_l=0\) is not required.
4. A bound on the total nonnegative epigraph multipliers by the cell volume
   controls the loss from replacing Phi by a finite tangent/chord envelope,
   provided the uniform envelope gap and all cell weights are recorded
   rationally.

Under these checks, weak duality supplies a genuine rational lower bound and
the repaired primal supplies a genuine upper bound. Without exact
post-quantization dual feasibility, a floating dual is numerical evidence
only.
