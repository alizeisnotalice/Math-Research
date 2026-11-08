# Finite grouped source-allocation program and certificates

## Scope

Fix finitely many groups \(j=1,\ldots,J\), each with finitely many receiving
rectangles \(R_{jr}\).  We allocate \(0\le G_j\le F\le2\), making every
assigned rectangle successful, while minimizing \(\sum_j\int\Phi_N(G_j)\).
This is a finite fixed-group problem.  It neither selects a grouping nor
finiteizes a continuum halo.

## 1. Cell program

Choose a bounded box containing all source rectangles and a finite partition
into cells \(C_\ell\), \(\ell=1,\ldots,L\), such that every listed receiver
rectangle is a union of cells.  Set
\[
 v_\ell=|C_\ell|,\qquad
 a_{jr,\ell}={|C_\ell\cap R_{jr}|\over|R_{jr}|}.
\tag{1}
\]
With cell capacities \(f_\ell\in[0,2]\) and strict levels \(b_{jr}>1\),
the allocation program is
\[
\begin{array}{ll}
\text{minimize}_{g_{j\ell}}&
 P(g)=\displaystyle\sum_{j,\ell}v_\ell\Phi_N(g_{j\ell})\\[2mm]
\text{subject to}&\displaystyle\sum_\ell a_{jr,\ell}g_{j\ell}\ge b_{jr},
 \quad1\le r\le m_j,\\
&0\le g_{j\ell}\le f_\ell.
\end{array}
\tag{P}
\]
The capacity is imposed separately for each group, \(g_{j\ell}\le f_\ell\),
not as \(\sum_jg_{j\ell}\le f_\ell\).  The latter is a different no-reuse
model.  The feasible set is compact, so a feasible \(P\) attains its minimum.
Rational endpoints and rational \(f_\ell,b_{jr}\) make all constraints
rational, though the logarithmic objective is not rational.

## 2. Fenchel dual and weak certificates

Define the restricted conjugate
\[
 \Phi_{N,[0,f]}^*(y)=\sup_{0\le t\le f}\{yt-\Phi_N(t)\}.
\tag{2}
\]
For \(\lambda_{jr}\ge0\), put
\[
 y_{j\ell}(\lambda)={\sum_r\lambda_{jr}a_{jr,\ell}\over v_\ell}.
\tag{3}
\]
The Lagrange dual is
\[
 \text{maximize}_{\lambda\ge0}\quad
 D(\lambda)=\sum_{j,r}\lambda_{jr}b_{jr}
 -\sum_{j,\ell}v_\ell\Phi_{N,[0,f_\ell]}^*(y_{j\ell}(\lambda)).
\tag{D}
\]
For every primal feasible \(g\) and nonnegative \(\lambda\),
\[
 D(\lambda)\le P(g).                                    \tag{4}
\]
This follows directly from the definition of the conjugate and the receiving
inequalities.  Thus a primal candidate and a dual candidate always certify
an interval for this finite program.  Strong duality can additionally be
used under relative Slater (for example, a relative-interior feasible point
with strict receiving inequalities); weak duality needs no qualification.

## 3. Explicit \(N=2\) conjugate

For \(\Phi_2(t)=t\) on \([0,1]\) and
\(\Phi_2(t)=t(1+\log t)\) on \([1,2]\),
\(\partial\Phi_2(1)=[1,2]\).  If \(0\le f\le1\),
\[
 \Phi_{2,[0,f]}^*(y)=
 \begin{cases}0,&y\le1,\\ f(y-1),&y\ge1.\end{cases}
\tag{5}
\]
If \(1<f\le2\),
\[
 \Phi_{2,[0,f]}^*(y)=
 \begin{cases}
 0,&y\le1,\\
 y-1,&1\le y\le2,\\
 e^{y-2},&2\le y\le2+\log f,\\
 f[y-1-\log f],&y\ge2+\log f.
 \end{cases}
\tag{6}
\]
The pieces agree at breakpoints; the linear branch is exactly the effect of
the subgradient interval at one.

## 4. Strict certificates from floating candidates

Floating solves only propose candidates.  With rational input data:

1. Round a primal candidate to rational \(g_{j\ell}\in[0,f_\ell]\) and
   verify every constraint exactly.  Strict margins allow sufficiently fine
   rounding; otherwise solve a rational correction LP or report no primal
   certificate.
2. Round multipliers to rational \(\lambda_{jr}\ge0\), which is exact dual
   feasibility.
3. Use outward rational intervals for logarithms to obtain \(U\ge P(g)\).
4. Since conjugates are subtracted in (D), use certified upper bounds for
   them to obtain \(L\le D(\lambda)\).

Then \(L\le P_*\le U\) rigorously.  Branch tests involving
\(2+\log f\) need interval certification; if ambiguous, maximize certified
upper bounds from adjacent branches.  Taylor enclosures or directed interval
arithmetic suffice.  Floating residuals alone do not.

## 5. Receiver-cell certificates beyond point clouds

Assume the source \(G\) is step-constant on an axis-parallel finite grid.
Fix an oriented anchor \(a\) and an open receiver cell \(P\) such that no
coordinate interval of \(P\) crosses a source-grid face and \(a_i\) is not
in the interior of that coordinate interval.  For a desired average
\(\sigma>1\), define
\[
 Q(x)=\int_{[a,x]}G(y)\,dy-\sigma\prod_i|x_i-a_i|.
\tag{7}
\]
On \(P\), the integral is multi-affine in \(x\): after splitting at the
fixed grid faces, every one-dimensional length is affine in its endpoint,
and tensor products preserve multi-affinity.  Hence \(Q\) is multi-affine.
Its minimum over the closed cell is attained at one of the \(2^N\) vertices.
Therefore exact verification
\[
 Q(v)\ge0\quad\text{at all vertices }v\text{ of }P
\tag{8}
\]
implies the average is at least \(\sigma\) throughout the nondegenerate
interior of \(P\).  Degenerate vertices with \(|v_i-a_i|=0\) are harmless:
the unnormalised inequality (7) is still meaningful, while all interior
rectangles have nonzero volume.  This turns each finite open receiver-cell
claim into finitely many exact inequalities, rather than a point sample.

## 6. Bounded measurable sources and finite rectangles

There are two distinct approximation statements.

**Replacing the source by cell averages.**  Given bounded measurable
\(0\le G_j\le F\le2\), average both functions on finite partitions refining
all listed rectangle faces.  Then \(0\le\bar G_j\le\bar F\le2\), every
listed rectangle average is preserved exactly, and
\[
 \int\Phi_N(\bar G_j)\longrightarrow\int\Phi_N(G_j).
\tag{9}
\]
This uses bounded martingale convergence and continuity of \(\Phi_N\) on
\([0,2]\).  Strict margins survive.  It approximates the source/allocation
pair, rather than retaining the original pointwise source.

**Keeping \(F\) fixed.**  Write \(G_j=\theta_jF\), with
\(0\le\theta_j\le1\) and \(\theta_j=0\) where \(F=0\).  Replace
\(\theta_j\) on each cell by its \(F\,dx\)-weighted average.  Then
\(G_{j,\mathcal P}=\theta_{j,\mathcal P}F\) remains between zero and \(F\),
preserves each cell and rectangle integral, and its objective converges by
bounded convergence.  This has finite allocation variables but cell costs
\(\int_{C_\ell}\Phi_N(\theta F)\), not the simple rational program \(P\).

Neither approximation proves a finite description of a full halo, chooses a
global grouping, or supplies the missing global source-budget inequality.

## 7. Exact finite reduction for a step-grid source and finite cell grouping

There is a stronger statement in the genuinely finite setting.  Suppose
\(F\) is already constant on a common finite axis-parallel refinement
\(C_\ell\), and that this refinement contains every source-grid face, every
assigned receiver-cell vertex, and every fixed-anchor coordinate.  For each
open receiver cell \(P\), include the finitely many vertex rectangles
\([a,v]\) in the constraint list.  It is convenient to write their
constraints without division by volume:
\[
 \int_{[a,v]}G_j\ge\sigma_{j,P}|[a,v]|.                  \tag{10}
\]
This also covers degenerate vertices, for which (10) is simply \(0\ge0\).
For nondegenerate vertices it is exactly the corresponding row of \(P\).

Let \(0\le G_j\le F\) be any measurable allocation which succeeds with
average at least \(\sigma_{j,P}>1\) throughout every assigned open cell.
Its rectangle integrals are continuous in the receiver endpoint, because
\(G_j\) is bounded.  Taking limits from the interior of each \(P\) gives
(10) at every vertex.  Now replace \(G_j\) by its ordinary Lebesgue cell
average \(\bar G_j\).  Since \(F=f_\ell\) is constant on \(C_\ell\),
\[
 0\le\bar G_j\le F\quad\hbox{on each cell}.             \tag{11}
\]
Every vertex rectangle is a union of full cells, so all integrals in (10)
are preserved exactly.  Jensen gives
\[
 \sum_\ell v_\ell\Phi_N(\bar g_{j\ell})
 \le\int\Phi_N(G_j).                                     \tag{12}
\]
The vertex test of Section 5 then restores success on the whole open
receiver cell for \(\bar G_j\).

Conversely, any feasible vector of cell values in the finite program lifts
to the step allocation taking that value on \(C_\ell\).  Its vertex
constraints imply whole-cell success by multi-affinity, and its objective
is exactly the finite objective.  Therefore, for this fixed anchor/finite
receiver-cell grouping and this step-grid source, the infimum over arbitrary
measurable allocations equals the optimum of the finite convex program:
\[
 \inf_{\substack{0\le G_j\le F\\\text{whole-cell successes}}}
 \sum_j\int\Phi_N(G_j)
 =\min\{P(g):g\text{ satisfies the vertex rows of the cell program}\}.\tag{13}
\]
This is an exact reduction, not an approximation and not a claim about a
full halo.

In particular, a refinement that includes all receiver faces/vertices (for
example the stated half-grid, “refinement 2”) can use (13).  A coarser grid
that omits some receiver faces (the stated unit-grid, “refinement 1”) does
not meet this hypothesis; it supplies only a restricted discretized
upper/lower model, not the exact equivalence (13).

## 8. Merging sources and the finite anchor-cover model

For one fixed anchor and a common charged threshold, retain at most one
source: if \(0\le G_q\le F\), the pointwise maximum \(G=\sup_qG_q\le F\)
preserves every previously successful rectangle and
\[
 J(G)\le\sum_qJ(G_q).
\tag{14}
\]
For finite families this is pointwise monotonicity of the nonnegative
modular; countable families follow by monotone convergence. The maximum
also preserves differently leveled constraints, but merging differently
weighted Hardy charges requires a separate cost comparison.

Let \(\mathcal A\) be a finite list of anchors, \(\mathcal P\) a finite
list of receiver cells, and choose one common level \(\sigma>1\).
Eligibility means the fixed anchor has one orientation on the open cell
and the partition is aligned as in Section 7. Define, at every closed-cell
vertex \(v\),
\[
 A_{aPv,\ell}=|C_\ell\cap[a,v]|,\qquad V_{aPv}=|[a,v]|.
\]
The exact finite cover model is
\[
\begin{array}{ll}
 \text{minimize}&\displaystyle\sum_{a,\ell}v_\ell\Phi_N(g_{a\ell})\\
 \text{subject to}&\displaystyle\sum_{a:(a,P)\ \mathrm{eligible}}z_{aP}
                     \ge1\quad(P\in\mathcal P),\\
 &\displaystyle\sum_\ell A_{aPv,\ell}g_{a\ell}
                     \ge\sigma V_{aPv}z_{aP},\\
 &0\le g_{a\ell}\le f_\ell,\qquad z_{aP}\in\{0,1\}.
\end{array}
\tag{MICP}
\]
There is no artificial big-M constant: when \(z_{aP}=0\), the corresponding
row is automatic by nonnegativity. Degenerate vertices give \(0\ge0\).
When \(z_{aP}=1\), the vertex test certifies the whole open cell.

This model is exact for the stated finite problem in which each receiver
cell is assigned to an eligible anchor succeeding on its entire interior
and one source is charged per anchor. From any measurable feasible
allocation, merge pieces sharing an anchor by (14), then cell-average as
in Section 7. No cost increases. Conversely every finite feasible vector
lifts to a step source and a valid whole-cell cover. This does not make a
continuum of anchors finite, nor provide a dimension-sharp bound for the
full halo. Barely relaxing the binary variables is not a valid physical
allocation; the integrality and perspective distinction is proved in
[integrality77.md](integrality77.md).
