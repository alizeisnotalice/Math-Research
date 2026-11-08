# Two fixed-dimension no-cost density reductions

Write
\[
 \Phi_N(t)=t(1+\log^+t)^{N-1},\qquad
 J(f)=\int_{\mathbb R^N}\Phi_N(f),\qquad
 \mathcal H(f)=\{M_{\rm S}f>1\}.                        \tag{1}
\]
The conclusions below are fixed-\(N\) density statements.  They provide no
uniform grid-size, anchor-count, or dimension-sharp budget estimate.

## A. Rational compact step inputs are dense for the cap-two supremum

Let
\[
 C_{N,2}=\sup_{\substack{0\le f\le2\\0<J(f)<\infty}}
 { |\mathcal H(f)|\over J(f)}.                           \tag{2}
\]
The same supremum results if it is restricted to nonnegative compactly
supported functions which are constant, with rational heights, on a finite
axis-parallel rational grid.  This includes the possibility that either
supremum is infinite.

To prove it, fix a cap-two \(f\) with finite modular.  First truncate to
\(f^{(R)}=f{\bf 1}_{[-R,R]^N}\).  Then
\(f^{(R)}\to f\) in \(L^1\), while \(J(f^{(R)})\le J(f)\).  On a rational
grid in \([-R,R]^N\), replace \(f^{(R)}\) by its ordinary cell average
\(u\).  It remains between zero and two, and Jensen gives
\[
 J(u)\le J(f^{(R)}).                                     \tag{3}
\]
Choose a common integer denominator \(d\), increasing sufficiently fast
along the refinement, and put on each cell
\[
 f_{R,d}=d^{-1}\lfloor d u\rfloor.                     \tag{4}
\]
Thus \(f_{R,d}\) has rational heights, satisfies \(0\le f_{R,d}\le2\),
and
\[
 J(f_{R,d})\le J(u)\le J(f).                             \tag{5}
\]
Taking a diagonal sequence of rational meshes, truncation radii, and
denominators (with \(d^{-1}|[-R,R]^N|\to0\)) gives
\[
 \|f_k-f\|_1\longrightarrow0,\qquad J(f_k)\longrightarrow J(f).
\tag{6}
\]
The second conclusion follows also from the fixed-\(N\) Lipschitz
continuity of \(\Phi_N\) on \([0,2]\).  Notice importantly that cell
averaging means \(f_k\) need **not** be pointwise bounded by the original
\(f\); only the cap and the modular inequality (5) are retained.

For every \(x\in\mathcal H(f)\), choose one positive-volume rectangle
\(R_x\ni x\) whose \(f\)-average is \(1+\delta_x>1\).  By (6),
\[
 \left|\int_{R_x}(f_k-f)\right|\le\|f_k-f\|_1
 <\delta_x|R_x|
\]
eventually.  The same fixed rectangle then witnesses \(x\in\mathcal
H(f_k)\).  Hence
\[
 {\bf 1}_{\mathcal H(f)}\le\liminf_k{\bf 1}_{\mathcal H(f_k)}
 \quad\hbox{a.e.},\qquad
 |\mathcal H(f)|\le\liminf_k|\mathcal H(f_k)|.          \tag{7}
\]
Together with (6), this shows that the rational compact-step supremum is at
least the ratio of every admissible \(f\); the reverse inequality is by
inclusion of classes.  No finite-corner description of the halos is needed
in this argument.

## B. Finite rational-anchor approximation of a real-anchor column cover

Here the original capacity \(F\) is fixed and is **not** replaced by a step
approximation.  Let \(E\subset\mathbb R^N\) be measurable with
\(|E|<\infty\).  Suppose a countable family of real-anchor true-success
columns is given:
\[
 (a_j,G_j,w_j),\qquad 0\le G_j\le F,\quad w_j\ge0,
 \quad \sum_jw_jJ(G_j)=B<\infty,                         \tag{8}
\]
and let
\[
 E_j=\{x:H_{a_j}G_j(x)>1\},\qquad
 \sum_jw_j{\bf 1}_{E_j}\ge1\quad\hbox{a.e. on }E.       \tag{9}
\]
The averages in (9) may use all signed orthants of \(a_j\).  For every
\(\epsilon>0\), there are finitely many columns
\((b_\ell,\widetilde G_\ell,\widetilde w_\ell)\), with
\(b_\ell\in\mathbb Q^N\) and \(0\le\widetilde G_\ell\le F\), such that
\[
 \sum_\ell\widetilde w_\ell J(\widetilde G_\ell)\le B,
 \qquad
 \int_E\left(1-\sum_\ell\widetilde w_\ell
   {\bf 1}_{\{H_{b_\ell}\widetilde G_\ell>1\}}\right)_+<\epsilon.
\tag{10}
\]
This is an \(L^1(E)\) coverage-deficit approximation, not an exact finite
rational cover of every receiver point.

### Proof of (10)

Discard zero-weight columns first. Then truncate the *columns*, not the sources.  For
\(W_m=\sum_{j\le m}w_j{\bf 1}_{E_j}\), the functions
\((1-W_m)_+\) decrease pointwise to zero on \(E\) and are bounded by one.
Dominated convergence and \(|E|<\infty\) give an \(m\) for which
\[
 \int_E(1-W_m)_+<\epsilon/3.                             \tag{11}
\]
The retained cost is at most \(B\).

For every retained \(j\), split the source into its anchor orthants,
\[
 G_{j,\eta}=G_j{\bf 1}_{\{\eta_i(y_i-a_{ji})>0\ \forall i\}},
 \qquad\eta\in\{-,+\}^N.                               \tag{12}
\]
Hyperplanes are null sets.  The pieces have disjoint supports, so
\[
 \sum_\eta J(G_{j,\eta})=J(G_j),                         \tag{13}
\]
with no \(2^N\) cost factor.  On the corresponding receiver orthant,
\(H_{a_j}G_j=H_{a_j}G_{j,\eta}\).  Let
\[
 E_{j,\eta}=E\cap E_j\cap
 \{x:\eta_i(x_i-a_{ji})>0\ \forall i\}.                 \tag{14}
\]

Each \(G_{j,\eta}\) is locally integrable (finite modular implies finite
\(L^1\) norm).  Thus its rectangle integral, and hence its anchored average
away from the anchor hyperplanes, is continuous in the receiver endpoint.
The success set \(E_j\cap\{\eta_i(x_i-a_{ji})>0\ \forall i\}\) is
consequently open inside its anchor orthant.  The additional intersection
with the measurable set \(E\) in (14) need not be open, but inner
regularity still gives compact
\(K_{j,\eta}\subset E_{j,\eta}\) with
\[
 \sum_{j\le m,\eta}w_j|E_{j,\eta}\setminus K_{j,\eta}|<\epsilon/3.
\tag{15}
\]
On every nonempty such compact set there are numbers \(d,\delta>0\) such
that
\[
 \min_{x\in K_{j,\eta},i}|x_i-a_{ji}|\ge d,\qquad
 H_{a_j}G_{j,\eta}(x)\ge1+\delta.                       \tag{16}
\]

Choose a rational \(b_{j,\eta}\) sufficiently close to \(a_j\), but
outward from this receiver orthant:
\[
 \eta_i=+\Rightarrow b_{j,\eta,i}<a_{ji},\qquad
 \eta_i=-\Rightarrow b_{j,\eta,i}>a_{ji}.                \tag{17}
\]
For \(x\in K_{j,\eta}\), the new rectangle contains the old one and,
by (16), its volume ratio can uniformly be made smaller than
\(1+\delta/2\).  Positivity gives
\[
 H_{b_{j,\eta}}G_{j,\eta}(x)
 \ge { |R(a_j,x)|\over|R(b_{j,\eta},x)|}
     H_{a_j}G_{j,\eta}(x)>1.                             \tag{18}
\]
Use the finitely many columns
\((b_{j,\eta},G_{j,\eta},w_j)\) with nonempty compact sets.

Let \(W_{\rm new}\) be their actual success coverage sum on \(E\), and put
\(C_K=\sum_{j\le m,\eta}w_j{\bf 1}_{K_{j,\eta}}\).  Then
\(C_K\le W_m\), while \(W_{\rm new}\ge C_K\).  Equations (11), (15), and
\[
 (1-C_K)_+\le(1-W_m)_++W_m-C_K
\]
yield
\[
 \int_E(1-W_{\rm new})_+
 \le\int_E(1-C_K)_+<2\epsilon/3<\epsilon.               \tag{19}
\]
Equation (13) gives the cost assertion in (10).  This proves the claimed
finite rational-anchor reduction while leaving every source piece under the
original capacity \(F\).

## Rational pricing versus real anchors

For a fixed locally integrable source and a real anchor \(a\), rational
anchors approaching \(a\) from a prescribed orthant give convergent
rectangle averages at every receiver not lying on an anchor coordinate
hyperplane.  Hence any strict real-anchor success persists for all
sufficiently close rational anchors.  If a receiver test measure is
absolutely continuous with respect to Lebesgue measure, Fatou therefore
permits rational-anchor pricing inequalities to pass to all real anchors.
This statement does **not** automatically hold for an atomic receiver
measure charging an anchor hyperplane: the required nondegenerate-sided
endpoint continuity then need not be available.
