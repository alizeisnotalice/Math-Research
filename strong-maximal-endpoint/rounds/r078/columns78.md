# Fractional source-column covers: an exact finite model

## Scope and standing finite interface

Fix a finite set of eligible oriented anchors \(\mathcal A\), a finite set
of open receiver cells \(\mathcal P\), and one common required average
\(\sigma>1\).  The source is a fixed step function \(0\le F\le2\) on a
finite axis-parallel grid.  The grid contains all coordinates of the
eligible anchors and receiver-cell vertices.  Each eligible anchor/cell
pair has one fixed orientation and satisfies the grid condition of
`r077/allocation77.md`, Section 5.  Thus its whole-cell success constraints
are exactly its finite list of vertex constraints.

Except for the explicitly marked same-anchor lemma, statements below concern
this finite, prescribed geometry.  They do not choose anchors, orientations,
or a finite partition for a continuum halo.

## 1. Columns and their exact costs

For an anchor \(a\) and a nonempty subset \(S\subseteq\mathcal P\) of cells
eligible for \(a\), call \((a,S,G)\) a *source column* if
\[
 0\le G\le F
 \quad\hbox{and}\quad
 {1\over |[a,x]|}\int_{[a,x]}G\ge\sigma
 \quad\hbox{for every }x\in P,\ P\in S.                 \tag{1}
\]
The rectangles in (1) have the prescribed orientation.  Define
\[
 c_a(S)=\min\left\{J(G):G\text{ satisfies (1)}\right\},
 \qquad J(G)=\int\Phi_N(G).                              \tag{2}
\]
Set \(c_a(S)=+\infty\) if there is no such \(G\).  This is a genuine
minimum, not merely an infimum: Section 7 of the cited finite reduction
replaces arbitrary measurable \(G\) by its cell averages without increasing
cost or losing (1), and the resulting finite convex program has a compact
feasible box.  Hence every finite \(c_a(S)\) is an exactly defined finite
convex-program value.

There are only finitely many relevant pairs \((a,S)\), at most
\(|\mathcal A|(2^{|\mathcal P|}-1)\).  A column may use a cost-minimizer in
(2); allowing a non-minimizing actual source produces no cheaper column.

## 2. The fractional column-cover LP

Let \(w_{aS}\ge0\) be a linearly charged weight of the column \((a,S)\).
The complete finite column LP is
\[
\begin{array}{ll}
\text{minimize}&\displaystyle \sum_{a,S}c_a(S)w_{aS}\\[1mm]
\text{subject to}&\displaystyle
 \sum_{a,S:\,P\in S}w_{aS}\ge1,\qquad P\in\mathcal P,\\
&w_{aS}\ge0,
\end{array}                                                \tag{CLP}
\]
where infeasible pairs (infinite cost) are omitted.  Its literal
``actual-column'' version permits a finite collection
\((a_i,S_i,G_i,w_i)\) with each \((a_i,S_i,G_i)\) satisfying (1), charges
\(\sum_iw_iJ(G_i)\), and has the same weighted coverage inequalities.
The two infima are equal.

Indeed, replace every actual \(G_i\) by a minimizer for \(c_{a_i}(S_i)\).
This retains its coverage and cannot increase its charge, proving
\(V_{\rm CLP}\le V_{\rm actual}\).  Conversely, choose a minimizer in (2)
for each positive LP variable.  This realizes that LP solution as actual
columns, proving the reverse inequality.  This equivalence is only for
the *linearly charged fractional-column model*.  A weight \(0<w<1\) is
not, by itself, a physical assertion that the scaled source \(wG\) succeeds.

### Same-anchor fractional consolidation

This lemma does **not** need the finite grid interface.  Let
\(\{H_{a,x}\}_{x\in X}\) be a common family of positive linear receiver
functionals for one fixed anchor \(a\), normalized so success means
\(H_{a,x}G>1\).  It applies, for example, to normalized rectangle averages;
the argument does not require a pointwise or Fourier representation.  Let
\((G_k,w_k)\) be a finite or countable family with \(0\le G_k\le F\), and
put \(E_k=\{x:H_{a,x}G_k>1\}\).  For a receiver \(x\), put
\[
 W(x)=\sum_{k:x\in E_k}w_k.                               \tag{2a}
\]
On \(E=\{x:W(x)\ge1\}\), define the common source
\[
 h=\min\left\{\Phi_N(F),\ \sum_kw_k\Phi_N(G_k)\right\},
 \qquad G=\Phi_N^{-1}(h),                                 \tag{2b}
\]
where the inverse is on \([0,\Phi_N(F)]\).  The present modular is strictly
increasing there; a right-continuous generalized inverse gives the same
conclusion in the non-strict version.  Thus \(0\le G\le F\) and
\[
 J(G)=\int h\le\sum_kw_kJ(G_k).                           \tag{2c}
\]

Fix \(x\in E\) with \(W(x)<\infty\), and set
\[
 A_x={1\over W(x)}\sum_{k:x\in E_k}w_kG_k.
\]
Convexity and \(W(x)\ge1\) imply, pointwise,
\[
 \Phi_N(A_x)
 \le {\sum_{k:x\in E_k}w_k\Phi_N(G_k)\over W(x)}
 \le \sum_kw_k\Phi_N(G_k),
 \qquad \Phi_N(A_x)\le\Phi_N(F).                         \tag{2d}
\]
Therefore \(G\ge A_x\), while positivity and linearity give
\(H_{a,x}A_x>1\).  Hence \(H_{a,x}G>1\).  If \(W(x)=\infty\), choose a
finite successful subfamily whose weights sum to at least one and apply the
same argument to that subfamily; its right side is bounded by the full
right side in (2b).  Thus one \(G\) succeeds on all of \(E\), with (2c).
This proves exact same-anchor consolidation for an arbitrary (also
non-finite) success set.

For the finite cell LP, take \(x\) to range over all receiver points in the
cells (or, after the grid reduction, their vertex constraints).  If the
required average is \(\sigma\ne1\), use \(H_{a,x}/\sigma\) in the lemma.

In particular, if there is only one eligible anchor, the fractional column
LP has no relaxation gap: all cells have \(W_P\ge1\), so (2b) gives one
source covering them all at no greater charge.  Conversely one such source
is a unit-weight column.  More generally, for one fixed anchor and a target
set \(S\),
\[
 c_a(S)\le\sum_kw_kc_a(S_k)
 \quad\hbox{if}\quad
 \sum_{k:P\in S_k}w_k\ge1\quad(P\in S).                 \tag{2e}
\]
This fractional subadditivity uses the common \(F\), common \(\sigma\),
and one common anchor.  Different anchors cannot be merged this way, since
their successful rectangles differ.

## 3. Dual, pricing, and what a certificate must verify

The LP dual is
\[
\begin{array}{ll}
\text{maximize}&\displaystyle\sum_{P\in\mathcal P}y_P\\[1mm]
\text{subject to}&\displaystyle
 \sum_{P\in S}y_P\le c_a(S)
       \quad\text{for every eligible finite-cost }(a,S),\\
&y_P\ge0.
\end{array}                                                \tag{D-CLP}
\]
Thus any \(y\) satisfying *all* displayed column inequalities gives the
weak-duality lower bound
\[
 \sum_Py_P\le V_{\rm CLP}.                                \tag{3}
\]
For a proposed dual vector, complete pricing at anchor \(a\) is the finite
global check
\[
 \min_{\varnothing\ne S\subseteq\mathcal P}
 \left[c_a(S)-\sum_{P\in S}y_P\right]\ge0.                \tag{4}
\]
Equivalently, because \(y\ge0\), no admissible source may have
\[
 J(G)<\sum_{P:\,G\text{ succeeds on }P}y_P.                \tag{5}
\]
The latter is only a reformulation; it is not licensed by a local search
over \(G\).  Under the standing grid interface, (4) can be checked by
enumerating all subsets and globally certifying the associated finite convex
programs.  If \(\underline c_a(S)\le c_a(S)\) is an inner Fenchel lower
certificate, then \(\sum_{P\in S}y_P\le\underline c_a(S)\) certifies that
column inequality.  A rational \(y\), rational grid data, and the interval
certificate mechanism of `allocation77.md` make this a rigorous finite
check (the logarithmic objective still requires directed bounds).

Likewise, an upper certificate consists of actual step sources satisfying
all selected vertex inequalities, nonnegative weights meeting the cover
rows, and an outward upper bound on their total charge.  A floating pricing
run, a floating inner solve, or checking only a selected column list is not
a full certificate.  In particular, restricting the primal to selected
columns gives a feasible-cover **upper** bound for the full LP, whereas its
restricted dual need not be feasible for (D-CLP) and gives no lower bound
for the full LP.  If every eligible subset is enumerated, standard finite
LP duality applies whenever the cover is feasible.

## 4. Integral cover and its relation to the finite MICP

The integral column cover replaces \(w_{aS}\ge0\) by
\(z_{aS}\in\{0,1\}\).  It may moreover impose
\[
 \sum_Sz_{aS}\le1\qquad(a\in\mathcal A).                 \tag{6}
\]
This restriction costs nothing.  If two selected columns sharing an anchor
use sources \(G_1,G_2\), their maximum is bounded by \(F\), succeeds on
\(S_1\cup S_2\), and has modular cost at most \(J(G_1)+J(G_2)\), since
\(\Phi_N\) is nonnegative and nondecreasing.  Iteration gives one source
per anchor.  The maximum preserves each original threshold; the common level specifies the present uniform pricing model. Differently weighted Hardy charges need a separate comparison.

After eliminating the subset choice, (6) is equivalent to the finite MICP
of `allocation77.md`: give each anchor one step source \(g_a\), select
binary assignments \(z_{aP}\), require each cell to have a selected anchor,
and impose the vertex rows only for its selected pairs.  The exact reduction
there proves equality of the two finite integral models.  The fractional LP
is a relaxation of this integral cover, not automatically a physical
no-reuse source allocation.

There is a simple finite-anchor comparison, but it is not
dimension-uniform.  Suppose \(|\mathcal A|=k\), and let \(w\) be any
feasible fractional cover.  Write
\[
 W_{aP}=\sum_{S\ni P}w_{aS}.
\]
For each cell \(P\), at least one anchor has \(W_{aP}\ge1/k\).  Assign
\(P\) to one such anchor.  At each anchor multiply every one of its column
weights by \(k\), retain its assigned cells, and apply the same-anchor
consolidation above.  The resulting single source at that anchor succeeds
on all of its assigned cells, and all anchor sources together have charge at
most
\[
 k\sum_{a,S}w_{aS}c_a(S).                                 \tag{7}
\]
Thus the finite integral allocation value is at most \(kV_{\rm CLP}\) for
an optimal fractional cover.  Standard set-cover arguments may give other
logarithmic dependencies on \(|\mathcal P|\), but neither estimate yields
a dimension-consistent conclusion when the required anchor count is
uncontrolled.

## 5. What would be needed beyond the finite model

An extension to a compact continuum of anchors and receiver points needs
more than column generation:

1. A compact parameter space after any tail truncation, a stable orientation
   convention, and continuity (or suitable lower semicontinuity) of
   rectangle averages and column costs under the chosen topology.
2. Uniform strict success margins and a common refining-grid argument, so
   that finitely many receiver cells approximate the intended success set
   without changing the capacity \(G\le F\) or the modular cost by an
   uncontrolled amount.
3. A measure-valued cover formulation with tightness and a valid dual
   pairing, plus a separation/pricing theorem that is global over the
   continuum.  Finite sampled pricing cannot replace this.
4. A quantitative source-budget inequality or an explicit global cover
   construction whose constants are dimensionally sharp.  The finite LP
   merely computes the best allocation for its prescribed columns; it gives
   neither a bound on how many anchors/cells are needed nor a way to charge
   a continuum of action-dependent sources without reuse.

Consequently, even an exact finite instance or an exhaustive finite-column
certificate does not by itself prove the desired dimension-consistent
global budget.  The missing statement is a uniform packing/transport or
allocation estimate linking efficient local columns to the full halo.
