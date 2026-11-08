# Step sources reduce the noncentred halo to finitely many corners

## Statement

Let \(G\ge0\) be compactly supported and constant on the cells of a finite
axis-parallel grid.  Enlarge the grid, if necessary, by the two faces of a
box containing \(\operatorname{supp}G\), and extend \(G\) by zero outside
that box.  Let \(B_i\) be the finite set of faces in coordinate \(i\), and
let \(\mathcal V=B_1\times\cdots\times B_N\).

For \(a\in\mathcal V\) and a sign vector \(\eta\in\{-,+\}^N\), define on
the open orthant
\[
 D_{a,\eta}=\{x:\ \eta_i(x_i-a_i)>0\text{ for every }i\}
\]
the signed fixed-corner average
\[
 H_{a,\eta}G(x)=
 {1\over\prod_i|x_i-a_i|}
 \int_{\prod_i[\min(a_i,x_i),\max(a_i,x_i)]}G(y)\,dy . \tag{1}
\]
For every threshold \(t\ge0\),
\[
 \{M_{\rm S}G>t\}
 =_{\rm a.e.}
 \bigcup_{a\in\mathcal V,\ \eta\in\{-,+\}^N}
 \{x\in D_{a,\eta}:H_{a,\eta}G(x)>t\}.                \tag{2}
\]
Here \(M_{\rm S}\) is the noncentred strong maximal operator over
positive-volume axis-parallel rectangles.  For \(t<0\), the analogous
statement is also true a.e. by nonnegativity; (2) is stated in its usual
nonnegative-threshold form.

## Proof

Discard the finite union of grid hyperplanes \(\{x_i\in B_i\}\); it has
Lebesgue measure zero.  The inclusion from right to left in (2) is immediate:
the rectangle in (1) contains \(x\).

For the converse, take an off-grid point \(x\) and a positive-volume
rectangle \(R\ni x\) with average greater than \(t\).  Splitting \(R\)
along all coordinate hyperplanes through \(x\) gives at most \(2^N\)
subrectangles, each having \(x\) as a corner.  The average on \(R\) is the
volume-weighted average of their averages, so one subrectangle \(R_0\) has
average greater than \(t\).  Write its endpoint opposite \(x\) as \(b\),
and retain its sign vector \(\eta\).

It remains to move \(b\) to grid faces without decreasing the average.
Fix all coordinates except \(b_i\).  On any grid interval not containing
\(x_i\), the section
\[
 u\longmapsto\int_{\prod_{r\ne i}[\min(b_r,x_r),\max(b_r,x_r)]}
                 G(y_1,\ldots,y_{i-1},u,y_{i+1},\ldots,y_N)\,dy_{-i}
\]
is constant.  Hence the numerator of the rectangle average is affine in
\(b_i\), while its denominator is affine and positive.  The resulting
ratio is fractional-linear, so on that interval its supremum is attained
at an endpoint (or it is constant).  If \(b_i\) and \(x_i\) lie in the same
grid interval, the displayed section is constant throughout the entire
segment between them; numerator and denominator acquire the same length
factor, and moving \(b_i\) to the grid face on its side of \(x_i\) leaves
the average unchanged.  No degeneracy occurs because \(x_i\) is off-grid.

If an endpoint lies beyond the outer support faces, slide it inward to the
nearest such face.  This leaves the numerator unchanged and decreases the
denominator, so it cannot lower the average.  Applying these one-coordinate
adjustments successively gives \(a\in\mathcal V\), with the same orthant
from \(x\), and
\[
 H_{a,\eta}G(x)\ge {1\over|R_0|}\int_{R_0}G>t.          \tag{3}
\]
This proves the reverse inclusion off the grid hyperplanes and hence (2).
Potential degenerate endpoint limits and receivers on source faces are
exactly why the statement is made a.e.; no zero-volume rectangle is used in
the definition of \(M_{\rm S}\).

## Finite count and its limitation

If coordinate \(i\) has \(m_i\) bounded grid cells, then it has at most
\(m_i+1\) relevant faces.  Thus (2) uses at most
\[
 2^N\prod_{i=1}^N(m_i+1)                                \tag{4}
\]
signed corner operators.  This is a finite exact description for a fixed
step source, but its cardinality depends on the grid and can be enormous.
There are only \(\prod_i(m_i+1)\) distinct anchor positions; the factor
\(2^N\) records their signed receiver directions.  If a later source-budget
argument splits a source into those directions, its orthant pieces are
spatially disjoint and their modular costs add to the original cost.  Thus
the operator count in (4) must not be misread as a \(2^N\) modular-cost
duplication.  In either interpretation it supplies no dimension-sharp
strong-maximal estimate by itself.

## Essential qualification for capped optimization

The hypothesis concerns the **actual source \(G\)**.  If a cap \(F\) is
step on a finite grid but an optimized allocation merely satisfies
\(0\le G\le F\), then \(G\) need not be step on that grid.  Its coordinate
sections need not be constant, and the fractional-linear endpoint argument
above fails.  Therefore (2) cannot be invoked merely from a step cap.

For a prescribed finite anchor/receiver-cell problem, `r077/allocation77.md`
gives a different exact cell-averaging reduction: one may replace a feasible
measurable \(G\) by a step allocation without increasing its cost while
preserving the finitely listed vertex constraints.  That result does not
finiteize the full halo, because the finite cells and anchors have already
been chosen.  Using it to justify (2) for an unconstrained global halo would
be circular.
