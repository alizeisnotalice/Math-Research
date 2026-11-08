# r077: a cap-two integrality barrier for fractional corner assignment

## Statement

Let \(N\ge2\), \(1<\sigma\le2\), and \(\eta>0\).  Put

\[
 P=(1/2,1)^N,\quad a_0=0,\quad a_1=(-\eta,0,\ldots,0),
\]

\[
 R_j(x)=[a_j,x],\qquad
 F=2\mathbf1_{[-\eta,1]\times[0,1]^{N-1}}.               \tag{1}
\]

Consider the finite mixed-integer convex allocation problem with
\(z_0,z_1\in\{0,1\}\), \(z_0+z_1\ge1\), \(0\le G_j\le F\), and the
closed-cell vertex constraints

\[
 \int_{R_j(v)}G_j\ge\sigma |R_j(v)|z_j
 \qquad (j=0,1;\ v\in\{1/2,1\}^N).                       \tag{2}
\]

Here vertices are those of the closure of P. For bounded sources, whole-open-cell averages at least σ extend to these vertex inequalities by continuity of unnormalised rectangle integrals. We write the vertex convention to obtain an explicit finite program.

With objective \(\sum_jJ(G_j)=\sum_j\int\Phi_N(G_j)\), the integer optimum
is exactly

\[
 {\cal I}_{\rm int}=\Phi_N(\sigma).                       \tag{3}
\]

Its continuous relaxation has a feasible point of cost
\(\sigma(1+\eta/2)\).  Consequently its integrality gap satisfies

\[
 {{\cal I}_{\rm int}\over{\cal I}_{\rm rel}}
 \ge{(1+\log\sigma)^{N-1}\over1+\eta/2}.                  \tag{4}
\]

The two anchors are genuinely different.  This is a barrier to using
fractional assignment variables as a physical source-budget upper bound;
it is not a counterexample to a strong-maximal upper bound.

## Integer optimum

If some \(z_j=1\), apply its constraint at \(v=(1,\ldots,1)\).  For j=0,
\(|R_0(v)|=1\), so Jensen gives

\[
 J(G_0)\ge\int_{[0,1]^N}\Phi_N(G_0)
 \ge\Phi_N\left(\int_{[0,1]^N}G_0\right)\ge\Phi_N(\sigma).
\]

For j=1, the rectangle has volume \(1+\eta\), and Jensen gives the stronger
bound

\[
 J(G_1)\ge(1+\eta)\Phi_N(\sigma)\ge\Phi_N(\sigma).         \tag{5}
\]

All other objective terms are nonnegative, proving the lower bound in (3).
It is attained by

\[
 z_0=1,\quad z_1=0,\qquad
 G_0=\sigma\mathbf1_{[0,1]^N},\quad G_1=0.                \tag{6}
\]

Indeed every R_0(x), \(x\in P\), lies in \([0,1]^N\), so its average is
\(\sigma>1\).  Thus (6) also gives strict actual success on the entire
receiver cell, not merely at the vertices.

## Fractional feasible point and relaxation lower bound

In the continuous relaxation take

\[
 z_0=z_1={1\over2},\qquad
 G_j={\sigma\over2}\mathbf1_{R_j(1)}.                     \tag{7}
\]

Since \(\sigma/2\le1\), these sources have
\(\Phi_N(\sigma/2)=\sigma/2\).  For every closed-cell vertex, and in fact
every x in P,

\[
 {1\over|R_j(x)|}\int_{R_j(x)}G_j={\sigma\over2}
 =\sigma z_j.                                             \tag{8}
\]

Thus (7) is feasible, including the whole-cell version of the constraints,
and has cost

\[
 J(G_0)+J(G_1)
 ={\sigma\over2}+\frac{\sigma}{2}(1+\eta)
 =\sigma(1+\eta/2).                                       \tag{9}
\]

Conversely, \(\Phi_N(t)\ge t\).  Applying the two top-vertex constraints
to any relaxed feasible point gives

\[
 \sum_jJ(G_j)
 \ge\sigma[z_0+(1+\eta)z_1]
 \ge\sigma(z_0+z_1)\ge\sigma.                              \tag{10}
\]

The exact relaxed optimum is not needed for (4): (9) is the required upper
bound on \({\cal I}_{\rm rel}\), while (3) is the exact integer value.
The lower bound (10) records that the relaxation is nevertheless a finite,
nondegenerate program.

## Dimension regimes

Write \(\delta=\sigma-1\).  Since

\[
 \log(1+\log(1+\delta))=\delta-\delta^2+O(\delta^3),       \tag{11}
\]

the logarithm of the right side of (4) is

\[
 (N-1)(\delta-\delta^2+O(\delta^3))
 -\log(1+\eta/2).                                         \tag{12}
\]

For the three indicated scales this gives, respectively,

\[
 \begin{array}{c|c}
 \delta & \log\hbox{ gap lower bound}\\ \hline
 N^{-1} & 1-\log(1+\eta/2)+o(1)\\
 (N\log N)^{-1/2} & \sqrt{N/\log N}+o(\sqrt{N/\log N})\\
 \sqrt{\log N/N} & \sqrt{N\log N}-\log N+o(\log N).
 \end{array}                                               \tag{13}
\]

The last two cases are valid for all sufficiently large N, when
\(\sigma\le2\).

## Interpretation

The fractional construction has local average only \(\sigma/2\), generally
below one.  Its feasibility comes solely from the artificial scaled
constraints \(\sigma z_j\).  It therefore cannot be read as a physical
successful-source allocation.  The example proves that any finite-grid
method which relaxes binary corner choice into this type of bare fractional
budget can underprice the true integer source allocation by the gap (4).

This does not rule out legitimate fractional-cover formulations.  In this
single-cell example, use the closed perspective convention at \(z_j=0\) and
the perspective cost

\[
 \sum_j z_jJ(G_j/z_j),\qquad 0\le G_j\le z_jF,\qquad
 \int_{R_j(v)}G_j\ge\sigma|R_j(v)|z_j.                    \tag{14}
\]

For every positive z_j, Jensen applied to \(G_j/z_j\) gives

\[
 z_jJ(G_j/z_j)\ge z_j|R_j(v)|\Phi_N(\sigma).
\]

Since both top-vertex volumes are at least one and
\(\sum_jz_j\ge1\), the perspective relaxation costs at least
\(\Phi_N(\sigma)\).  The integer solution (6) attains this value.  Thus the
perspective exactly repairs this toy example.

More generally, a legal fractional cover keeps each real source column
genuinely successful at its own threshold and charges it linearly, for
example by \(\sum_jw_jB_N(\tau_j)J(G_j)\), with coverage weights used only
in the union or integration step.  One may not lower that source's success
threshold to \(\sigma w_j\) and still charge \(J(w_jG_j)\) as if it were the
same physical column.  The present counterexample excludes that latter
bare scaling move, not all fractional source-allocation or cover methods.
