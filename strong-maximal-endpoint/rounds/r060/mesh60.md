# Non-uniform outer widths close the suffix escape cheaply

This is a counterexample within a *weighted* variant of the r059 geometric
mesh. It shows that r059's normalization price is specific to equal outer
slabs. It is not a lower bound for any endpoint constant.

Let \(N\) be even, \(M\ge3\), \(0<R\le1/(2M)\), and put

\[
 \alpha=2R,\qquad p_1=1-(M-1)\alpha,\qquad p_i=\alpha\quad(2\le i\le M),
 \qquad Q=1+{NR\over\alpha}=1+{N\over2}.
\]

The outer support has length one and consecutive slab widths \(p_i\).
Set \(r_0=0,r_j=Q^{j-1}\), \(h_j=r_j-r_{j-1}\), and, on outer slab \(i\),
start with \(1_{[0,r_M]}\), lower cell \(i\) by

\[
 {D_i\over h_i},\qquad D_i={Rr_i\over p_i},
\]

and, for \(i<M\), raise cell \(i+1\) by \(D_i/h_{i+1}\). For \(2\le i<M\) the recovery formula below holds; the downshift formula
holds for \(2\le i\le M\):

\[
 {D_i\over h_{i+1}}={R\over p_i(Q-1)}\le {1\over N},
 \qquad {D_i\over h_i}={R\over p_i}{Q\over Q-1}\le {1\over2}+{1\over N}.
\]

For \(i=1\), \(D_1/h_1=R/p_1\le1/2\). The condition \(2MR\le1\)
says precisely that \(p_1\ge\alpha\). Thus all values are nonnegative
(strictly positive for \(N>2\)), and every recovery height is at most
\(1+1/N\).

## Exact weighted rank identity

At inner prefix \(r_j\), direct cancellation gives

\[
 a_{ij}=\begin{cases}1-R/p_i,&i=j,\\1,&i\ne j.\end{cases}\tag{1}
\]

In particular, the weighted full outer average is exactly \(1-R\) at every
grid prefix. The same cancellation occurs inside each cell: its weighted
loss rate is \(R\), so the full-fibre inner average is \(1-R\) at every
anchored prefix, not only at grid points.

For the outer suffix beginning with slab \(l\), every fibre in it is
unchanged on any prefix ending before cell \(l\); hence its best inner
average is one if \(l\ge2\). Conversely, its aggregate loss rate is
nonnegative after that point. Therefore its exact inner maximum is one.
Fixed-inner-prefix outer averages are fractional-linear in the left outer
endpoint, so their maxima are attained at slab endpoints. With the outer
receiver at \(1+\eta\), the global adaptive and common-action values are

\[
 A^*(\eta)={1\over1+\eta},\qquad
 S^*(\eta)=\max\left\{{1-R\over1+\eta},\,{T\over T+\eta}\right\},
 \quad T=1-p_1=(M-1)\alpha.\tag{2}
\]

The first branch is the whole outer support; the second is the largest proper
suffix. This is an exhaustive outer reduction, including fractional slabs.
The branches cross at

\[
 \boxed{\eta_*={RT\over p_1-R}
 ={2(M-1)R^2\over1-(2M-1)R}.}\tag{3}
\]

This is \(O(MR^2)\), unlike the equal-width crossing \(\asymp MR\). At
\(b_0=1+\eta_*\), \(b_0A^*=1\) and \(b_0S^*=1-R\). With
\(b=b_0(1+N^{-2})\), strictly

\[
 bA^*=1+N^{-2}>1,\qquad bS^*=(1+N^{-2})(1-R)<1
\]

once \(R>N^{-2}/(1+N^{-2})\).

## Actual modular budget

Write \(\Phi_N(t)=t(1+\log^+t)^{N-1}\), and normalize the integral by
the area \(r_M\). Direct cell integration is

\[
 \mathcal B_N(bf)=\Phi_N(b)
 +{1\over r_M}\sum_i p_i h_i[\Phi_N(b\ell_i)-\Phi_N(b)]
 +{1\over r_M}\sum_{i<M}p_i h_{i+1}[\Phi_N(bu_i)-\Phi_N(b)],\tag{4}
\]

where \(\ell_i=1-D_i/h_i\) and \(u_i=1+D_i/h_{i+1}\).
The middle term is nonpositive. The weighted recovery-cell mass is

\[
 H={p_1h_2+\alpha\sum_{k=3}^M h_k\over r_M}
 ={p_1(Q-1)\over Q^{M-1}}+\alpha\left(1-{Q\over Q^{M-1}}\right)
 \le\alpha+Q^{-1}\quad(M\ge3).\tag{5}
\]

Thus

\[
 \mathcal B_N(bf)\le\Phi_N(b)+H\Phi_N(b(1+1/N)).\tag{6}
\]

More generally, take \(R=q_N/N\), where
\(q_N=\lceil(\log N)^2\rceil\), and any growing
\(M_N=o(N/q_N^2)\). Then \(2M_NR\to0\) and

\[
 N\eta_*=O\!\left({M_Nq_N^2\over N}\right)=o(1),
 \quad \Phi_N(b)=1+o(1),\quad H=O(R)+O(N^{-1}),
\]

while \(\Phi_N(b(1+1/N))=O(1)\). Hence

\[
 \boxed{\mathcal B_N(bf)=1+o(1)}.\tag{7}
\]

This permits, for example, \(M_N=\lfloor N/(\log N)^5\rfloor\), not only
\(M_N=\lfloor N^{1/4}\rfloor\). The matching lower estimate is elementary:
the unscaled average is \(1-R\), and \(\Phi_N(t)\ge t\) gives
\(\mathcal B_N(bf)\ge b(1-R)=1-o(1)\). Thus (7) is an actual normalized
modular-cost statement, not merely a one-cell lower estimate.

## Why this is not endpoint-dangerous

Embed the two-dimensional input in \(N\) dimensions by multiplying by the
indicator of the remaining \(N-2\) unit coordinates, then diagonally rescale
its support to the unit cube. The strictly normalized input obeys
\(bf\le h_N\mathbf1_{[0,1]^N}\), where

\[
 N\log h_N=N\log\{b(1+1/N)\}\longrightarrow1.
\]

The strong maximal halo is therefore dominated by the cube halo. The direct
orthant decomposition gives the following halo-volume bound (divide by (7) for the score)

\[
 h_N\sum_{k=0}^N{(2N\log h_N)^k\over(k!)^2}
 \longrightarrow I_0(2\sqrt2),\tag{8}
\]

which is bounded, whereas the fixed-action benchmark \(Q_N\) is exponential
in \(N\). This family disproves a *budget-only pointwise-regret* principle,
but supplies no near-endpoint witness. Any useful sharing theorem must retain
a joint efficiency/receiving-score hypothesis.

## Wolfram receipt

The script mesh60.wl directly integrates every cell for three rational finite
rows. It verifies (1), the full-prefix identity, the grid suffix upper bounds,
the crossing in (2)--(3), and the strict scaled separation. Its output is
mesh60.json. The endpoint reduction is analytic; finite enumeration is a
check of its algebra, not a substitute for it.
