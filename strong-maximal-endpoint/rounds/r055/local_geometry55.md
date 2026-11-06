# r055 local geometry: exchanged nearby receivers

Scope: regular cap-two two-equal-block micro budgets only. Take b>=4, zero-extend the kernels, and differentiate on a fixed quadrant large enough to contain all nearby supports; then the required two derivatives have no moving-boundary delta term. This is not a general
statement about \(C_N\).

Let \(N=2b\), \(d\xi=e^{-t_1-t_2}dt_1dt_2\), and
\(R_{b,c}(t)=e^tf_{b,c}(t)\).  For fixed \(d,a_0>0,h\), set
\[
u_\pm=(d\pm a_0\epsilon,d\mp a_0\epsilon),\qquad
c_\pm=(c\pm h\epsilon,c\mp h\epsilon),
\]
and use the r030 receiver kernels
\[
R_\pm(t)=e^{t_1+t_2}f_{b,c_{\pm,1}}(t_1-u_{\pm,1})
 f_{b,c_{\pm,2}}(t_2-u_{\pm,2}). \tag{1}
\]

## Exact pair reduction

For a strict target \(1<q<2\), minimize the cap-two Orlicz cost subject to
both receiver averages being at least \(q\).  If
\(\sigma(t_1,t_2)=(t_2,t_1)\), then \(R_-=R_+\circ\sigma\).  Replacing a
feasible \(H\) by \((H+H\circ\sigma)/2\) lowers convex cost and makes each
average equal to the average of the two old constraints.  Conversely a
symmetric function feasible for \(\bar R=(R_++R_-)/2\) is feasible for
both.  The KKT minimizer for \(\bar R\) is symmetric.  Hence exactly
\[
P_{\rm pair}(\epsilon)=\inf\{\int\Phi_N(H)d\xi:0\le H\le2,
\ \int H\bar R_\epsilon d\xi\ge q\}. \tag{2}
\]
This removes the duplicate-multiplier degeneracy at \(\epsilon=0\).

## Exact curvature comparison

At zero put
\[
\mathcal D_{a_0,h}=h(\partial_{c_1}-\partial_{c_2})
 -a_0(\partial_{t_1}-\partial_{t_2}),
\quad R_0=e^{2d}R_{b,c}(t_1-d)R_{b,c}(t_2-d).
\]
Then
\[
R_+'=\mathcal D_{a_0,h}R_0,\quad R_-'=-\mathcal D_{a_0,h}R_0,\quad
\bar R'_0=0,\quad\bar R''_0=\mathcal D_{a_0,h}^2R_0. \tag{3}
\]
Let \((H_0,\lambda)\) be the regular KKT pair.  Since the average kernel
has zero first derivative, the complete r050 response form does not occur
in its second derivative:
\[
P'_{\rm pair}(0)=0,\qquad
P''_{\rm pair}(0)=-\lambda T_{a_0,h},\qquad
T_{a_0,h}=\int H_0\mathcal D_{a_0,h}^2R_0d\xi. \tag{4}
\]
For the corresponding single moving receiver, translation covariance removes
the angular displacement (its total radial displacement is zero), and r050
gives
\[
P''_{\rm single}(0)=-\lambda\{T_{a_0,h}+B_{a_0,h}\}. \tag{5}
\]
Here, explicitly, \(g_{a_0,h}=\mathcal D_{a_0,h}R_0/R_0\) is the
logarithmic kernel score and
\(B_{a_0,h}=\mathfrak B(g_{a_0,h},g_{a_0,h})\ge0\), where
\(\mathfrak B\) is the r050 lower-coarea plus soft-response bilinear form
on scores (not a form evaluated on the unnormalised kernel derivative).
Thus
\[
\boxed{(\log P_{\rm pair})''_0=(\log P_{\rm single})''_0+
{\lambda\over P_0}B_{a_0,h}.} \tag{6}
\]
For pure translation \(h=0\), the single value is constant, so pair
curvature is nonnegative.  More generally, the nonnegative response term in
(6) cannot be used to cancel the pair's quadratic cost.  The stronger target-one asymptotic corollary from r054 is available only at target \(q=1\): there
\[
(\log P_{\rm pair})''_0\ge Np_\mu h^2+o(Nh^2).
\]
For a strict target \(q=1+\eta\), (2)--(8) remain exact under the stated
regularity, but this displayed r054 asymptotic is not being invoked.

The local symbolic identities are
\[
\ell=\mathcal D_{a_0,h}R_0/R_0,\quad
\mathcal D_{a_0,h}^2R_0/R_0=\ell^2+\mathcal D_{a_0,h}\ell,\quad
\ell\circ\sigma=-\ell.
\]
The chain-rule identities are analytic. The exploration script `local_geometry55.wl` checked the case a_0=1 in Wolfram 15.0; the final portable `code/pair_algebra55.wl` corroborates the general cone expansion and response bookkeeping.

## Pareto cusp: a real but limited two-cone gain

For \(\epsilon>0\), the actual union of the two upper cones has weight
\[
V(\epsilon)=e^{-2d}(2-e^{-2a_0\epsilon}),
\quad\log\frac{V(\epsilon)}{e^{-2d}}
=2a_0\epsilon-4a_0^2\epsilon^2+O(\epsilon^3). \tag{7}
\]
Together with (4),
\[
\log\frac{V(\epsilon)/P_{\rm pair}(\epsilon)}{e^{-2d}/P_0}
=2a_0\epsilon-
\frac12\{8a_0^2+(\log P_{\rm pair})''_0\}\epsilon^2+o(\epsilon^2). \tag{8}
\]
So every regular base point has a sufficiently close distinct symmetric pair
with a strictly better **two-cone** score than its one-cone score.  Taking a
strict target \(1+\eta<2\), then using the coordinatewise-monotone KKT
profile and r026 deep transfer, makes this an actual strict rectangle
witness.  Equivalently, an equality-target construction may be scaled by
\(1+\eta\) only while retaining its cap; target-\(1+\eta\) KKT is cleaner.
This does not compare with the radial scalar construction, which receives a
continuum of angular points.

## Receiver-complexity obstruction

r026 gives for each received cone, of weight \(w_j=e^{-u_{j1}-u_{j2}}\),
\[
J(H)\ge {w_j\over M_bQ_N},\qquad M_b\le C\sqrt{N/(b_1b_2)}. \tag{9}
\]
Indeed restrict to that cone and translate it: \(G(s)=H(u_j+s)\).  Its
translated budget satisfies \(J_G\le e^{u_{j1}+u_{j2}}J(H)\), while the
r026 one-receiver comparison gives \(J_G\ge(M_bQ_N)^{-1}\).  This proves
(9), rather than assuming the factor \(w_j\).
Thus a cone union of volume \(V\) has
\[
V/J(H)\le M_bQ_N E,\qquad E=V/\max_jw_j. \tag{10}
\]
For equal blocks \(M_b=O(N^{-1/2})\); the pair has \(E<2\), hence cannot
recover the half-log balanced-block loss relative to \(Q_N\).

For a Pareto curve \(u(v)=(d(v)+v,d(v)-v)\), \(v\in[-A,A]\),
\(|d'|\le1\), writing \(w=e^{-2d}\), vertical slicing gives
\[
V=\int_{-A}^{A}w(v)(1+d'(v))dv+w(A)
=\int_{-A}^{A}w(v)dv+\frac{w(A)+w(-A)}2.
\]
Hence \(E_{\rm curve}\le1+\Delta v=2A+1\).  For a finite ordered
subfamily, take the **polygonal Pareto interpolation** of its vertices;
on each affine segment the maximum of \(w=e^{-2d}\) occurs at an endpoint.
The curve estimate therefore uses the same \(\max_jw_j\), and the trivial
one-new-region-per-vertex bound gives
\[
E_{\rm finite}\le\min\{K,1+\Delta v\}. \tag{11}
\]
To exceed \(Q_N\) by \(N^\gamma\), a balanced two-block family therefore
needs angular span \(\Delta v\gtrsim N^{1/2+\gamma}\); a finite family needs
both that span and at least this many appropriately spread receivers.

This does not exclude growing-span, non-Pareto, or non-microblock mechanisms.
