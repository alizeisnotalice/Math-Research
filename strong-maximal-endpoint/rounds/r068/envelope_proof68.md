# Fixed-action conditional-band envelope: a Fubini bridge

## Statement

Fix one block dimension \(b\ge2\), one \(c>0\), and a compact
conditioning band
\[
I=[a_0,a_1]\Subset(0,\infty).
\]
Let \(K_c^s\) be the complete conditional joint-label convolution kernel
on \(H_b\times\{\pm1\}^b\).  Its exact conditional envelope mass is (proved by the nested-shell calculation in radial_envelope68.md)
\[
Z_I=\int\sup_{s\in I}k_{c,s}(y,J)\,dy\,d\#J
 =1+\log{H_{b,c}(a_1)\over H_{b,c}(a_0)},                 \tag{1}
\]
where
\[
H_{b,c}(s)=\sum_{k=0}^b{b\choose k}
 { (s-ck)_+^{b-1}\over (b-1)!}.                            \tag{2}
\]
The measure \(dy\) is the convolution-coordinate Lebesgue measure on
\(H_b\), and \(d\#J\) is counting measure.  With Euclidean Hausdorff
measure or normalized label probability, both factors change consistently
and the value of \(Z_I\) is unchanged.

Let \(G(w,X,\sigma)\ge0\) be measurable, put
\[
\mathcal J(G)=\int_{\mathbb R}e^{-w}
 \|G_w\|_{L^1(H_b\times\{\pm1\}^b)}\,dw<\infty,             \tag{3}
\]
and define, for a receiver radial coordinate \(u\),
\[
A_u(X,\sigma)=\int_{\mathbb R}{\bf1}_{\{w-u\in I\}}e^{-w}
 [K_c^{\,w-u}G_w](X,\sigma)\,dw.                           \tag{4}
\]
Then
\[
\boxed{\quad
\int\sup_{u\in\mathbb R}A_u(X,\sigma)\,dX\,d\#\sigma
 \le Z_I\,\mathcal J(G).\quad}                             \tag{5}
\]
Consequently the same estimate controls any measurable receiver selection
\(u=u(X,\sigma)\), since \(A_{u(X,\sigma)}\le\sup_uA_u\).

## Proof

For every fixed \(X,\sigma,u\), positivity and the condition \(w-u\in I\)
give
\[
\begin{aligned}
A_u(X,\sigma)
&\le\int_{\mathbb R}e^{-w}
 \sup_{s\in I}[K_c^sG_w](X,\sigma)\,dw.
\end{aligned}                                               \tag{6}
\]
The right side is independent of \(u\).  Take the supremum over \(u\),
integrate in \((X,\sigma)\), and use Tonelli:
\[
\begin{aligned}
\int\sup_uA_u
&\le\int_{\mathbb R}e^{-w}
 \left\|\sup_{s\in I}K_c^sG_w\right\|_1\,dw\\
&\le Z_I\int_{\mathbb R}e^{-w}\|G_w\|_1\,dw.
\end{aligned}
\]
The last step is convolution domination by the kernel envelope:
\[
\sup_{s\in I}K_c^sG_w
\le\left(\sup_{s\in I}k_{c,s}\right)*G_w,
\]
including the finite label convolution.  This proves (5).

The proof deliberately keeps \(w\) as the source variable through the
Fubini step.  It does not assume that \(G_w\) is independent of \(w-u\), or
that one source fibre can be reused unchanged for different receivers.

## Radial dilation covariance

Let the finite-depth radial dilation of a source family be
\[
G^D(w,X,\sigma)=G(w-bD,X,\sigma).
\]
Using \(v=w-bD\) in (4) gives the exact identity
\[
A_u^D(X,\sigma)=e^{-bD}A_{u-bD}(X,\sigma),
\qquad
\sup_uA_u^D=e^{-bD}\sup_uA_u.                              \tag{7}
\]
The budget in (3) has the same covariance:
\[
\mathcal J(G^D)=e^{-bD}\mathcal J(G).                      \tag{8}
\]
Thus (5) is compatible with the finite-depth radial dilation that
preserves the underlying source-to-receiver ratio.

## What this does not yet bridge

This is a fixed-\(c\), normalized conditional-kernel statement.  In the
actual physical averaging formula, the radial variable has its own
density \(f_{b,c}(s)\), coarea/Jacobian factors, and possibly a
receiver-dependent action.  Those factors are not represented merely by
the \(e^{-w}\) source weight in (4).  In particular, (5) does not by
itself:

- replace the actual \(f_{b,c}(s)\)-weighted radial average;
- permit \(c=c(u,X,\sigma)\) without a joint action envelope or another
  policy argument;
- build a common physical source profile or control reuse across multiple
  receiver cones;
- prove a general strong-maximal bound.

For a critical band of width \(O(\log N)\), (1) is only an
\(O(\log N)\) loss when \(B_b(c,s)=\partial_s\log H_{b,c}(s)\) is uniformly
\(O(1)\) on that band.  This is useful only after the missing physical
weight and allocation comparisons are supplied.

## Measurability of the receiver supremum

With finite J(G), almost every source fibre is L1. Set the exceptional fibres to zero. On a compact positive s-band the kernels are uniformly bounded and have bounded support, with boundaries of zero Hausdorff measure. Dominated convergence against each translated G_w proves joint continuity of K_c^s G_w in (s,X). The bound by a constant times ||G_w||1 permits dominated convergence in w, including the moving indicator endpoints (two null w-slices). Thus A_u(X) is jointly continuous in (u,X), for each finite label, and sup over u equals the supremum over rational u. In particular it is measurable and the displayed integrals and arbitrary measurable receiver selections are legitimate. Alternatively the explicit positive envelope convolution in the proof supplies the same measurable dominating density without selecting u.
