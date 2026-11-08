# Coarse uniform screen for critical vector policies

## Result and scope

Fix one global anchor chart and let every coordinate action satisfy
\[
 c_i(x)\in[c_-,c_+],\qquad c_\pm=c_r\pm K\rho,
 \qquad \rho=\sqrt{\log N/N}.
\tag{1}
\]
Let \(W=[N\mu-7\sqrt{N\log N},N\mu+7\sqrt{N\log N}]\).  For every
nonnegative source cost \(G\), the conditional vector kernels obey the
pointwise domination
\[
 K_{\mathbf c}^{s}G
 \le \exp\{O_K(\sqrt{N\log N})\}\,K_{c_-}^{s}G,
 \qquad s\in W.
\tag{2}
\]
Consequently fixed-anchor Fubini gives the global coarse screen
\[
 \int j_W^{\mathbf c(x)}G(x)\,d\nu(x)
 \le |W|\exp\{O_K(\sqrt{N\log N})\}\int G.
\tag{3}
\]
Here the action policy can be arbitrary measurable; neither a reanchor
Jacobian nor injectivity is used.  At the exponential threshold
\(J_{\rm ref}e^{aN}\), Markov gives
\[
 \nu\{x:j_W^{\mathbf c(x)}G(x)>J_{\rm ref}e^{aN}\}
 \le Q_{\rm ref}\exp\{-aN+o(N)\}\int G.
\tag{4}
\]

This is a fixed-global-anchor statement.  It does not handle arbitrary
moving anchors, and its exponential threshold is much too coarse to replace
the \(J_{\rm ref}N^\varepsilon\) screen needed for the efficient core.

## 1. Positive-kernel comparison

Write the unnormalised conditional slice density as
\[
 H_{\mathbf c}(s)={1\over\Gamma(N)}
 \sum_{J\subset[N]}\left(s-\sum_{i\in J}c_i\right)_+^{N-1},
\tag{5}
\]
and let \(D_{\mathbf c,J}(s)\) be the corresponding positive joint angular
slice.  The numerator of \(K_{\mathbf c}^sG\) is the sum over \(J\) of
the integral of \(G\) on \(D_{\mathbf c,J}(s)\).  Since \(c_i\ge c_-\),
\[
 D_{\mathbf c,J}(s)\subset D_{c_-,J}(s),
\tag{6}
\]
and positivity therefore bounds this numerator by the scalar \(c_-\)
numerator.  On the other hand \(c_i\le c_+\) gives
\[
 H_{\mathbf c}(s)\ge H_{c_+}(s).
\tag{7}
\]
Thus the directions of both comparisons are fixed, and
\[
 K_{\mathbf c}^{s}G
 \le{H_{c_-}(s)\over H_{c_+}(s)}K_{c_-}^{s}G.
\tag{8}
\]

Put \(\Delta=c_+-c_-=2K\rho\).  Term by term in (5), using \(|J|\le N\),
\[
 H_{c_+}(s)
 \ge H_{c_-}(s-N\Delta).                                 \tag{9}
\]
The lower endpoint of \(W-N\Delta\) is
\(N\mu-O_K(\sqrt{N\log N})>0\), hence this enlarged window remains in
the macrocompact region for all large \(N\).  The r068 logarithmic normalizer
estimate there gives
\[
 \left|\partial_s\log H_{c_-}(s)\right|\le C_K.
\tag{10}
\]
Integrating (10) over an interval of length \(N\Delta=O_K(\sqrt{N\log N})\)
and applying (9) proves
\[
 {H_{c_-}(s)\over H_{c_+}(s)}
 \le\exp\{O_K(\sqrt{N\log N})\},
\tag{11}
\]
which inserted in (8) is (2).

## 2. Fixed-anchor integration and the high-cost charge

For the fixed scalar action \(c_-\), the conditional kernels are positive
Markov kernels in the fixed global anchor chart.  Fubini therefore gives
the standard fixed-action \(L^1\) estimate after integration over a radial
window.  Combining it with (2) gives (3); the factor \(|W|\) is the coarse
window-length bound.  Applying Markov's inequality to (3) at level
\(J_{\rm ref}e^{aN}\), and using \(Q_{\rm ref}=J_{\rm ref}^{-1}\), gives
(4), since \(\log|W|+O_K(\sqrt{N\log N})=o(N)\).

The canonical optimizer source of `optimizer_review75.md` has, on its
critical rectangles,
\[
 {j_W^{\mathbf c}\Phi_N(F_N)\over J_{\rm ref}}
 \ge\exp\{\nu N-O(\sqrt{N\log N})\},\qquad \nu>0.
\tag{12}
\]
Choose any fixed \(0<a<\nu\).  Those receivers belong to the charged set
in (4), so this entire high-cost canonical family is globally negligible
even under arbitrary critical vector actions.  This is why the geometric
collapse example is not an endpoint-upper obstruction.

The conclusion is intentionally limited: (3) does not reach a polynomial
threshold \(J_{\rm ref}N^\varepsilon\), so it does not settle the remaining
low-cost efficient policy/anchor packing problem.
