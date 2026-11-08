# Single-corner Hardy frontier at reference precision

## Definition

Let
\[
 T_N=q_0+\alpha L,\qquad L=\log N,\qquad
 q_T={T_N\over N}-1.
\tag{1}
\]
The optimized cap-two Hardy frontier for a single near-corner charge is
defined by \(\log B_N(\tau_*)=T_N\).  For large \(N\), it is in the
interior branch of `nearcorner76.md`, because
\(q_T\to\mu-1>0\).  Therefore
\[
 \boxed{\ \tau_*=N\{q_T-\log(1+q_T)\}.\ }
\tag{2}
\]
This is a method/reference-scale frontier only.  It is not a finite-
\(N\) certificate for the actual strong-maximal constant.

## Expansion

Write \(p_\mu=(\mu-1)/\mu\) and
\[
 \beta_* =\mu-1-\log\mu.
\tag{3}
\]
The critical location has \(q_0=N\mu-d_N\), with
\[
 {d_N^2\over N}=\mu L+o(L).
\tag{4}
\]
Set \(h=-d_N/N+\alpha L/N\), so that
\(q_T=\mu-1+h\).  Taylor expansion of
\(f(q)=q-\log(1+q)\) at \(q=\mu-1\) gives
\[
 f'(\mu-1)=p_\mu,\qquad f''(\mu-1)=\mu^{-2}.
\tag{5}
\]
Because \(d_NL/N=o(L)\), \(L^2/N=o(L)\), and
\(d_N^3/N^2=o(L)\), (2)--(5) yield
\[
 \boxed{
 \tau_*=N\beta_*-p_\mu d_N+
 \left(p_\mu\alpha+{1\over2\mu}\right)L+o(L).
 }
\tag{6}
\]

The leading term agrees with the linear frontier in `nearcorner76.md`.
The \(d_N\) and \(L\) corrections record the reference target
\(q_0+\alpha L\), rather than a claim about an attained corner net or the
unknown finite-dimensional endpoint constant.

## A finite-reference Taylor remainder

For this defined reference target only, retaining its exact d_N avoids the slow replacement z_N² by log N. Put Delta_N=−d_N+alpha log N. Taylor's theorem also gives

\[
 \tau_*=N\beta_*+p_\mu\Delta_N+\frac{\Delta_N^2}{2\mu^2N}
 +O\!\left(\frac{(\log N)^{3/2}}{\sqrt N}\right).
\]

The error tends to zero, since the third derivative of q−log(1+q) is bounded near mu−1. This refines a method/reference calculation, not the unknown endpoint constant.
