# r025: angular deep-layer audit for genuine rectangle witnesses

Date: 2026-11-06.  This note concerns the positive orthant; signs only add
the already-present finite rectangle choices.  It supplies a scoped
obstruction for macroscopic angular perturbations and identifies the first
scaling where a genuinely correlated local model survives.  It does not claim
an upper bound for the full strong maximal operator.

Write \(X_i=-\log x_i\), \(s=\sum_iX_i\), and
\(\theta=X/s\in\Delta_{n-1}\).  With ordinary Lebesgue coordinates on the
simplex (of volume \(1/(n-1)!\)),
\[
 dx=e^{-s}s^{n-1}\,ds\,d\theta.\tag{1}
\]
The anchored rectangles used in the hyperbolic construction are genuine
rectangles.  At a receiver whose log vector is \(b\), their input log vector
has the form \(b+Z_{\bf r}\), where \(Z_{\bf r}\) is the familiar
coordinatewise increment vector and \(T_{\bf r}=\sum Z_{{\bf r},i}\).
Other fixed-ratio boxes also have angular freezing under a uniform log-tail,
but the comparison with \(Q_n\) below is asserted only for this anchored
increment family.

## Result 1: macroscopic angles are fibrewise and cannot improve \(Q_n\)

Let \(h_\theta(u)\) be a bounded finite angular partition of nonnegative,
nondecreasing-in-\(u\) profiles with compact \(u\)-support on the left.
For each cell choose bounded \(q(\theta)\) and a fixed rectangle parameter
\({\bf r}(\theta)\), and suppose with a strict margin that
\[
 \mathbb E h_\theta(T_{{\bf r}(\theta)}-q(\theta))>1.\tag{2}
\]
Use the actual input
\[
 f_L(e^{-X})=h_{X/s}(s-L).
\]
For a receiver with log vector \(b=(L-q(\theta))\theta\),
\[
 s(b+Z)-L=T_{{\bf r}}-q(\theta),\qquad
 \frac{b+Z}{|b+Z|_1}\longrightarrow\theta
\]
in probability, uniformly away from angular cell boundaries.  Dominated
convergence and (2) therefore make the receiver ray
\(s\ge L-q(\theta)\) a strict actual witness for all sufficiently large
\(L\), after deleting arbitrarily thin cell-boundary strips.

Formula (1) gives the two exact deep-layer factors
\[
 \int\Phi_n(f_L)dx\sim e^{-L}L^{n-1}
 \int_\Delta J(\theta)d\theta,\qquad
 J(\theta)=\int_\mathbb R e^{-u}\Phi_n(h_\theta(u))du,\tag{3}
\]
\[
 |\mathcal R_L|\sim e^{-L}L^{n-1}
 \int_\Delta e^{q(\theta)}d\theta.\tag{4}
\]
For every individual angular fibre, the defining scalar mechanism bound is
\(e^{q(\theta)}/J(\theta)\le Q_n\).  Integrating gives
\[
 \frac{|\mathcal R_L|}{\int\Phi_n(f_L)}\longrightarrow
 \frac{\int e^{q(\theta)}d\theta}{\int J(\theta)d\theta}\le Q_n.\tag{5}
\]
Thus finite block dependence, for example a profile depending on block sums
\(\sum_{i\in B}s_i/s\), is only a macroscopic angular partition and is
blocked by (5).  The cancellation is geometric: a patch of angular measure
\(p\) reduces both the cost in (3) and the actual receiver volume in (4) by
the same factor.  This applies to real rectangle witnesses, not to a formal
composition of one-dimensional operators.

The assertion is deliberately scoped.  It assumes a fixed angular scale,
bounded logarithmic rectangle displacements, finite angular partition, and
radial monotonicity used to fill the receiver ray.  It does not handle boxes
with exponentially changing aspect ratio in \(L\), nonmonotone radial
profiles, or angular windows that shrink with \(L\).

### Boundary strips and the correct weak-halo upper statement

Take finite cells \(C_1,\ldots,C_m\) in a compact interior simplex region,
with boundary strips of \((n-1)\)-measure \(O(\rho)\).  Suppose each cell
has a bounded profile \(h_j\), with a common left cutoff and finitely many
radial jumps.  At a receiver \(b=s\theta\), \(s\ge L-O(1)\), every anchored
increment law satisfies, uniformly in \({\bf r}\),
\[
 \left\|\frac{b+Z_{\bf r}}{|b+Z_{\bf r}|_1}-\theta\right\|
 \le C|Z_{\bf r}|/L
\]
on \(|Z_{\bf r}|\le L/4\).  Removing \(\rho_L=L^{-1/2}\) strips loses only
an \(O(L^{-1/2})\) fraction of the receiver volume.  A cell switch then
requires \(|Z_{\bf r}|\ge c\sqrt L\), whose probability is \(o(1)\) by the
uniform exponential tail.  Thus the strict version (2) gives the actual
strict receiver rays asserted above.

For a scope-limited *upper* theorem one must retain threshold plateaux.  Let
\(A_L\) be the supremum only over the actual anchored rectangles
\(\prod_i[-r_ix_i,x_i]\), \(0\le r_i\le1\), with receivers in the positive
cube.  On compact radial intervals the preceding estimate, a finite
\({\bf r}\)-net, and the continuous Gamma-mixture law of \(T_{\bf r}\) give
\[
 \limsup_L A_Lf_L(s,\theta)
 \le\sup_{\bf r}\mathbb Eh_j(u+T_{\bf r})
 \quad(\theta\in C_j^\circ,\ s=L+u).\tag{6a}
\]
Accordingly,
\[
 \limsup_{L\to\infty}\frac{|\{A_Lf_L>1\}|}{e^{-L}L^{n-1}}
 \le\sum_j|C_j|\int e^{-u}
 1_{\{\sup_{\bf r}\mathbb Eh_j(u+T_{\bf r})\ge1\}}\,du.\tag{6b}
\]
The non-strict \(\ge1\) is essential: it follows from
\(\limsup1_{\{a_L>1\}}\le1_{\{\limsup a_L\ge1\}}\), so no false
convergence of halo indicators at a plateau is used.  Uniform increment
tails remove very negative \(u\), and \(e^{-u}\) removes positive \(u\).

Finally, multiply each \(h_j\) by \(1+\delta\), use the scalar anchored
weak theorem \(W_n=Q_n\), and let \(\delta\downarrow0\).  Bounded-profile
dominated convergence gives that the right side of (6b) is at most
\(Q_n\sum_j|C_j|J_j\).  With (3), fixed-scale finite angular-cell inputs
cannot exceed \(Q_n\) for this genuine anchored family.  This says nothing
about all nonanchored rectangles of \(M_S\).

## Result 2: the critical \(L^{-1}\) angular scale leaves a correlated kernel

Fix \(\alpha\) in the interior of the simplex and put
\(P_\alpha z=z-\alpha\sum_i z_i\).  For a bounded compactly supported local
profile \(H(u,z)\), \(\sum_i z_i=0\), set
\[
 f_L(e^{-X})=H\!\left(s-L,\ L\left(\frac Xs-\alpha\right)\right).\tag{6}
\]
At a receiver
\(b=(L-q)\alpha+\xi\), \(\sum_i\xi_i=0\), the very same actual anchored
rectangle has the local limit
\[
 \left(s(b+Z)-L,\ L\left(\frac{b+Z}{|b+Z|_1}-\alpha\right)\right)
 \longrightarrow
 \left(T_{\bf r}-q,\ \xi+P_\alpha Z_{\bf r}\right).\tag{7}
\]
The two components on the right share \(Z_{\bf r}\), so this is genuinely a
correlated \(n\)-variable convolution model, not the scalar \(Q_n\) kernel.
Wolfram 15 independently simplified the scalar component identity
\[
 L\left(\frac{(L-q)\alpha+\xi+Z}{L-q+T}-\alpha\right)
 \to\xi+Z-\alpha T=\xi+P_\alpha Z,
\]
componentwise.

In the coordinates \(s=L+u\), \(\theta=\alpha+z/L\), (1) gives the local
cost limit
\[
 e^L\int\Phi_n(f_L)dx\longrightarrow
 \int e^{-u}\Phi_n(H(u,z))\,du\,dz.\tag{8}
\]
The factor \(L^{n-1}\) in the coarea Jacobian is exactly cancelled by the
angular-window volume \(L^{-(n-1)}\).  Local receiver windows have the same
scaling.  Hence a fixed-scale corner/shell perturbation in this model does
not wash out, but it also receives no free angular-volume gain.

There is no dimensional reduction hidden here.  Put \(w=X-L\alpha\).  Then
\(u=\sum_iw_i\), \(z=P_\alpha w\), and if
\(H(u,z)=F(w)\), the local window is just a diagonal log translation of the
general profile \(F\):
\[
 dx=e^{-L}e^{-\sum_iw_i}\,dw
\]
(with the usual \(2^n\) factor if all sign orthants are used).  The kernel
in (7) is exactly the general positive log-coordinate anchored model
\(F(w)\mapsto\mathbb EF(w+Z_{\bf r})\), written in sum/angular variables.
Thus it is a legitimate place to search for a correlated input, but it is
not a new cheap mechanism or an improvement over \(Q_n\) by itself.
