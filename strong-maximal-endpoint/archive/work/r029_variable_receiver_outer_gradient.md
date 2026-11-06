# r029: exact outer ratio gradient and a monotone finite-receiver update

Date: 2026-11-06.  This note concerns a finite collection of genuine
block-anchored receiver constraints.  It is an algorithmic/certification
interface, not a global upper bound and not evidence that the same-\(c\)
scalar benchmark bounds variable-ratio inputs.

Let \(f_{b,c}\) be the density of a sum of \(b\) copies of
\(Z_c=\operatorname{Exp}(1)+c\operatorname{Bern}((1+e^c)^{-1})\), and put
\(p_c=(1+e^c)^{-1}\).  With \(f_{0,c}=\delta_0\), the exact distributional
identity is
\[
 \partial_c f_{b,c}(t)=bp_c\{f_{b,c}(t)-f_{b-1,c}(t-c)\}.\tag{1}
\]
For \(b=1\), the second term is the required moving delta mass.  A quick
proof is to differentiate one convolution factor:
\(\partial_cf_{1,c}=p_c(f_{1,c}-\delta_c)\), then sum over the \(b\)
identical factors.  Formula (1) is valid weakly even across the truncated
power knots; for a discontinuous KKT profile it supplies an a.e. derivative,
while trial-step evaluation remains the definitive feasibility test.

For one receiver \(u=(u_1,u_2)\), block parameters \((c_1,c_2)\), and a
fixed bounded profile \(H\), define its actual anchored average
\[
 B_u(H;c_1,c_2)=\mathbb E H(u_1+T_{b,c_1},u_2+T_{b,c_2}).
\]
The two exact outer gradients are
\[
 \partial_{c_1}B_u=bp_{c_1}\left[B_u-
 \mathbb EH(u_1+c_1+T_{b-1,c_1},u_2+T_{b,c_2})\right],\tag{2}
\]
\[
 \partial_{c_2}B_u=bp_{c_2}\left[B_u-
 \mathbb EH(u_1+T_{b,c_1},u_2+c_2+T_{b-1,c_2})\right].\tag{3}
\]
For the physical ratio \(r=e^{-c}\), \(\partial_rB=-r^{-1}\partial_cB\).
Thus a variable-receiver search should optimize the two components in
\((c_{j,1},c_{j,2})\) separately (or a specified symmetry group jointly),
not prescribe the linear ansatz \(c\pm\theta v/n\).

## Feasibility-preserving alternating update

Fix finitely many receiver boxes \(U_j\), their anchored parameters
\(c_j\), a cap \(A\), and an input/profile \(H\) satisfying the strict
constraints
\[
 \inf_{u\in U_j}B_u(H;c_j)\ge1+\eta\quad(j=1,\ldots,M).\tag{4}
\]
With \(H\) held fixed, update each \(c_j\) by any bounded line-search step
whose *actual recomputed* minimum in (4) does not decrease.  Gradients
(2)--(3) are the ascent directions; evaluating (4), rather than a linear
prediction, handles knots and finite quadrature errors.

The old \(H\) is then feasible for the updated finite convex budget problem
\[
 \min_{0\le G\le A}\int e^{-\sum t}\Phi_n(G(t))dt
 \quad\text{subject to}\quad
 \inf_{u\in U_j}B_u(G;c_j^{\rm new})\ge1+\eta.
\]
After re-solving it, the new budget cannot exceed the old one.  If the
receiver boxes are held fixed, their union volume is unchanged, hence the
certified finite-box score is nondecreasing at every alternating step.  This
is the appropriate monotonicity lemma; moving receiver locations or changing
the union geometry needs a separate volume comparison and has no automatic
monotonicity.

## Actual-certificate interface

A list of receiver *points* has zero physical measure.  To turn a numerical
solution into a lower bound, retain boxes \(U_j\) with a positive margin in
(4), certify each block-kernel integral on the whole box (e.g. interval
quadrature plus a derivative/Lipschitz bound), and compute a rigorous lower
bound for \(|\bigcup_jU_j|\) in the reference \(e^{-u_1-u_2}du\) measure.
The block deep-layer reduction then transfers that strict finite-box witness
to actual reflected rectangles for each receiver point, with ratios
\(r_{j,\ell}=e^{-c_{j,\ell}}\).  The existing MATLAB grid/SQP output is a
useful initializer only: its SQP flags and a quadrature duality gap near
\(10^{-2}\) are not such a certificate.
