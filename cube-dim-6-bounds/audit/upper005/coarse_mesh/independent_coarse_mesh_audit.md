# Independent audit: fixed coarse-mesh obstruction

## Result

**Pass, with a strict scope restriction.** The numerical inequalities hold,
and their directions are correct only as statements about the specified
*upper-certificate budget*. They do not lower-bound the actual norm, actual
tail, \(I\), the model cost, or \(C_6\).

## Direction and interface audit

Let \(q=6/5\), and let \(B\) be the norm budget in the conditional cost
\[
6^q(B/g)^q.
\]
The target is below 64 exactly when \(B<16g/3\). The non-\(S_2\) exact point
has a continuous neighborhood, so its value gives an *essential* upper bound
\(g\le U=23\cdot10^{-6}\); hence \(B<16U/3\) is necessary for this route.
This uses the point in the right direction. A single point without the
non-\(S_2\)/continuity fact would not suffice for essential \(g\).

For \(\epsilon\le1/100\), \(H_\epsilon=\int_\epsilon^\infty p_tdt\)
contains the spectral range \(t\ge T\), so the prior uniform bound gives
\(H_\epsilon\ge L=6\cdot10^{-8}\). In the 120-leaf formula, every vertex
upper bound must be at least its actual value, hence at least \(L\). Since
every leaf has \(D_S\ge D_{\min}=25/192\), their volumes sum to \(1/120\),
and the 120 permutation copies reconstruct the face, the RHS of the summed
upper certificate, after the \(1/q\) power, is at least
\[
L\exp\{D_{\min}/(32\epsilon^2)\}.
\]
It is decreasing in \(\epsilon\), so its smallest value on this branch is at
\(1/100\). This is why `Dmin`, rather than `Dmax`, is sufficient: it
lower-bounds every term in the certificate RHS.

For \(\epsilon\ge1/100\), the prescribed whole-face tail allocation
\(K_{new}\epsilon^{7/6}\) is increasing, so it too is smallest at the split.
It is a floor on the *allocated upper budget*, not on the physical tail.

## Independent numerical checks

An independent 384-bit Arb script, with no kernel evaluation, gives
\[
16U/3=1.2266666666666667\cdot10^{-4},
\]
\[
L\exp(D_{\min}/(32(0.01)^2))
=2.816040899053389\cdot10^{10},
\]
\[
K_{new}(0.01)^{7/6}=1.952888275947603\cdot10^{-4}.
\]
Thus the tail branch exceeds the necessary ceiling by a factor about 1.59203,
and the coarse-mesh branch by about \(2.30\cdot10^{14}\). A separate MATLAB
double-precision scalar calculation reproduced these values; Arb is the
interval evidence.

The reproducible script is `independent_coarse_mesh_audit.py`, with its
machine-readable output in `independent_coarse_mesh_audit.json`.
