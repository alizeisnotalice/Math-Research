# r057: radial compensation for local common-action adaptation

Scope: freeze the symmetric target-one stationary KKT source \(H_0\) from
r056 in the two-block micro model. This is a local receiver-curve result,
not a result on the general strong maximal constant.

Let
\[
G(d,v,a)=E_{a,a}H_0(T_1+d+v,T_2+d-v). \tag{1}
\]
At \((d,v,a)=(0,0,c)\), write all derivatives without arguments.  Then
\[
G=1,\quad G_a=0,\quad G_d>0,\quad G_{aa}<0,\quad G(d,-v,a)=G(d,v,a). \tag{2}
\]
For \(b\ge6\), the zero-extended Gamma--binomial densities and bounded
\(H_0\) give the needed local \(C^4\) regularity.

## Base identities

Let \(k_\pm=k_1\pm k_2\), \(s_+=s_1+s_2\), and let \(B\) be the full
r050 lower-coarea plus soft-response form.  Integration by parts together
with the r054 Ward identity gives the exact identities
\[
G_d={2P\over\lambda}>0,\qquad
G_{vv}=-B(k_-,k_-)\le0,\qquad
G_{ad}=B(k_+,s_+). \tag{3}
\]
The last identity uses the frozen-profile derivative
\(\partial_dH_0=a h'(a)k_+\); all three occurrences of \(B\) include both
the lower jump and the soft part.  Consequently the required radial
compensation begins with
\[
d_2=-{G_{vv}\over2G_d}\ge0. \tag{4}
\]

## Coupled implicit equations and jets

Choose the common action to minimize the radial boundary by solving
\[
G(d(v),v,a(v))=1,\qquad G_a(d(v),v,a(v))=0. \tag{5}
\]
The \((d,a)\) Jacobian is lower triangular with rows \((G_d,0)\), \((G_{ad},G_{aa})\), and nonzero determinant \(G_dG_{aa}\), so
the implicit-function theorem gives unique local even functions.  Define
the radially corrected mixed coefficient
\[
K=G_{avv}-{G_{ad}G_{vv}\over G_d}
=G_{avv}+{B(k_+,s_+)B(k_-,k_-)\over G_d}. \tag{6}
\]
Then
\[
a(v)=c-{K\over2G_{aa}}v^2+o(v^2),\qquad
d(v)=d_2v^2+d_4v^4+o(v^4). \tag{7}
\]
Thus the raw fixed-radial coefficient \(G_{avv}\) is not decisive: radial
compensation can cancel it through the second term in (6).

Let \(d_{\rm fix}(v)\) solve \(G(d,v,c)=1\).  Its fourth jet is
\[
d_{4,\rm fix}=-{1\over G_d}\left({G_{dd}\over2}d_2^2+
{G_{dvv}\over2}d_2+{G_{vvvv}\over24}\right). \tag{8}
\]
The optimized common action has the exact fourth-order radial improvement
\[
\boxed{d_{\rm opt}(v)-d_{\rm fix}(v)
={K^2\over8G_{aa}G_d}v^4+o(v^4)\le0.} \tag{9}
\]
The displayed quartic coefficient is strictly negative iff \(K\ne0\); higher-order improvements are not excluded when it vanishes.  Wolfram 15.0 actually checked this Taylor
algebra in `geometry57.wl`, returning the equivalent numerator
\((G_{avv}G_d-G_{ad}G_{vv})^2/(8G_{aa}G_d^3)\).

## Actual Pareto-cone volume and strictness

Use anchors
\[
u(v)=(d(v)+v,d(v)-v),\qquad -A\le v\le A. \tag{10}
\]
For small \(A\), this is Pareto ordered.  The exact weighted volume of the
union of its upper cones is
\[
V[d;A]=\int_{-A}^{A}e^{-2d(v)}dv+
{e^{-2d(A)}+e^{-2d(-A)}\over2}. \tag{11}
\]
Since the two boundaries are even, (9) gives
\[
V[d_{\rm opt};A]-V[d_{\rm fix};A]
=-2{K^2\over8G_{aa}G_d}A^4+o(A^4). \tag{12}
\]
So nonzero \(K\) gives a strict fourth-order gain in actual cone-union
volume, with the frozen source's Orlicz denominator unchanged.

The equations above are equality reception.  For strict reception, shift
every anchor outward by a common small \(r>0\).  By \(G_d>0\), all averages
are then strictly above one on a compact arc; both volumes in (12) acquire
the same factor \(e^{-2r}\).  A fine finite subfamily retains strictness and
approximates the union from below, so the fixed-\(N\) deep-transfer argument
makes these actual strong rectangles.

## Remaining coefficient

No nonvanishing or rate for \(K\) is proved here.  If \(K=0\), the local
common-action mechanism has no quartic improvement; if \(K\ne0\), (12) is
a genuine micro-model gain over fixed common action.  Neither alternative
compares with the scalar radial certificate or proves a new general \(C_N\)
bound.


The improved profile-mass estimate and uniform common curvature are proved separately in [common_window57.md](common_window57.md).
