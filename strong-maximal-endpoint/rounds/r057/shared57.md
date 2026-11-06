# r057: pairwise and higher-order common-action regret

This note retains cross-height coupling by working with the full fiber loss,
rather than separate superlevel coverage. It supplies a valid pairwise lower
certificate, then shows by an actual step input that pairwise incompatibility
is not sufficient in general.

## 1. Pair and k-wise regret certificates

For one fixed outer interval `J`, set
\[
 L_y(I)=m_y-F_y(I)\ge0,\qquad
 R_J=A_J-S_J=\inf_I E_JL_y(I).
\]
For independent uniform `Y,Z` in `J`, define the pair incompatibility
\[
 d_2(Y,Z)=\inf_I\{L_Y(I)+L_Z(I)\}.
\]
For every common interval `I`,
\[
 E_{Y,Z}[L_Y(I)+L_Z(I)]=2E_YL_Y(I).
\]
Taking the infimum after the pointwise inequality `d_2<=L_Y+L_Z` gives
\[
 \boxed{\quad R_J\ge {1\over2}E_{Y,Z}d_2(Y,Z).\quad}\tag{1}
\]
This is tie-invariant and couples the entire payoff, so it is not vulnerable
to the cross-height defect of the r056 layer certificate.

More generally, for `k>=2`,
\[
 d_k(y_1,\ldots,y_k)=\inf_I\sum_{j=1}^kL_{y_j}(I)
\]
satisfies
\[
 \boxed{\quad R_J\ge {1\over k}E d_k(Y_1,\ldots,Y_k).\quad}\tag{2}
\]
Thus a k-fiber incompatibility certificate is a rigorous sufficient
condition for common-action loss.

For a finite action set of size `M`, with all losses in `[0,B]`, this
hierarchy asymptotically recovers the exact regret. If
\[
 \Delta_k={1\over k}E d_k(Y_1,\ldots,Y_k),
\]
then Hoeffding's exponential-moment bound and a union bound over the `M`
actions give
\[
 R_J-B\sqrt{{\log M\over2k}}\le\Delta_k\le R_J.\tag{3}
\]
Indeed `R_J-min_I k^{-1}sum_jL_{Y_j}(I)` is bounded by the maximum of the
`M` centered lower deviations, whose expected maximum has the displayed
bound. Equation (3) is a fractional-packing/sample dual for full regret:
unlike a layer test it preserves the fact that one interval must serve all
heights of all sampled fibers at once.

For the global two-level problem, it is enough to require, for every outer
interval with `A_J>A_*-eta`, a common `k` and lower bound
`Delta_k(J)>=delta`. The same near-outer split as in r054 then gives a global
pointwise loss at least `min(eta,delta)`. This is a real all-near-outer
condition, though no theorem here derives it from the Orlicz budget.

## 2. An actual three-fiber obstruction to pairwise certification

Pairwise loss can vanish identically while the full common-action regret is
positive. At inner receiving coordinate `x=0`, define
\[
 g_1=\mathbf1_{[0,2]},\quad
 g_2=2\mathbf1_{[1,2]}+\mathbf1_{[2,3]},\quad
 g_3=\mathbf1_{[0,1]}+2\mathbf1_{[2,3]}.
\]
Take three equal outer slabs and, on their union `J`, set
\[
 f(x,y)=g_j(x)\quad\hbox{on slab }j.\tag{4}
\]
For the inner intervals `[0,r]`, direct integration gives
\[
 q_1(r)=\begin{cases}1&0<r\le2,\\2/r&r>2,\end{cases}
\]
\[
 q_2(r)=\begin{cases}0&r\le1,\\2-2/r&1<r\le2,\\1&2<r\le3,\\3/r&r>3,
 \end{cases}
\quad
 q_3(r)=\begin{cases}1&r\le1,\\1/r&1<r\le2,\\2-3/r&2<r\le3,\\3/r&r>3.
 \end{cases}\tag{5}
\]
Each fiber maximum is one. Any pair has a common maximizing interval:
`r=2` for `(g_1,g_2)`, `r=1` for `(g_1,g_3)`, and `r=3` for `(g_2,g_3)`.
Consequently
\[
 d_2(g_i,g_j)=0\quad\hbox{for all }i,j.\tag{6}
\]
But no single interval maximizes all three, and
\[
 A_J=1,\qquad S_J=\sup_{r>0}{q_1(r)+q_2(r)+q_3(r)\over3}={8\over9},
 \qquad R_J={1\over9}.\tag{7}
\]
The triple loss is `d_3(g_1,g_2,g_3)=1/3`, so the order-three certificate
detects the obstruction whereas (1) does not. Under three independent
uniform slab samples, only the ordered triples containing all three types
contribute; hence `Delta_3=(1/3)(6/27)(1/3)=2/81>0`.

This is a physical finite-amplitude input, not an abstract action game. For
three unit outer slabs its exact two-dimensional Orlicz cost is
\[
 \int\Phi_2(f)=4+2\Phi_2(2)=8+4\log2.\tag{8}
\]
It is only a fixed-outer-interval statement: side slabs can change the
global outer optimization, so it makes no strong-maximal or `C_N` claim.

Geometrically, pairwise fiber-maximizer intersections need not imply a
global intersection. Here the action has been reduced to the ray parameter
`r`, but the exact maximizer set for `g_3` is disconnected (`(0,1] union {3}`), so
one-dimensional Helly intuition does not apply. A regularity theorem that
made all epsilon-maximizer sets genuine intervals (or convex endpoint sets)
could reduce the necessary order through Helly theory; arbitrary fibers do
not have that property.

## 3. No fixed arity is sufficient without regularity

The preceding pair failure extends to every prescribed arity. Fix `M>=2`
and `0<delta<1`, take `M` equal outer unit slabs, and on the inner interval
`[0,M]` define fiber `i` by
\[
 g_i=\mathbf1_{[0,M]}-\delta\mathbf1_{[i-1,i]}
       +\delta\mathbf1_{\{i<M\}}\mathbf1_{[i,i+1]}.
\tag{9}
\]
Set `f(x,y)=g_i(x)` on outer slab `i`, take the fixed outer interval
`J=[0,M]`, and receive at inner coordinate zero. Every prefix average is at
most one. At the integer prefix `r_j=j`, the exact table is
\[
 {1\over j}\int_0^j g_i=
 \begin{cases}1-\delta/j&i=j,\\1&i\ne j.\end{cases}\tag{10}
\]
Hence every fiber has maximum one, while every proper subset of the fibers
is jointly maximized by the prefix corresponding to a missing index.

The aggregate input has height `1-delta/M` on its first inner unit cell and
height one on every later cell. Therefore its best common prefix average is
\[
 S_J=1-{\delta\over M^2},\qquad A_J=1,\qquad R_J={\delta\over M^2}>0.
\tag{11}
\]
For any fixed `k<M`, every k-tuple of fibers has a proper distinct-index set,
so its joint loss is zero. Consequently
\[
 \Delta_k=0\quad(k<M),\tag{12}
\]
even though the full common-action regret is positive. This is an actual
bounded-amplitude input; taking `delta<=1/2` confines values on the support rectangle to
`[1/2,3/2]`. Its exact two-dimensional Orlicz cost is
\[
 \int\Phi_2(f)=M^2-\delta+(M-1)(1+\delta)\log(1+\delta).\tag{13}
\]

Thus a fixed-order pair/triple/... certificate cannot yield a universal
sharing theorem for arbitrary fibers. A useful endpoint theorem must either
let the incompatibility order grow with the action complexity, or impose a
convex/quasi-concave structure that bounds the Helly rank. As above, this is
only a fixed-outer-interval construction and makes no `C_N` assertion.

## 4. Scope for endpoint work

The certificates (1)--(3) retain the actual common interval and are immune
to arbitrary tie-breaking and the r056 height-by-height failure. They offer
a concrete finite-grid route: enumerate or bound k-wise losses while keeping
the exact denominator `sum_C |C|Phi_n(f|_C)`. What remains open is a
dangerous-input theorem giving a positive, feedback-scale lower bound for
some fixed or slowly growing `k` uniformly over all near-optimal outer
intervals. The flat/product and three-slab inputs show why that conclusion
cannot follow from nonnegativity, endpoint score, or pairwise variance alone.

`shared57.wl` was actually run in Wolfram 15.0.0. It verifies the exact
maxima, zero pair-loss matrix, positive triple loss, and cost in (7)--(8),
as well as the `M=5`, `delta=1/2` instance of (10)--(13). The inequalities
(1)--(3) are analytic consequences of the common-action regret definition,
not numerical inferences.

For completeness, the finite-class exponential bound needs no sampling experiment: for a centered loss in an interval of length B, the log-mgf has second derivative at most B²/4 under every exponential tilt. Integrating twice from zero gives log E exp(t(EL−L))<=t²B²/8. Independence, log-sum-exp over M actions, and optimization in t give the stated B sqrt(logM/(2k)) expected-deviation bound.
