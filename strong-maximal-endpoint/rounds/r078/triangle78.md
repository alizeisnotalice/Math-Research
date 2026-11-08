# A three-anchor, same-orientation fractional-cover test

## Scope

This note records a finite, explicitly rational certified gap between
the three-anchor integral allocation model and its fractional column-cover
relaxation.  It is a test of **this prescribed anchor list only**.  It is
not a counterexample to a general strong-maximal estimate, and it does not
license restricting the price search to these three anchors.

Let
\[
 Q=[-9/10,-1/10]^2,\qquad |Q|=16/25,\qquad
 F={25\over16}{\bf1}_Q,\qquad \sigma=1+10^{-8}.          \tag{1}
\]
Thus \(\int F=1\).  Set \(\varepsilon=10^{-7}\) and let
\[
 \begin{split}
 p_1&=(0,2/1000),&p_2&=(1/1000,21/20000),&p_3&=(2/1000,0),\\
 P_i&=p_i+(-\varepsilon,\varepsilon)^2.                 \tag{2}
 \end{split}
\]
The candidate anchors are \(a_j=(-h_{j1},-h_{j2})\), where
\[
 \begin{split}
 h_1&=\left({21\over20},{1-2077/10^6\over21/20}\right),\\
 h_2&=\left(1,1-2025/10^6\right),\\
 h_3&=\left({1-2075/10^6\over21/20},{21\over20}\right). \tag{3}
 \end{split}
\]
All endpoints are rational.  A common axis-parallel refinement containing
the faces of \(Q\), (2), and (3) therefore meets the exact finite reduction
interface of `r077/allocation77.md`.

## 1. Exact success and column cost formula

For \(x\in P_i\), the lower-left oriented rectangle \([a_j,x]\) contains
\(Q\): each \(h_{jr}>9/10\), while every coordinate of \(x\) is larger
than \(-10^{-7}>-1/10\).  (The first coordinate on \(P_1\) and the second
on \(P_3\) need not be positive.)
Consequently, for every \(0\le G\le F\),
\[
 \int_{[a_j,x]}G=\int_QG=:m(G),
 \qquad |[a_j,x]|=(h_{j1}+x_1)(h_{j2}+x_2).               \tag{4}
\]
The volume is coordinatewise increasing.  Put
\[
 V_{ji}=\max_{x\in\overline {P_i}}|[a_j,x]|
  =(h_{j1}+p_{i1}+\varepsilon)(h_{j2}+p_{i2}+\varepsilon).
\tag{5}
\]
Then \(F\) succeeds on all of \(P_i\) precisely when
\[
 \sigma V_{ji}\le1.                                    \tag{6}
\]
The word “precisely” here refers to this particular \(F\); for a general
subsource \(G\le F\), success on \(P_i\) is the scalar condition
\(m(G)\ge\sigma V_{ji}\).

For a nonempty eligible set \(S\) at anchor \(a_j\), write
\[
 M_{j,S}=\sigma\max_{i\in S}V_{ji}.                     \tag{7}
\]
Provided \(M_{j,S}\le1\), the exact minimum source charge is
\[
 c_{a_j}(S)=M_{j,S}\left[1+\log{M_{j,S}\over |Q|}\right].\tag{8}
\]
Indeed, (4) forces \(m(G)\ge M_{j,S}\).  Jensen on \(Q\) gives
\[
 \int_Q\Phi_2(G)\ge |Q|\Phi_2\left({m(G)\over|Q|}\right)
 \ge M_{j,S}\left[1+\log{M_{j,S}\over |Q|}\right],     \tag{9}
\]
because the relevant densities lie in \([1,2]\).  Equality is attained by
the constant subsource \(G=(M_{j,S}/|Q|){\bf1}_Q\), which is bounded by
\(F\) exactly when \(M_{j,S}\le1\).  Formula (8) is therefore an analytic
all-column price, not a numerical ansatz.

## 2. Certified three-anchor incidence and gap

Independent rational-volume and directed-log certificates establish
\[
 \sigma V_{jj}>1,\qquad \sigma V_{ji}<1\quad(i\ne j),   \tag{10}
\]
where the endpoint perturbation is already included in (5).  Thus anchor
\(a_j\) is eligible for exactly the two cells \(P_i\), \(i\ne j\), and the
complete list of nonempty finite-cost columns for the prescribed anchors has
nine elements: the two singleton subsets and their pair at each anchor.

One anchor cannot cover all three cells, whereas any two anchors can.
The integral model consequently needs at least two selected anchor sources.
The fractional cover taking weight \(1/2\) on each of the three pair
columns covers every \(P_i\) with total weight one.  Exhaustive certification
of all nine costs in (8), the integral dynamic program, and the full
nine-column LP and dual yields
\[
 1.333160990
 < {V_{\rm IP}^{\{a_1,a_2,a_3\}}\over
       V_{\rm CLP}^{\{a_1,a_2,a_3\}}}
 <1.333160992.                                            \tag{11}
\]
Thus the ratio is close to, but not asserted equal to, \(4/3\).

The completed finite certificate uses:

1. evaluate all nine rational volumes (5) outwardly and prove (10);
2. use directed intervals for the logarithms in (8);
3. enumeration of the finite integral choices (or the equivalent MICP); and
4. a nonnegative three-cell dual vector feasible for **all nine** columns.

More explicitly, if \(p_{12},p_{13},p_{23}\) denote the costs of the pair
columns covering the indicated two cells, then their half-weight cover has
cost \((p_{12}+p_{13}+p_{23})/2\).  The certified dual is
\[
 \lambda_1={p_{12}+p_{13}-p_{23}\over2},\quad
 \lambda_2={p_{12}+p_{23}-p_{13}\over2},\quad
 \lambda_3={p_{13}+p_{23}-p_{12}\over2}.                 \tag{12}
\]
Its interval verification gives \(\lambda_i\ge0\), equality in every pair
constraint, and \(\lambda_i\) no larger than every eligible singleton cost
covering cell \(i\).  Hence it is a full nine-column LP certificate, not
only a feasible fractional construction.

## 3. Why full anchor pricing remains indispensable

The three-anchor restriction is material.  The additional lower-left anchor
\[
 a_0=(-9/10,-9/10)                                       \tag{13}
\]
has rectangles \([a_0,x]\) containing \(Q\), and their volumes on the
three receiver cells are
\[
 (9/10+x_1)(9/10+x_2).                                   \tag{14}
\]
which are well below one throughout the displayed tiny cells.  Hence the
same scalar calculation produces an eligible single column covering all
three cells, with cost given by (8) after replacing \(M_{j,S}\) by
\(\sigma\max_i(9/10+p_{i1}+\varepsilon)(9/10+p_{i2}+\varepsilon)\).
The maximum in (14) occurs on \(P_2\).  Let
\(c_0=c_{a_0}(\{P_2\})=c_{a_0}(\{P_1,P_2,P_3\})\).  The rational volume
ordering shows that every other anchor column which covers \(P_2\) has
cost at least \(c_0\): its required mass is larger than the mass for the
natural-corner \(P_2\) column, and the function in (8) is increasing in
that mass.  Therefore the dual with price \(c_0\) on \(P_2\) and price zero
on \(P_1,P_3\) is feasible for all sixteen columns arising from
\(a_0,a_1,a_2,a_3\).  The single all-three-cell \(a_0\) column is a primal
solution of the same cost.  Consequently
\[
 V_{\rm CLP}^{\{a_0,a_1,a_2,a_3\}}
 =V_{\rm IP}^{\{a_0,a_1,a_2,a_3\}}=c_0.                 \tag{15}
\]
This analytic argument is stronger than the independently obtained interval
ratio \([1,1.000000000989421]\) for the same sixteen-column calculation.
It shows exactly why the apparent three-anchor gap is removed by this
natural corner.

Thus a dual checked only against \(a_1,a_2,a_3\) was not a full pricing
certificate even for this visibly enlarged finite list.  Any global
conclusion must price every eligible anchor or justify an independent
restriction of the anchor family.

Thus this construction can certify a finite same-orientation cross-anchor
integrality gap only under the stated eligible-anchor list.  It neither
contradicts the exact same-anchor consolidation lemma nor establishes a
general lower or upper bound for \(C_N\).
