# r080 independent audit: full-halo column budgets fail, in every fixed dimension

## 1. Audit of the two-dimensional smooth example

The construction in no_budget80.md is valid.  For
\[
 F_2(s,t)=2st\mathbf1_{[0,1]^2}(s,t),\qquad
 u(x)=x-\sqrt{x^2-1},
\]
the average over the lower-left rectangle from \((a,b)\in[0,1]^2\) to
\((x,y)\in(1,\infty)^2\) is
\[
 { (1-a^2)(1-b^2)\over2(x-a)(y-b)}.
\]
The one-dimensional factor is uniquely maximized at \(a=u(x)\).  Hence the
outer-corner maximal field is exactly \(2u(x)u(y)\), and the boundary
\(2u(x)u(y)=1\) has the asserted unique optimizing anchor.

On \(I=[25/24,13/12]\), one has \(u(I)=[2/3,3/4]\).  The curve used in the
outside proof,
\[
 (x,y)=(x,Y(x)-\delta),\qquad
 u(Y(x))={1\over2u(x)},\quad0<\delta<1/96,               \tag{1}
\]
lies in the strict halo and has unit area Jacobian in \((x,\delta)\).
The following two estimates check its budget argument.

* A successful column must have both first two anchor coordinates in
  \((0,1)\).  An anchor coordinate at least one has no positive-measure
  source intersection; an anchor coordinate at most zero has a
  one-dimensional factor at most \(1/(2x)\), which even combined with the
  optimal other factor is below one on this patch.
* Its successful rectangle has area at least
  \[
  (x-1)(y-1)\ge {1\over24}{1\over32}={1\over768}.
  \]
  Thus \(J_2(G)\ge\int G>1/768\).

For fixed anchor, the identity
\[
 {1-a^2\over2(x-a)}
 =u(x)-{(a-u(x))^2\over2(x-a)}
\]
and the depth estimates in the outside proof give
\[
 |a-u(x)|^2<10\delta.
\]
Since \(|u'(x)|\ge8/5\) on \(I\), the successful section of the depth
curve has \(x\)-length less than
\[
 {5\sqrt{10}\over4}\sqrt\delta<4\sqrt\delta.             \tag{2}
\]
These steps validate the stated exact-cover contradiction.  If
\(B=\sum_jw_jJ_2(G_j)\), Tonelli yields
\[
 \int_I W_\delta(x)\,dx\le3072B\sqrt\delta,              \tag{3}
\]
where \(W_\delta\) is weighted success coverage on (1).  Exact a.e.
coverage is impossible as \(\delta\downarrow0\).  Moreover its clipped
\(L^1\) deficit satisfies the quantitative bound already recorded there:
\[
 D\ge{|I|\over3}
 \min\left\{{1\over96},
       \left({|I|\over3072B}\right)^2\right\}.            \tag{4}
\]
Thus no sequence of uniformly bounded-budget column families can have
deficit tending to zero.

## 2. A proved all-\(N\) extension

For \(N\geq3\), the transverse-volume issue does not invalidate the
obstruction.  It can be removed with a fixed-anchor transverse Hardy
estimate.  Let \(D=N-2\) and
\[
 F_N(s,t,z)=2st\,\mathbf1_{[0,1]^N}(s,t,z),\qquad
 E_N=E_2\times(0,1)^D,                                  \tag{5}
\]
where \(E_2\) is the strict patch in (1).

Fix a column \(0\leq G\leq F_N\) and a fixed oriented anchor.  Since
\[
 H_aF_N=H_{a_1,a_2}F_2\,
             H_{a'}\mathbf1_{[0,1]^D}
 \leq H_{a_1,a_2}F_2,                                    \tag{6}
\]
success at a point of the depth-\(\delta\) slice forces its first
coordinate \(x\) into an interval of length at most \(4\sqrt\delta\), by
(2), uniformly in the transverse receiver \(z\).

Define the transverse marginal
\[
 M_G(z)=\int_{[0,1]^2}G(s,t,z)\,ds\,dt.
\]
Then
\[
 0\leq M_G\leq\tfrac12\mathbf1_{[0,1]^D},\qquad
 \int M_G=\int G\leq J_N(G).                             \tag{7}
\]
At any successful receiver, the first-two-coordinate rectangle has area at
least \(1/768\).  Integrating \(G\) first in these two coordinates therefore
shows that the transverse fixed-corner average obeys
\[
 H_{a'}M_G(z)>{1\over768}.                               \tag{8}
\]

For completeness, retain the real transverse corner \(a'\) and split
\(\mathbb R^D\) into its \(2^D\) orthants relative to \(a'\).  On each
orthant the fixed-corner average is a tensor one-sided Hardy operator with
\(L^2\) norm at most \(2^D\), applied to the matching orthant piece of
\(M_G\).  The output orthants and the corresponding input pieces are
disjoint, so the squared norms add without an extra orientation factor.
Chebyshev, (7), and \(M_G^2\leq\tfrac12M_G\) give
\[
 \left|\left\{z:H_{a'}M_G(z)>1/768\right\}\right|
 \leq 4^D768^2\int M_G^2
 \leq294912\,4^D J_N(G).                                 \tag{9}
\]
Combining the \(x\)-section bound with (9), one column covers at most
\[
 A_N\sqrt\delta\,J_N(G),\qquad
 A_N=4\cdot294912\cdot4^{N-2},                           \tag{10}
\]
of the \((x,z)\)-measure of a depth slice.

For a weighted family of total cost \(B\), its coverage \(\Theta_\delta\)
on \(I\times(0,1)^D\) satisfies
\[
 \int\Theta_\delta\leq A_NB\sqrt\delta.                  \tag{11}
\]
The slice itself has measure \(|I|\).  Therefore exact a.e. coverage of
\(E_N\) requires \(B=\infty\).  The same integration used for (4) gives
the explicit fixed-\(N\) deficit estimate
\[
 D_N\ge{|I|\over3}
 \min\left\{{1\over96},
       \left({|I|\over A_NB}\right)^2\right\}.            \tag{12}
\]
Together with Section 1 for \(N=2\), this closes the all-fixed-dimension
extension.  It does not claim a dimension-uniform constant because \(A_N\)
grows exponentially.

## 3. Exact compact-main packet transfer and near-optimal inputs

The packet argument can be stated without a tail assumption.  Fix any
compact cap-two main source \(H\), with \(0<J_N(H)<\infty\).  Choose a
translation \(H_R\) in its first two coordinates so that its support
projections lie in
\[
 [R,R+L_1]\quad\hbox{and}\quad[R,R+L_2].                 \tag{13}
\]
Let \(F_{{\rm bad},\theta}\) be the spatial \(\theta\)-dilation of (5),
supported in \([0,\theta]^N\).  Its bad receiver patch is
\(\theta E_N\), whose first two receiver coordinates are at most
\(13\theta/12\).  For \(R>13\theta/12\), the two sources have disjoint
supports, so their sum
\[
 F_{\theta,R}=F_{{\rm bad},\theta}+H_R
\]
is again cap-two.

Any oriented rectangle containing a receiver in \(\theta E_N\) and meeting
the translated main support has first-two side lengths at least
\[
 R-{13\theta\over12}.
\]
The union of the first-two source projections has respective lengths at
most \(\theta+L_1\) and \(\theta+L_2\).  Since the source height is at most
two, integration first in the transverse coordinates gives the uniform
spanning-rectangle estimate
\[
 {1\over|Q|}\int_QF_{\theta,R}
 \leq {2(\theta+L_1)(\theta+L_2)\over
                (R-13\theta/12)^2}.                    \tag{14}
\]
Choose \(R\) so that the right side is strictly below one.

This proves the key isolation statement.  If a source column
\(0\leq G\leq F_{\theta,R}\) succeeds at a receiver of \(\theta E_N\), its
successful rectangle cannot meet the main support.  Hence replacing it by
\[
 G_{\rm bad}=G\mathbf1_{\operatorname{supp}F_{{\rm bad},\theta}}
\]
leaves its successful average at every such receiver unchanged and can only
decrease \(J_N(G)\).  A finite-budget column cover of the full halo of
\(F_{\theta,R}\) would therefore restrict to a finite-budget cover of
\(\theta E_N\) under the bad packet.  Scaling back by \(\theta\) contradicts
the infinite-budget theorem of Section 2.  Thus:

\[
 \boxed{\text{Every compact cap-two main source can be augmented by an
 arbitrarily small packet to produce a cap-two source with infinite
 full-halo fixed-corner column budget.}}                 \tag{15}
\]

The added cost is exactly
\[
 J_N(F_{{\rm bad},\theta})=\theta^N J_N(F_N).
\]
For any \(\eta>0\), take \(\theta\) so this is at most \(\eta J_N(H)\).
Translation does not change either the main cost or its halo measure, and
maximality is monotone, so
\[
 {|\{MF_{\theta,R}>1\}|\over J_N(F_{\theta,R})}
 \geq {|\{MH>1\}|\over(1+\eta)J_N(H)}.                  \tag{16}
\]

The fixed-\(N\) rational compact-step density result of r079 says that
compact cap-two step sources have the same supremal halo/modular ratio as
the entire cap-two class.  Consequently, choose a compact \(H\) with
arbitrarily near-optimal score and then let \(\eta\downarrow0\) in (16).
There are therefore cap-two inputs of arbitrarily near-optimal score whose
full-halo fixed-corner column budget is infinite.  This conclusion does not
need an attachment to a prescribed noncompact near-extremizer.

It rules out the old r078 Section 6 style of a uniform global full-halo
source-column budget target.  It does not rule out a receiver-local residual
or another mechanism that is not required to furnish a single global
fixed-corner column cover of the entire halo.  Nor does it contradict a
general strong-maximal bound or itself prove a new \(C_N\) estimate.

For the following corollary set A_2=3072, and use (10) for N>=3.

## 4. Rational compact-step near-optimal corollary

The same obstruction is present inside rational compact step inputs, without
passing to a smooth source of literally infinite budget.  Let
\[
 F_m(s,t)=2{\lfloor ms\rfloor\over m}{\lfloor mt\rfloor\over m}
              \mathbf1_{[0,1]^2}(s,t),\qquad
 F_{N,m}(s,t,z)=F_m(s,t)\mathbf1_{[0,1]^D}(z).           \tag{17}
\]
For \(m>384\), the r080 step-source calculation gives
\[
 0\leq F_{N,m}\leq F_N,\qquad J_N(F_{N,m})\leq J_N(F_N),
 \qquad E_{4/m}\times(0,1)^D\subset\{MF_{N,m}>1\}.       \tag{18}
\]
Here \(E_{4/m}\) is the truncated strict outer-corner patch
\(\delta\geq4/m\).  The transverse factor in the witnessing rectangle is
one, so the two-dimensional inclusion persists for every transverse
receiver.

Apply (11) to a full column cover of the halo in (18).  It holds for almost
every depth \(\delta\geq4/m\); taking such depths down to \(4/m\) gives
\[
 B_m\geq {|I|\over A_N\sqrt{4/m}}
       ={ |I|\sqrt m\over2A_N}.                          \tag{19}
\]
Thus each individual rational step source has a finite possible column
budget, but these budgets cannot be uniformly bounded in \(m\).

Finally fix \(N\).  The r079 density theorem supplies a compact rational
step \(H\) whose cap-two halo/modular score is arbitrarily close to the
fixed-\(N\) cap-two supremum.  Choose a fixed rational \(\theta>0\) with
\[
 \theta^N J_N(F_N)\leq\eta J_N(H),
\]
and then choose a fixed rational separation \(R\) satisfying (14).  For
each \(m\), splice the translated \(H_R\) with the
\(\theta\)-scaled \(F_{N,m}\).  The result is again a compact rational step
source of height at most two.  Its score is at least
\[
 {|\{MH>1\}|\over(1+\eta)J_N(H)},                        \tag{20}
\]
while the bad-patch restriction and scaling of (19) give its full-halo
column-budget lower bound
\[
 B_{m,\theta,R}\geq
 {\theta^N|I|\sqrt m\over2A_N}\longrightarrow\infty.    \tag{21}
\]
Hence even the rational compact-step class contains inputs of arbitrarily
near-optimal score for which no uniform global full-halo fixed-corner column
budget can hold.  This is still a statement about the budget architecture,
not a claim that the strong-maximal weak constant diverges.
