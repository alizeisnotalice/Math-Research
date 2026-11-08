# r078: structure of legitimate same-anchor source-column costs

## Definition and basic facts

Fix an anchor a, a bounded cap-two source \(0\le F\le2\), a common genuine
success level \(\sigma>1\), and a finite family S of receiver cells whose
oriented fixed-corner rectangles are based at a.  Define

\[
 c_a(S)=\inf\left\{J(G):0\le G\le F,\quad
        {1\over|R|}\int_R G\ge\sigma
        \text{ for every assigned rectangle }R\in S\right\}. \tag{1}
\]

Strict success has the same infimum under the usual available strict margin;
the closed formulation is used for convexity and limiting arguments.

Two unconditional properties are:

\[
 S\subset T\Longrightarrow c_a(S)\le c_a(T),\qquad
 c_a(S\cup T)\le c_a(S)+c_a(T).                            \tag{2}
\]

For the second statement take near-minimizers \(G_S,G_T\) and use
\(G=\max(G_S,G_T)\).  It remains below F, succeeds on both families, and
\[
 J(G)\le J(G_S)+J(G_T)                                    \tag{3}
\]
pointwise.  The max construction alone is only ordinary subadditivity; the
stronger fractional statement below uses the particular Orlicz modular.

## Main result: fractional subadditivity at one anchor

Here write the cap-two modular as
\[
 J(G)=\int\Phi_N(G),\qquad
 \Phi_N(t)=t(1+\log^+t)^{N-1}\quad(t\geq0).               \tag{F1}
\]
The function \(\Phi_N\) is continuous, increasing and convex (with the
linear formula on \([0,1]\)); hence it has an increasing inverse on its
range.  This is the only convexity used below.

Let \(\{(S_k,w_k)\}_{k\in K}\) be a finite fractional cover of \(S\):
\(w_k\geq0\), and
\[
 W_i:=\sum_{k:i\in S_k}w_k\geq1\qquad(i\in S).           \tag{F2}
\]
Choose *real* columns \(0\leq G_k\leq F\), each successful on every
constraint in \(S_k\).  Define the capped Orlicz aggregate
\[
 H(y):=\min\left\{\Phi_N(F(y)),\ \sum_k w_k\Phi_N(G_k(y))\right\},
 \qquad G:=\Phi_N^{-1}(H).                                \tag{F3}
\]
It is important that the cap is imposed in modular coordinates.  It gives
measurable \(0\leq G\leq F\) and immediately
\[
 J(G)=\int H\leq\sum_k w_kJ(G_k).                         \tag{F4}
\]

For a fixed receiver constraint \(i\), all columns indexed by
\(k:i\in S_k\) are successful on the *same* rectangle \(R_i\).  Their
ordinary normalized average
\[
 A_i:=W_i^{-1}\sum_{k:i\in S_k}w_kG_k
\]
satisfies \(A_i\leq F\), and convexity gives
\[
 \Phi_N(A_i)
 \leq W_i^{-1}\sum_{k:i\in S_k}w_k\Phi_N(G_k)
 \leq\sum_k w_k\Phi_N(G_k).
\]
Together with \(\Phi_N(A_i)\leq\Phi_N(F)\), this says
\(\Phi_N(A_i)\leq H\), hence \(G\geq A_i\).  Averaging on \(R_i\),
\[
 |R_i|^{-1}\int_{R_i}G
 \geq |R_i|^{-1}\int_{R_i}A_i
 \geq\sigma.                                             \tag{F5}
\]
Thus \(G\) is a genuine common column, not a fractional-threshold column.
Taking near-minimizers in (F3)--(F5) proves the fractional-subadditivity
inequality
\[
 \boxed{\displaystyle c_a(S)\leq\sum_k w_kc_a(S_k).}    \tag{F6}
\]

The same proof covers a countable cover whenever the right side of (F6) is
finite.  The nonnegative series in (F3) is measurable; Tonelli proves (F4).
If \(W_i<\infty\), apply Jensen to a finite partial sum after adjoining
the missing weight \(1-W_{i,n}/W_i\) at the value zero, and then pass to
the full normalized sum by monotone convergence.  If \(W_i=\infty\),
select a finite subfamily with total cover weight greater than one, use that
subfamily in the definition of \(A_i\), and bound it by the full sum in
(F3).  Choose the near-minimizer errors so that their weighted sum is
summable, then let that sum tend to zero.  An uncountable cover needs an
explicit measurable kernel of columns and is not asserted here.

This is an XOS/fractional-subadditivity statement, not submodularity.  It
also applies unchanged to any family of common positive linear constraints.
The common-anchor/common-constraint qualification is essential: at different
anchors, a column may be successful on a different rectangle, so there is no
single \(A_i\) whose average can be used in (F5).

### A strict same-orientation cross-anchor failure

The preceding qualification cannot be dropped merely because all rectangles
have the same orientation.  Here is a cap-two example in every \(N\geq2\).
Put
\[
 \varepsilon={1\over8N},\qquad t=1+\varepsilon,\qquad
 \delta={\varepsilon\over4},\qquad h={\varepsilon^2\over64N},
\]
and let
\[
\begin{split}
 A&=[-1,-\delta]\times[-\delta,0]\times[-1,0]^{N-2},\\
 B&=[-\delta,0]\times[-1,-\delta]\times[-1,0]^{N-2},\\
 a_1&=(-1,-\delta,-1,\ldots,-1),\qquad
 a_2=(-\delta,-1,-1,\ldots,-1),\\
 E&=(0,h)^N.
\end{split}                                                   \tag{F7}
\]
The shared faces of \(A,B\) may be made half-open (or assigned arbitrarily),
so they are disjoint pointwise.  Set \(G_1=t\mathbf1_A\),
\(G_2=t\mathbf1_B\), and \(F=G_1+G_2\); then \(0\leq F\leq2\).
For every \(x\in E\), the two rectangles \([a_1,x]\) and \([a_2,x]\)
are both positively oriented and respectively contain \(A\) and \(B\).
Their averages obey
\[
 {1\over |[a_j,x]|}\int_{[a_j,x]}G_j
 \geq {t(1-\delta)\over
       (1+h)^{N-1}(1+h/\delta)}=:m.                         \tag{F8}
\]
Indeed \(h<\delta\), since \(h/\delta=\varepsilon/(16N)\), and
\[
 \log\!\left((1+h)^{N-1}(1+h/\delta)\right)
 \leq {\varepsilon^2\over64}+{\varepsilon\over16N}
 \leq {33\over1024}\varepsilon.
\]
Using \(e^u\leq1+2u\) here and
\(t(1-\delta)\geq1+47\varepsilon/64\) gives, in particular,
\[
 m\geq1+\varepsilon/3.                                    \tag{F9}
\]
Thus each \(G_j\) is strictly successful on the *whole* receiver cell \(E\)
at its own anchor, and weights \(w_1=w_2=1/2\) fractionally cover \(E\).

If one nevertheless applies the same Orlicz aggregation to these columns,
then, away from irrelevant faces,
\[
 H=\min\{\Phi_N(F),\tfrac12\Phi_N(G_1)+\tfrac12\Phi_N(G_2)\}
   =q\mathbf1_{A\cup B},\qquad q={\Phi_N(t)\over2}.       \tag{F10}
\]
Since \(1+\log(1+\varepsilon)\leq e^\varepsilon\),
\[
 q\leq\tfrac12e^{N\varepsilon}=\tfrac12e^{1/8}<1<t.
\]
Consequently \(\Phi_N^{-1}(H)=q\mathbf1_{A\cup B}\) has height below one;
every rectangle average of it, hence its entire strong-maximal field, is
strictly below one.  Its modular cost nevertheless has the exact weighted
identity
\[
 J\!\left(\Phi_N^{-1}(H)\right)
 =q(|A|+|B|)=\tfrac12J(G_1)+\tfrac12J(G_2).                \tag{F11}
\]
This is precisely where the proof of (F6) breaks across anchors: for a given
receiver \(x\), the two known averages are over \([a_1,x]\) and
\([a_2,x]\), not over one common positive linear constraint.  The example
only rules out transplanting that one-anchor aggregation formula.  It does
not assert an integrality gap for a legal cross-anchor column program, a
failure of a general strong-maximal estimate, or any statement about
critical efficiency.

## Why intersection/minimum is not a success operation

Even for one common constraint, a pointwise minimum need not preserve
success.  For common anchor zero take
\[
 R=[0,1]^2,\qquad F=2\mathbf1_R,
\]
and set
\[
 G_1=2\mathbf1_{[0,3/5]\times[0,1]},\qquad
 G_2=2\mathbf1_{[2/5,1]\times[0,1]}.
\]
Both have average \(6/5>1\) on the same rectangle \(R\), whereas their
pointwise minimum is \(2\mathbf1_{[2/5,3/5]\times[0,1]}\), of average
\(2/5<1\).  Thus the common intersection, or a pointwise minimum, cannot
be treated as a source column successful for the constraints that both
columns satisfy.  This only shows that a minimum/intersection proof of
submodularity is unavailable; it does not disprove submodularity itself.

The two-dimensional example from r077 also shows the opposite interaction:
for \(R_1=[0,1]^2\), \(R_2=[0,3/4]\times[0,1]\), and
\(F=2\mathbf1_{R_2}\), one grouped column
\((4/3)\mathbf1_{R_2}\) satisfies both closed constraints more cheaply than
the sum of two separately optimized columns.  Hence no unconditional
additivity claim is available in overlapping geometry.

## A proved block-additivity range

Let \(F=\sum_{k=1}^K F^{(k)}\) have disjoint support packets \(Q_k\), and
partition the receiver constraints into \(S_k\).  Suppose every rectangle
in \(S_k\) sees source only from its own packet:

\[
 R\cap\operatorname{supp}F\subset Q_k\qquad(R\in S_k).    \tag{4}
\]

Then any feasible column can be split into its restrictions to the disjoint
Q_k.  This preserves all constraints in the corresponding block and makes
the modular exactly additive over the blocks.  Therefore

\[
 c_a\left(\bigcup_kS_k\right)=\sum_k c_{a,k}(S_k),         \tag{5}
\]

where the right side uses the restricted packet \(F^{(k)}\).  Translated
isolated copies of the whole packet/anchor/receiver configuration satisfy
the analogous decomposition: both column cost and source modular scale by
K, so normalized cost is unchanged.

Condition (4) is essential.  A separated physical packet statement without
receiver-rectangle isolation does not imply (5), because a cross rectangle
may see several packets and a grouped source can lower cost.

## Legitimate fractional column covers

The r077 bare fractional relaxation is not a legal column cover: its
fractional columns have average \(\sigma z_j\), often below one.  In
contrast, let \(G_q\) be actual columns, each genuinely successful on a
set \(E_q\) at its own true threshold, and suppose weights \(w_q\ge0\)
satisfy

\[
 \sum_{q:x\in E_q}w_q\ge1\qquad(x\in E).                  \tag{6}
\]

Then the weighted union estimate is valid:

\[
 |E|\le\sum_qw_q|E_q|
 \le\sum_qw_q B_N(\tau_q)J(G_q).                          \tag{7}
\]

The weights are used only in the union/integration inequality.  They do
not scale down a source, its success threshold, or its modular.  Thus an
integer allocation gap does not invalidate legitimate fractional-cover
upper bounds; it only forbids identifying such a cover with the bare
scaled-threshold relaxation.

## Current structural boundary

For a common anchor, (F6) closes the fractional-subadditivity/XOS question:
legal fractional covers by real successful columns can be consolidated into
one capped real column at no greater weighted modular cost.  This is fully
consistent with the warning in the preceding section: the construction is a
capped Orlicz aggregate, not a pointwise minimum or common intersection.

What remains open here is submodularity, a quantitative approximation useful
for *different* anchors, and a source-budget allocation theorem that turns a
cross-anchor geometric cover into a sharp halo upper bound.  The legal
weighted-cover inequality (7) remains valid, while the r077 bare
scaled-threshold relaxation remains invalid; (F6) does not turn a column of
average \(\sigma z\) into a successful column when \(z<1\).
