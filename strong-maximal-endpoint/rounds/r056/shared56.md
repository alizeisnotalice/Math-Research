# r056: a tie-invariant, all-outer-interval sharing certificate

This note gives a usable replacement for fiber curvature/selector variance.
It is an exact regret identity plus a layer-cake lower certificate. It does
not prove that dangerous endpoint inputs satisfy the certificate.

## 1. Exact regret and a layer-cake lower bound

At a fixed receiving point and outer interval `J`, let
\[
 m_y=\sup_{I\ni x_1}F_y(I),\quad
 A_J={1\over|J|}\int_Jm_y,dy,\quad
 S_J=\sup_{I\ni x_1}{1\over|J|}\int_JF_y(I)\,dy.
\]
The exact common-action regret is
\[
 A_J-S_J=\inf_{I\ni x_1}{1\over|J|}\int_J[m_y-F_y(I)]\,dy.\tag{1}
\]
It is invariant under ties and needs no curvature.

For every height `t>0`, define
\[
 p_J(t)={|\{y\in J:m_y>t\}|\over|J|},\qquad
 p_{J,I}(t)={|\{y\in J:F_y(I)>t\}|\over|J|},
\]
\[
 D_J(t)=p_J(t)-\sup_{I\ni x_1}p_{J,I}(t)\ge0.\tag{2}
\]
Layer cake and `sup_I integral <= integral sup_I` give the rigorous,
curvature-free estimate
\[
 \boxed{\quad A_J-S_J\ge\int_0^\infty D_J(t)\,dt.\quad}\tag{3}
\]
Thus a positive-measure band of heights on which no single physical interval
covers all the adaptively active fibers gives a quantitative sharing loss.
Unlike the variance condition, (3) is insensitive to flat maxima and is
formulated directly in terms of actual interval averages.

The estimate is only sufficient, not an exact representation of regret. In
fact the layer deficit can vanish at every height while the regret is
strictly positive. At receiving point `(0,0)`, on the fixed outer interval
`J=[-1,1]`, take equal-height outer slabs and
\[
 g_1=\mathbf1_{[0,1]},\qquad
 g_2=\mathbf1_{[0,1]}+3\mathbf1_{[1,2]},
\]
\[
 f(x,y)=g_1(x)\mathbf1_{[-1,0]}(y)+g_2(x)\mathbf1_{[0,1]}(y).\tag{3a}
\]
For an inner interval `[0,r]`, the two averages are
\[
 q_1(r)=\begin{cases}1&0<r\le1,\\1/r&r>1,\end{cases}\quad
 q_2(r)=\begin{cases}1&0<r\le1,\\3-2/r&1<r\le2,\\4/r&r>2.\end{cases}
\]
Thus `m_1=1`, `m_2=2`, and
\[
 A_J={3\over2},\qquad S_J=\sup_{r>0}{q_1(r)+q_2(r)\over2}={5\over4}.
\tag{3b}
\]
So the exact regret is `1/4`. But for every height, a common interval covers
all fibers which are adaptively active at that height: use `[0,1]` below
height one, and `[0,2]` between heights one and two. Hence
\[
 D_J(t)=0\quad\hbox{for every }t>0.\tag{3c}
\]
The actual two-dimensional Orlicz cost is finite and exact,
\[
 \int\Phi_2(f)=2+3(1+\log3).
\]
This is deliberately a **fixed-outer-interval** limitation: the side outer
slab can have a larger common value, so no global strong-maximal gap is
claimed. It proves that height-by-height coverage cannot replace the full
regret (1), even for a finite-amplitude physical input.

An equivalent useful no-regret test is obtained from the epsilon-optimal
sets
\[
 K_y(\varepsilon)=\{I\ni x_1:F_y(I)\ge m_y-\varepsilon\}.
\]
If there is one interval in `K_y(epsilon)` for almost every `y` in `J`
(equivalently, the essential intersection is nonempty), then
`A_J-S_J<=epsilon`. Consequently a positive sharing theorem must exclude
large common intersections of these sets, not merely show that an arbitrary
selection of maximizers has variance. If a concrete regularity branch makes
the `K_y(epsilon)` convex in interval-endpoint space, Helly-type arguments
can test intersection; arbitrary nonnegative fibers do not have that
convexity, so it must be an explicit hypothesis rather than an inference.

## 2. Retaining every near-optimal outer interval

Let `A_* = sup_J A_J`. The following is a direct, all-outer-interval
certificate:
\[
 A_J>A_*-\eta\quad\Longrightarrow\quad
 \int_0^\infty D_J(t)dt\ge\delta.\tag{4}
\]
Then `M_2M_1f-M_Sf>=min(eta,delta)` at the receiving point. Indeed (3)
gives `S_J<=A_J-delta` for every near-optimal `J`, while outside that set
`S_J<=A_J<=A_*-eta`. This is the r054 two-level lemma with variance replaced
by a checkable, tie-invariant quantity.

More concretely, it suffices that for every near-optimal `J` there are a
height interval `T_J` and a number `rho_J>0` such that
\[
 \sup_{I\ni x_1}p_{J,I}(t)\le p_J(t)-\rho_J
 \quad(t\in T_J),\qquad |T_J|\rho_J\ge\delta.\tag{5}
\]
This is a genuine covering condition: it says that the adaptive superlevel
fibers cannot be served by one common interval over an entire value band.

## 3. A finite-grid certificate with the actual Orlicz denominator

For a step input on a finite rectangular cell complex and a fixed receiving
point, the common-average objective `S_J` may be reduced to cell breakpoints
(and the receiving coordinate): on a segment with constant fiber value, its
average is fractional linear in either endpoint, so an interior endpoint
cannot create a strict new maximum. Thus, for every grid outer interval `J`,
there is a finite payoff matrix
\[
 F_{J,y,i}=F_y(I_i),\qquad
 r_{J,i}=\sum_yw_y\{\max_\ell F_{J,y,\ell}-F_{J,y,i}\}.
\]
Exactly,
\[
 A_J-S_J=\min_i r_{J,i}.\tag{6}
\]
This is a finite common-action regret certificate, not a surrogate
independent-fiber optimization. For a step input, `M_1f(x_1,y)` is also a
step function of `y`, so the global suprema of both `A_J` and `S_J` occur at
outer grid endpoints. Hence a direct global finite certificate is obtained
by comparing all grid values of `A_J` and `S_J`. One useful finite sufficient
test is the grid disjunction
\[
 A_J\le A_*-\eta\quad\hbox{or}\quad r_{J,i}\ge\delta\ \text{for all }i.
\tag{7}
\]
This proves the corresponding global max gap directly through the finite
endpoint reduction. It does **not** by itself establish the continuous
quantifier in (4) for every fractional near-maximizing outer interval; that
would require an additional outer-interval argument.
The true denominator is retained exactly:
\[
 \mathcal B(f)=\sum_{\text{cells }C}|C|\,\Phi_n(f|_C).\tag{8}
\]
To turn a point certificate into an endpoint lower/upper comparison, one
must additionally verify the same inequalities on a receiving cell or
integrate them over a certified receiving region. Nothing in (6)--(8)
permits replacing the Orlicz cost by separate fiber budgets.

There is an important limitation: the layer coverage `sup_I p_{J,I}(t)`
need **not** be attained at an endpoint action, even when `S_J` is. A finite
layer-cake certificate therefore needs a separate arrangement analysis of
the inequalities `F_y(I)>t` (or a validated interval cover of the continuous
action set); it cannot simply reuse the endpoint matrix in (6).

For the two-slab example of r055, the endpoint actions `t=2,4` give the
matrix
\[
 \begin{pmatrix}1/2&0\\1/4&1/2\end{pmatrix}.
\]
Its adaptive value is `1/2` and its common value is `3/8`. But the actual
continuous action `t=3` covers both fibers through height `1/3`. Hence the
unrestricted layer deficit is only `D_J(t)=1/2` on
`t in (1/3,1/2)` (and zero otherwise), so
\[
 \int D_J=1/12<1/8=A_J-S_J.\tag{9}
\]
This strict inequality is a useful warning: (3) remains valid but is not
generally sharp, and a layer certificate has to retain its continuum of
common actions.

## 4. What this does and does not overcome

The r055 flat cube/product examples have `D_J=0` despite arbitrarily
variable maximizer selections, exactly as they should. The two-slab
construction has a positive `D_J` on its combined outer interval, but side
outer intervals have `D_J=0` and are already adaptive-optimal; this is why
(4), not merely one selected-outer-interval variance calculation, is needed.

The finite-grid test (7) is potentially useful for actual bounded-amplitude
searches because it enforces all outer alternatives and keeps (8). It does
not yet provide a feedback-scale theorem: an endpoint argument would need
to show that every input whose payoff competes at the missing `g log N`
scale satisfies a uniform version of (5), or else derive (5) from an
additional structural stability theorem. The r011 rearrangement bridge and
the global Orlicz budget alone do not yield such a condition.

`shared56.wl` was actually run in Wolfram 15.0.0. It verifies the exact
values in (3b) and (9), including the two different layer-cake limitations;
the theorem statements above are proved by layer cake and the stated finite
common-average reduction, not by the computation.
