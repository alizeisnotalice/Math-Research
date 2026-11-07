# r063 — identical angular responses and the missing excess hypothesis

## Outcome

The natural deep scalar packet has identical angular fibre responses, so its
exact constrained fibre argmax is constant whenever the capped scalar
response theorem applies. Thus labels alone cannot force variance. The
existing actual soft/ramp family proves physical strict reception for this
mechanism, but its cited results do not establish a fixed-cap family with
score \(\log Q_N-o(\log N)\). Accordingly, this note does **not** claim an
exact-Q counterexample to argmax-based P2. It identifies the hypothesis a
valid P2 must add: excess produced by non-common actions, or a proof that
score excess over \(Q_N\) entails such feedback excess.

There is nevertheless a genuine countertest at the broader critical budget
scale \(N\mu-d_N+O(\log N)\): the capped linear ramp in §5 has that budget,
strict physical reception after an arbitrarily small threshold decrease, and
the constant argmax in (2). Thus a P2 quantified merely over all critical
budget inputs is false; this is still not a \(Q_N\)-sharp statement.

Write \(L=\log N\), and let \(\mathsf{Sc}(f)\) denote the normalized actual
score of a strict physical witness.

## 1. The deep scalar packet has identical angular responses

Let \(a_i=-\log|x_i|\) in the deep cube and take the fixed radial input
\[
 f_D(x)=h_N\!\left(\sum_{i=1}^N a_i-D\right).
\]
At the radial boundary \(u=0\), disintegrate the angular variables into a
probability space \((Y,\nu)\). The depth-limit coarea weight is constant
on each angular simplex, while both the radial input and its anchored
rectangle kernel depend only on the radial coordinate. Hence the angularly
averaged fibre response is exactly
\[
 F_y(c)=E\,h_N(S_{N,c})=:F(c)
 \qquad\hbox{for \(\nu\)-a.e. }y. \tag{1}
\]
There is no angular dependence to hide in a chosen label. The radial
receiving-cone factor is \(\int_0^\infty e^{-u}\,du=1\) in the depth
limit, so the normalization does not change this conclusion.

Suppose now that \(h_N\) is capped and satisfies the r062 scalar
common-window theorem. If it has one successful action, strict concavity
makes the constrained argmax unique. Equation (1) then gives
\[
 a_y=a_*\quad\hbox{for \(\nu\)-a.e. }y,\qquad
 \operatorname{Var}_\nu(a_y)=0. \tag{2}
\]
When the corresponding vector theorem applies, the same identity gives a constant argmax vector and
zero trace covariance. Thus an argmax-based P2 cannot derive variance from
angular labels alone; it needs an additional non-common-action premise.

The actual soft/ramp construction in r017 proves strict physical reception
after bounded truncation and the same deep limit for a fixed scalar action.
r018 computes its score, but explicitly records that this is a lower-bound
family and must not be identified with \(Q_N\) or \(C_N\). Therefore the
present argument is an exact identical-response packet, not a proof that a
fixed-cap, \(Q_N\)-sharp family satisfies (2). This distinction closes the
labels-versus-argmax gap.

## 2. What the labels do and do not force

The point of (1)--(2) is physical, not a choice-of-maximizer artefact. The
rectangles that certify the anchored construction already carry constant
ratios. A label-preserving bridge cannot claim that their labels vary unless
it changes the witnessing rectangles or introduces a different source of
action dependence.

If P1 is formulated existentially, it must identify which *actual feedback
gain* forces the chosen disintegration away from this constant-label
realization. The natural packet does not by itself refute an existential
P1/P2: such a theorem could select a different packet. It must, however,
prove why the constant packet is discarded or why another packet carries
the excess. Strict concavity and unique fibre maximizers remove artificial
tie selections, but do not manufacture variation where underlying actions
are common.

The r026 single-cone reduction reinforces this diagnostic: conditional
averaging and monotone rearrangement return to one scalar \(Q_N\) terminal
with one fixed ratio. Its fixed-block loss concerns a particular correlated,
single-cone model; it is not evidence that all endpoint-scale physical
witnesses carry angular dispersion.

## 3. Corrected excess hypotheses

A necessary first repair is to restrict P2 to a quantified excess over the
fixed-action benchmark, for example
\[
 \Delta_Q(f):=\log{\mathsf{Sc}(f)\over Q_N}\ge\varepsilon L,
 \qquad\varepsilon>0. \tag{5}
\]
The constant-action packet has no feedback excess, and would be excluded if
it were \(Q_N\)-sharp. But (5) alone does **not** force variance: a
hypothetical new common-action physical mechanism might also beat the
presently defined anchored \(Q_N\) benchmark. A true upper-bound proof
therefore needs one of two additional, checkable alternatives.

1. **Common-mechanism reduction.** Prove that every physical
   label-preserving representation with negligible adaptive-common regret
   reduces, without logarithmic loss, to the scalar fixed-action class
   controlled by \(Q_N\). Then (5) forces non-negligible feedback/common
   discrepancy.

2. **Excess-to-regret transfer.** Directly prove that each witness
   satisfying (5) contains a normalized fibre packet of non-negligible
   endpoint mass for which
   \[
   \mathcal R_\nu(F):=
   \int F_y(a_y)\,d\nu(y)
   -\sup_{c\in I_N}\int F_y(c)\,d\nu(y)
   \ge\gamma(\varepsilon){L^2\over N}. \tag{6}
   \]
   The packet, its mass, and every discarded halo contribution must be
   stated explicitly. The scale is forced by bookkeeping: multiplying an
   \(L^2/N\) sharing loss by \(N/L\) produces a logarithmic score effect.

An actionable form of the second alternative is the bridge estimate
\[
 \log{\mathsf{Sc}(f)\over Q_N}
 \le o(L)+C{N\over L}\,\mathcal R_\nu(F), \tag{7}
\]
valid after all P1 physical accounting. Then (5) implies (6) with
\(\gamma(\varepsilon)=\varepsilon/(2C)\) for large \(N\). Formula (7) is
not proved here; it is the exact global comparison that converts a
logarithmic score excess into microscopic sharing regret.

The intrinsic condition is positive **feedback excess**, namely (6) itself.
It excludes the constant-action construction and does not confuse a new
common physical mechanism with a feedback gain.

## 4. When regret can force variance

The r062 curvature theorem gives only
\[
 \mathcal R_\nu(F)\ge {p_\mu-o(1)\over2}L\,
 \operatorname{Var}_\nu(a_y). \tag{8}
\]
This makes variance sufficient for regret; it cannot infer a lower variance
from a lower regret.

A converse needs more. A sufficient scoped assumption is that all maximizers
are uniformly interior to \(I_N\), and
\[
 -K L\le F_y''(c)\le-\kappa L
 \qquad(c\in I_N),\qquad 0<\kappa\le K. \tag{9}
\]
Then \(F_y'(a_y)=0\), and Taylor's theorem gives
\[
 {\kappa L\over2}\operatorname{Var}_\nu(a_y)
 \le\mathcal R_\nu(F)
 \le {K L\over2}\operatorname{Var}_\nu(a_y). \tag{10}
\]
Thus (6) implies
\[
 \operatorname{Var}_\nu(a_y)
 \ge {2\gamma(\varepsilon)\over K}{L\over N}. \tag{11}
\]

Boundary maximizers require another term. For
\(\bar a=\int a_y\,d\nu\), the upper Taylor bound contains
\[
 B_\partial:=-\int F_y'(a_y)(\bar a-a_y)\,d\nu(y)\ge0,
\]
so generally only
\[
 \mathcal R_\nu(F)\le B_\partial+
 {K L\over2}\operatorname{Var}_\nu(a_y) \tag{12}
\]
is available. A corrected variance conclusion must prove an interior action
buffer, control \(B_\partial\), or work with regret directly.

## 5. A capped critical-budget constant-argmax family

This section gives the promised broad-scale countertest. Take the common
action \(c=\mu\), write \(S=S_{N,\mu}\), and define the capped monotone
linear ramp
\[
 H_q(t)=
 \begin{cases}
 0,&t<q,\\
 \min\{2,1+(t-q)/N\},&t\ge q.
 \end{cases} \tag{13}
\]
It has expectation one exactly when
\[
 F(q)={E(S-q)_+-E(S-q-N)_+\over N}, \tag{14}
\]
where \(F(q)=P(S<q)\). Indeed,
\[
 EH_q=1-F(q)+{E(S-q)_+-E(S-q-N)_+\over N}.
\]

Let \(q_0=N\mu-d_N\). The refined critical Cramér calculation gives
\[
 F(q_0)={d_N\over N}+O\!\left({L^2\over N}\right),\qquad
 E(S-q_0)_+-E(S-q_0-N)_+=d_N+O(1), \tag{15}
\]
while differentiating the difference between the two sides of (14) gives
\[
 f_{S}(q)+{P(q<S<q+N)\over N}\asymp{L\over N}. \tag{16}
\]
Consequently there is a solution of (14) with
\[
 q=q_0+O(L). \tag{17}
\]
Only this broad error scale is needed below.

Its Orlicz cost is exact. With
\[
 g_N(u)=\log\Phi_N(1+u/N)-u,
\]
one has
\[
 J(H_q)=e^{-q}I_N,\qquad
 I_N=\int_0^N e^{g_N(u)}\,du+e^{g_N(N)}. \tag{18}
\]
Here \(g_N(0)=g_N'(0)=0\), and, for \(0\le u\le N\), writing
\(h=1+u/N\),
\[
 g_N''(u)=-{1\over N^2h^2}
 -{(N-1)(2+\log h)\over
 N^2h^2(1+\log h)^2}. \tag{19}
\]
For all sufficiently large \(N\), if
\(a=[4(1+\log2)^2]^{-1}\), then
\[
 -{2\over N}\le g_N''(u)\le-{a\over N}.
\]
The resulting Gaussian upper and lower integral bounds give
\[
 I_N\asymp\sqrt N,\qquad
 -\log J(H_q)=q-\tfrac12\log N+O(1)
 =N\mu-d_N+O(L). \tag{20}
\]
Thus \(H_q\) is a cap-two, critical-budget profile. Reducing q by any
positive amount makes its expectation strictly greater than one, with an
arbitrarily small change of score; the usual fixed-(N) truncation and
deep-cube limit produce an actual strict physical witness.

For the natural angular packet of this radial input, §1 gives identical
responses \(F_y(c)=F(c)\). The r062 capped common-window theorem makes its
constrained common-window maximizer unique, so the argmax is constant in y
and has zero variance. Therefore a P2 asserted for every cap-two profile
at the sole budget scale in (20) is false. The conclusion does not say
that this family has score \(\log Q_N-o(L)\), and hence does not settle the
stronger excess formulation in §3.

## 6. Necessary replacement for P2

The defensible replacement is a dichotomy, not unconditional variance:

> **P2-excess.** After P1 and all endpoint mass accounting, either the
> witness is common-action reducible to the \(Q_N\) scalar benchmark up to
> \(o(L)\) in log score, or every witness with
> \(\log(\mathsf{Sc}(f)/Q_N)\ge\varepsilon L\) has a labelled fibre packet
> satisfying the direct regret bound (6). A variance lower bound may be
> appended only with the two-sided/interior or boundary-slope control in
> (9)--(12).

This survives the constant-action packet, identifies the global inequality
still required, and keeps the work aimed at a true \(C_N\) upper bound rather
than an invalid universal fixed-model dispersion assertion.

## Scope

Nothing above proves (7), P1, or a \(C_N\) upper bound. The conclusion is
the necessary exclusion of bare endpoint efficiency from P2 and the
specific excess-to-regret comparison that would make a corrected P2
actionable.
