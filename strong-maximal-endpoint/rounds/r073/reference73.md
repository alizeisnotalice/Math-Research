# Action-to-reference transfer: a posterior-tail obstruction

## Verdict

There cannot be a uniform pointwise implication of the form

\[
 K_{c(x)}F(x)\ge1\quad\Longrightarrow\quad
 K_{c_{\rm ref}}F(x)\ge1-o(\log N/N)                    \tag{1}
\]

for all bounded compact physical inputs, even if the input is assembled
into a globally near-efficient construction by a disjoint small-scale
packet.  A one-coordinate rectangular example already produces a
`\Theta(\rho)` reference loss.  More seriously, the posterior-tail
mechanism below is a credible locally near-scalar-efficient obstruction:
assuming a standard uniform conditional saddlepoint/moderate-deviation
lemma, it produces a reference loss `N^{-1+o(1)}\gg\log N/N` while changing
the actual-action scalar cost by only `N^{o(1)}`.

The latter is not recorded as a finished theorem here, because the required
uniform posterior LAN/moderate-deviation comparison has not yet been proved
in the round files.  It identifies the exact missing analytic lemma and why
a deterministic common reference cannot be obtained from ordinary Orlicz
efficiency alone.

## A rigorous compact signed-box obstruction without local efficiency

Work at one receiver and put its log coordinate at zero.  The source height
below is the same on **every** source sign label.  For the resulting full
label increment Z, write

\[
 a(c)={2\over1+e^{-c}},\qquad
 \mathbb P_c(Z\le b)=q_c(b)=1-a(c)e^{-b}\quad(b>c),
 \qquad p_c={1\over1+e^c}.
\]

Thus `\partial_c\log q_c(b)=-p_c a(c)e^{-b}/q_c(b)`.
Fix `b>2`, `\lambda>0`, choose

\[
 B_N=\log {a(c_-)N\over\lambda},\qquad
 E_N=[0,b]\times[0,B_N]^{N-1},                            \tag{2}
\]

and let `c_-=c_0-\delta`, where `\delta\asymp\rho` remains inside the
action window.  Under `c_-`,

\[
 P_-:=\mathbb P_{c_-}(Z\in E_N)
 =q_{c_-}(b)(1-\lambda/N)^{N-1}.
\]

Choose `b,\lambda` so that `P_-^{-1}` is below a prescribed fixed cap, and
put the height `(1+\delta^2)P_-^{-1}` on this same increment box in every
source label.  Then

\[
 \mathbb E_{c_-}F_N=1+\delta^2>1,
 \qquad
 \mathbb E_{c_0}F_N=(1+\delta^2){P_{c_0}(E_N)\over P_-}
 =1-\kappa\delta+O(\delta^2+\delta/N).                  \tag{3}
\]

where

\[
 \kappa=p_{c_-}\left({a(c_-)e^{-b}\over q_{c_-}(b)}
                       +\lambda\right)+o(1)>0.
\]

The calculation follows by differentiating `\log P_c(E_N)`; the last
`N-1` coordinates contribute `-p_c\lambda+o(1)`.  Thus the deterministic
reference loss is `\Theta(\rho)`, far larger than `\log N/N`.

In physical variables the support is the compact signed product annulus

\[
 \prod_{i=1}^N\bigl([-1,-e^{-b_i}]\cup[e^{-b_i},1]\bigr),
 \qquad b_1=b,\quad b_i=B_N\ (i\ge2),
\]

and the height is fixed.  The kernel itself enforces the reachability
condition `z_i\ge cJ_i` for negative sources; no source-label restriction is
made.  It is therefore a genuine bounded compact-source counterexample to
an unrestricted pointwise transfer.  The strict actual inequality in (3)
and continuity of the box probability also give a positive neighboring
receiver halo.
It is *not* locally scalar-efficient: its conditional profile is nonzero at
small sums, so its full scalar cost has an order-one contribution.  Hence it
does not refute a theorem which assumes a genuinely near-optimal local
Orlicz budget.

## Global near efficiency alone still cannot give (1)

Let a compact near-efficient reference packet be supported at a macroscopic
scale.  Translate the preceding rectangle in every log coordinate by `D`.
Its physical modular and the physical measure of its receiver packet are
both multiplied by `e^{-ND}`.  For `D` large it is disjoint from the main
packet and has arbitrarily small relative modular.  At its translated
receiver the main packet is unreachable, because source-to-receiver
increments are one-sided nonnegative.  Use action `c_-` on that packet and
the reference action on the main packet.

The combined input remains globally near-efficient, whereas the translated
successful packet retains the `\Theta(\rho)` failure in (3).  Thus global
efficiency cannot yield an all-successful-receivers common-reference mean
statement.  At most it can support a theorem after an aggregate exceptional
set/packing argument.

## A sharper, locally efficient candidate: posterior tails

Let `c_a=c_{\rm ref}+\Delta c`, with

\[
 \Delta c=D_N\rho={\log(N)^{1/4}\over\sqrt N}.
\]

At the optimized shifted reference the first action derivative of the
scalar calibration vanishes, so the intended logarithmic action penalty is
`O(D_N^2\log N)=O(\sqrt{\log N})`; the weaker
`O(D_N\log N)=O(\log(N)^{3/4})` bound is already enough for `N^{o(1)}`.
Take a fixed-cap KKT scalar profile `h_{*,c_a}(s)`, with a
slightly smaller fixed cap so that the repairs below remain under the given
cap.  Truncate its very high radial tail at a polynomially negligible
probability level; this makes the physical source compact and costs `1+o(1)`.

For central sums `s=N\mu+O(\sqrt N)`, let `K` denote the jump count under
the posterior at `(c_a,s)`.  The exact posterior weights are

\[
 w_k(c,s)={N\choose k}(s-ck)_+^{N-1}.                     \tag{4}
\]

At `c=\mu`, the conditional normal approximation has

\[
 \operatorname {Var}(K\mid S=s)\sim {Np_\mu\over\mu^2},
 \qquad \partial_c\mathbb E(K\mid S=s)\sim-{Np_\mu\over\mu}. \tag{5}
\]

Thus the action displacement separates the two posterior means by
`\asymp\sqrt N\log(N)^{1/4}`, or `\asymp\log(N)^{1/4}` posterior standard
deviations.  On a central `s` slab of fixed positive sum probability, choose
an upper-tail event `A_s` with actual conditional probability

\[
 \delta_N={1\over N(\log N)^2};                           \tag{6}
\]

using a fractional boundary layer if needed.  Define the physical source,
in the receiver's increment coordinates and full label space, by deleting
that tail and repairing its height:

\[
 F(s,K)=h_{*,c_a}(s){\mathbf1_{K\notin A_s}\over1-\delta_N}               \tag{7}
\]

on the central slab, and use the unmodified `h_{*,c_a}` elsewhere.  Under
the actual action its conditional radial response is exactly
`h_{*,c_a}(s)`.  The Orlicz scaling cost is at most

\[
 (1-\delta_N)^{-N}=
 \exp(O((\log N)^{-2}))=1+o(1).                            \tag{8}
\]

The intended posterior moderate-deviation comparison is

\[
 \mathbb P_{c_{\rm ref}}(A_s\mid S=s)
 \ge\delta_N\exp\{c\log(N)^{3/4}+O(\sqrt{\log N})\}     \tag{9}
\]

uniformly on a positive-probability central slab, for some `c>0`.  If (9)
is proved, then the reference response of (7) loses

\[
 N^{-1}\exp\{c\log(N)^{3/4}+O(\sqrt{\log N})\}
 \gg {\log N\over N}.                                    \tag{10}
\]

The radial profile's own action-to-reference mean error is at most
polynomial-logarithmic over `N` in the intended calculation, and is
dominated by (10).  Hence this mechanism would disprove even a pointwise
local-near-efficiency transfer.

What remains to make this rigorous is specific, not conceptual:

1. a uniform conditional local CLT and moderate-deviation estimate for (4)
   through standardized deviations `\asymp\sqrt{\log N}`;
2. a uniform LAN likelihood-ratio expansion under
   `\Delta c=\log(N)^{1/4}/\sqrt N`, with error
   `o(\log(N)^{3/4})` on that tail; and
3. a direct bound on the reference mean of `h_{*,c_a}` strong enough that it
   is negligible compared with (10).

Without these three lemmas, (7)--(10) are an exploratory counterexample
mechanism rather than a proved theorem.

## What an alternative proof would have to pack

The ordinary global modular cannot prohibit either a small-scale packet or,
if (9) is established, the posterior-tail deletion: deleting source mass
reduces rather than increases the Orlicz budget.  Therefore a proof cannot
first demand a common deterministic reference mean at every successful
receiver.

The needed replacement is an action-indexed source packing.  At the actual
policy action one has the captured band quantity

\[
 A_{u,c}^{I}=e^{-u}\int_Ie^{-s}
 [K_c^s\Phi_N(F_{u+s})](X,\sigma)\,ds.                    \tag{11}
\]

An upper-bound argument must retain `c=\pi(u,X,\sigma)` in (11) and prove a
joint Carleson/stopping estimate that packs these action-labelled source
uses before replacing or integrating actions.  The required gain must be
better than the total-variation positive action envelope
`O(\sqrt{N\log N})`; otherwise it exceeds the allowed exponent.  Equivalently,
one needs either robust success over actions as an added hypothesis, or a
policy/source allocation theorem showing that action-sensitive posterior
tails such as (7) occupy a controlled aggregate receiver budget.  Fixed
reference (L^2) response and ordinary Bregman efficiency do not provide
that packing by themselves.

