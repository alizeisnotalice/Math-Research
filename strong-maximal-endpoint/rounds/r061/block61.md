# r061 — global angular-band test for the common-action block mechanism

**Outcome.** The desired global weak theorem is not proved here.  There is,
however, a rigorous obstruction to the tempting reduction that applies r046
at every receiver using only
\(
J_u\le e^{u_1+u_2}J(H)
\).  A far, long angular band has a receiving set whose weighted mass is
proportional to its actual global modular cost, while this *crude* relative
budget can be arbitrarily large.  Thus low-weight receivers cannot be
discarded by a global supremum bound of that form.  The example does **not**
contradict a `Q_N polylog(N)` theorem: its score is just a fixed scalar
one-dimensional score times an absolute factor.

This note concerns two equal blocks \(b=N/2\), a common fixed action
\(c=\mu\), and \(0\le H\le2\).  It is not a statement about general
rectangles or \(C_N\).

## 1. The band and its receiving strip

Let \(T_1,T_2\) be the two independent block sums, and put
\[
 S=T_1+T_2,\qquad D=T_1-T_2.
\]
At \(c=\mu\), \(ES=N\mu\), \(ED=0\), and
\(
\operatorname {Var}D=N\mu
\).  Fix
\[
 a=3/2,\qquad P(S\ge q_N)=3/4,
 \qquad W=8\sqrt{N\log N}.
\]
For any \(L>\max\{W,q_N\}\), set \(R=L+q_N\), and on \(\mathbb R_+^2\) define
\[
 H_{R,L}(t_1,t_2)
   =a\,\mathbf1_{\{t_1+t_2\ge R,\ |t_1-t_2|\le L\}}. \tag{1}
\]
This is an actual cap-two input.  Write \(x=u_1+u_2\) and
\(y=u_1-u_2\).  Every point in
\[
 E_{R,L}=\{R-q_N\le x\le R-q_N+1,\quad |y|\le L-W\} \tag{2}
\]
is received with a uniform strict margin.  Indeed, on \(S\ge q_N\) and
\(|D|\le W\), (1) is hit (the extra `+1` in (2) only lowers the required
threshold for \(S\)).  Therefore
\[
 T_\mu H_{R,L}(u)
 \ge a\{3/4-P(|D|>W)\}
 \ge {3\over2}\left({3\over4}-{\mu\over64\log N}\right)>1 \tag{3}
\]
for all sufficiently large \(N\).  The final inequality is the elementary
Chebyshev bound; an exponential tail would only improve it.  Notice that
the strip lies in the positive receiver quadrant: the least coordinate is
at least \(W/2\), since \(R=L+q_N\).

No independence between \(S\) and \(D\) was used in (3).

## 2. Exact modular and receiving-volume identities

Use \(x=t_1+t_2\), \(y=t_1-t_2\), so
\(dt_1dt_2=dx\,dy/2\).  Because \(R\ge L\), the angular constraint in
(1) is not clipped by the positive quadrant for any \(x\ge R\).  Hence,
with \(\Phi_N(h)=h(1+\log^+h)^{N-1}\),
\[
 \boxed{J(H_{R,L})=L\Phi_N(a)e^{-R}.} \tag{4}
\]
The same change of variables gives the exact target mass
\[
 \boxed{\nu(E_{R,L})=(L-W)e^{-(R-q_N)}(1-e^{-1}).} \tag{5}
\]
Thus
\[
 {\nu(E_{R,L})\over J(H_{R,L})}
 =(1-W/L)(1-e^{-1}){e^{q_N}\over\Phi_N(a)}. \tag{6}
\]
It has no dependence on how far out the band is placed.  In particular this
family is a calibration, not an above-\(Q_N\) counterexample: (6) is the
score of the elementary scalar terminal with amplitude \(a\), up to the
displayed fixed factor.

## 3. Exact failure of the *crude* relative-budget localization

For a receiver \(u\), let its actual upper-cone budget be
\[
 J_u=\int_{\mathbb R_+^2}e^{-s_1-s_2}\Phi_N(H(u+s))\,ds.
\]
The universally valid global estimate is
\(
J_u\le e^{u_1+u_2}J(H)
\).  At the central receiver
\(u_1=u_2=(R-q_N)/2\), it gives
\[
 e^{u_1+u_2}J(H_{R,L})=L\Phi_N(a)e^{-q_N}. \tag{7}
\]
But direct integration in the translated positive cone is
\[
\boxed{J_u=\Phi_N(a)\big((q_N+1)e^{-q_N}-e^{-L}\big).} \tag{8}
\]
The imposed \(L\ge q_N\) is essential for this displayed branch.  In the
other branch \(L<q_N\), the exact integral is instead
\(J_u=L\Phi_N(a)e^{-q_N}\), so the crude bound happens to be sharp.
Consequently the ratio between the available global upper bound and the
actual local budget is asymptotic to \(L/(q_N+1)\), which is unbounded.
Taking, for example, \(L=N^K\) or much larger is allowed after moving
\(R=L+q_N\); (4)--(6) remain unchanged apart from the harmless
\(1-W/L\).

This pinpoints the gap in the proposed angular partition proof.  r046 can
be invoked if an **actual** local budget is at the endpoint scale.  It
cannot be invoked merely from the global input budget through (7), and the
weighted size of these remote receiver strips is not negligible compared to
the global cost by (6).  This example does not show that (8) itself is endpoint-critical; the
fixed-amplitude band is not near the sharp endpoint benchmark. It shows
only that the global-to-local budget ratio is unbounded. A successful
argument must combine efficiency selection with local truncation or an
integrated substitute, rather than assume that the crude bound is sharp.

## 4. A rigorous reduction that survives arbitrary angular geometry

For any nonnegative \(H\), any common action, and a received set
\(E=\{T_cH>1\}\), write \(r=(u_1+u_2)/2\), \(v=(u_1-u_2)/2\), and let
\(
 \ell_E(r)=|\{v:(r,v)\in E\}|
\).  Define the angular fibre mass
\(
 h(r)=\int_\mathbb R H(r,v)\,dv
\), with the convention that both sides may be infinite.  Positivity and
Fubini yield the exact inequality
\[
 \ell_E(r)
 \le\int_\mathbb R T_cH(r,v)\,dv
 =E\,h\!\left(r+{S\over2}\right). \tag{9}
\]
The last identity uses only translation in the angular coordinate; no
independence of \(S,D\) is required.  Hence
\[
 \nu(E)\le2\int_\mathbb R e^{-2r}E h(r+S/2)\,dr. \tag{10}
\]
This is an exact, globally valid fibre reduction.  It also exposes the
remaining essential analytic task: the actual modular controls
\(
\int e^{-2r}\int\Phi_N(H(r,v))\,dv\,dr
\), whereas (10) contains the *linear* fibre mass \(h\).  Convexity gives
no reverse estimate from that modular to \(\Phi_N(h)\) on an infinite
angular fibre (a long constant band is precisely the obstruction).  Moreover, the raw right side of (10) can be infinite: exchanging its
nonnegative integrals produces the factor `E exp(S)=infinity` for the
fixed Gamma--Bernoulli law. Thus (9)--(10) alone cannot close a
`Q_N polylog(N)` weak bound.

## 5. Local Wolfram check

[`block61.wl`](../../code/block61.wl) uses the exact finite Gamma--Binomial mixture
for \(S\), then evaluates (4)--(8) for \(N=64,128\) and very long bands.
The saved [`block61.json`](block61.json) records the numerical quantile,
the rigorous Chebyshev lower average, and the large ratio between (7) and
(8).  These are numerical diagnostics of the exact formulas, not interval
certificates and not evidence for a global theorem.

## Status

The requested global weak estimate remains open.  What has been ruled out
is a proof that deletes all angularly remote receivers by applying r046
through only \(e^{u_1+u_2}J(H)\).  The viable next step is an **integrated
local-budget** decomposition, retaining the local cone restriction in (8),
or a vector-valued modular inequality that upgrades (9) while charging
angular fibre length linearly rather than through \(\Phi_N(\int H\,dv)\).

## 6. Budget-efficient inputs are universally close to one

The preceding angular issue does not prevent a useful, independent
stability statement.  This section records a proof that uses no KKT
equation, monotonicity, or product ansatz.  Put
\[
 \rho_N=\sqrt{\log N/N},\qquad
 d_N=\sqrt{\mu N}\,z_N,\qquad
 \bar\Phi(z_N)=\sqrt{\mu/N}\,z_N.
\]
Thus \(d_N=O(\sqrt{N\log N})\).  Constants below may depend on a fixed
cap \(A\), a fixed number of balanced blocks, and the fixed `O(log N)`
slack, but not on \(H\) or the individual near-\(\mu\) actions.

### Fixed balanced blocks

Let \(T=(T_1,\ldots,T_m)\) be fixed balanced independent blocks with
\(|c_j-\mu|\le C\rho_N\), and let
\[
 J(H)=\int_{\mathbb R_+^m}e^{-\sum t_j}\Phi_N(H(t))\,dt
 \le\exp\{-N\mu+d_N+O(\log N)\}.
\]
If \(H\ge0\) and \(EH(T)\ge1\), then
\[
 \boxed{E|H(T)-1|=O(\rho_N).} \tag{11}
\]
If in addition \(H\le A\) for a fixed cap, then also
\(
\|H(T)-1\|_2=O(\rho_N^{1/2})
\).

Here is the proof.  Take the r046 central box
\[
 C_a=\{ |T_j-b_j\mu|\le a\sqrt{b_j\log N}\quad(1\le j\le m)\}.
\]
For arbitrarily large fixed \(K\), r046 (with its actual modular budget)
gives
\[
 E[H(T)\mathbf1_{C_a^c}]≤ N^{-K}. \tag{12}
\]
The stable block-density bound used in r046 gives, on this box,
\[
 R(t):=e^{\sum t_j}\prod_jf_{b_j,c_j}(t_j)
 \le\exp\{N\mu+C_a\sqrt{N\log N}\}. \tag{13}
\]
Polynomial density factors have simply been absorbed into the displayed
constant.  With \(\delta=B\rho_N\), the elementary pointwise implication
\[
 h>1+\delta\quad\Longrightarrow\quad
 h\le {\Phi_N(h)\over[1+\log(1+\delta)]^{N-1}} \tag{14}
\]
therefore shows
\[
 \begin{split}
 E[H(T)\mathbf1_{C_a\cap\{H>1+\delta\}}]
 &\le {e^{N\mu+C_a\sqrt{N\log N}}J(H)
       \over[1+\log(1+\delta)]^{N-1}}\\
 &\le\exp\{d_N+C_a\sqrt{N\log N}+O(\log N)
 -(N-1)\log[1+\log(1+B\rho_N)]\}. \tag{15}
 \end{split}
\]
The elementary inequality
\[
 \log[1+\log(1+\delta)]\ge {\delta\over1+\delta} \tag{16}
\]
is enough here: its difference has derivative with the sign of
\(\delta-\log(1+\delta)\) and vanishes at zero.  Thus choose fixed \(B\)
larger than the coefficient in the first two terms of (15).
The right side is then \(N^{-K}\) (indeed exponentially small on the
\(\sqrt{N\log N}\) scale).  Combining (12), (15), and the contribution
at most \(\delta\) from \(1<H\le1+\delta\) gives
\[
 E(H-1)_+\le B\rho_N+o(\rho_N). \tag{17}
\]
As \(EH\ge1\),
\(
E(1-H)_+\le E(H-1)_+
\), proving (11).  The cap gives
\(
|H-1|^2\le A|H-1|
\), proving the second.

The same estimates also prove that \(EH(T)<\infty\): below the threshold
the central contribution is at most \(1+\delta\), (15) controls the
overshoot, and (12) controls the outer box.  Thus no finiteness assumption
has been smuggled into the argument.

### Feedback terminal extension (unbounded \(L^1\))

There is a one-dimensional version that applies to every predictable
feedback policy with actions \(c_i\in[0,\infty]\) (the endpoint is the
natural limiting law).  Let \(S_N\) be its terminal sum, and assume
\[
 h\ge0,\quad
 \int_0^\infty e^{-t}\Phi_N(h(t))dt
 \le e^{-N\mu+d_N+O(\log N)},\quad Eh(S_N)\ge1. \tag{18}
\]
Then, without a height cap,
\[
 \boxed{E|h(S_N)-1|=O(\rho_N).} \tag{19}
\]
If \(h\le A\) for a fixed cap, it also follows that
\(\|h(S_N)-1\|_2=O(\rho_N^{1/2})\).  The cap cannot be removed from
this \(L^2\) corollary.

The one-step law obeys both \(\sup_zq_c(z)\le2\) and
\(q_c(z)\le2e^{-z}\).  The latter, iterated conditionally over the
simplex, gives the all-policy inequality
\[
 f_{S_N}(t)\le {2^Ne^{-t}t^{N-1}\over\Gamma(N)},\qquad
 R_\pi(t):=e^tf_{S_N}(t)\le {2^Nt^{N-1}\over\Gamma(N)}. \tag{20}
\]
Conditioning only before the last step and using \(\sup_zq_c(z)\le2\)
gives the softer but useful \(f_{S_N}(t)\le2\).  For clarity, the uniform
MGF calculation is explicit.  For finite \(c\),
\[
 \ell_c(\theta):=\log Ee^{\theta Z_c}
 =-\log(1-\theta)+\log(1+e^{-(1-\theta)c})-\log(1+e^{-c}).
\]
Its derivative at zero is \(m(c)\le\mu\), and, for
\(0\le\theta\le1/2\),
\[
 \ell_c''(\theta)={1\over(1-\theta)^2}
 +{c^2e^{(1-\theta)c}\over(1+e^{(1-\theta)c})^2}\le8. \tag{21}
\]
Indeed, with \(y=(1-\theta)c\), the second term is at most
\(4y^2e^{-y}\le16/e^2<3\), while the first is at most four.  The same
bound holds at \(c=\infty\) by its direct limit.  Taylor's formula gives
\[
 \log E(e^{\theta Z_c}\mid\mathcal F)
 \le\mu\theta+4\theta^2. \tag{22}
\]
Iterating it and optimizing the Chernoff parameter yields, uniformly over
all predictable policies,
\[
 P(S_N>N\mu+x)\le e^{-x^2/(16N)}\qquad(0\le x\le4N). \tag{23}
\]
There is deliberately no lower-tail centralization: a policy may have mean
strictly below \(N\mu\).

The cap is unnecessary for the upper moderate-deviation tail.  Fix the
constant \(C_0\) in the budget assumption and set
\(s=N\mu-d_N-C_0\log N\), so that \(J(h)\le e^{-s}\).  Choose
\(\lambda=e^sN^{-(K+2)}\).  By (20), Stirling, and \(s=N\mu+o(N)\),
\[
 \left({R_\pi(t)\over\lambda}\right)^{1/(N-1)}
 \le C_K(1+t/N), \tag{24}
\]
The uncapped conjugate bound from r046,
\[
 \Phi_N^*(y)\le y\exp\{y^{1/(N-1)}\}\qquad(y\ge0), \tag{25}
\]
combined with Young--Fenchel **once** on
\(E_a=\{S_N>N\mu+a\sqrt{N\log N}\}\), gives
\[
 \begin{split}
 E[h(S_N)\mathbf1_{E_a}]
 &\le\lambda J(h)+\lambda\int_{E_a}
    \Phi_N^*(R_\pi/\lambda)e^{-t}\,dt\\
 &\le N^{-K-2}+E\left[e^{C_K(1+S_N/N)}\mathbf1_{E_a}\right]. \tag{26}
 \end{split}
\]
For fixed \(p>1\), take \(N\) sufficiently large that
\(pC_K/N\le1/2\); then (22) gives
\(Ee^{pC_KS_N/N}=O_{p,K}(1)\).  Hölder and (23) therefore make the last
term in (26)
\[
 O_{p,K}\!\left(N^{-(1-1/p)a^2/16}\right). \tag{27}
\]
Choose \(a\), after \(p,K\), to make it \(O(N^{-K})\).

On the remaining upper-central event
\(S_N\le N\mu+a\sqrt{N\log N}\), the likelihood ratio to
\(e^{-t}dt\) is at most \(2e^{N\mu+a\sqrt{N\log N}}\).  The same
threshold calculation (14)--(17), with \(\delta=B\rho_N\), now uses
(26)--(27) for the tail rather than a cap.  First choose \(a\), then
choose \(B\) large relative to \(a\) and the coefficient of \(d_N\).
This proves (19).  Under a cap, \(|h-1|^2\le A|h-1|\) gives the \(L^2\)
claim.

This is genuine all-input/profile stability at the endpoint modular scale,
but it is not score-sensitive enough to control an \(O(1/N)\) action or
kernel perturbation.  Its \(L^2\) error is only \(O(\rho_N^{1/2})\) even
when a cap is imposed.  No improvement of the unresolved logarithmic
coefficient follows from (11) or (19) alone.


## 7. Uniform amplitude localization and a simultaneous cap reduction

The proof contains a stronger conclusion before `Eh>=1` is used. For every
fixed `K>0` there is a fixed `B_K` such that, uniformly over all policies and
all nonnegative profiles with the stated budget,
\[
 E_\pi[h(S_N)\mathbf1_{\{h(S_N)>1+B_K\rho_N\}}]\le N^{-K}        \tag{28}
\]
for sufficiently large `N`. Choose the weighted outer-tail exponent larger
than `K`, then choose the central threshold constant in (15) to make its
contribution smaller than any power. Thus (28) does not assume that the
policy receives the input. The fixed-balanced-block proof gives the same
statement uniformly in its stated near-mu action window.

Consequently let `epsilon=N^(-K)`, `K>2`, and define one profile
\[
 \widehat h=\frac{\min\{h,1+B_K\rho_N\}}{1-\epsilon}.
\]
For **every** policy which had `E_pi h>=1`, (28) gives
`E_pi widehat h>=1`; the same single capped profile works simultaneously
for that whole set of policies. Its height is `1+O_K(rho_N)` and its actual
modular budget obeys
\[
 \log\frac{J(\widehat h)}{J(h)}\le O_K(N^{1-K}).                 \tag{29}
\]
For completeness, for `a>=1,x>=0` one has
`Phi_N(ax)<=a(1+log a)^(N-1)Phi_N(x)`. If `x>=1` this follows by comparing
`1+log(ax)` with `(1+log a)(1+log x)`; if `x<1`, use either linearity below1
or `log(ax)<=log a` above1. Therefore the logarithm of the scaling factor
is at most `N log a`. Apply this with `a=(1-epsilon)^(-1)` after truncating,
which can only decrease `J`. This proves (29). Zero-budget inputs do not
receive, so the displayed ratio is used only when `J(h)>0`.

This is a reduction of the abstract terminal/block profiles, not a
pointwise amplitude truncation theorem for arbitrary physical spatial
inputs. It preserves the stated receiver averages, without asserting a
uniform strict margin when the original average was exactly1.
