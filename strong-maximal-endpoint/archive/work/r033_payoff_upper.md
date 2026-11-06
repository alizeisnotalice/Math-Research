# r033: the correct single-payoff Bellman interface and its local obstacle

Date: 2026-11-06.  This note uses the r032 terminal payoff, not a pure
CDF:
\[
\phi_N(r)=
\begin{cases}e^{-r},&r>0,\\1-r/N,&r\le0,\end{cases}\qquad
v_t(r)=E_\mu\phi_N(r-S_t).
\]
The Bellman quantity relevant for a sharp \(R_N\) upper is
\(\sup_cK_cv_{t-1}\), or its corrected supersolution analogue.

## Exact spatial identities

Put
\[
R_t(r)=E[e^{S_t-r}1_{\{S_t<r\}}],\quad
Q_t(r)=P(S_t>r),\quad f_t=F_t'.
\]
Direct differentiation gives
\[
v_t'=-R_t-Q_t/N,\qquad
v_t''=R_t-(1-1/N)f_t. \tag{P1}
\]
Thus \(-v_t'\ge Q_t/N\), the success-side linear slope floor.  In a
region where \(R_t\) is negligible, the formal quadratic action gain is
suppressed:
\[
 {\mu^2p\over8}{(v_t'')^2\over-v_t'}
\lesssim Nf_t^2. \tag{P2}
\]
This is the genuine advantage of the single payoff over a pure CDF
source in far low-density states.

At \(c=\mu\), the exact parameter derivatives begin with
\[
\partial_cK_cv_{t-1}|_\mu
=p[K_\mu v_{t-1}(r)-v_{t-1}(r-\mu)],
\]
\[
\partial_c^2K_cv_{t-1}|_\mu
=p(2p-1)[K_\mu v_{t-1}(r)-v_{t-1}(r-\mu)]
+p\,v_{t-1}'(r-\mu). \tag{P3}
\]
If a state-Taylor expansion across a one-step displacement had a
controlled third remainder, it would suggest the local maximization
candidate
\[
 {\mu^2p\over8}{(v_t'')^2\over-v_t'}. \tag{P4}
\]

## Correction: critical-tube derivatives are Laplace-cancelled

The first draft of this note incorrectly used a pure Gamma calculation at
\(r/t=\mu>1\), which is a *large deviation for pure Gamma* and is not the
critical regime of the actual Exp--Bernoulli mixture.  The correct
identity is obtained from the positive-polynomial formulas.  If
\(S_t=G+cK\), then
\[
 {R_t(r)\over f_t(r)}
 ={1\over t}E[G\mid S_t=r]. \tag{P5}
\]
There is a sharper direct proof, requiring no posterior mean estimate.
With \(g_t^\eta\) the centered tilted density and
\(m_\eta=r/t\), untilting gives for \(u\ge0\)
\[
 {f_t(r-u)\over f_t(r)}
=e^{\eta u}{g_t^\eta(-u)\over g_t^\eta(0)}.
\]
The Fourier bounds \(g_t^\eta(0)\asymp t^{-1/2}\),
\((g_t^\eta)'(0)/g_t^\eta(0)=O(t^{-1})\), and
\(\|(g_t^\eta)''\|_\infty/g_t^\eta(0)=O(t^{-1})\) imply, for all \(u\),
\[
\left|{g_t^\eta(-u)\over g_t^\eta(0)}-1\right|
\le C{u+u^2\over t}. \tag{P6}
\]
For \(u\le\sqrt t\) this is Taylor's theorem; for \(u>\sqrt t\) it
follows by the global density upper bound.  Since
\[
R_t(r)=\int_0^\infty e^{-u}f_t(r-u)\,du
\]
and \(|\eta|<1/2\), integration of (P6) yields the uniform expansion
\[
 {R_t(r)\over f_t(r)}={1\over1-\eta}+O(t^{-1}). \tag{P7}
\]
Also \(f_t'/f_t=-\eta+O(t^{-1})\). In the lower-tail endpoint tube,
\(\eta\to0\), not a fixed large-deviation parameter. Together with
\[
 R_t'=f_t-R_t,\qquad
 v_t'''=f_t-R_t-(1-1/N)f_t', \tag{P8}
\]
this gives
\[
 D_t:=-v_t'=R_t+Q_t/N,\quad
 {v_t''\over f_t}={\eta\over1-\eta}+{1\over N}+O(t^{-1}),
\]
\[
 {v_t'''\over f_t}=-{\eta^2\over1-\eta}-{\eta\over N}+O(t^{-1})
=O(\eta^2+t^{-1})\quad(t\le N). \tag{P9}
\]
Thus the nominal optimizer in the critical tube is in fact
\[
\delta_*={\mu\over2}{v_t''\over-v_t'}
=O\left({(|\eta|+t^{-1})f_t\over D_t}\right),
\]
not order one.  There is, however, a necessary additive finite-time
remainder.  At \(\eta=0\) the relative version of the following bound is
false (the third moment survives).  The local maximum-gain lemma has the
form
\[
\sup_{|c-\mu|\le c_0}
\{K_cv_{t-1}(r)-K_\mu v_{t-1}(r)\}
\le(1+\varepsilon){\mu^2p\over8}
 {\eta^2f_t(r)^2\over R_t(r)+Q_t(r)/N}
+C_\varepsilon {f_t(r)\over t^2}. \tag{P10}
\]
The \(f_t/t^2\) term handles the centre strip and is summable in the
Poisson accounting.  The pure-Gamma identity \(R/f=r/t\) is valid, but it
must be read at that model's own mean \(r/t\simeq1\); it is not a
counterexample to (P7).

For a direct uniform third-derivative proof, use the convolution
cancellation
\[
R_t(r)=f_t(r)-f_t'(r)+
\int_0^\infty {u^2\over2}e^{-u}f_t''(r-\theta_u u)\,du,\qquad
0<\theta_u<1.
\]
Hence \(v_t'''=f_t'/N-\) the displayed remainder.  The common-tilt global
bound
\[
|f_t''(r+z)|\le C(\eta^2+t^{-1})e^{|\eta||z|}f_t(r)
\]
makes the integral finite uniformly for \(|\eta|<1/2\), proving (P9).
It also gives the parameter estimates
\[
A_t={p\mu\over2}\eta f_t+O((\eta^2+t^{-1})f_t),\qquad
B_t=-pD_t(1+o(1)),
\]
from which Young's inequality yields (P10).

## Actual remaining obstacle

The single payoff removes the r032 global-CDF linear-corrector
counterexample in the success-side far tail through (P2).  For the
uncorrected reference payoff, the all-action maximum itself admits an
exact global root characterization:
\[
D(c)=Gv(r)+e^{-c}Gv(r-c)-(1+e^{-c})v(r-c),\qquad
D'(c)=(1+e^{-c})v'(r-c)<0.
\]
The endpoint signs \(D(0)>0\), \(D(r+\mu)\le0\) give a unique global
maximizer, so far actions are not an additional issue at this reference
stage; the local root estimate identifies it near \(\mu\) in the
small-\(\eta\) tube.

What remains unproved is preservation of this structure **after adding
the corrector potential**.  For a Bellman MAX upper majorant the correct
orientation is
\[
W_t=v_t+U_t,\qquad U_t=K_\mu U_{t-1}+g_t,
\]
not the \(v-U\) direction used for the r029 CDF-minimizing construction.
Indeed the required pointwise condition is
\[
g_t(r)\ge\sup_c\{K_cv_{t-1}(r)-K_\mu v_{t-1}(r)
 +(K_c-K_\mu)U_{t-1}(r)\}. \tag{P11}
\]
The unproved state-localized all-control comparison is therefore for
\[
\sup_c\{K_c(v_{t-1}+U_{t-1})(r)-K_\mu(v_{t-1}+U_{t-1})(r)\}
\]
must be bounded by the source in (P10) plus a summable global remainder.
The bare global-root monotonicity uses \(v'<0\), which need not survive
\(W'=v'+U'\).  This, rather than far actions for the uncorrected payoff,
is the current concrete Bellman obstacle.

Equivalently, the still-open comparison is
\[
\sup_c\{K_c(v_{t-1}+U_{t-1})-K_\mu(v_{t-1}+U_{t-1})\},
\]
not a local CDF approximation.  This requires a nonlinear or
state-localized Bellman argument; it is not supplied by the presently
known corrector estimates.

Accordingly no sharp \(R_N\) upper, and no statement about \(C_N\), follows
yet.

## Final audit of formal r033 §§4--5

Date: 2026-11-06.  The common-tilt calculation in formal §4 has the right
orientation and retains the indispensable
\(d_t=R_t+Q_t/N\) denominator.  The all-\(u\) step in (15) is justified by
the same split used above: Taylor for \(u\le\sqrt t\), and the global
tilted density bound for \(u>\sqrt t\), followed by integration against
\(e^{-(1-\eta)u}\).  The displayed convolution cancellation then gives
the third-derivative estimate and explains why an additive
\(Cf_t/t^2\) centre-strip term is necessary in (18).

One proof-compression should be expanded to make (18) fully rigorous.
The unique global root from \(D_c<0\) does not by itself prove that the
root lies in the tight action interval on which the third derivative is
controlled.  Add the following localization sentence: on the moderate
small-\(\eta\) band, (15) implies only the uniform lower bound
\(d_t\ge c f_t\) (not \(d_t\asymp f_t\), which fails in sufficiently deep
lower tails), and
\[
h_t=C(|\eta|+t^{-1})f_t/d_t=o(1).
\]
Taylor's formula at \(c=\mu\pm h_t\), with \(C\) large, gives opposite
signs for \(\partial_cK_cv_{t-1}\).  Strict decrease of the global
residual then puts its unique root in \([\mu-h_t,\mu+h_t]\).  With this
one line, the parameter Taylor and Young estimate prove (18), including
the \(\eta=0\) additive term.

Formal §5 is correctly labelled rough.  The Gamma-shift domination (20)
has the right stochastic direction for the decreasing remaining-distance
payoff, and the exponential bound (21) follows from
\(\phi_N(r)\le e^{-\theta r}\) for \(\theta\ge1/N\).  Neither helper is
used as, or claimed to be, a sharp \(O(\log N)\) upper bound.
