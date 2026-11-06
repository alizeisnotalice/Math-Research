# r032: sharp feedback upper — candidate dual corrector and first blocker

Date: 2026-11-06.  Goal: determine whether the r029 lower-corrector can
be reversed to prove
\[
\log\mathcal R_N\le N\mu-d_N+\frac{3+\mu}{12}\log N+o(\log N).
\]
The analysis below finds the correct local coefficient, but also identifies
the first presently unclosed all-control estimate.  It is not a proof of
the upper bound.

## The necessary Bellman inequality

Fix a baseline \(P=K_a\), with \(m=E Z_a\), \(b=a-m\), and
\(V_t=P^t1_{[0,\infty)}=F_t\).  If
\[
U_t=PU_{t-1}+g_t,\qquad W_t=V_t-U_t,
\]
then for an arbitrary action \(c\),
\[
K_cW_{t-1}(r)-W_t(r)
=g_t(r)-\{PV_{t-1}(r)-K_cV_{t-1}(r)\}
-(K_c-P)U_{t-1}(r). \tag{U1}
\]
Therefore the all-policy CDF lower bound required for an \(R_N\) upper
would follow from the pointwise, all-state condition
\[
g_t(r)\ge\sup_c\left[
 PV_{t-1}(r)-K_cV_{t-1}(r)+(K_c-P)U_{t-1}(r)\right]. \tag{U2}
\]
It would imply \(K_cW_{t-1}\ge W_t\) for every action, hence
\[
P_\pi(S_N\le q)\ge F_N(q)-U_N(q)
\]
for every predictable policy.  This is the exact reverse of r029's
single-policy inequality; simply reversing its signs without (U2) is
invalid.

## Local calculation and the only possible sharp source

Take \(a=a_N=\mu-\frac{\sigma}{2}\sqrt{L/N}\), or more cleanly the
common-negative-tilt Laplace optimizer \(c_\eta\).  At remaining time
\(t\), write \(x=r-tm\).  The r029 derivative computation gives locally
\[
 {A_t(r)\over p f_t(r)}
=b-\frac{x}{2t}+o(t^{-1/2}),\qquad B_t(r)=pf_t(r)(1+o(1)). \tag{U3}
\]
Consequently the largest local quadratic reduction of the failure CDF is
\[
\sup_{\delta\ {\rm local}}\{PV_{t-1}-K_{a+\delta}V_{t-1}\}
=\frac p8\,\frac{(x-2tb)^2}{t^2}f_t(r)(1+o(1)). \tag{U4}
\]
The square, including its cross term, is mandatory.  Under the
baseline bridge at the endpoint, \(E(X_t/t)\simeq2b\), so replacing it by
\(b^2+x^2/t^2\) creates a spurious order-\(L\) source.

The corresponding candidate source for a dual corrector is therefore
\[
 g_t(r)=(p/8+\varepsilon) f_t(r)
 \min\left\{\frac{(r-tm-2tb)^2}{t^2},1\right\}
\varepsilon\,{f_t(r)\over t}, \tag{U5}
\]
with a smooth version of the minimum.  The final small positive term is
needed near \(x=2tb\), where the leading square vanishes but corrector
perturbation errors do not.  A common-tilt bridge would give
\[
 U_N(q)=(\mu p/8+o(1))f_N(q)\sum_{t}t^{-1}
\]
on the late time window, exactly the coefficient
\(g=(\mu-1)/8\) needed for a matching \(R_N\) upper.

For the exact Laplace baseline \(c_\eta\), the condition
\(\partial_cM_c(\eta)=0\) yields
\(e^{\eta c_\eta}/M_{c_\eta}(\eta)=1\).  This is the natural formulation
of the bridge-bias cancellation; it is stronger and cleaner than
approximating \(a_N\) alone.

## What can be proved locally

On a fixed standardized band, and for
\(|c-a|=O(t^{-1/2})\), r029's common-tilt estimates give
\[
 |(K_c-P)U_{t-1}(r)|
\le o(1)\,f_t(r)/t
\]
provided the source is restricted to a fixed band.  Then the added
\(\varepsilon f_t/t\) in (U5) absorbs the error, and (U2) follows for
local actions.  The r022 global one-step functional
\[
\mu-m(c)+\lambda\{v(c)-\mu\}/2
\]
has its unique minimizer near \(c_\eta\), so it is a plausible mechanism
to exclude actions separated from \(a\); but that theorem concerns an
independent sum and does not yet control the extra
\((K_c-P)U_{t-1}\) term in (U2).

## First explicit blocker: loss of a uniform relative derivative bound

For the source (U5), write \(g_j=f_j\rho_j\).  Unlike r029's fixed-band
source, \(\rho_j\) is not uniformly \(O(1/j)\): it is order one outside
the bridge-centered tube.  In the early part of the Poisson convolution
\[
U_t=\sum_{j\le t}P^{t-j}g_j,
\]
the r029 proof of relative derivatives uses
\(\int g_j^\eta\rho_j\lesssim1/j\) after tilting at the *current* state.
For (U5), such a bound is available only when the current tilt is the
terminal bridge tilt, because then
\[
\int f_j^\eta\rho_j
\lesssim {1\over j}+(m_\eta-m-2b)^2. \tag{U6}
\]
At an arbitrary state/tilt the second term can be order one and the
source is saturated.  Summing the early terms can then be order \(t\),
not \(O(\log t)\).  Thus the global estimates
\[
 |U_t^{(\ell)}(r+u)|
\lesssim L\,t^{-\ell/2}e^{C|u|/\sqrt t}f_t(r)
\]
used in r029 are not presently true for this global sharp source.
Without them, \((K_c-P)U_{t-1}\) cannot be made a lower-order error in
(U2), especially for actions or states outside the fixed central band.

This is a genuine all-state issue, not a controlled-occupation issue:
r029 only needed the derivative estimate where its selected policy
operates, whereas an upper requires (U2) at every \(r\) and every \(c\).
The r021 Gaussian barrier controls all states, but its initialization and
outer-band patch cost an unspecified \(\Theta(L)\) shift.  Grafting it
onto (U5) would lose the coefficient being sought.

## Other required, still open patches

1. **Far actions.**  A global action-exclusion lemma must combine the
r022 mean/variance deficit with the corrector perturbation, proving that
far \(c\) make the right side of (U1) nonnegative.  Local Taylor and the
quadratic completion alone do not do this.

2. **Small remaining times.**  Starting the sharp source at
\(t_0=L^4\) loses only \(O(\log L)\) of its intended harmonic bridge
mass, but setting \(U_t=0\) before \(t_0\) does not satisfy (U2): the
baseline action need not minimize \(K_cV_{t-1}\) there.  A sharp finite
time Bellman/block initialization with only \(o(L)\) terminal effect is
needed.  r031's same-\(c\) block Fisher regularity does not by itself
control an arbitrary adaptive block.

Accordingly, the r029 Poisson construction determines the candidate
sharp coefficient and a local dual comparison, but it does **not** yet
close a matching all-policy upper bound.  The first concrete missing
mathematics is a globally relative Poisson-derivative (or nonlinear
Bellman supersolution) estimate that preserves the tilted
bridge-cancellation in (U4) while controlling every state and action.

## Exact global quadratic-source bridge formula

The preceding obstruction is not merely an artifact of a crude derivative
bound.  For a fixed baseline \(c\), write
\[
S_t=G_t+cK,\qquad G_t\sim\Gamma(t,1),\quad
K\sim{\rm Bin}(t,p_c).
\]
Conditional on \(S_t=q\), the posterior weights of \(K\) are proportional
to
\[
 {t\choose K}p_c^K(1-p_c)^{t-K}
 \frac{e^{-(q-cK)}(q-cK)_+^{t-1}}{\Gamma(t)}. \tag{E1}
\]
For \(t>1\), define
\[
C_t(q)=E\left[
 { (q-cK)^2\over t+1}
 {c^2K(t-K)\over t-1}\ \middle|\ S_t=q\right]. \tag{E2}
\]
If \(S_j\) is the partial sum of the first \(j\) increments, exchangeable
conditioning, the Dirichlet law for the Gamma allocation, and the
hypergeometric law for the Bernoulli allocation give the exact identities
\[
E(S_j\mid S_t=q)={jq\over t},\qquad
 {\rm Var}(S_j\mid S_t=q)={j(t-j)\over t^2}C_t(q). \tag{E3}
\]
There is no omitted posterior-variance term in (E3): conditional on \(K\),
the conditional mean of \(G_j+cK_j\) is already \(jq/t\), independent of
\(K\).

Consider the uncapped all-state quadratic source
\[
g_j(s)=f_j(s){p\over8}
 { (s-jm-2jb)^2\over j^2},\qquad b=c-m. \tag{E4}
\]
For its Poisson sum \(U_t=\sum_{j\le t}P^{t-j}g_j\), Bayes and (E3) give
the exact formula
\[
{U_t(q)\over f_t(q)}
={p\over8}\left[
t\alpha^2+{C_t(q)\over t}(H_t-1)\right],
\qquad
\alpha={q\over t}-m-2b. \tag{E5}
\]
Indeed, the mean square contributes \(j^2\alpha^2/j^2\) at every level,
and
\(\sum_{j=1}^t(t-j)/j=t(H_t-1)\).

If sources are activated only for \(m_0\le j\le t\), the equally exact
version is
\[
{U_t(q)\over f_t(q)}
={p\over8}\left[
(t-m_0+1)\alpha^2+
{C_t(q)\over t^2}
\{t(H_t-H_{m_0-1})-(t-m_0+1)\}\right]. \tag{E6}
\]
In particular Jensen already gives the global envelope
\[
{U_t(q)\over f_t(q)}\ge {p\over8}\,t\alpha^2. \tag{E7}
\]

Formula (E7) is a precise version of the all-state failure above: the
linear corrector \(V-U\) cannot be a uniformly small perturbation of the
baseline CDF away from the tilted bridge tube \(\alpha=O(t^{-1/2})\).
It also suggests that a viable sharp upper may need to exponentiate or
shift the CDF by a nonlinear state potential involving \(t\alpha^2\),
rather than use a linear Poisson subtraction.  Establishing a Bellman
inequality for such a nonlinear shift remains open; (E5)--(E7) do not
prove the sharp upper.

## Exact terminal-centering reference

For a specified terminal \(q\) and horizon \(N\), choose the reference
action \(c\) by
\[
 2c-m(c)={q\over N}. \tag{E8}
\]
This has a unique solution for every target \(q/N>-1\). Indeed
\(m(0)=1\), \(2c-m(c)\to\infty\), and
\[
{d\over dc}\{2c-m(c)\}=2-m'(c)\ge3/2,
\]
since \(m'(c)\le p_c\le1/2\). At the relevant lower-tail target,
\[
c=\mu+\frac12(q/N-\mu)+O((q/N-\mu)^2). \tag{E9}
\]

With \(b=c-m(c)\), the quadratic-source centre is
\[
m(c)+2b=2c-m(c)=q/N. \tag{E10}
\]
Therefore the exact formula (E5), evaluated at terminal \((N,q)\), has
\(\alpha=0\) identically:
\[
 {U_N(q)\over f_N(q)}
 ={p_c\over8}{C_N(q)\over N}(H_N-1). \tag{E11}
\]
This removes the problematic \(N\alpha^2\) term from the terminal bridge
accounting without deleting early times.

Let \(\eta\) be the common tilt satisfying \(m_{c,\eta}=q/N\). In the
present regime \(c\to\mu\), \(\eta\to0\). Under this tilt the posterior
of \(K\) given \(S_N=q\) has Gamma-density weights around its tilted
binomial centre. Standard local lower bounds for the total tilted Gamma
density and its global upper bound give posterior variance \(O(N)\).
Applying this to (E2) yields
\[
 {C_N(q)\over N}\longrightarrow
 {1\over(1-\eta)^2}+c^2p_{c,\eta}(1-p_{c,\eta})
=v_{c,\eta}\longrightarrow\mu. \tag{E12}
\]
Thus the whole-time quadratic source has terminal mass
\[
 U_N(q)/f_N(q)=(\mu p_\mu/8+o(1))\log N,
\]
precisely the desired \(g\log N\) scale.

This is a strong simplification of the candidate terminal accounting, not
a completed upper proof. The unresolved requirement is still the
all-state/all-action Bellman inequality (U2); terminal centering does not
restore the global relative derivative estimates that fail in (E7).

## Audit of formal r032 §§4--6

Date: 2026-11-06.  §4 correctly states a **sufficient but unproved**
single-terminal Bellman MAX assertion.  The terminal
\(\phi_N(r)=e^{-r}\) for \(r>0\), \(1-r/N\) for \(r\le0\) has the right
remaining-distance orientation for (9).  The text does not reverse the
r029 CDF construction or assert the sharp upper as proved.

§5 is correct.  In particular the posterior in (12) has no residual
exponential/Bernoulli factor: \(p_c^k(1-p_c)^{t-k}\) and the
\(e^{-(q-ck)}\) Gamma factor cancel their \(e^{\mp ck}\) parts.  Dirichlet
and hypergeometric conditioning give (14), and summing
\((t-j)/j\) gives (15)--(16).  The all-state envelope (19) is consequently
a valid counterexample to a uniform logarithmic linear-corrector bound.

There is one small proof-compression issue in §6, not a mathematical
counterexample.  To justify
\(C_N/N=v_{c,\eta}+O(N^{-1/2})\), explicitly record the posterior
domination.  Under the common tilt, let \(K_0\sim{\rm Bin}(N,p_{c,\eta})\).
The tilted total Gamma density at its mean is at least \(c_0/\sqrt N\),
while its global supremum is at most \(C_0/\sqrt N\).  Hence
\[
 P(K=k\mid S_N=q)\le C\,P(K_0=k),\qquad
 E[(K-Np_{c,\eta})^2\mid S_N=q]\le CN. \tag{A1}
\]
Writing \(A=(1-\eta)^{-1}\) and
\(\Delta=K-Np_{c,\eta}\), one has
\[
q-cK=NA-c\Delta,\quad
E K(N-K)=N^2p_{c,\eta}(1-p_{c,\eta})+O(N^{3/2}). \tag{A2}
\]
Substitution in (13) gives
\[
C_N/N=A^2+c^2p_{c,\eta}(1-p_{c,\eta})+O(N^{-1/2})
=v_{c,\eta}+O(N^{-1/2}). \tag{A3}
\]
Thus the asserted conclusion is correct once (A1)--(A3) are added; the
bare statement “conditional variance \(O(N)\), directly gives” is too
compressed on its own.
