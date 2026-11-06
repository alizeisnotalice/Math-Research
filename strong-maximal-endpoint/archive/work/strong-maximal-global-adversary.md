# Adversarial audit: global rearrangement / feedback bridge

Date: 2026-11-06. This audits the proposed general-input bridge only, not the
separate entropy estimate for the feedback problem.

## Verdict

I found no counterexample to the bridge. There is one real qualification which
must be explicit: after the first coordinate maximal operation an intermediate
function is generally not globally \(L^1\) (already \(M1_{[0,1]}\) has a
\(1/|x|\) tail). Thus an \(L^1\)-only rising-sun/rearrangement lemma cannot be
iterated as written. A spatial-cutoff extension repairs this completely. With
it, the reduction is valid for all nonnegative inputs with finite
\(\Phi_n\)-modular.

## Safe one-step lemma

For nonnegative \(q\in L^1(\mathbb R^n)\), fibrewise rising-sun gives
\[
sD\le F_q(D)+A_s(q),\quad
D=|\{M_iq>s\}|,\quad
F_q(d)=\int_0^d q^*,\quad A_s(q)=\int(q-s)_+ . \tag{1}
\]
The two one-sided halos \(H^\pm\), with \(I=H^+\cap H^-\), satisfy
\[
s|H^+\cup H^-|=\int_{H^+\cup H^-}q+\int_I(q-s).
\]
The terms at right are bounded by \(F_q(D)+A_s(q)\).

For decreasing \(g\), define
\[
Tg(x)=\sup_{0\le r\le1}{\int_0^xg+\int_0^{rx}g\over(1+r)x}. \tag{2}
\]
Its even realization \(g^\#(y)=g(2|y|)\) obeys
\[
(Mg^\#)(x/2)=(Mg^\#)^*(x)=Tg(x). \tag{3}
\]
At \(y>0\), right intervals have average at most \(g^\#(y)\); a left interval
may be changed to \([-ry,y]\), \(0\le r\le1\), by adding missing inner mass or
removing outer lower mass. For \(g^\#\), (1) is an equality. Concavity of
\(F_g(d)+A_s(g)-sd\), together with \(D\ge|\{q>s\}|\), gives
\[
(M_iq)^*\le T(q^*)\quad\text{a.e.} \tag{4}
\]
Use the rightmost zero if there is a plateau.

The required extension is:

> **Cutoff lemma.** Equation (4) holds for every nonnegative locally integrable
> \(q\), with extended values allowed on the right.

Apply (4) to \(q_R=q1_{[-R,R]^n}\). As \(R\uparrow\infty\), \(q_R\uparrow q\),
\(M_iq_R\uparrow M_iq\), and distribution functions and decreasing
rearrangements converge at continuity points. Monotone convergence in the
integrals in (2) gives (4). This handles non-\(L^1\) intermediates.

For a bounded compact seed, the r011 report can instead establish every needed
finiteness directly. If \(0\le f\le H1_{[-R,R]^n}\), then after coordinates
\(1,\ldots,k\),
\[
M_k\cdots M_1f\le
H\prod_{j\le k}M_1 1_{[-R,R]}(x_j)
\prod_{j>k}1_{[-R,R]}(x_j).
\]
Each displayed one-dimensional factor is at most one and tends to zero at
infinity; every positive level set of the product is therefore bounded. The
next coordinate is still unused and has compact fibre support. Thus the
rising-sun identity applies fibrewise, while all halo and excess integrals in
(1) are finite. This is an equally valid alternative to the cutoff lemma.

## Iteration and approximation

For bounded compact \(f\), (4) and order preservation of \(T\) give
\[
(M_n\cdots M_1f)^*\le T^nf^*. \tag{5}
\]
The uncentred strong rectangular maximal operator is at most this iterated
operator; equality is not true in general because an inner coordinate interval
in the iterated operator may vary with the outer integration variable. Thus
\[
C_n\le R_n:=
\sup_{g\downarrow}{|\{T^ng>1\}|\over\int_0^\infty\Phi_n(g)},\qquad
\Phi_n(t)=t(1+\log^+t)^{n-1}. \tag{6}
\]

For arbitrary nonnegative finite-modular \(f\), use
\(f_N=\min(f,N)1_{[-N,N]^n}\). Then \(f_N\uparrow f\), all finite-rectangle
averages (and hence coordinate and strong maximal functions) increase to their
limits, strict superlevel sets increase, and
\(\int\Phi_n(f_N)\uparrow\int\Phi_n(f)\). This proves (6) without applying an
\(L^1\) statement directly to a nonintegrable intermediate.

## Feedback and boundary normalization

Put \(h(u)=g(e^{-u})\). Then
\[
Tg(e^{-u})=\sup_{0\le r\le1}E h(u+Z_r),\quad
Z_r=Y+(-\log r)J, \tag{7}
\]
where \(Y\sim{\rm Exp}(1)\), \(J\sim{\rm Bernoulli}(r/(1+r))\), independently;
at \(r=0\), \(Z_0=Y\). Finite-horizon Bellman recursion gives
\[
T^ng(1)=\sup_{\pi}E_\pi h(S_n),\qquad S_n=\sum_{k=1}^nZ_{r_k}, \tag{8}
\]
with predictable controls. This is not a fixed-ratio independent-sum model.
For bounded Borel \(h\), rational actions plus measurable epsilon-greedy
selection prove (8); unbounded finite-cost \(h\) follows by monotone
truncation. No exchange of \(\sup\) and expectation is assumed.

Define
\[
J(h)=\int_0^\infty e^{-u}\Phi_n(h(u))\,du,\qquad
V_n(h)=\sup_\pi E_\pi h(S_n).
\]
Then
\[
R_n=\sup_{\substack{h:[0,\infty)\to[0,\infty)\ {\rm nondecreasing}\\
                     J(h)<\infty,\ V_n(h)\ge1}}J(h)^{-1}. \tag{9}
\]
If \(D=|\{T^ng>1\}|\), choose \(a<D\), dilate by \(a\), and erase input on
\(x>1\). All stages evaluated at receiver \(x=1\) are unchanged because their
arguments only decrease. Its cost is
\[
a^{-1}\int_0^a\Phi_n(g)\le a^{-1}\int_0^\infty\Phi_n(g).
\]
Let \(a\uparrow D\). Conversely set \(g(e^{-u})=h(u)\) on \(0<x\le1\), and
\(g=0\) above one. If \(V_n(h)>1\), the strict level set has measure at least
one. If \(V_n(h)=1\), use \((1+\eta)h\) then let \(\eta\downarrow0\). Bounded
profiles suffice: scale first, truncate next, then use monotone convergence.

## Numerical check

MATLAB R2026a: a random 160-cell nonnegative step input, exhaustive discrete
uncentred intervals, and the cell analogue
\(Tg(m)=\max_{0\le k\le m}(P_m+P_k)/(m+k)\) gave
\[
\min_m\{Tg(m)-(Mf)^*(m)\}=0,\qquad\#\{m:Tg(m)<(Mf)^*(m)\}=0.
\]
This is only a sanity check; (1)--(4) are the proof.

## Scope

The bridge is safe after the cutoff qualification. It does not itself control
the adaptive supremum in (9); a fixed-control tail calculation cannot replace
that subsequent entropy argument.

# Independent audit of r012: adaptive entropy upper bound

Date: 2026-11-06. I independently audited the r012 report, including the
conditional change of measure, path likelihood sign, profile decomposition,
and r011 promotion. Verdict: **no mathematical blocker found** for
\[
\log C_n=n\mu-\sqrt\mu\,\sqrt{n\log n}+O(\sqrt n).
\]
The proof does not use its numerical calibration as evidence.

## Conditional tilt and likelihood ratio

For a predictable action \(r_i=e^{-c_i}\), let \(P_i\) be the original
conditional kernel. Its down-tilt is
\[
{dQ_i\over dP_i}(z)={e^{-\lambda_i z}\over E_{P_i}e^{-\lambda_iZ_i}},
\qquad m_{Q,i}=E_{Q_i}Z_i .
\]
Its entropy and the path density are exactly
\[
k_i=-\lambda_i m_{Q,i}-\log E_{P_i}e^{-\lambda_iZ_i},\qquad
\log{dQ\over dP}=\sum_i k_i-
\sum_i\lambda_i(Z_i-m_{Q,i}). \tag{A}
\]
Thus the r012 sign in (7) is correct. The last sum is a \(Q\)-martingale, and
the proof uses it under \(Q\), as required.

If \(m(r)>\mu-\delta\), put
\(\Delta=m(r)-(\mu-\delta)\in(0,\delta]\). Since
\[
\Delta=\int_0^\lambda v_r(a)\,da,\qquad v_r(a)\ge(1+a)^{-2},
\]
\(\lambda\le\Delta/(1-\Delta)\le2\delta\). Combining the variance-defect
bound with the Lipschitz bound gives \(v_r(a)\ge\mu-\varepsilon_\delta\).
Parameterizing entropy by mean loss then yields
\[
k_i=\int_0^\Delta\lambda(b)\,db
\le{\Delta^2\over2(\mu-\varepsilon_\delta)}. \tag{B}
\]
This remains valid when the action is adaptive.

With \(\delta=(d+2\log n+4\sqrt V\sqrt n)/n\), tilted conditional means are
at most \(\mu-\delta\), so the \(Q\)-martingale fluctuation puts
\(S_n\le n\mu-d-2\log n\) with probability at least \(15/16\). Also,
\[
\operatorname{Var}_Q\sum_i\lambda_i(Z_i-m_{Q,i})\le4Vn\delta^2.
\]
Controlling its absolute value by Chebyshev and intersecting the events gives
the stated \(7/8\) event. On it, (A) upper-bounds \(dQ/dP\), hence lower-bounds
\(P\). This verifies the otherwise easy-to-reverse likelihood-ratio step.

The expansion from (8) to (9) is uniform for
\(0\le d\le2\sqrt{\mu n\log n}\): here
\(n\delta^2\varepsilon_\delta=O((\log n)^{5/4}n^{-1/4})=o(1)\), while the
remaining cross terms are \(O(d/\sqrt n+1)\).

## Profile forcing and adaptive density

From \(J(h)\le e^{-s}\) and monotonicity,
\[
e^{-t}\Phi_n(h(t))\le J(h), \tag{C}
\]
because the tail integral from \(t\) onward is at least
\(e^{-t}\Phi_n(h(t))\). This proves the pointwise profile bounds in r012:
\(h\le n^{-2}\) below \(s-2\log n\), \(h\le1\) below \(s\), and the displayed
linear surplus on the next interval of length \(n/4\). The latter uses
\(\log(1+y)\ge y/(1+y)\) for \(y=\log h\); constant 4 has slack.

Under the original adaptive law,
\[
S_n=n\mu-\sum_i(\mu-m(r_i))+M_P,\qquad E M_P^2\le nV.
\]
The predictable drift is nonnegative, so
\[
E(S_n-s)_+\le d+\sqrt{nV}. \tag{D}
\]
This proves the middle-profile estimate without asserting that the adaptive
sum has unconditional variance \(O(n)\).

The original one-step density is
\[
f_r(z)={e^{-z}\over1+r}\{1+1_{\{z\ge-\log r\}}\}\le2e^{-z}. \tag{E}
\]
Though \(r_i\) depends on the past, conditional multiplication gives the
joint-density envelope \(2^ne^{-\sum z_i}\); simplex integration gives
\[
f_{S_n}(t)\le {2^ne^{-t}t^{n-1}\over(n-1)!}. \tag{F}
\]
This is feedback-valid. On \(h(t)\ge e^{4t/n}\), (F) and
\(\Phi_n(h(t))\ge h(t)(4t/n)^{n-1}\) bound the high-amplitude tail by
\[
J(h){2^n(n/4)^{n-1}\over(n-1)!}=e^{-cn+o(n)};
\]
its exponential base is \((e/2)/B<1\). The complementary tail follows from
the conditional positive-mgf estimate using a fixed small Chernoff parameter
on \(e^{4S_n/n}1_{\{S_n>s+n/4\}}\). Both are uniform in policy and profile.

Together with \(E_\pi h(S_n)>1\), (C)--(D) force the claimed
\(P(S_n\le s-2\log n)\lesssim\sqrt{\log n/n}\). The baseline-one decomposition
is legitimate: subtract the low-tail deficit, add middle surplus, and bound
the complete far-tail expectation.

## Final comparison and promotion

For \(d=\sqrt{\mu n\log n}-K\sqrt n\), the lower-tail exponent in (9) is
\[
-\tfrac12\log n+
(K/\sqrt\mu-A_0\sqrt\mu)\sqrt{\log n}-O_K(1).
\]
For large fixed \(K\), it exceeds the profile-forced
\(O(\sqrt{\log n/n})\), a contradiction. This proves the \(\mathcal R_n\)
upper bound. The r011 bounded-profile definition is enough; strict feasibility
survives scaling then truncation, and bounded compact approximation promotes
\(C_n\le\mathcal R_n\) to all finite-modular inputs.

The r003 lower estimate is stronger than needed:
\[
n\mu-\sqrt{\mu n(\log n-2\log\log n)}-O(\log n)
=n\mu-\sqrt\mu\sqrt{n\log n}+o(\sqrt n).
\]
It therefore supplies the other half of the r012 \(O(\sqrt n)\) statement.

## Local numerical checks (non-proof)

MATLAB R2026a gave \(\mu=1.278464542761\), \(B=3.591121476669\), and at
\(r_*=B^{-1}\), \(m(r_*)=v(r_*)=\mu\). It also checked the density/envelope
ratio in (E) on a representative action grid (maximum \(1.999800019998\)).
These only detect sign or normalization errors; the proof above is analytic.

## Exposition requests, not blockers

For publication-quality presentation, display the one-line fixed-parameter
Chernoff calculation for the low-amplitude far tail and explicitly state the
baseline-one decomposition leading to (10). The cited inequalities imply both,
but displaying them makes adaptive uniformity immediate.

# Independent audit of r013: stopped policy and third term

Date: 2026-11-06. I audited the final r013 report as requested, focusing on
the stopped-policy coupling, Fan--Grama--Liu (FGL) hypotheses, far tails, and
the strict-profile promotion. **Verdict: accept; no blocking flaw found.**

## Stopping/coupling direction

Let \(\tau\) be the first step at which the original deficit
\(D_i=\sum_{j\le i}(\mu-m(r_j))\) exceeds \(D_*\), and let \(\pi_K\) follow
\(\pi\) through \(\tau\), then use \(r_*\). Coupling by the same primitive
randomness makes the paths identical on \(\{\tau>n\}\), and the stopping event
itself is identical for both policies. On that event,
\[
M_K\le-d\quad\Longrightarrow\quad S_n=n\mu-D_n+M_K\le n\mu-d.
\]
On the stopped event, original success implies
\[
M_\pi=S_n-n\mu+D_n>D_*-d\ge8\sqrt{Vn\log n}
\]
over the entire FGL window. The centered conditional positive-mgf estimate,
iterated predictably, with
\(\eta=u/(4Vn)\), gives
\[
P_\pi(M_\pi>u)\le\exp[-u^2/(8Vn)]=n^{-8}.
\]

The correct set algebra is
\[
P_\pi(\mathrm{fail})
\ge P_K(M_K\le-d)-P_\pi(\mathrm{stopped\ and\ successful}).
\]
It follows by cancelling the common stopping probability:
\[
P_K(M_K\le-d)\le P_\pi(\mathrm{fail,\ un-stopped})+P(\mathrm{stopped}).
\]
Thus r013 §5 subtracts the right event, not the usually fatal full stopping
probability. This validates (7)'s direction.

## FGL A1/A2

For
\(\xi_i=-(Z_i-m(r_i^K))/(\sigma\sqrt n)\), the uniform moment envelope gives
\[
|E(\xi_i^k\mid{\cal F}_{i-1})|
\le {k!\over2}\left({32\over\sigma\sqrt n}\right)^{k-2}
E(\xi_i^2\mid{\cal F}_{i-1}),\qquad k\ge3.
\]
This is exactly FGL A1. The stopped deficit has
\(D_n^K\le D_*+\mu-1\); Cauchy--Schwarz and the variance-defect estimate give
\[
|\langle X\rangle_n-1|
\le Cn^{-1/4}(\log n)^{1/4}=:\delta_n^2,
\]
which is FGL A2 (the theorem's parameter is \(\delta_n\), not
\(\delta_n^2\)). I checked the cited FGL v3 Theorems 2.1--2.2 and its
log-relative estimate (14): their allowed range is
\(x\le\alpha_0\epsilon_n^{-1}\), \(\delta_n\le\alpha_0\), both satisfied here.

At \(x\le2\sqrt{\log n}\), the FGL logarithmic error is
\[
x^3\epsilon_n+x^2\delta_n^2+
(1+x)(\epsilon_n|\log\epsilon_n|+\delta_n)
=O(n^{-1/8}(\log n)^{5/8})=o(1).
\]
Hence the relative, uniform normal-tail statement (6), including its lower
tail, is justified. This is importantly stronger than an absolute CLT bound.

## Improved profile and far tails

The length \(b=n^{3/4}\) is compatible with the near-threshold inversion:
for \(0\le t-s\le b\),
\[
(h(t)-1)_+\le(1+O(n^{-1/4}))(t-s)/n.
\]
The original adaptive sum has the decomposition
\(S_n=n\mu-D_n+M_\pi\), with \(D_n\ge0\) and
\(E|M_\pi|\le\sqrt{nV}\), proving the middle contribution estimate.

For the low-amplitude far tail, setting
\(\theta=cn^{-1/4}\) and \(\eta=\theta+4/n\), predictable conditional mgf
iteration gives exponent
\[
-\theta(b-d)+2Vn(\theta+4/n)^2+4\mu.
\]
Since \(\theta b=c\sqrt n\), \(\theta d=o(\sqrt n)\), choose a fixed small
\(c\) so the exponent is \(-c_1\sqrt n\). This is uniform in the policy.
The high-amplitude branch again uses the feedback-valid density envelope
\(f_{S_n}(t)\le2^ne^{-t}t^{n-1}/(n-1)!\), and its modular bound remains
\(e^{-c_2n}\). I find no hidden independence assumption in either branch.

## Third-term promotion

With
\[
d_A=n\mu-s+2\log n,
\quad
s=n\mu-\sigma\sqrt{nL}
\sigma\log L\sqrt{n/L}+C\sqrt{n/L},
\]
one has
\[
(d_A/(\sigma\sqrt n))^2
=L-2\log L-2C/\sigma+o(1).
\]
Mills plus the uniform relative lower tail yields the coefficient
\(e^{C/\sigma}/\sqrt{2\pi}\) in (10). The profile upper coefficient is
\(\sigma+o(1)\); every fixed
\(C>\sigma\log(\sigma\sqrt{2\pi})\) is excluded. Strict feasibility
\(E_\pi h(S_n)>1\) is exactly what makes the baseline-one argument a strict
contradiction. The bounded \(\mathcal R_n\) definition and r011 monotone
approximation then promote the result to all finite-modular inputs. The r008
\(Q_n\) lower construction supplies the matching liminf, so the claimed
general third term follows by \(Q_n\le C_n\le\mathcal R_n\).

Minor editorial suggestion only: retain the cancellation line in §5 (it is
the crucial coupling direction), and explicitly say that FGL's \(X_n\) is
\(M_K/(\sigma\sqrt n)\).

# Feedback log-scale mechanism: independent heuristic audit

Date: 2026-11-06. Scope: the proposed log-scale feedback mechanism, not a
new theorem about \(\mathcal R_n\) or \(C_n\).  Verdict: the Duhamel algebra
and coefficient are internally consistent; the \(\log n\) term is a credible
feedback mechanism.  A rigorous Bellman comparison or adaptive local limit
theorem is still missing, so it must not yet be claimed as a general upper or
lower bound.

## Local Bellman expansion and sign

Use the survival-value convention \(W(t,x)=P(S_t\ge t\mu+x)\), so
\(W_x<0\) in the lower-tail region.  Put \(c=-\log r=\mu+a\).  At the
mean-maximising action,
\[
m(c)=\mu+\tfrac12m''(\mu)a^2+O(a^3),\quad
v(c)=\mu+v'(\mu)a+O(a^2),
\]
with
\[
m''(\mu)=-(\mu-1)/\mu,\qquad v'(\mu)=\mu-1.
\]
In the moving frame, the action-dependent second-order term is
\[
(\mu-m(c))W_x+\tfrac12(v(c)-\mu)W_{xx}
=-\tfrac12m''(\mu)a^2W_x+\tfrac12v'(\mu)aW_{xx}+O(a^2W_{xx}+a^3W_x).
\]
In the lower tail both \(W_x,W_{xx}<0\), and its local maximizer is near
\[
a_*={v'(\mu)W_{xx}\over2m''(\mu)W_x}<0.
\]
Substitution yields
\[
\partial_tW={\mu\over2}W_{xx}
-D{W_{xx}^2\over W_x}+\hbox{higher-order terms},\qquad
D={\mu(\mu-1)\over8}. \tag{H1}
\]
Because \(W_x<0\), the nonlinear correction in (H1) is positive.  This
resolves the apparent sign ambiguity: treating \(W\) as an increasing CDF
would reverse the interpretation.

## Duhamel calculation

Let \(W_0(t,x)=\overline\Phi(x/\sqrt{\mu t})\).  Then
\[
-D{(W_{0,xx})^2\over W_{0,x}}
=D\,{x^2\over\mu^2t^2}\,p_t(x),
\]
where \(p_t\) is the \(N(0,\mu t)\) density.  Duhamel convolution and the
Brownian-bridge identity give, at terminal time \(T\),
\[
W_1(T,x)=Dp_T(x)\int_0^T
 E\!\left[{B_t^2\over\mu^2t^2}\,\middle|\,B_T=x\right]dt.
\]
Since
\[
E(B_t^2\mid B_T=x)={t^2\over T^2}x^2+\mu t(1-t/T),
\]
the induced positive quantile shift is
\[
\Delta q=-{W_1(T,x)\over W_{0,x}(T,x)}
={D\over\mu}\left[z^2+\log T+O(1)\right],\qquad
x=-z\sqrt{\mu T}, \tag{H2}
\]
after a microscopic lower-time cutoff.  Hence
\[
{D\over\mu}={\mu-1\over8}=:g.
\]
The \(gz^2\) piece agrees with the fixed-action log gain from r016/r017;
the \(g\log n\) piece is specifically the Brownian-bridge variance contribution
and is the proposed feedback surplus.

MATLAB R2026a arithmetic check: \(\mu=1.278464542761\),
\(D=0.044500880542\), \(D/\mu=g=0.034808067845\).  With
\(\beta=\sigma\sqrt{\log n/n}\), (H2) predicts
\(g(z^2+\log n)=0.225662920318\) at \(n=256\) and
\(0.302136431813\) at \(n=1024\), compared with the supplied Bellman gains
\(0.29745\) and \(0.36723\).  The residuals of about \(0.07\) are compatible
with the unspecified \(O(1)\), but are not asymptotic verification.

## What can safely be claimed

1. The formal leading feedback correction, **conditional on a valid
   parabolic Bellman expansion**, is \(g\log n\), and its coefficient is not
   an algebraic sign error.
2. The logarithm comes from mesoscopic bridge times
   \(\int dt/t\), not from the singular microscopic endpoint.  On
   \(K\le t\le n/K\), the candidate control is \(a_*=O(t^{-1/2})\) on
   bridge-typical paths; cubic Taylor errors should be summable while the
   main term remains \(g\log n+O(\log K)\).  This identifies a plausible
   route for a clipped-feedback Bellman subsolution.
3. It is not yet a verified sub/supersolution.  One needs uniform discrete
   derivative bounds for a value function on a parabolic lower-tail strip,
   control of off-strip paths, a comparison principle that preserves the
   inequality through the exact nonlocal Bellman operator, and accumulated
   remainder \(o(\log n)\).

## Current rigorous obstacle

The r013 FGL control has relative-tail error
\(\kappa_n=n^{-1/8}(\log n)^{5/8}\).  At the relevant tail its conversion to
a quantile uncertainty is roughly
\[
{\sqrt n\over\sqrt{\log n}}\kappa_n
=n^{3/8}(\log n)^{1/8},
\]
which dwarfs \(\log n\).  Thus FGL, even with the stopping trick, cannot
resolve this scale.  A proof of either a \(g\log n\) feedback gain or a
matching log-scale obstruction requires an adaptive local limit/parametrix
or a verified Bellman barrier; merely sharpening the existing relative
moderate-deviation argument will not suffice.

# Horizon-smoothed Gaussian barrier: audit

Date: 2026-11-06. The proposed \(O((\log N)^2)\) feedback upper route has no
sign-level flaw in this review. It can be a genuine theorem if the following
uniform one-step lemma is written and proved:
\[
E_rW_{k-1}(q-Z_r)\ge W_k(q)-\delta
\quad\text{for all }r,q,k. \tag{B1}
\]
Here \(W_k(q)=\Phi((q-\mu k-H_k)/(\sigma\sqrt{k+\tau_0}))\),
\(\tau_0=C_0\log N\), and \(\Delta H_k=K\log N/(k+\tau_0)\).

The orientation is correct. With \(G\sim N(\mu,\mu)\),
\[
EW_{k-1}(q-G)=W_k(q+\Delta H_k).
\]
Thus \(V_k=W_k-\epsilon-k\delta\) is a Bellman subsolution under (B1), and
\[
P_\pi(S_N\le q)\ge W_N(q)-2N^{-p}
\]
holds for every adaptive policy. The stated \(H_0\) makes
\(V_0\le1_{\{q\ge0\}}\).

In the padded moderate band, Taylor expansion at \(\mu\), with
\(D=\mu-m(r)\), gives exactly
\[
E W_{k-1}(q-Z_r)-E W_{k-1}(q-G)
=Df+\tfrac12(v(r)+D^2-\mu)f'+R. \tag{B2}
\]
The mean-loss sign is positive. The inequality
\(|v-\mu|\le C\sqrt D\), small
\(|f'/f|=O(\sqrt{\log N/(k+\tau_0)})\), and the uniform exponential moments
from the density bound \(f_r(z)\le2e^{-z}\), imply
\[
\text{right side of (B2)}\ge-C{\log N\over k+\tau_0}f.
\]
This validates the proposed variance-stability argument.

One line still must be explicit: in the right half of the band, density
decreases, so the positive semigroup shift must be lower-bounded as
\[
W_k(q+\Delta H_k)-W_k(q)
\ge c\,\Delta H_k\,W_k'(q). \tag{B3}
\]
It follows from the exact Gaussian density ratio after choosing \(C_0\) large
relative to \(K,A\). Use nested moderate strips to avoid a boundary gap. The
proper constant order is: choose \(A,p,C_0\) for small tilt exponent; choose
\(K\) against the now-uniform Taylor constant; enlarge \(C_0\) relative to
\(K,A\). This avoids a circular choice.

Far-left is valid once \(W_k\le N^{-p-2}<\delta\). Far-right is valid after
splitting at \(Z_r-\mu=(A/2)\sigma\sqrt{(k+\tau_0)\log N}\): the complementary
tail is \(\le2e^{-cA\sqrt{(k+\tau_0)\log N}}\), and otherwise the preceding
barrier is still in an \((A/2)\sqrt{\log N}\) Gaussian upper tail. Choose
\(A,C_0\) to make both \(O(N^{-p-2})\). This supplies the global glue without
independence.

Given domination for all \(q\), integrate \(P(S_N>t)\le1-W_N(t)+2N^{-p}\)
only to \(CN\), and integrate a fixed-parameter conditional-mgf Chernoff
tail afterwards. This legitimately gives
\[
E(S_N-s)_+\le E(G_N-s)_++O(N\epsilon)+e^{-cN}=d+O((\log N)^2).
\]
With \(b=\sqrt N(\log N)^2\), the profile error is
\(O((\log N)^{5/2}/N)\), while the shifted Gaussian lower CDF at
\(s-2\log N\) exceeds \(\beta=d_N/N\) by \(c(\log N)^3/N\) if the terminal
shift coefficient is dominated by \(K_0\). Since
\((\log N)^3/N\gg(\log N)^{5/2}/N\), the stated contradiction and
\(O((\log N)^2)\) upper remainder follow.

The supplied exGaussian one-step probe is directionally consistent but is
not a uniform certificate. Safe status: strong near-complete proof blueprint;
the central missing deliverable is a fully quantified proof of (B1), including
(B3), nested bands, and parameter choices.

# Final audit of r020 Gaussian barrier proof

Date: 2026-11-06. I read the completed r020 report. **Accept: no blocker
found.** The earlier blueprint's missing pieces are now supplied at the
required level for the claimed \(O((\log n)^2)\) upper remainder.

* Section 3 has the correct one-step orientation. For
  \(Y=\mu-Z_r\), the Taylor difference from the Gaussian increment is
  \(Df+\frac12(v+D^2-\mu)f'+R\). The \(Df\) sign is positive. The
  \(\sqrt D\) variance defect and the \(D^2\) term are correctly absorbed
  using small \(|f'/f|\), leaving \(-C_*Lf/\tau\). The exponential density
  envelope makes the unbounded-increment Taylor remainder uniform.
* The constant order is noncircular: first take \(C_0\) large only to make
  the exponential-moment tilt uniformly small, obtaining a \(C_*\) stable
  under further increases; choose \(K\) next; then enlarge \(C_0\) relative
  to \(K,A\) for the explicit density-ratio lower bound. The later \(K_0\)
  depends only on the resulting finite terminal-shift coefficient.
* Section 4's exact ratio compares
  \(f_{\tau+1}(x-u)\) to \(f_\tau(x)\), and its lower bound yields the
  essential right-shift gain
  \(W_k(q+h_k)-W_k(q)\ge(K/4)(L/\tau)f_\tau(x)\). This correctly covers the
  previously missing decreasing-density side. The central-to-\(2A\) padding
  and the explicit ratio avoid a transition-band gap.
* The far-left bound is immediate from \(W_k\le\delta\). In the far-right
  split, the stated two events cover the complement exactly; the exponential
  \(Z_r\) tail and normal tail are uniformly below \(\delta\). The base
  \(V_0\le1_{\{q\ge0\}}\) also checks. Conditional backward iteration is
  valid for arbitrary history-dependent (and, by conditioning, randomized)
  actions, so (7) really is all-policy domination.
* Integrating the additive CDF error only to \(T\), then using a separate
  conditional-mgf Chernoff tail, is the necessary and correct cure for the
  otherwise divergent Kolmogorov error. The \(b_N=\sqrt N L^2\) profile
  bookkeeping gives \(O(L^{5/2}/N)\), while the Gaussian \(L^2\) displacement
  changes the calibrated tail by \(cL^3/N\); the final contradiction has the
  correct direction and scale.

The MATLAB/Wolfram checks are appropriately labelled diagnostics rather than
proof. No corrective action requested.

# r021 candidate \(O(\log N)\) barrier: blocking scaling flaw

Date: 2026-11-06. The heat identity is correct:
\[
(\partial_t-\tfrac\mu2\partial_{xx})(Hf_t)
=K\,{x^2/(\mu t)+1\over t}f_t
\]
for \(H=K[x^2/(\mu t)+2\log(t/t_0)]\). Its formal sign is also desired for
the lower CDF correction \(W=\Phi-Hf\). It does **not** yield a uniform
one-step barrier with \(t_0=C_0\log N\).

Put \(x=-\sigma\sqrt t\,z\), \(z>0\), and take \(H=Kz^2\) at \(t=t_0\).
Direct differentiation gives
\[
{W_x\over f_t}
=1+{2Kz\over\sigma\sqrt t}
-{Kz^3\over\sigma\sqrt t}. \tag{C1}
\]
For \(z\asymp(\sqrt t/K)^{1/3}\), this is negative. At
\(t=C_0\log N\), \(z\asymp(\log N)^{1/6}\) is well inside a fixed
\(A\sqrt{\log N}\) moderate band, and its tail is vastly larger than
\(\delta=N^{-p}\); a cap at \(4A^2\log N\) does not affect it.

This blocks the actual-versus-Gaussian step: an action with fixed mean deficit
has first Taylor term \(D W_x\), negative of order \(f_t\) at (C1), while the
proposed heat gain is merely \(K(1+z^2)f_t/t=O(t^{-2/3}f_t)\). Thus the
advertised one-step lower inequality cannot be obtained and is expected to
fail in this regime.

The shifted alternative \(\Phi((x-H(x,t))/(\sigma\sqrt t))\) does not repair
the proposal automatically: at the same states
\[
H|x|/(\mu t)\asymp Kz^3/(\sigma\sqrt t)=O(1),
\]
then grows deeper in the band. Its nonlinear and discrete-Taylor terms are
therefore not small. A fixed \(C_0\gg K^2\) cannot fix this as
\(\log N\to\infty\); making \(C_0\) grow with \(N\) loses the advertised
\(O(\log N)\) initialization.

MATLAB R2026a illustration, with \(C_0=10^4,K=100,L=10^6,t=C_0L\), gave
\(W_x/f_t=0.1333,-1.9584,-6.0399\) at \(z=10,15,20\), while the corresponding
log left tails were \(-53.2,-116.1,-203.9\), versus
\(\log\delta=-5,000,000\).

Conclusion: r021 remains a formal identity, not an \(O(\log N)\) theorem. A
redesign needs a state/time cap before \(z\sim t^{1/6}\), plus a new argument
for the intervening tails, or a nonlinear barrier with an independently
proved one-step comparison.

## Correction: exact nonlinear shift remains viable

The linear-barrier obstruction does not apply to the exact shifted CDF
\(W(t,x)=F(t,x-H(t,x))\). For
\(H=K[x^2/(\mu t)+2\log(t/t_0)]\), direct chain rule gives
\[
(\partial_t-\tfrac\mu2\partial_{xx})W
=\left[-(H_x-\tfrac12H_x^2){y\over t}
+{\mu\over2}H_{xx}-H_t\right]F_y
\]
\[
=\left[-{K\over t}-{Kx^2\over\mu t^2}
+{2KxH\over\mu t^2}
+{2K^2x^2(x-H)\over\mu^2t^3}\right]F_y,\quad y=x-H. \tag{S1}
\]
The leading sign is coercive, \(-K(1+z^2)F_y/t\). On
\(|z|\le A\sqrt{\log N}\), \(t\ge C_0\log N\), its positive terms are a
small relative fraction after \(C_0\) is chosen sufficiently large after
\(K,A\) (terms \(KA/\sqrt{C_0}\), \(K^2A^2/C_0\), and
\(K\log(t/t_0)/\sqrt t\)). Also
\[
W_x=(1-H_x)F_y>0
\]
on this strip. Hence the earlier \(W_x<0\) objection is specific to
\(\Phi-Hf\), not to the exact shift.

The exact heat-Duhamel sign is
\[
P_1W(t-1,x)-W(t,x)
=-\int_0^1P_s[(\partial_t-\tfrac\mu2\partial_{xx})W(t-s,\cdot)](x)\,ds.
\]
Since (S1) is negative on an inner strip, this gives positive reference
gain. With a smooth nonnegative cap
\(H=K[\psi_L(z)+2\log(t/t_0)]\), equal to \(Kz^2\) on the inner strip and
constant outside a wider strip, Gaussian one-step escape has probability
\(\exp[-cC_0(\log N)^2]\ll\delta\); the capped operator is globally bounded.
This supplies the appropriate Duhamel patch mechanism.

The remaining unproved but plausible lemma is the discrete moment comparison
measured relative to the new density \(F_y\). On the inner strip, derivative
bounds for \(W_x,W_{xx},W_{xxx}\), the mean-deficit square completion, and
the uniform exponential moments should give loss
\(-C(1+z^2)F_y/t\), while Duhamel gives \(+cK(1+z^2)F_y/t\). The cap makes
the far bands manageable; in particular \(H\le C_K\log N\) and, with large
\(C_0\), the shifted far-right argument remains near one. Thus no conceptual
blocker was found for the exact nonlinear candidate, though it remains
unproved until these quantified lemmas are written.

The explicit proposed cap
\(\psi_L(z)=L\psi(z/\sqrt L)\), equal to \(z^2\) for
\(|z|\le2A\sqrt L\), with the supplied degree-seven transition polynomial
and a constant beyond \(3A\sqrt L\), has the required \(C^3\) matching:
at the inner endpoint its first three \(z\)-derivatives match \(z^2\), and
at the outer endpoint they vanish. Its derivative bounds give
\[
|H_x|=O(KA/\sqrt{C_0}),\quad |H_{xx}|=O(K/t),\quad
|H_{xxx}|=O(K/(t^{3/2}\sqrt L)).
\]
For a full Duhamel lower bound, write
\(R=(\partial_t-\mu\partial_{xx}/2)W\) exactly. On the central strip
\(R\le-cK(1+z^2)F_y/t\); in
\[
P_1W(t-1)-W(t)=-\int_0^1P_sR(t-s)\,ds,
\]
restricting to a fixed \(|U|\le1\) Gaussian subevent supplies a positive
constant fraction with density and curvature comparable to their value at
\((t,x)\). The rest of the central strip has the same favorable sign; only
escape to the cap exterior needs the global bounded-operator estimate, and
has probability \(e^{-cC_0L^2}\ll\delta\). This is the precise local-tail
decomposition needed in a formal proof.

## Final audit of r021 report

Date: 2026-11-06. **Accept: no blocker found.** The completed r021 proof uses
the exact nonlinear shift and successfully resolves the linear-barrier defect.

The capped shift has the required \(C^3\) derivative scales. Its exact
residual is negative on the inner strip, and the signed Duhamel identity gives
positive Gaussian-semigroup gain. The shifted-density Taylor comparison keeps
the mean-loss term favorable because \(W_x=(1-H_x)f>0\); square completion
and the third-order remainder have the claimed \(C(1+z^2)f/t\) scale after
the stated choice of constants.

The Duhamel proof correctly uses a fixed \(|U|\le1\) positive-probability
subevent for density comparability, retains the favorable remainder of the
wide inner strip, and charges only cap escape to the exponentially small
outer error. Far bands, base, predictable-policy backward induction, finite
CDF-error integration, and the final \(L^2/N\) versus \(L^{3/2}/N\)
contradiction also check. Consequently the terminal \(O(\log N)\) shift
supports the reported general \(O(\log n)\) upper remainder.

# Sharp-coefficient base obstruction and positive-support alternatives

Date: 2026-11-06. Any direct Gaussian-smoothed base of the r021 type has an
unavoidable order-\(\log N\) initialization cost. If
\[
W_0(q)=\Phi_G((q-b)/(\sigma\sqrt{t_0}))
\]
must satisfy \(W_0(q\le0)\le\epsilon=N^{-p}\), then
\[
b\ge\sigma\sqrt{t_0}\,\Phi_G^{-1}(1-\epsilon)
=\sigma\sqrt{2p\,t_0\log N}\,(1+o(1)). \tag{G1}
\]
Since full-strip comparison needs \(t_0\gtrsim\log N\), this forces
\(b=\Omega(\log N)\) with a coefficient containing the nonsharp smoothing
constant. Sharp local compensation alone cannot identify the terminal
\(\log N\) coefficient.

For \(q\ge0\), the exact one-step minimum CDF is
\(\inf_rP(Z_r\le q)=\tanh(q/2)\). The Gamma\((2,1)\) CDF
\[
G_2(q)=1-e^{-q}(1+q)
\]
is below it: this reduces to
\(1+q\ge2/(1+e^{-q})\), whose derivative is at least \(1/2\).
Thus Gamma shape two is a legitimate positive-support lower base at one step.

It cannot simply be joined to the ordinary Gaussian heat semigroup. For
\(q<0\) and every \(Z_r\ge0\), \(E\,G_2(q-Z_r)=0\), while
\(E\,G_2(q-G)>0\) for \(G\sim N(\mu,\mu)\). MATLAB quadrature at \(q=-1\)
gave the latter as \(0.001895398038067\). The earlier value
\(0.159659918\) was erroneous: it used the wrong sign for the mean of
\(u=q-G\), which is \(q-\mu\), not \(\mu-q\). Reproducible MATLAB:

    mu = fzero(@(x) x - 1 - exp(-x), 1.3); q = -1;
    value = integral(@(u) gammainc(u,2,'lower') .* ...
      exp(-(u - (q-mu)).^2/(2*mu)) / sqrt(2*pi*mu), ...
      0, Inf, 'AbsTol', 1e-15, 'RelTol', 1e-12);

The strictly positive corrected value still means a Gamma base plus Gaussian heat
reference violates the one-step lower comparison on the negative half-line.

A serious alternative is a positive-support Gamma reference semigroup:
increments \(\Gamma(\mu,1)\) have mean and variance \(\mu\), and an initial
\(\Gamma(\alpha,1)\) adds only \(O(1)\) terminal mean/variance. But its
generator is nonlocal, so r021 heat Duhamel and the nonlinear barrier do not
transfer automatically. An exact short-horizon Bellman base has the same
unresolved splice requirement: it must dominate a later smooth barrier with
\(o(\log N)\) location loss. These are research directions, not
sharp-coefficient proofs.

# r025 sharp-feedback local coefficient audit

Date: 2026-11-06. Let \(c=\mu+\delta\), with
\(m(c)=1+c/(1+e^c)\) and
\(v(c)=1+c^2p_c(1-p_c)\). At \(c=\mu\),
\[
m''(\mu)=-(\mu-1)/\mu=:-a,\qquad v'(\mu)=\mu-1=:b.
\]
For a lower-tail CDF \(F\), \(F_x,F_{xx}>0\). Relative to the
mean-\(\mu\), variance-\(\mu\) Gaussian reference, the leading action term is
\[
{a\delta^2\over2}F_x+{b\delta\over2}F_{xx}.
\]
Minimizing it gives
\[
\delta_*=-{b\over2a}{F_{xx}\over F_x}
=-{\,\mu\over2}{F_{xx}\over F_x}. \tag{Q1}
\]
For a normal lower tail, \(F_{xx}/F_x=-x/(\mu t)\), hence
\[
\delta_*(t,x)={x\over2t};
\]
it is negative in the lower tail and must be clipped near the short-time
boundary because \(c\ge0\).

The optimized local term is
\[
-D_{\rm fb}{F_{xx}^2\over F_x},\qquad
D_{\rm fb}={b^2\over8a}={\mu(\mu-1)\over8}. \tag{Q2}
\]
For a survival function the signs of both first two derivatives reverse in
the lower tail, the Bellman extremum reverses, and the same positive
quantile gain results. Wolfram actual symbolic differentiation gave
\(\mu=1.2784645427610737\), \(m''=-0.2178117057198001\),
\(v'=0.2784645427610738\), and
\(D_{\rm fb}/\mu=(\mu-1)/8=0.03480806784513421\).

Heat Duhamel applied to a normal terminal lower tail
\(x=-z\sqrt{\mu T}\) yields the formal quantile gain
\[
{D_{\rm fb}\over\mu}\,[z^2+\log(T/t_{\rm cut})+O(1)].
\]
The \(z^2\) portion is the fixed-action gain. The Brownian-bridge variance
integral is the genuinely feedback-specific harmonic surplus. At the
endpoint scale \(z_N^2\sim\log N\), it predicts an *additional*
\[
{\,\mu-1\over8}\log N
\]
over the fixed-action \(g z_N^2\) gain, so the formal total leading
logarithmic gain is \(2g\log N\), \(g=(\mu-1)/8\).

This is a local mechanism, not a theorem for \(\mathcal R_N\). The harmonic
integral is singular at short horizons, the control (Q1) leaves the legal
action range without clipping there, and the Gaussian initialization in r021
already costs an unspecified \(O(\log N)\) shift. Thus neither the formal
\(g\log N\) feedback surplus nor its apparent coefficient can presently be
registered as a sharp general upper coefficient. A rigorous program needs a
short-time/base splice with \(o(\log N)\) loss plus a uniform discrete
Bellman comparison for the clipped path-dependent control.

For the sharp formal coefficient, write a shifted normal ansatz as
\(H=bz^2+a\log t\). Its leading residual is
\[
(\partial_t-\mu\partial_{xx}/2)W/f
\simeq (b-a)/t-bz^2/t.
\]
The local optimized Hamiltonian is \(-gz^2f/t\), where
\(g=(\mu-1)/8\), so the sharp formal choice is \(a=b=g\), not r021's
safe \(a=2b\)-type logarithmic overcompensation. The Brownian-bridge
variance then produces the \(g\log n\) quantile gain even though the local
forcing has no constant \(1/t\) term. This gives the formal
\(g[z_n^2+\log n]\sim2g\log n\) feedback gain versus the central fixed
action, while the fixed optimized-action gain is \(gz_n^2\sim g\log n\).
It remains a formal coefficient calculation only.

## Final audit of r025 auxiliary heat equation

Date: 2026-11-06. **Accept.** The local Bellman algebra and the stated
auxiliary heat-equation theorem are correct, and their scope is accurately
limited.

For a failure CDF, the optimized local Hamiltonian is
\(-\mu g\,w_{xx}^2/w_x\), \(g=(\mu-1)/8\), with
\(\delta_*=x/(2t)\) for the normal lower tail. Thus its forcing is
\(-g z^2f/t\). Under the heat bridge from \(s\) to \(t\),
\[
E(X_s^2\mid X_t=x)={s^2x^2\over t^2}+{\mu s(t-s)\over t},
\]
so the normalized Duhamel integrand is \(z_t^2/t+1/s-1/t\). Integrating
from \(t_0\) gives
\[
U=-g[(1-t_0/t)(z^2-1)+\log(t/t_0)]f,
\]
with zero initial data and the claimed forcing. The simple-shift residual
has the stated \((b-a)/t-bz^2/t\) sign.

The report explicitly retains the necessary limitations: no uniform
discrete-feedback approximation, no short-time/clipping control, no sharp
\(\mathcal R_n\) expansion, and no \(C_n\) lower bound. The finite-grid
Bellman observations are labelled diagnostics. No correction requested.

## Exact-kernel perturbation route (scoped proposal)

Let \(P=T_\mu\) be the exact one-step operator of the fixed \(c=\mu\)
Gamma--Bernoulli increment,
\[
T_c\phi(x)=E\,\phi(x+\mu-Z_c).
\]
The exact Laplace identity
\(\partial_cM_c(\lambda)=p_c(M_c(\lambda)-e^{-\lambda c})\) gives
\[
\partial_cT_c\phi=p_c[T_c\phi-\phi(x+\mu-c)].
\]
At \(c=\mu\), \(p=(\mu-1)/\mu\), hence
\[
A\phi=p(P\phi-\phi),\qquad
B\phi=p(2p-1)(P\phi-\phi)+p\phi'. \tag{P1}
\]
Here, exactly,
\[
A=\left.\partial_cT_c\right|_{c=\mu},\qquad
B=\left.\partial_c^2T_c\right|_{c=\mu};
\]
thus (P1) is for the entire **centred backward transition** \(T_c\), not
merely for a distributional part of the kernel. The uncentred operator
\(K_c\phi(x)=E\phi(x-Z_c)\) instead obeys
\[
\left.\partial_cK_c\right|_\mu\phi
=p[K_\mu\phi-\phi(x-\mu)],
\]
\[
\left.\partial_c^2K_c\right|_\mu\phi
=p(2p-1)[K_\mu\phi-\phi(x-\mu)]+p\phi'(x-\mu).
\]
They are equivalent through \(T_c\phi(x)=K_c\phi(x+\mu)\).
Wolfram symbolic differentiation of the exact Laplace transform returned
zero residuals for both identities. Their leading terms reproduce the local
variance response and quadratic mean-loss Hamiltonian.

A prospective bounded-feedback construction is dyadic after
\(t_0=n^\alpha\): at \(t_j=2^jt_0\), freeze on the next block
\[
c_j=\mu+\lambda\,\mathrm{clip}\!\left({X_{t_j}\over2t_j}\right),
\qquad X_t=S_t-\mu t.
\]
A fixed number of blocks gives only \(O(1)\) bridge-variance gain; the
\(\Theta(\log(n/t_0))\) dyadic scales could accumulate a positive multiple
of \(\log n\). A \(\lambda\)-fraction of the locally optimal action formally
has gain fraction \(2\lambda-\lambda^2\).

This is not a theorem. The required proof uses uniform local-limit and
density-derivative estimates for powers of the fixed exact \(P\), a Duhamel
expansion based on (P1) with summable third-order remainder for \(t\ge t_0\),
and a clipped-path tail estimate. It supplies neither a \(C_n\) lower bound
nor a sharp \(\mathcal R_n\) coefficient at present.

## Exact nonanticipative first-variation identity

Let \(F_m^0(x)=P_\mu(Z_1+\cdots+Z_m\le x)\), and let a policy be perturbed
from the fixed \(c=\mu\) baseline by
\[
c_i=\mu+\varepsilon u_i(S_{i-1}),\qquad |u_i|\le U.
\]
For a fixed failure threshold \(q\), differentiation at \(\varepsilon=0\)
gives the exact identity
\[
{d\over d\varepsilon}P_\varepsilon(S_n\le q)\big|_0
=\sum_{i=1}^nE_0[u_i(S_{i-1})G_i(S_{i-1})], \tag{P2}
\]
where
\[
G_i(s)=p\{F_{n-i+1}^0(q-s)-F_{n-i}^0(q-s-\mu)\},\qquad
p=(\mu-1)/\mu. \tag{P3}
\]
This follows directly from the centred/un-centred derivative formulas above,
conditioning on the past, and then running the unchanged future baseline.
It is a genuine nonanticipative gradient; no independent-action assumption is
present.

For a block \(B_j=\{a_j+1,\ldots,b_j\}\), constrain the action to be frozen
at \(c=\mu+\varepsilon u_j(S_{a_j})\). Then (P2) becomes
\[
\sum_jE_0[u_j(S_{a_j})\Gamma_j(S_{a_j})],\qquad
\Gamma_j(s)=E_0\!\left[\sum_{i\in B_j}G_i(S_{i-1})
\mid S_{a_j}=s\right].
\]
The admissible choice \(u_j=-\operatorname{sign}\Gamma_j\) makes this
\(-\sum_jE|\Gamma_j(S_{a_j})|\), strictly negative whenever a block gradient
is nonzero. This is a clean exact finite-dimensional proof mechanism for
strict feedback improvement at a generic lower-tail threshold.

It is not yet a logarithmic theorem: controlling the second variation for
\(\varepsilon\) on the natural \(t_j^{-1/2}\) scale and proving a sharp
lower asymptotic for the dyadic sum of block gradients require the fixed-kernel
local-limit/Hessian estimates described above. A crude operator-norm Taylor
bound forces \(\varepsilon\) too small and loses the \(\log n\) scale.

## Exact greedy one-step improvement lemma (conditional form)

For the centred full transition \(T_c\), retain \(P=T_\mu\) and the
derivatives \(A,B\) from (P1). Its exact third derivative at \(\mu\) is
\[
C\phi=p(1-6p+6p^2)(P\phi-\phi)
+p(3p-2)\phi'-p\phi''. \tag{P4}
\]
Wolfram differentiation of the exact Laplace transform verified (P4) with
zero symbolic residual.

For an exact forward-policy telescope, use instead the uncentred operator
\[
K_c\phi(r)=E\phi(r-Z_c),\qquad P^K=K_\mu,\qquad
V_m=(P^K)^m{\bf1}_{[0,\infty)}.
\]
Thus \(V_m(r)=P_\mu(Z_1+\cdots+Z_m\le r)\) literally. Its derivatives are
\[
A_K\phi=p(P^K\phi-\phi(\,\cdot-\mu)),
\]
\[
B_K\phi=p(2p-1)(P^K\phi-\phi(\,\cdot-\mu))
+p\phi'(\,\cdot-\mu),
\]
\[
C_K\phi=p(1-6p+6p^2)(P^K\phi-\phi(\,\cdot-\mu))
+p(3p-2)\phi'(\,\cdot-\mu)-p\phi''(\,\cdot-\mu). \tag{P4K}
\]
These are exactly the translations of (P1),(P4).

On a state-time region where
\[
B_K V_m>0,\qquad
|\delta|=\lambda\,|A_KV_m|/(B_KV_m)\le\eta,
\]
Taylor's theorem gives, for
\(\delta=-\lambda A_KV_m/(B_KV_m)\),
\[
K_{\mu+\delta}V_m
\le P^KV_m-(\lambda-\lambda^2/2){(A_KV_m)^2\over B_KV_m}
+{M_m\over6}|\delta|^3, \tag{P5}
\]
where \(M_m\) bounds the third parameter derivative on the action interval.
If \(M_m|A_KV_m|/(B_KV_m)^2\) is suitably small, (P5) gives pointwise strict
gain.

It has an exact nonanticipative telescope: at forward step \(i\), use the
remaining threshold \(q-S_{i-1}\), apply (P5) to baseline future
\(V_{n-i}\), and choose the displayed clipped action. Conditional
expectation lowers the baseline future-failure conditional expectation by
the gain, and summation gives
\[
P_{\rm greedy}(S_n\le q)
\le P_\mu(S_n\le q)-E_{\rm greedy}\sum_i\mathrm{gain}_i.
\]
This is a rigorous finite state/time lemma once its local derivative
conditions are verified. The still-missing asymptotic input is an occupation
lower bound for the controlled process on \(BV_m>0\) regions plus a uniform
third-derivative ratio at natural dyadic action scales; without it, no
\(\kappa\log n\) conclusion follows.

## r028 pointwise greedy conditions on a moderate lower-tail band

Let \(f_t\) be the density of the fixed-\(\mu\) sum of \(t\) increments and
put \(r=t\mu+x\). Fix \(0<c_0<C_0<\infty\), \(L=\log n\), and restrict to
\[
t\ge L^4,\qquad -C_0\sqrt{tL}\le x\le-c_0\sqrt t. \tag{R1}
\]
A standard negative-saddlepoint/Fourier local-density lemma for the fixed
Gamma--Bernoulli kernel gives, uniformly on (R1),
\[
{f_{t-1}'((t-1)\mu+x)\over f_{t-1}((t-1)\mu+x)}
=-{x\over\mu t}\{1+O(L/\sqrt t)\},
\]
\[
|f_{t-1}''((t-1)\mu+x)|
\le C{1+x^2/t\over t}f_{t-1}((t-1)\mu+x), \tag{R2}
\]
and fixed bounded shifts of the density are comparable. One proof tilts by
the unique \(\lambda=-x/(\mu t)+O(L/t)\), applies the existing uniform
Fourier density/derivative estimates to the centred tilted convolution, and
then untilts; \(t\ge L^4\) makes every displayed relative error small.

Apply r027's uncentred exact derivatives to
\(\phi=V_{t-1}=K_\mu^{t-1}1_{[0,\infty)}\). With
\(p=(\mu-1)/\mu\), the first derivative satisfies
\[
a=K_\mu'\phi(r)
=p\{V_t(r)-V_{t-1}(r-\mu)\}
={p\mu\over2}f_{t-1}'((t-1)\mu+x)
+O(L/t)f_{t-1}, \tag{R3}
\]
while
\[
b=K_\mu''\phi(r)
=p f_{t-1}((t-1)\mu+x)\{1+O(\sqrt{L/t})\}>0. \tag{R4}
\]
For \(c\in[\mu/2,3\mu/2]\), the exact third-derivative identity and
bounded-shift comparison give
\[
M:=\sup_c|K_c'''\phi(r)|\le C f_{t-1}((t-1)\mu+x). \tag{R5}
\]
Consequently
\[
{|a|\over b}\le C\sqrt{L/t}\le\mu/2,\qquad
{M|a|\over b^2}\le C\sqrt{L/t}\le1
\]
for sufficiently large \(n\). r027's greedy rule is therefore activated
throughout (R1), has
\[
\delta=-a/b={x\over2t}\{1+O(L/\sqrt t)\},
\]
and supplies the pointwise gain
\[
{a^2\over3b}\ge \kappa\,{x^2\over t^2}
f_{t-1}((t-1)\mu+x), \tag{R6}
\]
for a positive constant \(\kappa=\kappa(c_0,C_0,\mu)\).

This is a precise target lemma; to call it proved in the project, the
negative-saddlepoint Fourier estimates in (R2) must be written with the
displayed uniform constants and their bounded-shift variant. Even granted,
(R6) is only pointwise. A \(\kappa\log n\) feedback theorem additionally
needs a lower bound on controlled-path occupation of (R1) across dyadic
times; none is claimed here.

## Final audit of r027 exact greedy feedback section

Date: 2026-11-06. **Accept.** In r027 §6 the full uncentred kernel
\(K_c\), rather than a distribution-only derivative, is used consistently.
The three displayed derivatives agree with the exact Laplace transform and
with the translated \(T_c\) formulas.

For \(a=K_\mu'\phi(r)\), \(b=K_\mu''\phi(r)>0\), and
\(M=\sup_{c\in[\mu/2,3\mu/2]}|K_c'''\phi(r)|\), the hypotheses
\(|a|/b\le\mu/2\), \(M|a|\le b^2\) keep
\(\mu-a/b\) in the stated interval and make the Taylor remainder at most
\(a^2/(6b)\). Hence the one-step decrease
\(-a^2/(3b)\) is correct.

Taking \(V_t=K_\mu^t1_{[0,\infty)}\), and at forward step \(i\) applying
the rule to \(V_{n-i}(q-S_{i-1})\), gives an exact conditional telescope:
the terminal potential is \(1_{\{S_n\le q\}}\), while the baseline recursion
is \(K_\mu V_{n-i}=V_{n-i+1}\). Thus (11) correctly has the controlled-path
expectation of the sum of gains. Requiring at least three future convolution
steps is a sufficient smoothness guard for the \(C^2\) third-derivative
formula.

The report properly leaves the occupation lower bound and natural-scale
third-derivative estimates open, so it does not falsely claim a
\(\kappa\log n\), \(\mathcal R_n\), or \(C_n\) theorem.

## Fixed-\(\mu\) saddlepoint Fourier lemma: proof record

Let \(Z=Z_\mu\), \(m_\lambda=E_\lambda Z\) under the negative tilt
\(dP_\lambda=e^{-\lambda Z}dP/Ee^{-\lambda Z}\), and let \(f_k\) be the
density of \(S_k=Z_1+\cdots+Z_k\). Fix \(0<c_0<C_0\), \(L=\log n\), and
\[
k\ge L^4,\qquad -C_0\sqrt{kL}\le x\le-c_0\sqrt k,\qquad s=k\mu+x.
\]
There is a unique \(\lambda>0\) with \(m_\lambda=s/k\), and
\[
\lambda=-{x\over\mu k}+O(L/k),\qquad \lambda=O(\sqrt{L/k}).
\]

Under this tilt, \(Z\) is an independent sum of an
\(\operatorname{Exp}(1+\lambda)\) variable and a bounded Bernoulli shift.
For its centred characteristic function \(\varphi_\lambda\),
\[
|\varphi_\lambda(u)|\le(1+u^2/(1+\lambda)^2)^{-1/2}.
\]
All tilted cumulants through order five are uniform. Let \(g_{k,\lambda}\)
be the density of \(S_k-km_\lambda\) under the tilt. Fourier inversion and
the local expansion \(\log\varphi_\lambda(u)=-v_\lambda u^2/2+O(|u|^3)\)
give
\[
g(0)=(2\pi kv_\lambda)^{-1/2}+O(k^{-1}),\quad
g'(0)=O(k^{-3/2}),
\]
\[
g''(0)=-{g(0)\over kv_\lambda}+O(k^{-2}),\qquad
\sup_y|g^{(j)}(y)|\le Ck^{-(j+1)/2}\quad(0\le j\le3).
\]
For the rates, compare \(\varphi_\lambda(u)^k\) on a small Fourier interval
to \(\exp(-kv_\lambda u^2/2)\) by
\(Ck|u|^3e^{-bku^2}\). The \(j=0,1,2\) integrals are respectively
\(O(k^{-1}),O(k^{-3/2}),O(k^{-2})\); the \(j=1\) Gaussian integral vanishes
by oddness. The complement is exponentially small by the displayed global
characteristic-function bound.

Untilting gives
\[
f_k(s+u)=M_\mu(\lambda)^k e^{\lambda(s+u)}g(u),
\]
and hence
\[
{f_k'(s)\over f_k(s)}=\lambda+O(k^{-1}),\qquad
{f_k''(s)\over f_k(s)}
=\lambda^2-{1\over kv_\lambda}+O(\lambda/k+k^{-3/2}).
\]
Moreover, for all real \(u\),
\[
f_k(s+u)\le C e^{\lambda|u|}f_k(s),\qquad
|f_k''(s+u)|\le C{1+L\over k}e^{\lambda|u|}f_k(s).
\]
The exponential weight is integrable against every fixed-\(c\) kernel shift
because its density is bounded by \(2e^{-z}\).

With \(k=t-1\), \(s=k\mu+x\), and
\(\phi=V_{t-1}=K_\mu^{t-1}1_{[0,\infty)}\), Taylor expansion of
\(V_{t-1}(s-(Z-\mu))\), using the last inequality as integral remainder,
proves
\[
a={p\mu\over2}f_{t-1}'(s)+O((1+L)t^{-1}f_{t-1}(s)),
\]
\[
b=p f_{t-1}(s)\{1+O(\sqrt{L/t})\},\qquad
M\le C f_{t-1}(s).
\]
Here \(a=K_\mu'V_{t-1}(t\mu+x)\), \(b=K_\mu''V_{t-1}(t\mu+x)\), and
\(M\) is r027's compact-\(c\) third-derivative supremum. Thus the r027
greedy conditions hold and its pointwise gain is
\(\ge\kappa x^2f_{t-1}(s)/t^2\) on this band. The remaining issue for a
logarithmic theorem is controlled-path occupation, not the one-step
Fourier/Taylor estimate.

## r028 sharpened single-step coefficient audit

Date: 2026-11-06. Let \(F_k(r)=\Pr_\mu(S_k\le r)\), with density
\(f=f_k\), put \(s=k\mu+x\), \(r=s+\mu\), where
\[
 -C\sqrt{kL}\le x\le-c\sqrt k,\qquad k\ge L^4.
\]
Let \(\lambda>0\) be the negative-saddlepoint parameter whose tilted
increment mean is \(s/k\). The claimed sharpened coefficient is
algebraically correct, subject to the fixed-\(\mu\) Fourier lemma.

The tilted mean has the exact form
\[
m_\lambda={1\over1+\lambda}+
 {\mu p e^{-\mu\lambda}\over1-p+p e^{-\mu\lambda}},
\qquad p={\mu-1\over\mu}.
\]
Thus \(m_\lambda=\mu-\mu\lambda+O(\lambda^2)\), and inversion gives
\[
\lambda=-{x\over\mu k}+O(x^2/k^2),\qquad
c'k^{-1/2}\le\lambda\le C'\sqrt{L/k}. \tag{S1}
\]
Wolfram 15.0 actually returned the series
\(m-m\lambda-\tfrac12m(m-3)\lambda^2+O(\lambda^3)\).

The Fourier proof record sharpens its global bounded-shift derivative
bound to
\[
 |f''(s+u)|\le C(\lambda^2+k^{-1})e^{\lambda|u|}f(s)
 \quad(u\in{\mathbb R}). \tag{S2}
\]
After untilting,
\(f''=M(\lambda)^ke^{\lambda(s+u)}
(\lambda^2g+2\lambda g'+g'')\). The global Fourier bounds
\(\|g\|_\infty\ll k^{-1/2}\), \(\|g'\|_\infty\ll k^{-1}\),
\(\|g''\|_\infty\ll k^{-3/2}\), and \(g(0)\gg k^{-1/2}\), give
(S2), since \(\lambda k^{-1/2}\le(\lambda^2+k^{-1})/2\).

Write
\[
D=K_\mu F_k(r)-F_k(s)
 =E\,F_k(s-(Z_\mu-\mu))-F_k(s).
\]
A centred second-order Taylor expansion, using (S2) and the finite
exponential moment of \(|Z_\mu-\mu|\) at \(\lambda<1/2\), gives
\[
D={\mu\over2}f'(s)+O((\lambda^2+k^{-1})f(s))
 ={\mu\over2}\lambda f(s)
   \{1+O(\lambda+(k\lambda)^{-1})\}. \tag{S3}
\]
The sign is positive: \(x<0\) makes \(\lambda>0\), hence \(f'(s)>0\).

For the full uncentred kernel derivatives,
\[
a=K'_\mu F_k(r)=pD,\qquad
b=K''_\mu F_k(r)=p f(s)+p(2p-1)D.
\]
Consequently
\[
a={p\mu\over2}\lambda f(s)
 [1+O(\lambda+(k\lambda)^{-1})],\qquad
b=p f(s)[1+O(\lambda)]. \tag{S4}
\]
The exact third-derivative formula from r027, (S2) for bounded translates,
and the exponentially decaying kernel tail give
\[
\sup_{c\in[\mu/2,3\mu/2]}|K_c'''F_k(r)|\le C f(s). \tag{S5}
\]
Hence r027's admissibility and third-remainder conditions hold uniformly
for sufficiently large \(k\).

With \(\delta=-a/b\), (S1)--(S4) give
\[
\delta={x\over2k}[1+O(\sqrt{L/k})]. \tag{S6}
\]
Taylor's theorem, retaining the leading quadratic term, yields
\[
K_\mu F_k(r)-K_{\mu+\delta}F_k(r)
 ={a^2\over2b}[1+O(\lambda)]
 ={p\over8}{x^2\over k^2}f_k(s)
   [1+O(\sqrt{L/k})]. \tag{S7}
\]
The last coefficient is
\[
{(p\mu\lambda f/2)^2\over2pf}
 ={p\mu^2\lambda^2f\over8}
 ={p x^2f\over8k^2}(1+O(\sqrt{L/k})).
\]
A Wolfram 15.0 exact-algebra check returned
\((m-1)x^2/(8mk^2)=p x^2/(8k^2)\).

This accepts only the uniform one-step lemma (improving r027's safe
\(p/12\) lower coefficient to asymptotic \(p/8\)). It does not establish a
controlled-path occupation bound, a cumulative \(\log n\) gain, or any
claim about \(\mathcal R_n\) or \(C_n\).

## Final audit of r029 §7

Date: 2026-11-06. **Accept.** The soft lift correctly uses
\(\omega_q(r)=u_q(q-r)\), which is decreasing in the remaining distance.
Thus \(P^N\omega_q(q)=E_Pu_q(S_N)\); using the increasing \(u_q\) itself
as a remaining-distance terminal would reverse the argument.

The capped \(u_q\) is Lipschitz with \(\|u'_q\|_\infty=2/N\) and has
derivative-measure total variation \(4/N\).  Therefore the three
derivative estimates (22) are valid in the distributional/BV sense:
convolve \(d\omega'_q\) with \(f_j\) and \(f'_j\), rather than assuming a
classical second derivative at its two corners.  The exact kernel
derivative identity and the all-state bounds make the soft telescope
uniform over the active policy states.  Its \(O(L/N)\) loss is
\(o(f_N(q_0)L)\), since \(f_N(q_0)L\asymp L^2/N\).

At the r019 *uncapped plateau* root \(q_0\), shifting to
\(q'=q_0+(1-6\varepsilon)gL\) loses
\((1-6\varepsilon)g f_N(q_0)L(1+o(1))\) in the baseline average.
The CDF corrector returns at least
\((1-4\varepsilon)g f_N(q')L\), while cap loss is exponentially small.
This leaves the strict positive margin displayed in (25).  The no-atom
statement is valid for the adaptive sum because every step retains an
independent continuous exponential component, so the plateau/CDF
threshold convention is immaterial.

Finally the cost inequality has the correct direction:
\(h_{q'}\le H_w(\cdot-q')\) implies
\(J(h_{q'})\le e^{-q'}K_N(w_N)\).  It yields a bounded admissible
profile for \(\mathcal R_N\).  The report explicitly and correctly does
not infer a new lower bound for \(C_N\).

## Final audit: r028 §§5--7

Date: 2026-11-06. **Accept; no mathematical blocker in the Fourier and
one-step portion.** The formal text consistently takes \(k\) to be the
number of *future* \(Z_\mu\) increments, \(F_k\) their lower CDF,
\(s=k\mu+x\), and evaluates the uncentred kernel at \(r=s+\mu\). Thus
\[
K_\mu F_k(r)-F_k(s)
=E[F_k(s-(Z_\mu-\mu))-F_k(s)]
\]
is exactly the centred Taylor difference used in §6; there is no
off-by-one or translation sign error.

In §5 the Fourier local error after multiplying by \(|u|^j\) is
\(O(k^{-(j+2)/2})\): it gives respectively \(k^{-1},k^{-3/2},k^{-2}\)
for \(j=0,1,2\), and the first Gaussian derivative vanishes by parity.
Together with the global tilted characteristic-function estimate this
supports (10), including the \(g'(0)=O(k^{-3/2})\) improvement. Untilting
then gives precisely the sharp \(\lambda^2+k^{-1}\) form of (11), rather
than an absolute Berry--Esseen error. The exponentially weighted Taylor
remainder is uniform because \(\lambda=O(\sqrt{L/k})<1/2\) and the positive
tail of the action kernel is exponential.

For §6, the exact full-kernel identities give
\(a=pD\) and \(b=pf(s)+p(2p-1)D\). The global shifted bounds control both
the kernel difference and the two density terms in \(K_c'''F_k(r)\);
integrating the exponential shift factor against the uniform exponential
action tail proves \(M\le Cf(s)\). The cubic parameter remainder relative
to \(a^2/(2b)\) is \(O(M|a|/b^2)=O(\lambda)\), so (14)--(15) and the
\(p/8\) coefficient follow exactly as claimed.

The MATLAB grid bound for \(M\) is correctly labelled a diagnostic rather
than a certification, and no numerical computation is used to establish
the lemma. §6 explicitly leaves the controlled occupation problem open;
therefore it makes no unjustified cumulative logarithmic-gain or
\(\mathcal R_n,C_n\) assertion. Cosmetic suggestion only: render the few
literal sqrt(L/k) terms as \(\sqrt{L/k}\) for typesetting clarity.

## r029 fixed-baseline Poisson corrector: audit and proof skeleton

Date: 2026-11-06.  This is a genuinely viable way to remove the
controlled-occupation circularity.  No sign obstruction was found, but a
formal proof must fix the following exact time convention.

Let \(P=K_{c_0}\), \(m=m(c_0)\), \(V_t=P^tH=F_t\), where
\(H=1_{[0,\infty)}\).  For a decision immediately preceding \(t-1\)
baseline future increments, set
\[
G_t(r)=PV_{t-1}(r)-K_{\alpha_t(r)}V_{t-1}(r).
\]
The correct source is indexed by the *output* level:
\[
g_t(r)=\kappa\,{(r-tm)^2\over t^2}\,
 \psi\!\left({r-tm\over\sqrt t}\right)f_t(r), \tag{P1}
\]
with \(g_t=0\) outside the chosen time window.  Here \(\psi\) is a smooth
two-sided fixed-band cutoff (eventually \(\theta(2-\theta)\)).  Define
\[
U_0=0,\quad U_t=PU_{t-1}+g_t,\quad W_t=V_t-U_t. \tag{P2}
\]
Then, exactly,
\[
K_{\alpha_t}W_{t-1}
=W_t+g_t-G_t-(K_{\alpha_t}-P)U_{t-1}. \tag{P3}
\]
Thus \(G_t\ge g_t+|(K_{\alpha_t}-P)U_{t-1}|\) proves the global
supersolution inequality and, by iteration,
\[
\Pr_\alpha(S_N\le q)\le F_N(q)-U_N(q). \tag{P4}
\]
There is no requirement that \(W_t\) itself be nonnegative.

This indexing is essential.  If a draft calls \(t\) the number of
remaining increments, then the source in (P2) is \(g_{t+1}\), not \(g_t\).
The \(r=s+\mu\) notation from r028 must similarly be translated before
using it here.

The choice (P1) has an important exact terminal identity:
\[
 {P^{N-t}g_t(q)\over f_N(q)}
 =\kappa E\!\left[
 {X_t^2\over t^2}\psi(X_t/\sqrt t)\ \middle|\ S_N=q\right],
 \qquad X_t=S_t-tm. \tag{P5}
\]
It follows directly by convolving \(f_{N-t}(q-y)f_t(y)\).  Therefore the
terminal lower bound is a standard fixed-baseline bridge calculation, not
a claim about the controlled path.  A source written with the density
\(f_{t-1}\) instead can work only after a uniform bounded-shift density
comparison; (P1) avoids that unnecessary mismatch.

### Corrector derivatives

For \(r=tm+O(\sqrt t)\), \(t_0\le t\), and \(0\le\ell\le3\), the desired
bound
\[
 |U_t^{(\ell)}(r+u)|
 \le C L t^{-\ell/2}e^{C|u|/\sqrt t}f_t(r) \tag{P6}
\]
has a direct proof.  Expand
\(U_t=\sum_{j\le t}P^{t-j}g_j\).
For \(j\le t/2\), take the \(\ell\) derivatives on \(f_{t-j}\).  Tilt all
three densities at the common parameter giving tilted mean \(r/t\).  The
global tilted Fourier derivative bound for the prefix is
\(O((t-j)^{-(\ell+1)/2})\); after division by the central tilted
\(f_t(r)\asymp t^{-1/2}\), this is
\(O(t^{-\ell/2}e^{C|u|/\sqrt t})\).  Since
\(\int g_j\le C/j\), summation contributes \(L\).

For \(j>t/2\), smoothness of the cutoff and the fixed standardized band
give
\[
 |g_j^{(\ell)}(v)|\le
 Cj^{-1-\ell/2} f_j(v)\,1_{\rm expanded\ band}.
\]
Put derivatives on \(g_j\), use
\(f_{t-j}*f_j=f_t\), and sum \(O(t)\) terms.  This gives the same
\(t^{-\ell/2}\) scale without an extra logarithm.  Untilting at the common
parameter supplies the displayed exponential factor in \(u\).  This
establishes the advertised derivative mechanism, provided the existing
Fourier estimates are recorded uniformly for small tilts of either sign.

### Parameter error

The exact parameter derivative can be Taylor-expanded around the mean of
one increment.  At the distinguished \(\mu\),
\[
m'(\mu)=0
\]
(independently checked in Wolfram 15.0 from
\(m(c)=1+c/(1+e^c)\) and \(e^\mu=1/(\mu-1)\)).  Hence for
\(c_0=\mu+O(\sqrt{L/N})\), uniformly over the small action interval,
\[
 |K'_cU|\le C\{ |c_0-\mu|\,|U'|+|U''|+|U'''|\},\qquad
 |K''_cU|\le C(|U'|+|U''|+|U'''|). \tag{P7}
\]
With \(|\delta|\le C\theta/\sqrt t\), (P6) gives
\[
 |(K_{c_0+\delta}-P)U_{t-1}|
\le C\theta {f_t(r)\over t}
\left\{|c_0-\mu|L+{L\over\sqrt t}\right\}. \tag{P8}
\]
This is \(o(\theta f_t(r)/t)\) for
\(t\ge L^4\), since \(|c_0-\mu|=O(\sqrt{L/N})\).  The same conclusion
holds with \(|c_0-m|\), which is comparable at this scale.

On a two-sided band bounded away from zero, the pointwise greedy
calculation around \(c_0\) gives, after reserving a fixed small amount,
\[
G_t(r)\ge (p/8-\varepsilon)\theta(2-\theta)
 {x^2\over t^2}f_t(r),\qquad x=r-tm. \tag{P9}
\]
The mean bias is negligible on \(t\le N/(C_1L)\) once \(C_1\) is fixed
large.  Equations (P8)--(P9), with \(\kappa<p/8-\varepsilon\), establish
the comparison condition in (P3).

Finally, for \(q=Nm-d_N+O(L)\), the fixed-baseline bridge in (P5) has,
for \(t\le N/(C_1L)\), standardized mean \(O(C_1^{-1/2})\) and variance
\(\mu(1-t/N)+o(1)\).  A tilted local limit theorem therefore gives
\[
U_N(q)\ge(\mu p/8-o(1))f_N(q)\sum_{L^K\le t\le N/(C_1L)}{1\over t}
=(\mu p/8-o(1))f_N(q)L, \tag{P10}
\]
after taking the cutoff to cover both signs, then its inner radius to zero
and outer radius to infinity.  This establishes the *corrector/CDF*
mechanism conditional on writing the uniform small-tilt Fourier lemma and
the near-\(\mu\) version of (P9).  It still does not by itself prove the
proposed \(\mathcal R_N\) statement: the soft-terminal/profile lift and
its Orlicz accounting remain separate required steps.

## r029 soft-terminal/profile lift: closed conditional on the CDF corrector

Date: 2026-11-06.  **Accept.** With the remaining-state convention
\(r=q-\) (past sum), the correct soft terminal is decreasing:
\[
 w_q(r)=u_q(q-r),\qquad
 u_q(s)=\min\!\left\{(e^{(s-q-w)/N}-1)_+,1\right\}. \tag{Sft1}
\]
Then \(P^Nw_q(q)=E_Pu_q(S_N)\), as required.  Using \(P^ju_q\) instead
would have the wrong orientation.

Set
\[
h_q(s)=1_{\{s\ge q\}}+u_q(s).
\]
It is nondecreasing, \(0\le h_q\le2\), with a unit plateau of width \(w\).
The actual Wolfram 15 check gives
\[
\|u'_q\|_\infty=2/N,\qquad {\rm TV}(du'_q)=4/N; \tag{Sft2}
\]
the three contributions to the latter are the initial \(1/N\) corner, the
rise from \(1/N\) to \(2/N\), and the \(2/N\) cap corner.

Let \(A_j=P^jw_q\), \(P=K_{c_0}\).  Uniformly for \(c_0\) in its
\(O(\sqrt{L/N})\) window about \(\mu\), the density Fourier bounds
\(\|f_j\|_\infty\ll j^{-1/2}\) and \(\|f'_j\|_\infty\ll j^{-1}\), together
with (Sft2), give, for \(j\ge1\),
\[
\|A'_j\|_\infty\le2/N,\qquad
\|A''_j\|_\infty\le C/(N\sqrt j),\qquad
\|A'''_j\|_\infty\le C/(Nj). \tag{Sft3}
\]
Indeed, \(A''_j\) is \(f_j\) convolved with the signed measure \(dw'_q\),
and \(A'''_j\) is \(f'_j\) convolved with that same measure.  This handles
the plateau and cap corners without assuming a nonexistent classical
second derivative.

Write \(m_0=m(c_0)\), \(b_0=c_0-m_0\).  The exact all-\(c\) identity
\[
K'_c\phi(r)=p_c[K_c\phi(r)-\phi(r-c)] \tag{Sft4}
\]
and a centred Taylor expansion about \(m_0\) show
\[
\|K'_{c_0}A_{t-1}\|_\infty
\le C\{|b_0|/N+1/(N\sqrt t)\}. \tag{Sft5}
\]
The exact second and third parameter derivatives, together with (Sft3),
are \(O(N^{-1})\).  Hence an active action
\(|\delta_t|\le C\theta/\sqrt t\) changes the soft potential by at most
\[
C\left\{{|b_0|\over N\sqrt t}+{1\over Nt}\right\}. \tag{Sft6}
\]
For inactive states the action is \(c_0\), so the change is exactly zero.
The standard nonhomogeneous Markov-operator telescope (each kernel is a
positive sup-norm contraction) yields
\[
|E_\pi u_q(S_N)-E_Pu_q(S_N)|
\le C\sum_{t_0\le t\le N/(C_1L)}
\left\{{|b_0|\over N\sqrt t}+{1\over Nt}\right\}
=O(L/N). \tag{Sft7}
\]
Here \(|b_0|=O(\sqrt{L/N})\); its summed term is only \(O(N^{-1})\).

Let \(q_0\) be the r019 *uncapped plateau* root for \(c_0\), not the soft
ramp root.  Its results give \(f_N(q_0)\asymp L/N\), density flatness on
\(O(L)\) shifts, and
\[
 B_{\rm unc}(q_0+\eta L)
=1-\eta L f_N(q_0)(1+o(1)). \tag{Sft8}
\]
For the capped profile the baseline cap loss is exponentially small:
the cap only changes the unbounded tail where
\(S_N>q+w+N\log2\), a linear positive deviation.  A fixed positive
exponential-moment Chernoff bound for the baseline kernel gives
\[
B_{\rm cap}(q)=B_{\rm unc}(q)+O(e^{-cN}) \tag{Sft9}
\]
uniformly on this \(O(L)\) root window.  No adaptive tail estimate is
needed: the policy comparison is applied directly to the bounded \(u_q\).

At \(q'=q_0+\eta L\), the CDF corrector supplies
\[
F_P(q')-F_\pi(q')\ge(g-\varepsilon)f_N(q')L,\qquad
g={\mu-1\over8}. \tag{Sft10}
\]
Combining the CDF plateau gain, (Sft7), and (Sft8)--(Sft9), and choosing
\(\eta=g-2\varepsilon\), gives \(E_\pi h_{q'}(S_N)>1\) for all large
\(N\).  The inequality is strict because \(f_N(q')L\asymp L^2/N\),
whereas the soft loss is \(O(L/N)\).

Finally \(J(h_{q'})\le J(h^{\rm unc}_{q'})=e^{-q'}K_N(w)\), with
\(\log K_N(w)=o(1)\).  Thus this bounded, admissible profile gives
\[
\log R_N\ge q'-\log K_N(w)
=\log Q_N+(g-2\varepsilon+o(1))L. \tag{Sft11}
\]
The r022 equality for \(\log Q_N\) and the r019 plateau lower profile
have the same logarithmic coefficient, justifying the last replacement.
Letting \(\varepsilon\downarrow0\) yields the proposed
\(\liminf(\log R_N-\log Q_N)/\log N\ge(\mu-1)/8\), conditional only on
the r029 CDF-corrector theorem.  It is a statement about the rearrangement
feedback model \(R_N\), not a direct \(C_N\) lower bound.
