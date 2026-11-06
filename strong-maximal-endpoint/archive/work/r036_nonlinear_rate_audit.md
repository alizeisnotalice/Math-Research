# r036: nonlinear exponential correction versus the right-LD obstruction

Let `eta in (0,1)`, `y=Lambda_mu'(eta)`, and write
\[
 D(\eta)={p_\mu\over8}(y-\mu)^2\eta(1-\eta),\qquad
 J_\alpha(y)=I_\mu(y)-\alpha D(\eta). \tag{1}
\]
At the fixed positive large-deviation point `r=ty`, the nonlinear candidate
`v exp(alpha Jquad/v + beta H_{t,2} f/v)` has logarithmic rate `J_alpha`:
indeed `v/f -> [eta(1-eta)]^{-1}` and the exact quadratic bridge formula
gives `Jquad/f ~ (p_mu/8)t(y-mu)^2`.  Fixed `beta` only contributes lower
order terms to this rate.

## The correct one-step rate test

The relevant condition is stronger than merely comparing the value rate to
the envelope. If a smooth candidate has rate `J`, then at `r=ty` its
one-step MAX log-residual tends to
\[
 H_J(y)=J(y)-yJ'(y)+\Lambda_*(J'(y)),\tag{2}
\]
provided `0<J'(y)<1`. This follows by the same bounded-shift tilted local
limit as r035: the ratio at shift `z` is
`exp(J-yJ'+zJ')`, then the action supremum integrates to
`M_*(J')`. A necessary rate-level supersolution condition is `H_J<=0`.

There is a useful sufficient criterion. If `J` is convex and `J<=I_*`
pointwise, then at a differentiability point, with `theta=J'(y)`, convex
duality gives
\[
 \Lambda_*(\theta)\leq J^*(\theta)=\theta y-J(y),
\]
hence `H_J(y)<=0`. This is a rate-function statement only; it does not
control finite-time prefactors, the central tube, or the terminal condition.

## Small positive tilts

The following coefficients have an exact logistic-derivative derivation
(using only `exp(-mu)=mu-1`):
\[
 A={\mu(\mu-1)\over8},\quad
 B_\Lambda={\mu(\mu-1)(9-2\mu)\over48},\quad
 B_D=A(2-\mu).
\]
More precisely, envelope optimization gives
\(\Lambda_*-\Lambda_\mu=A\theta^3+B_\Lambda\theta^4+O(\theta^5)\).
Since the dual optimizer for `I_*` is
\(\theta_*=\eta-3A\eta^2/\mu+O(\eta^3)\), Legendre transformation gives
\[
 I_\mu-I_*=A\eta^3+B_{\rm gap}\eta^4+O(\eta^5),\qquad
 B_{\rm gap}=B_\Lambda-{9A^2\over2\mu}
 ={\mu(\mu-1)(99-43\mu)\over384}.
\]
Together with `D=A eta^3+B_D eta^4+O(eta^5)`, this yields the exact strict
difference
\[
 K:=B_{\rm gap}-B_D={\mu(\mu-1)(3+5\mu)\over384}>0.\tag{3}
\]
Thus the previous decimal expansions can be read rigorously as
\[
\begin{aligned}
 A&=.0445008805420\ldots,\\
 B_{\rm gap}&=.0408166013373\ldots,\quad
 B_D=.0321089631894\ldots,\quad K=.00870763814790\ldots .
\end{aligned}
\]
For reference, one efficient derivation differentiates the logistic law with
`D_c-p(1-p)D_p`, uses `m=1+cp`,
`v=1+c^2p(1-p)`, and the third cumulant
`k_3=2+c^3p(1-p)(1-2p)`, then inserts `c_1=mu/2` in the envelope condition.
Wolfram independently evaluated all exact expressions and their equality
under `exp(-mu)=mu-1`.
\]
Thus the linear correction with `alpha=1` is still short by
`K eta^4+O(eta^5)`. The direct one-step computation gives the
same leading obstruction
\[
 H_{J_1}=K\eta^4+O(\eta^5)>0. \tag{4}
\]
More generally, for fixed `alpha>1`,
\[
 H_{J_\alpha}=-(\alpha-1)A\eta^3+O(\eta^4),\tag{5}
\]
so a fixed excess over one repairs the r035 rate obstruction for all
sufficiently small *fixed* positive tilts. This is a concrete scoped
advantage of exponentiating the correction; it is not a full Bellman proof.

## The eta-to-one end is not the first rate obstruction

Put `eps=1-eta`. Direct expansion gives
\[
 y=\epsilon^{-1}+\mu/2+O(\epsilon),\quad
 I_\mu=\epsilon^{-1}+\log\epsilon+O(1),\quad
 D={p_\mu\over8}\epsilon^{-1}+O(1).
\]
For `0<alpha<8/p_mu`, set `kappa=1-alpha p_mu/8`. Then
\[
 J_\alpha=\kappa/\epsilon+\log\epsilon+O(1),\qquad
 J_\alpha'=\kappa+O(\epsilon),
\]
and hence
\[
 H_{J_\alpha}=\log\epsilon+O(1)\longrightarrow-\infty.\tag{6}
\]
Thus the `eta -> 1` end does not itself defeat a modest fixed `alpha>1` at
the rate level. It is instead a useful warning that a very large alpha can
force `J'` out of `(0,1)`.

## Numerical reconnaissance only

The actual scan in `r036_nonlinear_rate_scan.wl` found
`(I_mu-I_*)/D` approximately `1.00002, 1.01827, 1.04534, 1.00467,
0.55382, 0.08295` at eta `=.0001,.1,.4,.6,.9,.99`, respectively. Thus
`alpha=1` fails near zero and in a moderate interval, while the scan did not
find a value-rate obstruction for alpha about `1.05`. It is not an interval
certificate. The same script directly evaluates (2), and should be used to
test convexity of `J_alpha` before any global rate claim.

The terminal critical inflation is harmless at the desired endpoint:
there `Jquad=O(f_N log N)` while `v=1+o(1)`, so exponentiation changes the
value by `O(f_N log N)` at first order and its higher powers are smaller.
What remains open is a uniform finite-time all-action comparison, including
the transition from this LD rate extension to the critical tube.
