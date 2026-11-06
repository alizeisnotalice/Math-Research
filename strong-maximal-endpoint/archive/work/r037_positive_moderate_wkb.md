# r037: positive moderate tilts for the exponentiated potential

Set `alpha=1+eps`, with `eps>0` fixed, and work first with the fixed
reference `mu`.  Put
\[
 y(\eta)=\Lambda_\mu'(\eta),\qquad
 D(\eta)={p_\mu\over8}(y(\eta)-\mu)^2\eta(1-\eta),
 \qquad J_\alpha=I_\mu-\alpha D. \tag{1}
\]
This note records a precise WKB transport lemma which reduces the proposed
positive-moderate extension to uniform tilted local-limit estimates. It does
not claim the complete global Bellman construction.

## 1. Rate margin

From the exact r036 coefficients,
\[
 H_{J_\alpha}(y(\eta))
 =-(\alpha-1)A_3\eta^3+O(\eta^4),
 \qquad A_3={\mu(\mu-1)\over8}>0. \tag{2}
\]
Consequently there are absolute `c_0,C_0>0` such that, after choosing
\(\delta\le c_0\varepsilon\),
\[
 H_{J_\alpha}(y(\eta))\le-c_0\varepsilon\eta^3,
 \qquad0<\eta\le\delta. \tag{3}
\]
This is the margin that a finite-time calculation must preserve. A mere
comparison of `J_alpha` with an envelope rate is not enough; the Hamiltonian
`H_J=J-yJ'+Lambda_*(J')` is the relevant object.

## 2. Uniform positive-moderate transport lemma

Let `eta_0=A sqrt(log(t)/t)` and suppose
\(\eta_0\le\eta\le\delta\), `r=t y(eta)`. Let
\[
 U_t=\alpha J_t+s_t f_t,\qquad V_t=v_t\exp(U_t/v_t),
\]
where `J_t` is the exact quadratic Poisson potential and `s_t` is the
summable r036 source. The following is the right target statement.

**Lemma (uniform WKB transport; sufficient form).**  If the fixed-tilt
Gamma--Bernoulli local limit expansion and its first three differentiated
forms are used uniformly for `eta in [eta_0,delta]`, then for all sufficiently
large `t`, uniformly in every action `c`,
\[
 {K_cV_{t-1}(r)\over V_t(r)}
 \le \exp\{H_{J_\alpha}(y(\eta))
       +C(\eta/t+\eta^4)\}. \tag{4}
\]
Thus, after taking `delta` a sufficiently small multiple of `eps` and then
`A,t_0` sufficiently large,
\[
 \sup_cK_cV_{t-1}(r)<V_t(r) .\tag{5}
\]
Indeed, `eta^4<=delta eta^3`, while
\[
 {\eta/t\over\varepsilon\eta^3}
 ={1\over\varepsilon t\eta^2}
 \le {1\over\varepsilon A^2\log t}.
\]
This shows exactly why the lower boundary must be at least the
`sqrt(log t/t)` scale.

The lemma includes the source increment in `V_t`; harmless r036 summable
terms can be appended to the error in (4). For the matching reference `a_N`,
the displacement `a_N-mu` is in general only `O(eta_0)` near `t=N`, not
little-oh. One must retain the explicit matched-centre/affine drift terms of
r036; after that recentering the same lemma is the appropriate target.

## 3. Proof architecture for (4), and the non-small-w issue

Here `w=U/v` can be as large as `tD(eta)`, so a Taylor expansion in `w` is
invalid. The needed estimates are instead logarithmic. Uniform saddlepoint
algebra gives
\[
 -\partial_r\log V_t(ty)=J_\alpha'(y)+O(\eta/t+\eta^4), \tag{6}
\]
and, for bounded shifts,
\[
 \log{V_{t-1}(ty-z)\over V_t(ty)}
 =J_\alpha-yJ_\alpha'+zJ_\alpha'
   +O\big((\eta/t+\eta^4)(1+z^2)\big). \tag{7}
\]
The term `-alpha D_y` in `J_alpha'` is essential. It cancels the apparent
large derivative from `w`; replacing (6) by a bound on `|w|` loses the result.

For large shifts one needs a one-sided rather than two-sided oscillation
bound. The explicit leading term gives, for `z>=0`,
\[
 w_{t-1}(ty-z)\le w_t(ty)+C(\eta^2+t^{-1}). \tag{8}
\]
At leading order this is just monotonicity of `D(y)` on the positive side:
the time decrement can increase it only by `D_y y-D=O(eta^2)`, whereas the
left shift contributes `-zD_y<=0`.  The exact bridge/posterior derivatives
need to be used to make the remainder in (8) uniform.

Combining (8) with the reference log-concavity/Doob estimate supplies an
integrable envelope of the form
\[
 q_c(z){V_{t-1}(ty-z)\over V_t(ty)}
 \le C e^{-(1-\theta_0)z},\qquad \theta_0<1, \tag{9}
\]
uniformly in `c`.  Split at a slowly increasing fixed-moment cutoff, use (7)
on the first part, and use (9) on the remainder. Since
\(\int q_c(z)e^{J_\alpha'z}\,dz=M_c(J_\alpha')\), this produces (4), with
the action supremum becoming `Lambda_*` precisely as in the Hamiltonian.
The bound `q_c<=2e^-z` makes all polynomial weighted moments uniform in `c`.

## 4. What is still genuinely unproved

Equations (6)--(8) are not automatic consequences of r036, whose proof used
`w=o(1)` on the natural tube. To promote this note to a proof one must write
the uniform tilted expansions for the exact `C_t` posterior term through the
range `eta_0<=eta<=delta`, including its time difference and first derivative.
The most delicate requested cancellation is that the error in (7) be
`O(eta/t+eta^4)`, not a generic `O(1/t)`: a bare `O(1/t)` would overwhelm the
margin at the lower boundary. This is a concrete analytic sublemma, not an
occupation or splice issue.

No numerical grid is used here. The r036/r035 Wolfram computations establish
the exact small-tilt coefficients and diagnose the rate function, but they do
not certify (4).
