# r035: rigorous obstruction to the naive global quadratic minimum

This note checks the following specific candidate only:
\[
 W_t=O_t^-\wedge O_t^+\wedge I_t,\qquad
 I_t=v_t+\alpha U_t^{\rm quad}+\beta H_{t,2}f_t,
\]
where the quadratic correction is the all-time centred fixed-`mu` Poisson
potential, and `alpha>0` is fixed. It does not rule out a different nonlinear
or tilted global extension.

Choose a fixed sufficiently small `eta>0` and put
\(y=m_\mu(\eta)\in(\mu,2)\). Let `N=t` and `r_t=ty`. The exact fixed-kernel
saddlepoint relation is, for each fixed `z>=0`, 
\[
 {f_{t-1}(ty-z)\over f_t(ty)}
 \longrightarrow {e^{\eta z}\over M_\mu(\eta)}. \tag{1}
\]
For completeness, tilt both densities by `eta`. The exponential factors
leave `e^{eta z}/M_mu(eta)`, while the tilted arguments are `y-z` and `0`,
respectively; their centred local-density ratio tends to one. This is the
ordinary fixed-tilt local limit theorem, at a bounded displacement.

The exact bridge formula of r032 gives
\[
 {U_t^{\rm quad}(ty)\over f_t(ty)}
 ={p_\mu\over8}\left[t(y-\mu)^2+{C_t(ty)\over t}(H_t-1)\right]
 \sim {p_\mu\over8}t(y-\mu)^2, \tag{2}
\]
because `C_t/t` has a positive finite tilted limit. The analogous expression
at `(t-1,ty-z)` has the same leading factor. The reference terminal value is
only `O(f_t(ty))`: this follows directly by fixed-tilt integration of the
two terminal pieces (`R/f -> 1/(1-eta)` and the upper-tail linear piece is
also `O(f)`). The `beta H_{t,2}f_t` term is `o(t f_t)` as well. Therefore
\[
 {I_{t-1}(ty-z)\over I_t(ty)}
 \longrightarrow {e^{\eta z}\over M_\mu(\eta)}. \tag{3}
\]

The coarse left screen does not clip this region. Indeed, since `y<2`, its
Gamma term is exponentially small at both current and bounded-shifted
previous states, whereas
\[
 O_t^-(ty)\longrightarrow1+\mu-y>0,
 \qquad I_t(ty)=O(t f_t(ty))\longrightarrow0. \tag{4}
\]
The same holds with `t-1,ty-z` for each fixed z. Thus `W_t(ty)=I_t(ty)` and
`W_{t-1}(ty-z)=I_{t-1}(ty-z)` eventually on every fixed bounded z interval.
The right exponential screen also does not change this conclusion. Write
\(\Lambda_*(\theta)=\log M_*(\theta)\) and
\(\Lambda_\mu(\theta)=\log M_\mu(\theta)\). Then
`Lambda_* >= Lambda_mu`, while the inequality is strict at the chosen
`eta`, since its strict maximizer is `c_eta != mu`. The fixed-`mu` Cramer
rate at `y` has the unique maximizer `eta`. Hence the envelope rate
\[
 I_*(y)=\sup_{0<\theta<1}\{\theta y-\Lambda_*(\theta)\}
\]
is strictly smaller than
\(I_\mu(y)=\eta y-\Lambda_\mu(\eta)\): otherwise a maximizing sequence
for `I_*` must also maximize the strictly concave fixed-`mu` expression and
therefore converge to `eta`, contradicting the strict gap there. Consequently
`O_t^+(ty)=exp(-t I_*(y))` is exponentially *larger* than
`I_t(ty)=exp(-tI_mu(y)+o(t))`; the same statement holds at the bounded-shifted
previous states. Thus all three branches select `I` where used above.

Let `c_eta` be the strict maximizer of `M_c(eta)`. It differs from `mu` for
`eta>0`, and
\(M_{c_\eta}(\eta)>M_\mu(\eta)\). Since the integrand is nonnegative,
Fatou on each bounded interval, followed by monotone convergence in its
endpoint, yields from (3)
\[
\begin{aligned}
 \liminf_{t\to\infty}{K_{c_\eta}W_{t-1}(ty)\over I_t(ty)}
 &\ge {1\over M_\mu(\eta)}\int_0^\infty q_{c_\eta}(z)e^{\eta z}\,dz\\
 &= {M_{c_\eta}(\eta)\over M_\mu(\eta)}>1. \tag{5}
\end{aligned}
\]
No domination of the unbounded `z` tail is needed for this lower bound.
Consequently, for all large `t`,
\[
 {\cal T}W_{t-1}(ty)>W_t(ty),
\]
so the naive global `min` is not a MAX supersolution. This failure happens
in a fixed positive large-deviation region, not at the left coarse-cut
interface and not because of the `f/t^2` central remainder.

Numerical calibration only: WolframScript at `eta=0.05` found
`mu=1.278464542761...`, `c_eta=1.311493133190...`,
`y=1.345270200195...`, and
`M_{c_eta}(eta)/M_mu(eta)=1.00000587213355...>1`. This confirms the strict
sign but is not used in (1)--(5). Wolfram exited; MATLAB was not used.
