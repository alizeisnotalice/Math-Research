# r040: rate uplift that removes the r039 small-positive-tilt obstruction

Let
\[
 E(\eta)=I_\mu(y(\eta))-I_*(y(\eta)),\qquad
 F(\eta)=E(\eta)-D(\eta),\qquad
 D={p_\mu\over8}(y-\mu)^2\eta(1-\eta). \tag{1}
\]
The r036 exact coefficients give
\[
 F(\eta)=K_4\eta^4+O(\eta^5),\qquad
 K_4={\mu(\mu-1)(3+5\mu)\over384}>0. \tag{2}
\]
The proposed positive-side uplift adds `t F(eta)` to the logarithm of the
old nonlinear inner potential (with a smooth positive-side cutoff discussed
below). At the rate level this changes
\[
 J_{\rm old}=I_\mu-(1+\varepsilon)D
 \quad\hbox{to}\quad
 J_{\rm new}=I_* -\varepsilon D. \tag{3}
\]

## 1. Hamiltonian sign, not merely a value-rate comparison

At a smooth point, put `theta_*=I_*'(y)`. Since `H_{I_*}=0`, Taylor
expansion of the Hamiltonian yields
\[
\begin{aligned}
H_{J_{\rm new}}
&=-\varepsilon D+
 {\varepsilon^2\over2}\Lambda_*''(\theta_*)(D')^2
 +O(\varepsilon^3|D'|^3)\\
&=-\varepsilon D+O(\varepsilon^2\eta^4)<0 \tag{4}
\end{aligned}
\]
uniformly for `0<eta<=delta`, once `eps*delta` is sufficiently small. Here
`D~A_3 eta^3` and `D'=O(eta^2)` as a function of `y`. This is strictly
stronger than observing `J_new<I_*`.

Thus the fixed-small-eta r039 obstruction disappears: its former positive
rate excess `F` is exactly removed. The rate calculation does not prove a
finite-time supersolution.

## 2. Why the right coarse screen can now clip

For a fixed `eta>0`, the inner/right-screen ratio has exponential factor
\[
 \exp\{t[I_*(y)-J_{\rm new}(y)]\}
 =\exp\{t\varepsilon D(\eta)\}.
\]
Hence the inner is larger than `O^+` and the pointwise minimum selects the
global right screen at sufficiently large fixed eta. Including polynomial
prefactors, the crossover solves to first order
\[
 \eta_{\rm cross}\asymp
 \left({\log t\over\varepsilon t}\right)^{1/3}. \tag{5}
\]
This lies far outside the natural positive boundary
\(\eta_0=A\sqrt{\log t/t}\): their ratio tends to infinity for fixed
epsilon. Therefore the intended arrangement is internally consistent:
the exponentiated local branch remains active at the natural boundary, and
the right screen takes over before any fixed positive tilt can recreate the
r039 Fatou counterexample.

## 3. A finite-time uplift candidate and its actual analytic requirement

Let `eta_t(r)` solve `Lambda_mu'(eta)=r/t`. A candidate positive-side
extension is
\[
 \widetilde V_t(r)=V_t(r)
 \exp\{t\,\chi_t(\eta_t(r))F(\eta_t(r))\}, \tag{6}
\]
where `V_t` is the r036 potential and `chi_t` is a smooth cutoff, zero on the
negative side and one for `eta>=eta_0`. It should transition on a scale
comparable to but below `eta_0`; since `F=O(eta^4)`, the cutoff derivatives
are of the same order as the ordinary `F_y=O(eta^3)` terms and do not impose
a new leading source. At the natural boundary the uplift itself is
\(tF(\eta_0)=O((\log t)^2/t)=o(1)\), so it does not alter the critical
endpoint accounting.

For the deterministic phase `R_t(r)=tF(eta_t(r))`, its one-step Taylor
remainder is favorable:
\[
 \partial_{rr}R_t=O(\eta^2/t),\qquad
 \partial_{tr}R_t=O(\eta^2/t),\qquad
 \partial_{tt}R_t=O(\eta^2/t). \tag{7}
\]
Thus its bounded-shift defect is `O(eta^2/t)`, in particular no worse than
`O(eta/t)`. The original `tD` phase has the larger but still desired
`O(eta/t)` second-order defect.

The unresolved point is the *density prefactor*, not (7). A bare
saddlepoint estimate leaves a `1/t` term, too large relative to the
Hamiltonian margin at `eta_0`. One must normalize at the reference action
using the exact identity `P f_{t-1}=f_t`; then the action-window
(`c-mu=O(eta)`) remainder is expected to improve to `O(eta/t)`. A rigorous
two-parameter Edgeworth statement for the Gamma--Bernoulli density and the
auxiliary posterior densities is still required. Without that statement,
(6) is a well-motivated candidate, not a global Bellman supersolution.

## 4. Actual Wolfram diagnostic

`r040_rate_check.wl` was run under WolframScript (no MATLAB). For
`eps=.01,.005,.002,.001` and `eta=6 eps`, it found positive `F` and negative
`H_new`; for example at `eps=.01`,
`F=1.13537716399e-7` and `H_new=-1.00487981093e-7`.
All listed `J_new'` values lie strictly in `(0,1)`. These values support the
signs in (2)--(4) but are not an interval certificate.
