# r039: strict intermediate-time obstruction to a fixed global exponentiation

This note concerns the particular global candidate obtained by taking the
minimum of the r034 coarse barriers with the r036 fixed-reference nonlinear
potential. It does not rule out a time-dependent right-tail extension which
changes the rate function.

Let `alpha=1+eps`, with `eps>0` sufficiently small, and write
\[
 J_\alpha(y(\eta))=I_\mu(y(\eta))-\alpha D(\eta),\qquad
 D=A_3\eta^3+B_D\eta^4+O(\eta^5),
\]
where `y=Lambda_mu'(eta)`.  The exact r036 coefficients give
\[
 I_\mu-I_*=A_3\eta^3+B_{\rm gap}\eta^4+O(\eta^5),
\quad K_4:=B_{\rm gap}-B_D
 ={\mu(\mu-1)(3+5\mu)\over384}>0. \tag{1}
\]
Thus the difference between the candidate and the true all-action rate is
\[
 \Delta:=J_\alpha-I_*
 =-\varepsilon A_3\eta^3+K_4\eta^4
 +O(\varepsilon\eta^4+\eta^5). \tag{2}
\]
This confirms that the fourth-order coefficient is exactly `K_4`, not merely
a numerical fit.

Choose a fixed `c>A_3/K_4` and set `eta=c eps`. Then explicitly
\[
 \Delta=\varepsilon^4c^3(cK_4-A_3)+O(\varepsilon^5)>0. \tag{3}
\]
For all `alpha in (1,1+eps_0)` with `eps=alpha-1` and a sufficiently small
absolute `eps_0`, this applies. The corresponding Hamiltonian has the
same leading sign. Indeed, if `theta_*=I_*'(y)`, then Taylor expansion of
\[
 H_{J_\alpha}=J_\alpha-yJ_\alpha'+\Lambda_*(J_\alpha')
\]
around `I_*` gives
\[
 H_{J_\alpha}=\Delta+{1\over2}\Lambda_*''(\theta_*)
 (\Delta')^2+o((\Delta')^2)
 =-\varepsilon A_3\eta^3+K_4\eta^4+o(\varepsilon^4)>0. \tag{4}
\]
Here `J_alpha'` stays in `(0,1)` for small epsilon.  In particular this is
the opposite sign to the small-tilt margin used in r037, because r037 chose
`eta<=delta eps` with `delta<A_3/K_4`, while this note takes `c` above the
threshold.

## Intermediate-time Fatou counterexample

Let
\[
 t_N=\lfloor\sqrt{\log N}\rfloor,\qquad r_N=t_Ny(c\varepsilon).
\]
Keep the terminal parameter at the global value `N` while comparing times
`t_N` and `t_N-1`.  The matched reference displacement obeys
\(t_N|a_N-\mu|\to0\).  For the actual r036 source
\(s_t=s_0+B\sum_{j>t_0}(\log j)^6/j^{3/2}\), `s_t` converges and
`s_t-s_{t-1}->0`.  Hence its density contribution is lower order than the
quadratic term `J_t\asymp t_N f_t` and cancels in the bounded-shift one-step
ratio (use `f/v -> eta(1-eta)`).  A merely bounded, arbitrarily oscillatory
sequence `s_t` would not justify this conclusion, since its one-step
difference can survive in (4).  The additive correction `Lambda t_N b_N^2/N`
is polynomially small in `N`, hence negligible relative to
`exp(-t_NJ_alpha(y))`.

For the nonlinear inner branch `I_{t,N}`, fixed-tilt local limits therefore
give, for every fixed `z>=0`, 
\[
 {I_{t_N-1,N}(r_N-z)\over I_{t_N,N}(r_N)}
 \longrightarrow
\exp\{J_\alpha-yJ_\alpha'+zJ_\alpha'\}. \tag{5}
\]
This ratio needs more than a logarithmic large-deviation principle. A
sufficient rigorous input is the following uniform fixed-tilt expansion on a
small compact positive-eta interval:
\[
 I_{t,N}(ty)=A(\eta)t^{\beta(\eta)}e^{-tJ_\alpha(y)}(1+o(1)), \tag{5a}
\]
with the same expansion after one bounded state shift and with `A,beta,J`
having one uniform `y` derivative. It implies (5) directly by expanding the
`t-1` expression at `y+(y-z)/(t-1)`.

Here (5a) is available in principle from the explicit Gamma--binomial
formula, not from a bare LDP. Under a fixed compact tilt, Fourier inversion
gives a uniform first-order Edgeworth expansion for the density and its
first derivative: the exponential component makes the law nonlattice, the
mgf is analytic for `eta<1`, and the noncentral-frequency integral decays
uniformly.  The exact identity
\[
 f_tC_t=tg_{t+2,t}+\mu^2t p_\mu(1-p_\mu)g_{t,t-2}(\cdot-\mu)
\]
then gives the required posterior expansion. More explicitly,
\[
 {v_t\over f_t}=h_0(\eta)+{h_1(\eta)\over t}+O(t^{-2}),\qquad
 {J_t\over f_t}=tq_0(\eta)+q_1(\eta)\log t+q_2(\eta)
 +O(\log t/t).
\]
Thus `U_t/v_t=tD+E log t+F+o(1)`. For the r036 source,
`s_t->s_infty` and `s_t-s_{t-1}->0`, so it only changes `F` and cancels in
the one-step ratio. This produces (5a), including the polynomial exponent
`beta(eta)`. The same Edgeworth estimates control the `R,Q,call` pieces of
`v_t/f_t` by their fixed-tilt tail integrals.

For clarity, the required prefactor calculation can be made fully explicit.
Let `nu` and `kappa_3` be the variance and third cumulant of one increment
under the `eta` tilt. Uniformly for `|z|<=C log t`,
\[
 {f_t(ty-z)\over f_t(ty)}=e^{\eta z}
 \left[1+{\kappa_3z/(2\nu^2)-z^2/(2\nu)\over t}
 +O\!\left({\operatorname{polylog}t\over t^2}\right)\right]. \tag{5b}
\]
Exponential-tail integration of (5b) against the three terminal pieces gives
\[
 {v_t\over f_t}=h_0+{h_1\over t}+o(t^{-1}),\quad
 h_0={1\over\eta(1-\eta)},\quad
 h_1={\kappa_3h_0'\over2\nu^2}-{h_0''\over2\nu}. \tag{5c}
\]
The exact auxiliary-density identity for `C_t` yields
`C_t/t=nu+O(t^-1)`.  Consequently, with
\(q_0=p_\mu(y-\mu)^2/8\) and `B=p_mu/8`,
\[
 {J_t\over f_t}=q_0t+B\nu\log t+B\nu(\gamma_E-1)+o(1). \tag{5d}
\]
Combining (5c)--(5d) gives the stronger version of (5a)
\[
 \log I_{t,N}(ty)=-tJ_\alpha(y)
 +\left[-\tfrac12+{\alpha B\nu\over h_0}\right]\log t
 +F(\eta)+o(1), \tag{5e}
\]
where `F` is smooth on the chosen compact eta interval (it includes the
`h_1` correction and the convergent source limit). This proves the
bounded-shift ratio (5), rather than merely its logarithmic rate.  The
Edgeworth remainder follows by expanding the analytic tilted characteristic
function through second order; Gamma smoothing supplies the uniform Cramer
tail. The `z>C log t` portion is exponentially integrable, but the Fatou
argument itself needs only fixed `z` and does not rely on this tail bound.
The left coarse barrier is asymptotic to one in this regime (`r_N=o(N)`),
whereas the inner value tends to zero. The right exponential barrier has
rate `I_*<J_alpha`, hence is exponentially larger than the inner branch.
The pointwise minimum thus selects the inner branch at the current state and
at every bounded-shifted preceding state.

Let `c_*` maximize `M_c(J_alpha')`.  Nonnegativity and Fatou on each bounded
`z` interval, then monotone convergence of those intervals, give
\[
 \liminf_{N\to\infty}
 {K_{c_*}W_{t_N-1,N}(r_N)\over W_{t_N,N}(r_N)}
\geq \exp\{H_{J_\alpha}(y)\}>1. \tag{6}
\]
No domination of the unbounded left-shift tail is used.  Consequently the
fixed-reference global-min candidate with a convergent density source fails
its one-step MAX supersolution inequality along this intermediate time scale,
despite a finite-time initialization and despite any late-time tube estimate.

## Actual numerical check

`r039_global_inflation_check.wl` was run with WolframScript and wrote
`/tmp/r039_global_inflation_check.json`. It found
\[
 A_3=.0445008805420\ldots,\quad K_4=.00870763814790\ldots,
 \quad A_3/K_4=5.1105569370\ldots.
\]
For `c=6`, both `Delta` and `H` are positive after division by `eps^4` for
`eps=.01,.005,.002,.001`; for `c=5` both are negative. These are diagnostics
only; strictness in (2)--(5) is analytic.

This does not contradict r037: the r037 moderate-tilt lemma must restrict
its upper constant to a value below `A_3/K_4` (and is only a local result).
