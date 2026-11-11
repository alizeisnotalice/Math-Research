# A Gaussian comparison through every sublinear number of free coordinates

Fix $B\in\mathbb R$ and an integer sequence $0\le k=k_N=o(N)$, $n=N-k$. Use $L=\log N$, $\rho=\sqrt{L/N}$, $\mu=1+e^{-\mu}$, $\sigma=\sqrt\mu$, $\bar\Phi(z_N)=\sigma z_N/\sqrt N=\tau_N$, $d_N=\sigma\sqrt N z_N$, $q_{0,N}=\mu N-d_N$, $c_Q=(9-\mu)/24$, $p_N=1-\tau_N-BL/N$. Here $\Phi$ is the standard Gaussian cdf. Let $E\subset\mathbb R^k\times\mathbb R^n$ have positive finite volume. For almost every prefix position $z$, require that the true suffix section $E_z$ be an unconditional downset with convex logarithmic image as in r094, or its translate, invertible diagonal scaling, or the empty set, up to null modification. Centers, shapes and scales can vary freely with $z$. Then
\[
 \log\sup_E \frac{|\{M_{S,N}1_E>p_N\}|}{|E|}=q_{0,N}+c_QL+o_B(L).                    \tag{1}
\]
The error is uniform across inputs along each stipulated sequence $k=o(N)$. There is no claim of an error constant uniform over all possible rates $k/N\to0$, or a conclusion for a fixed positive fraction of free coordinates. The input is a true indicator with its actual volume; the prefix feedback used below enters only as an upper comparison.

## 1. Suffix CDF with an $O(1)$ error in Gaussian rank throughout a relative moderate-deviation window

For $y\in[\frac12\sqrt L,2\sqrt L]$ put $q=\mu n-\sigma\sqrt n\,y$. Let $S=\sum_{i=1}^n\lambda_i Z_{c_i}$ with independent $Z_c=Y+cJ$, where $Y$ is unit exponential and $J$ independent Bernoulli of probability $(1+e^c)^{-1}$; $c=\infty$ denotes no shift. Allow every deterministic action vector and arbitrary $\lambda_i\ge0$, $\sum\lambda_i=n$. Uniformly in this whole window,
\[
 \Pr(S\le q)\ge \bar\Phi\!\left(y+\frac{c_Qy^2+C}{\sigma\sqrt n}\right)              \tag{2}
\]
with an absolute finite $C$. The relative comparison matters: an absolute $O(1/n)$ cdf error alone would not give a bounded rank error everywhere in this window.

**Matched iid reference.** Use r094 §6's $F(s,c)=\log\mathbb Ee^{-sZ_c}$, $\ell(s)=\inf_c F(s,c)$ and minimizing $c(s)$. At zero $c(0)=\mu$, $m(\mu)=v(\mu)=\mu$, $m'(\mu)=0$, $m''(\mu)=-(\mu-1)/\mu$, $v'(\mu)=\mu-1$, $\kappa_3(\mu)=\mu(3-\mu)$. The stationarity equation to first order gives $c(s)=\mu-\mu s/2+O(s^2)$. In the expansion of $F(s,\mu+h)$ the $h$-dependent cubic terms are $s[(\mu-1)/\mu]h^2/2+s^2(\mu-1)h/2$. At the optimum these contribute $-\mu(\mu-1)s^3/8$. The ordinary third cumulant contributes $-\mu(3-\mu)s^3/6$; the respective coefficients sum to $-\mu c_Q$. Thus
\[
 \ell(s)=-\mu s+\mu s^2/2-\mu c_Q s^3+O(s^4).                                        \tag{3}
\]
Take the unique small saddle $-n\ell'(s_*)=q$, let $s_0=y/(\sigma\sqrt n)$, $c_*=c(s_*)$. Then $s_*=s_0+O(s_0^2)$ and the iid action $c_*$ has this same exact saddle. Taylor at $s_0$ in the objective, including the stationary quadratic minimum, gives
\[
 I_*=n\ell(s_*)+s_*q=-y^2/2-c_Q y^3/(\sigma\sqrt n)+O(y^4/n),\qquad
 b_*=s_*\sqrt{nF_{ss}(s_*,c_*)}=y[1+O(\rho)].                                      \tag{4}
\]
The exact tilted Mills factor is $M(b)=\int_{-\infty}^0e^{bt}\phi(t)dt=e^{b^2/2}\bar\Phi(b)$. Its logarithmic derivative is $O(1/b)$ on $b\ge1$. For the iid reference $F_*(q)$, r094's tilted absolute standardized density error $O(n^{-1/2})$ integrates against $e^{b_*t}1_{t\le0}$ to a *relative* $O(n^{-1/2})$ error on $M(b_*)$. Equations (3)–(4) therefore imply
\[
 \log F_*(q)=\log\bar\Phi(y)-c_Q y^3/(\sigma\sqrt n)
                  +O(\rho+y^4/n+n^{-1/2}).                                             \tag{5}
\]
In particular $\log F_*(q)=-y^2/2-\log y+O(1)$ uniformly. This centered-at-$y$ computation does not identify the Gaussian tail at $y$ with the tail at $z_N$; the two may differ by an unbounded factor.

**Why the arbitrary-weight comparison still applies here.** It suffices to consider vectors with $\Pr(S\le q)\le2F_*(q)$, all other vectors already meeting a stronger lower comparison. In r094 §6.2–6.4, keep the true saddle $s$ and its upper action and Jensen deficits $A,P\ge0$ in precisely the notation there. For clarity $P=n\log(1+s)-\sum\log(1+s\lambda_i)$ and $A=\sum(F(s\lambda_i,c_i)-\ell(s\lambda_i))$. The uniform fourth-moment estimate excludes $q\ge ES$ in this competitive range, so $s>0$ exists. The r094 global curvature, saddle minimum at most zero, and $q=\mu n-O(\sqrt{nL})$ first give $s=O(\rho)$. The uniform lower prefactor $c/(1+s\sqrt{V_s})$, $s^2V_s\le C(ns^2+P)$, and strongly curved baseline give
\[
 I_*+c n(s-s_*)^2+A+P/2-\log[C(1+\sqrt{ns^2+P})]\le\log(2F_*(q)).
\]
Use (4)–(5) on **both sides**, so that the $-y^2/2$ terms cancel, and the $-\log y$ offsets the $O(\sqrt L)$ prefactor. Thus the three deficits $n(s-s_*)^2,A,P$ are $O(1)$, with constants uniform as $y/\sqrt L\in[1/2,2]$. Also $s\asymp s_*\asymp\rho$. These are exactly the inputs to r094 §6.4: writing $h_n=ns^3$, $b=s\sqrt{V_s}$, its derivative estimates give $|b^2-b_*^2|\le C(P+A+\sqrt{Ah_n})$, and the absolute standardized tilted density error is $O(n^{-1/2}+P^{3/2}/L^{3/2})$. Comparing the exact Mills factors and retaining the positive rate deficits, r094 (15) gives
\[
 \log[\Pr(S\le q)/F_*(q)]\ge A+P/2-C[(P+A+\sqrt{Ah_n})/L+n^{-1/2}+P^{3/2}/L^{3/2}]
                           \ge-C'/\sqrt n .                                              \tag{6}
\]
The same constants work here: the scale used in these inequalities is $s^2n\asymp y^2\asymp L$, $h_n=O(L^{3/2}/\sqrt n)$, and no centering at the special implicit root enters after the three deficits have been bounded. For noncompetitors (6)'s last inequality is automatic. Finally the Gaussian hazard at $u\ge1$ is $u+O(1/u)$. Shifting $y$ by $c_Qy^2/(\sigma\sqrt n)$ changes the log Gaussian tail by $-c_Qy^3/(\sigma\sqrt n)+O(y/\sqrt n+y^4/n)$. In this window $y^4/n=o(y/\sqrt n)$ and $\rho\asymp y/\sqrt n$. An additional fixed $C/(\sigma\sqrt n)$ in the Gaussian argument reduces its log tail by at least a large fixed multiple of $y/\sqrt n$. Equations (5)–(6) then give (2). In particular the conversion to a bounded *unstandardized rank* is uniform, without an absolute Edgeworth remainder at tiny tail densities.

**From (2) to original sections.** Write $g(z,w)=M_{S,n}1_{E_z}(w)$, $g^*$ its decreasing rearrangement on the whole $N$-dimensional volume space. At each $x=q$ in the above window, (2) bounds the upper hit probability at rank $x$ for every deterministic normalized weighted sum in r094's supporting-plane inequality. The laws are continuous; it follows that their strict hit quantile at level $1-\bar\Phi(y+(c_Qy^2+C)/(\sigma\sqrt n))$ is at most $x$. Integrating the actual source-angle bound r094 (5) on each section gives its halo volume at that level at most $e^x|E_z|$. The same is true for translated, diagonally rescaled or empty sections. Integrating over $z$ and using the right-continuous generalized inverse yields
\[
      1-g^*(|E|e^x)\ge \bar\Phi(y+(c_Qy^2+C)/(\sigma\sqrt n)).                         \tag{7}
\]
Take a consistent rearrangement representative; any null exceptional ranks can alternatively be removed by continuity of the minorant. Measurability uses rational box suprema and Fubini, not a choice of shape parameters for the individual sections.

## 2. Propagating a lower Gaussian cdf through a predictable prefix

Fix $\epsilon>0$, $q=q_{0,N}+(c_Q+\epsilon)L$. If the true strict halo volume exceeds $|E|e^q$, r103 §3 (r011 rising-sun comparison iterated only on the prefix axes and applied to this whole-space $g^*$) gives a predictable policy with
\[
              \mathbb E g^*(|E|e^{q-S_k})>p_N.                                         \tag{8}
\]
Indeed one true physical prefix average of $g$ dominates the real box mean, and an interior normal volume rank and an approximately optimal policy preserve the strict sign. Since this is an upper argument we need only exclude (8) for all policies. The sharp lower will use a true source separately.

For the given policy set $m_i=\mathbb E_{i-1}Z_i$, $D_i=\mu-m_i\ge0$, $v_i=\operatorname{Var}_{i-1} Z_i$, $M_j=\sum_{i\le j}(Z_i-m_i)$. By r022 there are absolute constants with $0\le D_i\le C$ and $|v_i-\mu|\le C\sqrt{D_i}$. On the *same* filtered experiment put
\[
 e_i=\min(D_i,K\rho\sqrt{D_i}),\quad \widetilde Z_i=Z_i-m_i+\mu-e_i,\quad
 x_j=q-\sum_{i\le j}\widetilde Z_i=q-\mu j-M_j+\sum_{i\le j}e_i .                       \tag{9}
\]
Choose the absolute $K$ as below. Pathwise $\widetilde Z_i\ge Z_i$ and thus $x_k\le q-S_k$. Since $1-g^*(|E|e^x)$ is nondecreasing in $x$, a lower for its expectation at $x_k$ is sufficient. We do not rerun the feedback rule on the lifted path and do not assert that lifted increments belong to the original increment family.

Let $W=C_0(\sqrt{kL}+L)$. Uniformly over all policies the probability of an exit $\max_{j\le k}|M_j|>W$ is $O(N^{-8})$, taking fixed $C_0$ sufficiently large. This follows by the conditional bound $\mathbb E_{i-1} e^{t(Z_i-m_i)}\le e^{C_1 t^2}$ for either sign of small $t$, the two signed exponential supermartingales, and their maximal inequality. The bound is trivial for $k=0$. Prior to an exit, with $t=N-j$ the residual standardized depth obeys
\[
 \frac{\mu t-x_j}{\sigma\sqrt t}=z_N(1+o(1))
\]
uniformly in $j$ and all policies: the numerator is $d_N+O_{B,\epsilon}(L+W+Kk\rho)$ and $t/N=1+O(k/N)$. The deterministic relative error $\gamma_N$ here can be taken $C_{K,\epsilon,C_0}(\sqrt{k/N}+k/N+\sqrt{L/N})\to0$. In particular for the terminal dimension $n$ choose
\[
 D_0=c_Q z_N^2(1+\gamma_N)^2+C=c_Q L+o(L)
\]
with $C$ from (2). Equation (7) bounds $1-g^*(|E|e^{x_k})$ below by $H_n(x_k)$ on the no-exit event, where
\[
 A_t=D_0+C_2(t-n)\rho^2,\qquad H_t(x)=\Phi\!\left(\frac{x-\mu t-A_t}{\sigma\sqrt t}\right),\quad n\le t\le N.                    \tag{10}
\]
The linear drift will pay the one-step error. Since $k=o(N)$, $C_2 k\rho^2=o(L)$ for each fixed $C_2$.

At a pre-exit state $(t,x)=(N-j,x_j)$ with $j<k$, put $u=x-\mu$ and $H=H_{t-1}$. Its Gaussian density $H'(u)>0$ has log slope $O(\rho)$ in absolute value. More generally for every real $v$,
\[
 |H'''(u+v)|\le C_3 H'(u)e^{C_3\rho|v|}(\rho^2+1/N+v^2/N^2),
 \qquad |H''(u)|\le C_3\rho H'(u).                                                        \tag{11}
\]
Indeed the ratio $H'(u+v)/H'(u)$ equals the exponential of a linear term of slope $O(\rho)$ minus $v^2/[2\mu(t-1)]$. Differentiating the density twice gives (11). The slope constant may be taken absolute for every fixed $K,C_2$ and large enough $N$ since the state depth divided by $z_N$ tends to one and $A_t=O(L)$. This avoids any circular choice of the constants.

Conditional Taylor through order two, using (11) with $v=e_{j+1}-(Z_{j+1}-m_{j+1})$, gives
\[
 \mathbb E_j H(x-\widetilde Z_{j+1})
 =H(u)+e_{j+1}H'(u)+\tfrac12(v_{j+1}+e_{j+1}^2)H''(u)+O(\rho^2 H'(u)).
\]
The same expansion for $G\sim\mathcal N(\mu,\mu)$ uses centered variance $\mu$ and zero linear term. Remainders are uniform: the centered true increments have a uniform conditional exponential absolute moment at a fixed positive argument, $e_i=O_K(\rho)$, and the Gaussian has the same required integrability. Thus the $|v|^3$ and $|v|^5$ remainders in (11) integrate with the stated cost, even for one step which *leaves* the tube. The $e_i^2 H''$ term can be kept as nonnegative since $u$ is below the Gaussian mean for all these states; simply bounding it after fixing $K$ also gives a smaller order. Consequently the difference from $\mathbb E H(x-G)$ is at least
\[
 [e_{j+1}-C_4\rho\sqrt{D_{j+1}}-C_5\rho^2]H'(u)\ge -C_6\rho^2 H'(u).                  \tag{12}
\]
Choose fixed $K\ge C_4$. For $e_i=K\rho\sqrt{D_i}$ the linear loss is nonnegative. For $e_i=D_i$ use $D_i-C_4\rho\sqrt{D_i}\ge-C_4^2\rho^2/4$. The constants in the Taylor remainder can be bounded independently of fixed $K$ for large $N$ since $e_i\le1$ eventually.

Gaussian convolution is exact: $\mathbb EH(x-G)=H_t(x+C_2\rho^2)$. For $0\le h\le C_2\rho^2$ the density $H_t'(x+h)/H_{t-1}'(u)=1+o(1)$ uniformly at pre-exit states. This follows by taking logs in the explicit densities: changing the variance by $\mu$ changes the exponent by $O(L/N)$ and changing the center relative to $u$ by $O(C_2\rho^2)$ costs $O(C_2\rho^3)$. Thus fixed $C_2>2 C_6$ pays (12), showing
\[
              \mathbb E_j H_{t-1}(x-\widetilde Z_{j+1})\ge H_t(x).                      \tag{13}
\]
Stop this bounded submartingale at the first exit or at $k$. On paths with no exit use the terminal comparison (7), (10). On the other paths the true lower-test value is nonnegative and the stopped $H$ is at most one. Optional stopping at this bounded time therefore gives
\[
 \mathbb E[1-g^*(|E|e^{q-S_k})]\ge H_N(q)-O(N^{-8}).                                   \tag{14}
\]
Here $A_N=c_Q L+o(L)$ along the given $k$ sequence, hence $H_N(q)=\bar\Phi(z_N-(\epsilon L+o(L))/(\sigma\sqrt N))$. Its excess over $\tau_N$ is $(\epsilon+o(1))L^2/N$ by the root and Mills' bounds. This dominates both the fixed $B L/N$ and the exit error. Inequality (14) contradicts (8) and gives the upper of (1) after $\epsilon\downarrow0$. For the lower use the full $N$-coordinate equal-weight bounded hyperbolic cross of r093: every prefix section is a log-halfspace cut by a box, a full box, or empty; the entire physical halo achieves the fixed-vector coefficient. No prefix or Gaussian model is substituted for the lower example.

**Local algebra check.** Single-core Wolfram expanded the matched action's cubic cumulant expression with $h=-\mu s/2$ and returned zero for its difference from (3). It also returned zero for the square-completion residual used in (12). The probabilistic estimates and the positive-volume statement in (1) follow from the analytic argument.
