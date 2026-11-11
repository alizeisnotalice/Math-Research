# A Gaussian heat comparison for arbitrary prefix coordinates

Write $L=\log N$, $k=N-n$, $p_N=1-\tau_N-BL/N$ with $B$ fixed, and use the r100 section class. Suppose
\[
                         nL/N\longrightarrow\infty .                            \tag{1}
\]
In particular any fixed positive suffix proportion is allowed, and any $n\gg N/\log N$ is allowed. For this class the **whole physical** indicator score has sharp logarithm $q_{0,N}+c_QL+o_B(L)$. Condition (1) is strictly stronger than just $n=N^{1-o(1)}$. The obstruction to dropping it in this proof is identified below. Put $g_0=(\mu-1)/8$, $a_3=(3-\mu)/6$, $c_Q=a_3+g_0$, $\sigma^2=\mu$.

We apply the r011 global-volume rearrangement only along the prefix axes, as justified in r103. For the suffix maximal output $g(z,y)=M_{S,n}1_{E_z}(y)$ write
\[
 D_n(b)=1-g^*(|E|e^{\mu n+b}).
\]
If the full halo exceeds volume $|E|e^{\mu N+b}$, some predictable prefix policy has $E_\pi D_n(b+\mu k-S_k)<1-p_N$. Thus we seek an all-state *lower* comparison for this deficit. The stochastic comparison is only an upper device for the true halo. Put
\[
 f_t(b)=\frac{\phi(u)}{\sigma\sqrt t},\quad u=b/(\sigma\sqrt t),\quad
 P_s h(b)=E h(b+\sigma\sqrt s\,G),\quad G\sim N(0,1).
\]

## 1. Initial boundary at suffix time

Fix a large absolute $A$ (for example $24$). There exist fixed $U_0,C>0$ such that, for $U_0\le U=-b/(\sigma\sqrt n)\le A\sqrt L$,
\[
 D_n(b)\ge \Phi(-U)-(c_QU^2+C)f_n(b).                                          \tag{2}
\]
This includes a genuine interval of growing moderate-deviation levels, not just the one special $z_n$ in r094. Here is the change to the r094 proof needed for uniformity. For each log support plane consider any independent weighted sum $S=\sum\lambda_iZ_{c_i}$, $\sum\lambda_i=n$, at $q=\mu n-\sigma\sqrt n\,U$. If its actual CDF exceeds $2\Phi(-U)$ it already suffices for (2). Otherwise the r094 §6.2 crude lower prefactor and §6.3 localization apply with this test value. In the notation of those sections $I_*=\min_s[n\ell(s)+sq]=-U^2/2+O_A((1+U^3)/\sqrt n)$ and $s_*\sim U/(\sigma\sqrt n)$ uniformly. The true saddle $s$ exists since at and above the mean a fixed positive left probability follows from the uniform fourth-moment estimate; take $U_0$ large. From $\Psi(s)+sq\le0$ obtain $s\le C U/\sqrt n$, by the same baseline curvature argument. Comparison to $\log[2\Phi(-U)]=-U^2/2-\log U+O(1)$ then gives
\[
 c n(s-s_*)^2+A_{\rm def}+P_{\rm def}/2
    \le C+\tfrac12\log(1+P_{\rm def}/U^2).
\]
Here the two nonnegative defects are exactly those of r094. Constants can be fixed before raising $U_0$, since $U\ge1$ and $(1+U^3)/\sqrt n=o(1)$ uniformly. Thus all three defects are $O(1)$; for $U_0$ sufficiently large, $s\asymp U/\sqrt n$, $\sum(\lambda_i-1)^2\le Cn/U^2$, and at least a fixed fraction of the coordinates regularize the characteristic function.

For explicit error accounting, let $h=ns^3=O(U^3/\sqrt n)$, $R=P_{\rm def}+A_{\rm def}+\sqrt{A_{\rm def}h}$. The compact-action derivative estimates r094 (12) are unchanged because the compact bound on every $s\lambda_i$ came from $P_{\rm def}=O(1)$. The summed centered Taylor estimates in r094 (13) give $|b_s^2-b_*^2|\le CR$, $b_s,b_*\asymp U$ for the true and iid tilted standard-deviation parameters. The density argument (14) gives $O(n^{-1/2}+P_{\rm def}^{3/2}/U^3)$, and the exact two Gaussian Mills factors lose $O(R/U^2)$ relatively. Consequently the log-CDF ratio versus the actual common iid action $c_*=c(s_*)$ is bounded below by
\[
 A_{\rm def}+P_{\rm def}/2
 -C[R/U^2+n^{-1/2}+P_{\rm def}^{3/2}/U^3]\ge -C'n^{-1/2}.
\]
Indeed $U_0$ absorbs all terms linear in the bounded defects, and completing the square in the mixed term leaves $O(h/U^4)=O(1/(U\sqrt n))$. This extension of the refined comparison does not use a Gaussian additive CDF error against a deep tail.

To evaluate the matched iid CDF, the exact optimizing Laplace function in r094 has
\[
 \ell(s)=-\mu s+\tfrac12\mu s^2-\mu c_Q s^3+O(s^4).
\]
In fact $F(s,c)=-m(c)s+v(c)s^2/2-\kappa_3(c)s^3/6+O(s^4)$ locally; the optimized $c(s)=\mu-\mu s/2+O(s^2)$, $m''(\mu)=-(\mu-1)/\mu$, $v'(\mu)=\mu-1$, and $\kappa_3(\mu)=\mu(3-\mu)$. Substitution adds $-\mu g_0 s^3$ to the third cumulant term. Scalar minimization gives
\[
 I_*=-U^2/2-(c_Q/\sigma)U^3/\sqrt n+O(U^4/n),\qquad
 b_*=U+O(U^2/\sqrt n).
\]
The iid tilted density estimate in r094 (14) makes its exact CDF $e^{I_*}M(b_*)[1+O(n^{-1/2})]$, $M(x)=e^{x^2/2}\bar\Phi(x)$. Since $|(\log M)'|\le C/U$ here, division by $\Phi(-U)$ gives the log lower bound $-(c_Q/\sigma)U^3/\sqrt n-O(U^4/n+U/\sqrt n)$. Using $e^x\ge1+x$ and $U\Phi(-U)\le\phi(U)$ yields (2) for the scalar CDF with additive $C f_n$: the fourth-order term in rank units is at most $C U^3/\sqrt n=o(1)$. Fubini and the true angular support-plane integral r094 (5) give (2) for the physical suffix deficit. For each tested rank choose the threshold equal to the complement of the proved CDF lower value, so every sectional volume is at most its source mass times $\exp(\mu n+b)$. Generalized quantiles can be taken right-continuous.

For *all* real $b$, including the center and the right of the mean, an independent, coarse **physical** initial bound is
\[
 D_n(b)\ge\Phi\left(\frac{b-C_H\log n}{\sigma\sqrt{n+t_0}}\right)-2n^{-4},
 \qquad 0\le t_0=O(\log n).                                                       \tag{3}
\]
This follows from r091 (6), integrating over suffix fibers. A right side $\le0$ can simply be discarded; elsewhere the corresponding positive comparison threshold is in the allowed quantile range. On $|u|\le A\sqrt L$ the lower bound is $\Phi(u)-C_1 L f_n(b)-2n^{-4}$ by the mean value theorem: in rank units the shift is $O(L+|u|L/\sqrt n)$ and the intervening densities are comparable. An accurate unweighted iid CDF at each central or positive state was not assumed for arbitrary supporting normals.

## 2. A comparison at every state

Choose a smooth $0\le\chi\le1$ equal to $1$ on $u\ge -U_0$ and $0$ on $u\le-U_0-1$, with fixed bounded derivatives through four. For integers $n\le t\le N$ define
\[
 V_t(b)=\Phi(u)-[c_Q u^2+C_2+K\log(t/n)]f_t(b)
       -C_1L P_{t-n}[\chi(\cdot/(\sigma\sqrt n))f_n]\,(b).                       \tag{4}
\]
Here $C_2\ge C$ in (2) and $K>g_0+1$ are absolute constants. Bounds (2),(3) show $D_n\ge V_n-2n^{-4}$ on $|u|\le A\sqrt L$, since the two initial estimates overlap throughout the transition of $\chi$. On the far negative side $V_n\le\Phi(u)\le N^{-20}$; on the far positive side (3) gives $D_n\ge1-2n^{-4}-N^{-20}$, whereas $V_n\le1$. Thus there is an all-state initial comparison with error $2n^{-4}+N^{-20}$.

The last term in (4) solves the heat equation exactly. Its ratio to $C_1 L f_t(b)$ is
\[
 H_t(b)=E\chi\left(\sqrt{n/t}\,u+\sqrt{1-n/t}\,G\right),                        \tag{5}
\]
by the product of Gaussian densities, including $t=n$. Hence $|\partial_b^j H_t|\le C_j t^{-j/2}$ for $0\le j\le4$, with no singularity at the initial time. The exact one-step Gaussian residual is
\[
 \frac{P_1V_{t-1}(b)-V_t(b)}{f_t(b)}
    =c_Q(u^2-1)/t+K\log(t/(t-1)).
                                                                              \tag{6}
\]
Indeed the conditional second moment of the standardized Gaussian variable at time $t-1$ given terminal $b$ equals $((t-1)/t)u^2+1/t$. In particular a support correction costing order $L$ at the initially uncontrolled states has no heat residual at later times.

For each deterministic action $c$ put $Y=\mu-Z_c$, $\Delta=EY=\mu-m(c)\ge0$. For $|u|\le A\sqrt L$ we have, uniformly for $n<t\le N$,
\[
 \frac{\inf_c E_cV_{t-1}(b+Y)-P_1V_{t-1}(b)}{f_t(b)}
 =-g_0u^2/t-a_3(u^2-1)/t+O_A((1+L)^C/t^{3/2}).                                  \tag{7}
\]
Here a fixed polynomial degree $C$ suffices. Details follow to justify uniform minimization. Let $v=v(c)$, $\kappa=\kappa_3(c)$, and $G_1$ be centered normal of variance $\mu$. Taylor-expand both $V_{t-1}(b+Y)$ and $V_{t-1}(b+G_1)$ through degree three. After normalization by $V_{t-1}'(b)>0$ the action-dependent expression is
\[
 \Delta+(v-\mu+\Delta^2)\lambda'/2
    +(-\kappa+3\Delta v+\Delta^3)w'/6,qquad
 \lambda'=V_{t-1}''/V_{t-1}',\quad w'=V_{t-1}'''/V_{t-1}' .
\]
Write $P_L=(1+L)^{C'}$ with a fixed sufficiently large degree. Differentiating (4),(5) and the Gaussian gives
\[
 V'_{t-1}/f_t=1+O(P_L/\sqrt t),\quad
 \lambda'=-u/(\sigma\sqrt t)+O(P_L/t),\quad
 w'=(u^2-1)/(\mu t)+O(P_L/t^{3/2}).
\]
For the fourth-order Taylor remainder, at $b+v'$ the Gaussian density divided by $f_t(b)$ is bounded by $C\exp(C_A\sqrt{L/t}\,|v'|)$ (dropping the negative quadratic); derivatives add polynomial factors in $\sqrt L+|v'|/\sqrt t$. The resulting fourth derivative is bounded by $C_A f_t(b)P_L t^{-3/2}(1+|v'|^C)e^{|v'|/8}$ for sufficiently large $n$. Uniform exponential moments of $Y$ and of $G_1$ integrate the remainder. This also justifies comparing the infimum of the full expression to the infimum of its Taylor polynomial.

Compared to the action $c=\mu$, any minimizing/approximately minimizing action (within $O(P_L/t^{3/2})$) has $\Delta\le P_L/t$ after increasing the degree. Indeed $|v-\mu|\le C\sqrt\Delta$, $\kappa$ is uniformly bounded, and comparison first yields $\Delta\le C\sqrt\Delta|\lambda'|+O(|w'|+P_L/t^{3/2})$; the $\Delta^2|\lambda'|$ term is absorbed since $\Delta$ is bounded. There $|\kappa-\mu(3-\mu)|\le C\sqrt\Delta$ by r022. Thus replacing the last term by $-a_3(u^2-1)/t$, and replacing $\lambda'$ by $\lambda=-u/(\sigma\sqrt t)$ in the variance term, costs $O(P_L/t^{3/2})$ (increase the polynomial degree as needed). The same comparison applies to minimizers of the simpler $\Delta+(v-\mu)\lambda/2$. Its infimum is $-\mu g_0\lambda^2+O(|\lambda|^3)$ by r022 (15). This expansion holds for either sign of small $\lambda$: the unique, nondegenerate maximum of $m$ is interior and the minimizing smooth branch persists on both sides by the same implicit function argument. Finally $V'_{t-1}/f_t=1+O(P_L/\sqrt t)$ affects its already $O(P_L/t)$ infimum only at remainder scale. This proves (7).

Adding (6) and (7) leaves $-g_0/t+K\log(t/(t-1))+O_A(P_L/t^{3/2})$. The signed quadratic coefficient cancels because $c_Q=a_3+g_0$. Since $n\gg(\log N)^C$ for every fixed $C$ in (1), the combined residual is nonnegative throughout this window, uniformly for every action at large $N$.

Outside the window the following absolute comparison is sufficient. Formula (4) gives $V_t\le1$, $\inf_b V_t(b)\ge-N^{-10}$ uniformly for all large $N$: to see the second assertion use positivity on $|u|\le A\sqrt L$ (the subtracted terms have ratio $O_A(L(1+|u|)/\sqrt t)=o(1)$ to $\Phi(u)$) and Gaussian polynomial tails to the left and right. For $u<-A\sqrt L$, $V_t\le N^{-20}$; hence any conditional expected previous value exceeds $V_t-O(N^{-10})$. For $u>A\sqrt L$, on $Y\ge-b/2$ the preceding standardized state is $\ge A\sqrt L/2$, where (4) gives $V_{t-1}\ge1-O(N^{-20})$; the remaining probability is $O(e^{-cb})$ by the uniform exponential increment bound. For $A=24$ the Gaussian tail even at half the cutoff is much smaller than $N^{-20}$. Thus at **all** real states $E_c V_{t-1}(b+\mu-Z_c)\ge V_t(b)-O(N^{-10})$. Conditional iteration using the actual predictable prefix controls yields
\[
 E_\pi D_n(b+\mu k-S_k)\ge V_N(b)-O(n^{-4}+N^{-9}).                              \tag{8}
\]
No assumption on the path staying in a standardization window is required, including when the cumulative conditional drift is not negligible.

## 3. Terminal rank, new range and limit

Test $b=-d_N+(c_Q+\epsilon)L$ for fixed $\epsilon>0$. Writing $u_N=b/(\sigma\sqrt N)$ gives
\[
 \Phi(u_N)=\tau_N+f_N(b)(c_Q+\epsilon)L[1+o(1)],\qquad u_N^2=L+o(L).
\]
Here $z_N\sim\sqrt L$ and the shift of the standardized coordinate is $O(L/\sqrt N)$, negligible on the local density scale. Under (1), $\sqrt{n/N}\,u_N\to-\infty$, so the bounded Gaussian average (5) tends to zero. Also $\log(N/n)\le\log L+o(L)$, and the constant $C_2$ costs $o(L)$ in the rank. Therefore
\[
 V_N(b)\ge\tau_N+f_N(b)[\epsilon L-o(L)]>1-p_N+O(n^{-4}+N^{-9}),
\]
because $f_N(b)\sim L/N$, and $1-p_N=\tau_N+BL/N$. Equations (8) and the strict physical comparison exclude this rank. Let $\epsilon\downarrow0$. The full $N$-dimensional equal-weight cross from r093 lies in the section class and gives the matching lower.

If $nL/N$ stays bounded, the Gaussian average of the coarse central/right initial correction in (5) need not tend to zero. Even $n=N^{1-o(1)}$ can have $nL/N\to0$. At those dimensions the same one-step cancellation still holds whenever $n$ dominates fixed powers of $L$, but the initial term may cost order $L$ at the terminal rank. This is a boundary-data gap for the present proof, not a counterexample to a sharper physical theorem.

**Independent local check.** A Python standard-library calculation took (4) without its constant, log and initial heat terms, integrated the actual exponential/Bernoulli action at a single step by composite Simpson on $[0,40]$ with spacing $0.025$, then optimized the action scalar numerically. At $t=4096$ and $u=-2,0,2$, the residual times $t/f_t$ was approximately $-0.050099,-0.034717,-0.018710$; at $t=16384$ the values were $-0.042422,-0.034652,-0.026726$. The predicted common leading residual is $-g_0=-0.0348081\ldots$. This floating check validates signs and scales only; (6) is exact and (7) has the analytic uniform error stated above.
