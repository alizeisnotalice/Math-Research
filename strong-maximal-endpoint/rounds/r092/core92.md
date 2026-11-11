# A one-indicator upper reduction at the critical density and its feedback obstruction

Write $L=\log N$, $\rho=\sqrt{L/N}$, $\mu=1+e^{-\mu}$, $\sigma=\sqrt\mu$, $\bar\Phi(z)=\sigma z/\sqrt N$, $d=\sigma\sqrt N z$, $\tau=d/N$, $q_0=\mu N-d$, $c_Q=(9-\mu)/24$, $g=(\mu-1)/8$, $c_R=c_Q+g$. Use the physical indicator constant $H_N(p)$ of r091 and the actual modular $J$ of r089–091. No Gaussian or feedback indicator constant is identified with $H_N$.

## 1. A single physical near-one indicator suffices for an upper bound

Fix $P>0$. There is a finite $B_P>0$ such that for all sufficiently large $N$ and every measurable $0\le F\le 1+2\sigma\rho$ of positive finite original cost, put $e=\min(F,1)$, $S_0=\int e$, $\beta_N=1-\tau-B_PL/N$. There is a measurable $U$ of volume $\le Q_N N^{-P}J(F)$ such that every successful rectangle $R$ through every point outside $U$ satisfies
\[
 m_R(F-1)_+\le\tau+B_PL/N,\qquad m_R e>\beta_N.                  \tag{1}
\]
The true general bound consequently is
\[
 |\{M_SF>1\}|\le Q_NN^{-P}J(F)+H_N(\beta_N)S_0,
 \qquad C_N\le(1+O(L^{-1/4}))[Q_NN^{-P}+H_N(\beta_N)].          \tag{2}
\]
The first error in (1) is $L/N$ on the upper side of the excess, a full logarithm smaller than the previous two-sided $O_P(L^2/N)$. A receiver may have many witnesses; they all obey (1) outside the same physical exceptional set. The second inequality of (2) applies r089's critical-cap reduction after the first physical inequality, using $S_0\le J(F)$. In particular a true estimate $H_N(\beta_N)\le Q_N N^{g-\delta+o(1)}$ for one fixed $P>0$ and $0<\delta\le g$ gives the same improvement for $C_N$ (with closure when $\delta=g$). No estimate at exactly $1-\tau$ has been presumed interchangeable with one at $\beta_N$ for the actual $H_N$.

**Proof.** Choose $s=q_0-(P+5)L$. On the r090 envelope use the increasing continuous transform $T(v)=(v-1)_+$, which removes the unit jump entirely. For $t\ge s$ its value satisfies
\[
 (W_s(t)-1)_+\le (t-s)/N+C\rho^2.
\]
This follows from the uniform quadratic inverse estimate up to the envelope cap, and the right side continues increasing beyond the cap. At $t<s$ the transform is zero. With r021's comparison $G'$ used in r090, $\mathbb E(s-G')_+=O_P(1)$ and $\mathbb EG'-s=d+O_P(L)$. Thus $\mathbb ET(W_s(G'))\le(d+O_P(L))/N+C\rho^2=\tau+O_P(L/N)$. This bound for the smooth excess avoids the $O_P(L^2/N)$ error in the occupation probability of the baseline jump. Take $B_P$ large enough for the $6N^{-4}$ comparison error as well. r090 (4) bounds the actual strict halo $U=\{M_S((F-1)_+)>\tau+B_PL/N\}$ by $e^sJ(F)\le Q_NN^{-P}J(F)$ for purified inputs. For any original gray input its purification has *identically the same entire excess function* and precisely the same $J$, so the same ordinary physical $U$ works directly. The identity $F=e+(F-1)_+$ and the strict mean on $R$ prove the second part of (1).

For every $0\le e\le1$ of finite mass, $|\{M_S e>p\}|\le H_N(p)\int e$. Indeed subdivide into fine regular cubes, in each keep a measurable indicator portion of measure $\int_I e$. The indicator corollary of r090 §1 preserves total mass, and on a fixed successful rectangle the integral converges by boundary-cell cancellation. Fatou and the *physical* indicator definition prove the estimate for each $0<p<1$. Applying it to (1) proves (2). The two terms pay the actual exceptional volume and the entire base mass, without recharging every upper height separately.

## 2. Physical price and incidence consequences

Let $K>0$ be fixed and $H=\{M_SF>1\}$ satisfy $|H|\ge Q_NN^{-K}J$. Take $P=K+3$. The r091 indicator quantile estimate (6) gives
\[
 \log H_N(\beta_N)\le q_0+C_H L+O_K(1).                            \tag{3}
\]
Indeed $1-\beta_N+2N^{-4}=\tau+O_K(L/N)$ and $\phi(z)/(\sigma\sqrt N)\sim L/N$. Gaussian inverse differentiation on this interval changes the *unstandardized Gaussian rank* by $O_K(1)$: to justify local uniformity, perturb $-z$ first by a fixed sufficiently large multiple of $1/\sqrt N$, where the density ratio is $1+o(1)$. The change of standard deviation in r021 contributes $O(L^{3/2}/\sqrt N)=o(1)$. Equation (3) is one-sided. It gives no continuity estimate of this order for the true $H_N$.

By the first bound of (2) and $\log Q_N=q_0+c_Q L+o(L)$,
\[
 1\le J/S_0\le N^{K+C_H-c_Q+o(1)}.                                 \tag{4}
\]
Here and below the source varies freely with $N$. Even the coarser r090 deficit and r091 (6) give some polynomial in (4); eliminating the jump identifies the dependence on the coarse r021 shift rather than an additional $K$-dependent quantile coefficient. In particular (4) improves the previous $\exp(d+O_K(L^2))$ bound on $J/S_0$. For all levels $\lambda_y=1+\tau+y/\sqrt N$ within the cap, $y\ge0$, the exact price bound $|E_y|\Phi_N(\lambda_y)\le J$ gives uniformly
\[
 |E_y|/S_0\le\exp(-d-y\sqrt N+O_K(L)),\qquad E_y=\{F>\lambda_y\}.  \tag{5}
\]
This includes noncompact gray sources; $S_0$ is finite. For the r091 near-top level $y_N=\sigma\sqrt{L-8\log L}$ and its original good receiver set $G$, let $w$ denote the load of any measurable successful-rectangle selection there. Its lower tail mass $\int_{E_{y_N}}w\ge cN^{-1/2}L^{7/2}(1-O(N^{-2}))|G|$ comes from the r091 one-test theorem, including the gray transfer. Using $\int we\le|G|$ gives
\[
 \log\frac{|E_{y_N}|^{-1}\int_{E_{y_N}}w}{S_0^{-1}\int we}
 \ge d+y_N\sqrt N-O_K(L),\qquad \liminf\frac{1}{d}\log\frac{|E_{y_N}|^{-1}\int_{E_{y_N}}w}{S_0^{-1}\int we}\ge2.            \tag{6}
\]
No common choice of $U$ with the r091 exceptional set is needed for the proof: (4)–(5) concern the total source and the r091 load concerns its own $G$. This is a necessary incidence growth for efficient sources, not a universal upper control of that growth.

**Matching near-cap load and volume.** Put $W=\sqrt{N/L}$ and $t_N=\sqrt{L-8\log L}$. On r091's further retained $G_{\rm tail}$ a strict indicator threshold $p_N'=cN^{-1/2}L^{7/2}/2$ works for $E_{y_N}$. Its volume is at least $(1-o(1))Q_NN^{-K}J$. Mills' bounds give $\bar\Phi(t_N)\asymp N^{-1/2}L^{7/2}$; multiplication of this tail by any fixed positive constant changes its quantile by $O(1/\sqrt L)$, as follows by integrating the normal hazard across $t_N\pm C/\sqrt L$. The additive $2N^{-4}$ is negligible relative to $p_N'$. Thus r091 (6) for this *small* indicator threshold implies $\log H_N(p_N')\le\mu N+y_N\sqrt N+O_K(W)$. Since $G_{\rm tail}$ lies inside its strict indicator halo, this gives the opposite bound on the layer volume. Together with its exact original price bound,
\[
 \log(|E_{y_N}|/J)=-d-y_N\sqrt N+O_K(W).                            \tag{6a}
\]
For all the original $G$'s witness selections, $\int we\ge (1-o(1))|G|$, $S_0\le J$, and $\int_{E_{y_N}}w\le |G|$. Equations (6) and (6a) now give
\[
 \log\frac{|E_{y_N}|^{-1}\int_{E_{y_N}}w}{S_0^{-1}\int we}
 =d+y_N\sqrt N+O_K(W)=2d-3\sigma W\log L+O_K(W).                  \tag{6b}
\]
In particular the ratio of this logarithm to $d$ tends to $2$. Indeed Gaussian Mills in the defining root yields $z^2=L-2\log L-\log(2\pi\mu)+o(1)$, since $\bar\Phi(z)\sim\phi(z)/z$ and $z^2/L\to1$. Expanding $t_N$ and $z$ gives the last equality. The $O_K(W)$ uncertainty is larger than $L$.

## 3. Fixed and feedback calibrations of the critical indicator threshold

For a terminal law $S$, write $s_p(S)=\sup\{s:\mathbb P(S>s)>p\}$. For any fixed $B\in\mathbb R$, $p=1-\tau-BL/N$ is in $(0,1)$ eventually. Uniformly over the deterministic independent increment vectors of r022 with sum $S$, mean $M$, variance $V$, and third cumulant $K_3$,
\[
 s_p(S)=M-\sqrt V z+\frac{K_3}{6V}(z^2-1)+O_B(1).                \tag{7}
\]
To verify it, use r022 (10), the absolute one-term Edgeworth formula uniform over vectors on the whole real line. At the uncorrected quantile $M-\sqrt V z$, the density scale $f_*=\phi(z)/\sqrt V\asymp L/N$ since $V\asymp N$. The Edgeworth correction in the cdf is $-f_*K_3(z^2-1)/(6V)$. Across any displacement $O(L)$ the Gaussian Taylor remainder divided by $f_*$ is $O(L^{5/2}/\sqrt N)=o(1)$; the change of the Edgeworth correction divided by $f_*$ is also $O(L^{5/2}/\sqrt N)$ by differentiation of $(1-x^2)\phi(x)$. The absolute $O(1/N)$ error divided by $f_*$ is $O(1/L)$. Hence displacing the stated corrected location by either sign of a large constant depending on $B$ brackets $\tau+BL/N$. Continuous positive density resolves strictness.

Optimizing (7) gives
\[
 \sup_{\text{fixed vectors}}s_p(S)=q_0+c_Q L+o(L).                 \tag{8}
\]
Here are the localization details. Write $D=\mu N-M$ and $\lambda=d/(\mu N)$. By r022 (2), $|V-\mu N|+|K_3-N\mu(3-\mu)|\le C\sqrt{ND}$. Since $K_3/V=O(1)$, any vector within the relevant region $s_p(S)\ge q_0-O(L)$ must satisfy $D\le C L+C z\sqrt D$, hence $D=O(L)$. A vector at $c_i=\mu$ reaches this region. There $\sqrt V z=d+\lambda(V-\mu N)/2+o(1)$ by Taylor, and $K_3/(6V)=(3-\mu)/6+O(\sqrt{L/N})$. The coordinatewise minimization r022 (15) shows the least value of $D+\lambda(V-\mu N)/2$ is $-g z^2+O(N\lambda^3)=-gL+o(L)$, attained by identical minimizers. This proves (8), using $(3-\mu)/6+g=c_Q$. The shifted iid reference $a_N=\mu-\sigma\rho/2$ attains the same coefficient, by Taylor of $m,v$ or (7). r028's actual physical lower construction applied to a step profile $1_{t>s}$ gives true indicators at threshold $p$ (equivalently use height $1/p$ on that set), hence $\log H_N(p)\ge q_0+c_Q L+o(L)$ in the liminf sense. Its lower direction uses common-source geometric construction, not feedback.

For the full all-predictable *pure hit* model the verified statement here is
\[
 c_R\le\liminf\frac{\sup_\pi s_p(S_\pi)-q_0}{L}
 \le\limsup\frac{\sup_\pi s_p(S_\pi)-q_0}{L}\le C_H.             \tag{9}
\]
The last bound is r021 directly, as in (3). For the lower bound let $q_b$ be the r029 (20) plateau root and $t=s_p(S_{a_N})$ the actual shifted-iid hit quantile for $a_N=\mu-\sigma\rho/2$. The continuous positive density gives $F_a(t)=\tau+BL/N$ exactly. By (7)–(8) for this reference and r029 (20), $t=q_0+c_Q L+o(L)=q_b+o(L)$. Thus for fixed small $\epsilon>0$, $s=t+(1-6\epsilon)gL$ lies in the entire $q_b+O(L)$ window of r029 (17),(19). The uniform local density in r029 (20) gives $F_a(s)-F_a(t)=(1-6\epsilon)g f_N(q_b)L(1+o(1))$. At this very $s$, the r029 controlled CDF reduction exceeds that increment by a positive $\epsilon$-multiple of $f_N(q_b)L$ (use the same slack in r029 (25), or take its internal loss parameter smaller first). Hence $F_\pi(s)<\tau+BL/N$ for all large $N$. Let $\epsilon\downarrow0$, yielding the lower in (9) without an assumption about the soft surplus under arbitrary hit policies.

The sharp r044 upper is for the **soft endpoint objective**; it does not on its own bound the pure hit quantile in (9) by $c_R$. An arbitrary hitting strategy can change its terminal surplus on paths that already crossed the target. The exact coefficient for the standalone controlled indicator upper is not proved here. Its proved lower already precludes any improvement over $c_R$ using arbitrary feedback as the only upper surrogate for (2). The $O(L/N)$ change of threshold affects a fixed-vector rank by $O_B(1)$ in (7); no pointwise transfer of that regularity to actual $H_N$ or all policies is used.

## 4. Integrating all layer prices and the limit of the rankwise feedback tests

There are two distinct relaxations, depending on which geometry remains available. The following all-layer inequality uses the **actual** indicator constant. For $v=|\{M_SF>1\}|\in(0,\infty)$ and $E_t=\{F>t\}$, set
\[
 K_N(b)=\int_0^1\min\{1,bH_N(p)\}\,dp.
\]
Any bounded nonnegative source of finite positive cost satisfies
\[
 1<\int_0^{\|F\|_\infty} K_N(|E_t|/v)\,dt,\qquad
 J(F)=\int_0^{\|F\|_\infty}\Phi_N'(t)|E_t|\,dt.                    \tag{10}
\]
The a.e. derivative in the second formula suffices because the modular is absolutely continuous on bounded intervals. For each receiver success one common rectangle has mean greater than one, and $m_RF=\int m_R1_{E_t}dt\le\int M_S1_{E_t}(x)dt$. These functions are jointly measurable by using a countable rational family of rectangles. Integrate over the halo, then integrate each indicator maximal output over its levels $p\in(0,1)$. The intersection with the same halo has size at most $\min(v,H_N(p)|E_t|)$. This proves (10) and keeps only one genuine J, though the upper side forgets compatibility of optimal layer rectangles. Monotonicity of $H_N$ guarantees the needed measurability. In probability terms $K_N(e^{-s})=\mathbb P(X_N+Y>s)$ except immaterial zero-endpoint conventions, with $X_N=\log H_N(U)$, $U$ uniform on $(0,1)$, and $Y$ an independent unit exponential. This follows by conditioning on $U$; it accounts exactly for averaging over receiver ranks. Even replacing the pointwise maximum over policies below by this integrated bound cannot remove a model lower bound.

Here is a sharp statement for a clearly specified strongest *common-policy rank model*. At critical cap $A_N=1+2\sigma\rho$, take a nondecreasing measurable $h:[0,\infty)\to D_{A_N}$, $D_A=\{0\}\cup[1,A]$, and a probability $\nu$ on $D_{A_N}$. Require
\[
 \int v\,d\nu(v)>1,\qquad
 \int T(v)\,d\nu(v)\le\sup_\pi\mathbb E_\pi T(h(S_N))             \tag{11}
\]
for every bounded nondecreasing continuous $T$ on $D_{A_N}$ with $T(0)=0$. Use the entire original price $j(h)=\int e^{-s}\Phi_N(h(s))ds$. The supremum of $1/j(h)$ in this relaxation equals the critical-cap purified feedback norm exactly. For the upper direction $T(v)=v$ gives $\sup_\pi\mathbb E h(S_N)>1$, so some single policy is strictly feasible. For the converse take $\nu=(h(S_N))_\#\mathbb P_\pi$ under that **same policy** for all tests. In fact this choice simultaneously respects tests of all threshold indicators and all bounded measurable transforms when each test is compared to its own supremum, and all their joint statistics if compared to the genuinely chosen single policy. The crucial calibration is
\[
 \log\sup_{(h,\nu)\ \mathrm{in}\ (11)} 1/j(h)=q_0+c_R L+o(L).       \tag{12}
\]
Thus the exact overestimate relative to $Q_N$ of this specified rank-only model is $gL+o(L)$ even with the r089 sharp cap, the purified height gap and simultaneous test constraints. No conclusion of this kind is drawn for genuinely sharper physical indicator estimates in (2) or (10).

**Calibration at exactly the critical cap.** The upper of (12) is r044 without clipping. For the lower take the r029 §7 plateau $H_w(t-q)$ at $q=q_b+(1-6\epsilon)gL$, $w=L/2+\log L$, with its actual policy $\pi$. Replace the height-two cap of r029 by $A_N$. The replacement keeps the step $1_{t\ge q}$ and still has a soft remainder with derivative bounded by $2/N$ and derivative total variation at most $4/N$. Indeed it is initially zero, then grows exponentially with derivative between $1/N$ and $A_N/N$, then is constant. The r029 (22)–(23) convolution/parameter-derivative proof uses exactly these two derivative bounds for the soft part and the *unchanged* policy; its $|\mathbb E_\pi u-\mathbb E_{a_N}u|=O(L/N)$ remains uniform for the newly clipped soft part.

Under the reference iid law at $a_N=\mu-\sigma\rho/2$, the additional clipping loss relative to the unbounded plateau is $O(1/(N\sqrt L))$. To check this, write $t_*=q+w+N\log A_N\ge \mu N+\sigma\sqrt{NL}-O(L)$ (as $z\le\sqrt L$) and $M=Nm(a_N)\le\mu N$. With $\theta=\sqrt L/(\sigma\sqrt N)$, the common fixed-neighborhood cumulant expansion gives $\log\mathbb E e^{\theta(S-t_*)}\le-L/2+o(1)$: the variance per step is $\mu+O(\rho)$, $N\theta^3=o(1)$ and $\theta L=o(1)$. The removed mean is
\[
 A_N\int_0^\infty N^{-1}e^{u/N}\mathbb P(S>t_*+u)\,du
 \le A_N\frac{e^{-L/2+o(1)}}{N\theta-1}=O(1/(N\sqrt L)).
\]
The r029 pure-step CDF gain is unchanged, and its soft part differs from reference by $O(L/N)$ before or after clipping. Both errors are $o(f_N L)$ against the strict slack of r029 (25), where $f_N L\sim L^2/N$. The new profile therefore remains strictly feasible. Its original price only decreases, $e^{-q}\le j(h)\le e^{-q}K_N(w)$ with the r029 plateau cost factor $K_N(w)\to1$ (this cost factor is distinct from the indicator function $K_N(b)$ in (10)). Since the original positive-set rank mass is $e^{-q}$ exactly, the calibrated example even has genuine price-to-base ratio $1+o(1)$. Let $\epsilon\downarrow0$ after the dimension limit, giving (12). This uses r029 for the actual CDF improvement, without attributing $c_R$ to any actual physical weak norm.

If $H_N$ in (10) is replaced by the rank surrogate $H_N^{\rm hit}(p)=\exp(\sup_\pi s_p(S_\pi))$, the corresponding $K_N(e^{-s})$ dominates $\mathbb P_\pi(S_N>s)$ for every fixed policy, even before the independent extra exponential smoothing. Indeed each $p<\mathbb P_\pi(S_N>s)$ has $H_N^{\rm hit}(p)\ge e^s$. Use normalized global layer volumes $b_t=\int_0^\infty e^{-r}1_{h(r)>t}dr$ of the same strictly feasible example: for a nondecreasing $h$ this equals $e^{-s_t}$ for its layer cutoff, including empty layers by the limit convention. Then $\int K_N(b_t)dt\ge\mathbb E_\pi h(S_N)>1$. Its price is exactly $\int\Phi'_N(t)b_tdt=j(h)$. Thus the integrated-layer/feedback-hit relaxation fails to exclude the explicit $c_R$ examples even with one consistent cost, while its own upper could be higher than $c_R$ because it may optimize the policy separately at each threshold. Formula (12) gives the exact obstruction for (11), not an unproved sharp pure-hit asymptotic.

**A single real rectangle can carry the same joint price and height data.** For the preceding capped profile, let $\lambda$ be the law of $h(T)$ with $T$ a unit exponential, and $\nu$ its law under the successful feedback sum. Both are probability measures on the same disconnected height range, and $\int\Phi_N d\lambda=j(h)$. For every finite $N$, $\nu\le M_N\lambda$ for some finite $M_N$. On $0\le t\le t_*$, the last conditional increment density is $\le2$, so the terminal density divided by $e^{-t}$ is $\le2e^{t_*}$. Above $t_*$ all mass is pushed to the same capped value, with $\lambda$-atom at least $e^{-t_*}$. These facts give the domination on every Borel set. Inside the unit physical cube, choose one subcube $R$ of volume $m<\min(1/2,1/M_N)$. Put the law $\nu$ in $R$, and the probability law $(\lambda-m\nu)/(1-m)$ in its complement using atomless uniform coordinates. Extend the resulting source by zero off the unit cube. Its entire price is exactly $j(h)$ and its occupied global layer volumes are $\lambda((t,\infty))$ for $t>0$. Every receiver in $R$ has *the same actual witnessing rectangle* with mean $\int v d\nu>1$ and the same joint law for every transformed test. The demonstrated physical volume is only $m$, not one. This construction verifies one-rectangle consistency and makes explicit the missing macroscopic volume comparison; it does not establish the model score physically.

## 5. Scale check and open comparison

The final [experiment92.py](../../code/experiment92.py) uses local mpmath at 80/160 digit targets with adaptive guard digits to evaluate thirty-six critical Gaussian tail-rank windows and their positive-part moment diagnostics. Independent [verify92.wl](../../code/verify92.wl) computes the same roots with a different initialization and higher precision guard, and checks twelve symbolic identities. Both retain the small rank displacement before multiplying by $\sqrt N$, rather than subtracting two ranks of order $N$. At large $N$ the serialized 65-digit quantiles may coincide; the separately computed rank difference and its strict monotonicity remain the checked observables. Full replay and scope are in [certification92.md](certification92.md) and [validation92.json](validation92.json). These are scalar checks, not a finite physical onset or a calibration of $C_H$ or $B_P$.

The residual geometric question for (2) is now one true indicator near-one upper at the stated location $\beta_N$. A fixed-ratio lower matches $c_Q$ there, whereas admitting all predictable hit policies already admits at least $c_R$. Establishing an improvement requires bounding the *physical* restricted weak indicator norm there, or controlling the fraction $S_0/J$ together with physical indicator extremality. An argument using only (11), or price plus the feedback surrogate in (10), cannot yield that improvement. The coefficient of general $\log C_N$ remains bracketed between $c_Q$ and $c_R$.
