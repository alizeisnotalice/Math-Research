# Common-rectangle height tests and the indicator halo

Write $H_N(\beta)=\sup_{0<|E|<\infty}|\{M_S1_E>\beta\}|/|E|$, $0<\beta<1$, for the actual axis-parallel noncentred operator on $\mathbb R^N$. The first calculation determines a joint-covering constant and realizes the every-witness qualifications on physical indicator shells. The second derives the global sizes of nested height layers and their incidence with selected actual rectangles for every source of polynomially near-optimal efficiency. The logarithmic gap for unrestricted $C_N$ remains open.

## 1. Exact common-rectangle constant at proportional layer thresholds

Fix $1=p_0>p_1>\cdots>p_k>0$. Take genuine nested measurable sets $E=E_0\supset E_1\supset\cdots\supset E_k$ with $|E_i|=p_i|E|$. Let $A$ consist of receivers admitting one common rectangle $R$ with
\[
 m_R1_{E_i}>\beta p_i\quad(0\le i\le k).
\]
Then, exactly in each dimension and for any finite number of layers,
\[
 \sup_{E,(E_i)}\frac{|A|}{|E|}=H_N(\beta).                                    \tag{1}
\]
The same constant results if a prescribed $\eta>0$ imposes all the further restrictions
\[
 \beta<m_R1_E<\beta+\eta,\qquad
 |m_R1_{E_i}/m_R1_E-p_i|<\eta\quad (1\le i\le k).                             \tag{2}
\]
For arbitrary nonnegative prices $w_i$ with $\sum_iw_ip_i>0$, replacing $|E|$ by $\sum_iw_i|E_i|$ divides the exact constant by $\sum_iw_ip_i$. This includes the actual increments of a step source's physical modular when its nested heights are fixed.

**Proof.** The base condition puts $A$ inside $\{M_S1_E>\beta\}$. For the reverse comparison, fix any finite-volume $E$. On $I\cap E$ for each half-open cube $I$ of a refining regular grid, use an atomless uniform coordinate and take nested portions of fractions $p_i$. Each resulting $E_i^{(j)}$ has exactly the stipulated global volume. For each fixed bounded rectangle $R$,
\[
 \sup_i |m_R1_{E_i^{(j)}}-p_i m_R1_E|\longrightarrow0.
\]
Only boundary cubes contribute to the discrepancy; the bound by their volume inside $R$ is independent of the tested fraction. Given any $x$ in the strict indicator halo, take a rectangle with density greater than $\beta$. If needed, expand it continuously about its center while retaining $x$; its density is continuous by boundary-nullness and eventually tends to zero since $|E|$ is finite. Thus there is a fixed witness with density strictly between $\beta$ and $\beta+\eta$. The lower thresholds have a positive margin since $p_k>0$. This same rectangle eventually meets all the restrictions. Fatou gives (1) and (2). Strict finite witness conditions allow slightly enlarged rational rectangles and measurable halo sets.

More generally, the atomless coordinate can assign any prescribed bounded occupied height law $K$ on $[1,A_*]$. Its conditional CDF and all its tails converge uniformly on each fixed rectangle of positive occupied volume, by the same boundary estimate applied to every threshold. For exact full global law $K$, the original price is exactly $|E|\int\Phi_N(v)\,dK(v)$. If $\int v\,dK(v)=\beta^{-1}$, the witness chosen above is also strictly successful for the resulting nonconstant-height physical source for all sufficiently fine grids. Thus the supremum remains $H_N(\beta)/(\int\Phi_N\,dK)$ for the gated joint halo requiring strict physical success, the base band (2) and any fixed positive conditional CDF tolerance relative to $K$. To see the measurability for this full-law version one can take strict tolerance in Kolmogorov distance: rectangle integrals of normalized threshold indicators vary uniformly under changing a rectangle by a small symmetric difference with its occupied volume bounded below. This argument does not bound the entire unfiltered physical halo by the gated one.

### All successful rectangles on an actual indicator shell

Let $E$ be a finite union of $K_0$ bounded closed rational axis-parallel boxes with disjoint interiors. Let $g:[0,1)\to[1,A_*]$ be a bounded quantile of any probability law $K$, with mean $a>1$ and $\beta=1/a$. Choose a sequence $\varepsilon\downarrow0$ dividing all box side lengths into integers. On $E$ set
\[
 f_\varepsilon(v)=g\left(\left\{\textstyle\sum_{i=1}^N v_i/\varepsilon\right\}\right),
 \qquad f_\varepsilon=0\ \,\text{off }E .
\]
Integration along one coordinate through complete periods gives the global law $K$ on $E$ exactly, and $J(f_\varepsilon)=|E|\int\Phi_N\,dK$ for every mesh. Fix any $x\notin E$ and $b_x=\operatorname{dist}_\infty(x,E)>0$. For every rectangle $R\ni x$ meeting $E$, its longest side has length at least $b_x$. For any $0\le w\le1$ on the phase circle, with mean $\bar w$, one has
\[
 \left|m_R\left[1_E(w(\{\textstyle\sum v_i/\varepsilon\})-\bar w)\right]\right|
 \le 2K_0\varepsilon/b_x .                                                         \tag{3}
\]
Indeed, integrate the longest coordinate on the intersection of $R$ with each constituent box. Along each transverse fiber all complete periods cancel and two endpoint fragments have total length at most $2\varepsilon$. The transverse volume is at most the cross-section of $R$, giving (3) after summing. Rectangles disjoint from $E$ have zero discrepancy. The bound works for arbitrarily eccentric rectangles, every phase threshold and every witness, simultaneously. Using stripes along just one fixed axis would not yield this uniform estimate.

In particular
\[
 \sup_{R\ni x}|m_R f_\varepsilon-a m_R1_E|\le 2A_* K_0\varepsilon/b_x.             \tag{4}
\]
Let $\eta>0$ and $\beta<M_S1_E(x)\le\beta+\eta$. For all sufficiently small meshes this receiver has an actual strict success by (4). **Every** strict witness through it then has occupied density between $\beta-O(\varepsilon/b_x)$ and $\beta+\eta$, mean between $1$ and $1+a\eta+O(\varepsilon/b_x)$, and occupied conditional Kolmogorov distance from $K$ at most $2K_0\varepsilon/(b_x m_R1_E)$. These estimates follow by applying (3) to each threshold of $g$ and using (4). The assertions are uniform on shell subsets with a common positive lower bound on $b_x$ and on $M_S1_E(x)-\beta$, with $N,E$ fixed. No dimension-uniform mesh size or lower bound on the volume of an arbitrary shell is assumed.

To see the particular r090 scales, write $\tau=d/N=\bar\Phi(z)$, $L=\log N$ and take the occupied height
\[
 V=1+\tau+\frac{\sigma}{\sqrt N}\max(-z,\min(Z,z)),\qquad Z\sim N(0,1).
\]
Here $a=1+\tau$, $1\le V\le1+2\tau\le1+2\sigma\sqrt{L/N}$ eventually, and $\sqrt N(V-1)-d/\sqrt N$ has Kolmogorov distance at most $\tau$ to $N(0,\mu)$. Put $\eta=L^2/N$, $\beta=(1+\tau)^{-1}$. Since $1-\beta=\tau+O(L/N)$ and $\tau=o(LN^{-1/4})$, for every shell receiver all its successful rectangles eventually satisfy the three r090 qualifications with dimension-independent slack constants, taking the mesh fine enough for that receiver and dimension. This provides actual common-side incidences for every successful rectangle on the shell, including spatially arbitrary finite rational-box unions $E$. For example a single unit cube has a positive-measure shell for each such $\beta$ and sufficiently small positive $\eta$: take receivers in the open transverse cube, first coordinate $u>1$, for which $M_S1_E=1/u$.

The original price is important. The clipped normal has an atom $\tau$ at $1+2\tau$. Thus
\[
 \frac{\int\Phi_N(v)\,dK(v)}{\Phi_N(1+\tau)}
 \ge \tau\frac{\Phi_N(1+2\tau)}{\Phi_N(1+\tau)}
 =\exp\{d+O(L)\}.
\]
The equality on the right concerns the explicit lower-bound expression, by the uniform quadratic expansion and $\log\tau=O(L)$. The fine mixing demonstrates that common sides, an almost fixed support deficit, strict average, and the full central Gaussian height law can coexist with the incidence pattern of an indicator shell. It does not transfer that shell to a near-extremal source at the original cost. A bound improving the general coefficient must still control how the expensive levels are spatially reused across receivers relative to their true volumes. Equation (1) rules out any uniform multiplicative gain over the indicator constant for the indicated proportional tests, even with arbitrarily many tested layers or arbitrarily narrow positive tolerances.

## 2. Global layer volumes and common-rectangle incidence in every efficient physical source

Fix $K>0$, write $h=L N^{-1/4}$, and consider a critical-cap source of finite positive original cost $J$ whose physical strict halo $H$ has $|H|\ge Q_N N^{-K}J$. For fixed $A<\infty$, let
\[
 \lambda_y=1+d/N+y/\sqrt N,\qquad E_y=\{F>\lambda_y\},\qquad |y|\le A .
\]
The following conclusions hold for purified and gray sources. For sufficiently large $N$ depending on $K,A$, uniformly on the entire compact interval,
\[
 \boxed{\log\frac{|E_y|}{J}=-d-y\sqrt N+O_{K,A}(N^{1/4}L).}                         \tag{5}
\]
Thus two fixed central layers have exponentially different global sizes on the $\sqrt N$ scale despite each carrying positive volume in every retained witness.

**Indicator estimate used in the proof.** For $0<p<1$, $0<1-p+2N^{-4}<1$, r011 and r021 give
\[
 \log H_N(p)\le\mu N+C_H L+\sigma\sqrt{N+t_0}\,
           \Phi^{-1}(1-p+2N^{-4}),\qquad t_0=O(L).                                  \tag{6}
\]
Indeed, if a strictly successful physical indicator halo contains decreasing-volume rank $a$ in the sense of r011, the normalized profile from $1_D^*$ at that rank is $1_{S_N>\log(a/|D|)}$ under the feedback law. For general finite-volume $D$, truncation only decreases $|D|$. The last increment has a density, so the endpoint choice makes no difference; if $a<|D|$ the profile is identically one on the terminal domain. Consequently there is a policy for which the CDF at $\log(a/|D|)$ is less than $1-p$. Apply r021's all-policy lower bound on the CDF and let the rank approach the strict halo volume. No multidimensional indicator lower estimate is assumed.

**Proof of (5).** Put $P=K+2$ and use r090 with exponent $P+1$, reserving room for the single additional tail test below. Write $G=H\setminus B$ with $|B|\le Q_NN^{-(P+1)}J$, so $|G|\ge(1-N^{-2})Q_NN^{-K}J$. On every successful rectangle through $G$, the original whole-rectangle CDF of $\sqrt N(F-1)-d/\sqrt N$ is within $C_P h$ (absorbing the dependence on $P+1$ in $C_P$) of $\Phi(y/\sigma)$, including gray inputs by r090's physical purification transfer. Thus for any $|y|\le A$,
\[
 G\subset\{M_S1_{E_y}>\bar\Phi(y/\sigma)-(C_P+1)h\}.
\]
This is an actual common receiver set for the full continuum of layers. In (6), for this threshold the inverse quantile is $y/\sigma+O_{P,A}(h)$, since the derivative of $\Phi^{-1}$ is bounded on a slightly enlarged compact probability interval. Hence $\log H_N(p)\le\mu N+y\sqrt N+O_{K,A}(N^{1/4}L)$. Subtract this from the physical lower estimate on $|G|/J$ to obtain the lower side of (5). For the upper side the original cost gives $\Phi_N(\lambda_y)|E_y|\le J$. The critical-window expansion in r090 yields $\log\Phi_N(\lambda_y)=d+y\sqrt N+O_A(L)$. Positivity of $|E_y|$ follows from $|G|>0$ and the strict positive indicator threshold. The sources may vary with $N$.

**Integrated incidence.** Enumerate bounded rational rectangles and, for $x\in G$, choose measurably the first containing $x$ in its interior on which $m_RF>1$. Outward perturbation of a strict rectangle makes this possible. More generally take any measurable such selection $R(x)$, and define its load on the actual source space by
\[
 w(v)=\int_G\frac{1_{R(x)}(v)}{|R(x)|}\,dx .
\]
Tonelli gives $\int w=|G|$ and, simultaneously for every $|y|\le A$,
\[
 \frac1{|G|}\int_{E_y}w(v)\,dv=\bar\Phi(y/\sigma)+O_K(h).                         \tag{7}
\]
Each identity integrates the **same** physical rectangles over the **same** retained receivers. In particular for $-A\le u<v\le A$, (5) and (7) give
\[
 \boxed{\log\frac{|E_v|^{-1}\int_{E_v}w}{|E_u|^{-1}\int_{E_u}w}
              =(v-u)\sqrt N+O_{K,A}(N^{1/4}L).}                                      \tag{8}
\]
This holds for every measurable successful-witness selection on $G$. It is a necessary incidence constraint on all polynomially near-optimal inputs. An asymptotic positive separation in (8) requires $(v-u)\sqrt N\gg N^{1/4}L$.

**Finite base measure for gray sources and a moving upper layer.** Put $e=\min(F,1)$ and $S_0=\int e$. Then $0<S_0\le J$ even when $\{F>0\}$ has infinite measure. On every retained witness $m_Re=1-m_R(1-F)_+>b$ with $b=1-d/N-(C_P+1)L^2/N$, by the gray deficit law. Since $\Phi_N(t)/t$ is nondecreasing,
\[
 |G|\le|\{M_S(e/b)>1\}|\le C_N\!\int\Phi_N(e/b)
      \le C_N\Phi_N(1/b)S_0 .
\]
Using $\log C_N=\log Q_N+O(L)$ and $\log\Phi_N(1/b)=d+O_K(L^2)$ gives
\[
 J/S_0\le\exp(d+O_K(L^2)),\qquad
 |E_y|/S_0\le\exp(-y\sqrt N+O_K(L^2)).                                      \tag{9}
\]
The second inequality holds uniformly for $0\le y$ with $\lambda_y\le1+2\sigma\sqrt{L/N}$: the original cost controls $\Phi_N(\lambda_y)|E_y|$, and the quadratic expansion gives $\log\Phi_N(\lambda_y)=d+y\sqrt N+O(L)$ on this range. On $E_y$ the density $e$ equals one. Also $b|G|\le\int w e\le|G|$ by Tonelli. For purified inputs $S_0$ is the volume of the positive set and $\int we$ is its ordinary witness load.

**One additional tail test near the top of the cap.** Specify
\[
 t_N=\sqrt{L-8\log L},\qquad y_N=\sigma t_N,\qquad
 \Delta=\sigma/\sqrt L,\qquad a=\Delta/(2\sqrt N).
\]
For large $N$, $\lambda_{y_N}$ is inside the physical cap since $t_N\le\sqrt L$ and $d\le\sigma\sqrt{NL}$. We prove that, after deleting one additional small set, **every** successful rectangle has $m_R1_{E_{y_N}}\ge cN^{-1/2}L^{7/2}$. The ordinary $O_K(h)$ Kolmogorov error is larger than this probability.

First let $F$ be purified. For $Y(v)=\sqrt N(v-1)-d/\sqrt N$, take the increasing ramp $\psi(w)=\max(0,\min(1,(w-y_N)/\Delta))$. On $D=\{0\}\cup[1,1+4\sigma\sqrt{L/N}]$ define $T(0)=0$, $T(v)=v-a\psi(Y(v))$ for $v\ge1$. For two occupied heights $v_2>v_1$, the subtracted difference lies between zero and $a\sqrt N(v_2-v_1)/\Delta=(v_2-v_1)/2$. Hence $T$ is increasing, continuous on $D$, has $T(1)=1$, and takes values in $[0,3]$ for large $N$. Write $D_*=P+5$ and $s=q_0-D_*L$. Use the Gaussian $G'$, envelope $W_s$ and occupied subprobability $\nu_s$ of r090 (4)–(6). Set $b_s=\int\psi\,d\nu_s$ and $\xi=C_{D_*}L^2/N+6N^{-4}$ with the constant of the upper bound on $\mathbb E W_s(G')-1$. Then
\[
 \mathbb E T(W_s(G'))+6N^{-4}\le 1-a b_s+\xi.
\]
The r090 (4) bridge bounds the strict physical halo $B^-:=\{M_ST(F)>1-a b_s+\xi\}$ by $e^sJ(F)$. Specifically, the normalized rearranged terminal of the transformed source is $T(F^*(re^{-t}))$ at volume rank $r$, and $r^{-1}\int_0^r\Phi_N(F^*(v))\,dv\le J(F)/r$ bounds the **original** terminal envelope. Monotonicity of $T$ preserves this domination, allowing the r021 increasing-test comparison. The price here is the actual $J(F)$, including the low part, throughout.

If $m_RF>1$, then $m_RT(F)>1-a r_R$ for $r_R=m_R[1_{F\ge1}\psi(Y(F))]$. Consequently every such $R$ through $x\notin B^-$ satisfies $r_R\ge b_s-\xi/a$, and $m_R1_{E_{y_N}}\ge r_R$. By r090 (5)–(6), $\nu_s(\mathbb R)=1-O_{D_*}(\sqrt{L/N})$ and its CDF differs from $\Phi(\cdot/\sigma)$ by $O_{D_*}(L/\sqrt N)$. Using the mass of $(y_N+\Delta,\infty)$,
\[
 b_s\ge \bar\Phi(t_N+1/\sqrt L)-O_{D_*}(L/\sqrt N).
\]
Since $(t_N+1/\sqrt L)^2=L-8\log L+2+o(1)$, Mills' formula gives $\bar\Phi(t_N+1/\sqrt L)\sim \exp(-1)N^{-1/2}L^{7/2}/\sqrt{2\pi}$. The perturbation error is $\xi/a=O_{D_*}(N^{-1/2}L^{5/2})$. Thus $m_R1_{E_{y_N}}\ge cN^{-1/2}L^{7/2}$ for all the asserted rectangles. The positive constant $c$ may be chosen absolute once $D_*$ is fixed and $N$ is sufficiently large depending on $D_*$.

For gray $F$, choose the r090 sequence of purified $F_j$ having exactly the same $J$ and unchanged values on $\{F>1\}$. Every fixed strict original rectangle is eventually strict for $F_j$. The tested high layer $\{F_j>\lambda_{y_N}\}=E_{y_N}$ is unchanged for every $j$; the inserted heights zero and one do not contribute. For the numerical threshold $cN^{-1/2}L^{7/2}$ one obtains precisely the same lower bound outside $B^-=\liminf_jB_j^-$: every point of its complement is outside infinitely many $B_j^-$, and Fatou gives $|B^-|\le e^sJ$. No assertion that the gray source's zero set survived purification is required.

There is just **one** additional physical exceptional halo per purified input. Since $D_*=P+5$ and $\log Q_N=q_0+\alpha L+o(L)$, eventually $e^s\le Q_N N^{-(P+1)}$. Combining with the original $B$ gives $|B\cup B^-|\le 2Q_N N^{-(P+1)}J\le Q_NN^{-P}J$. Define $G_{\rm tail}=G\setminus B^-$. Then $|G_{\rm tail}|/|G|\ge1-O(N^{-2})$ by $K$-efficiency. The stronger upper-layer property holds pointwise on $G_{\rm tail}$. For the original load $w$ of **any** chosen successful rectangles on all of $G$, positivity and restriction to $G_{\rm tail}$ give $\int_{E_{y_N}}w/|G|\ge (|G_{\rm tail}|/|G|)cN^{-1/2}L^{7/2}$. Combining (9), $\int we\le|G|$, and the $O(L)$ logarithm of that retained layer fraction yields
\[
 \boxed{\displaystyle
 \log\frac{ |E_{y_N}|^{-1}\int_{E_{y_N}}w}{ S_0^{-1}\int we}
 \ge y_N\sqrt N-O_K(L^2),\qquad
 \liminf_{N\to\infty}\frac1d\log
 \frac{ |E_{y_N}|^{-1}\int_{E_{y_N}}w}{ S_0^{-1}\int we}
 \ge1.}                                                                       \tag{10}
\]
The liminf applies along every sequence of the stipulated $K$-efficient sources because $t_N/z\to1$ and $L^2=o(d)$.

For comparison, the single uniform CDF conclusion of r090 already suffices at the lower layer $\sigma\sqrt{L/2-4\log L}$: its unperturbed Gaussian tail is asymptotic to $N^{-1/4}L^{3/2}/\sqrt\pi\gg h$. This gives the corresponding $1/\sqrt2$ lower coefficient on $G$ without the extra exceptional halo. For the new test, local Python evaluation at $L=128,512,2048,8192$ bounded $\bar\Phi(t_N+1/\sqrt L)/(N^{-1/2}L^{5/2})$ by Mills' elementary inequalities in the intervals $[25.91,26.21]$, $[82.72,82.90]$, $[309.37,309.54]$ and $[1212.56,1212.73]$ (rounded outward). These are scale checks; no finite onset uniform in $K$ is inferred.

The global tail estimates and (8) distinguish homogeneous mixing from efficient spatial use. Their $N^{1/4}L$ uncertainty exceeds the target log-volume accuracy $o(L)$. They do not furnish an upper bound on how much of the required load can in fact be carried by those high layers.

## 3. Independent finite verification

The final [experiment91.py](../../code/experiment91.py) uses exact Fraction integration of phase strips by a positive-part second primitive and a four-corner difference. [verify91.wl](../../code/verify91.wl) independently clips rational polygons by the two boundaries of each diagonal strip. The sixteen two-box cases include four off-grid receiver rectangles per case, and four further cases reproduce the independent unit-square audit exactly. Original $\Phi_2$ prices, strict mean margins and nonzero boundary discrepancies are evaluated separately.

Twelve scalar dimensions up to $2^{4096}$ compare both moving Gaussian tails, their actual error scales and the top-atom lower bound for the clipped-normal price. The latter is explicitly a lower-bound expression, not the full expected modular. Python uses stable `log1p`; Wolfram uses the original nested logarithm with additional guard digits. The 80/160 digit targets, exact arithmetic qualifications, source hashes and resource reports are recorded in [certification91.md](certification91.md) and [validation91.json](validation91.json). No finite starting dimension for the physical theorem follows from these scalar samples.
