# Deep-core obstruction to optimized power certificates, local dual prices, and true slice comparisons

Use the strong rectangle basis and the physical periodic density/cover dual of r095–096. Write $N=n+1$, $L=\log N$, $d_N=\sigma\sqrt N z_N$, $\tau=d_N/N=\bar\Phi(z_N)$, $p=p_N=1-\tau-BL/N$ for fixed $B$, $q_{0,N}=\mu N-d_N$, $\log Q_N=q_{0,N}+c_QL+o(L)$. When a periodic suffix weak constant is used, the ordinary physical $H_n$ is a valid upper by large-box transference. The results below do not presume that the upper coefficient of unrestricted $H_N$ equals $c_Q$.

## 1. Gluing local source prices into one true dual certificate

Fix a finite family of physical rectangles $R_i$ of positive volume. Its original price $P_p(\mathcal R)$ is the infimum of $\int f$ needed to have $m_{R_i}f\ge p$ with $0\le f\le1$; r096 (7) gives its sharp dual and the strict physical infimum by purification/mixing. Assign each row $i$ to a group $l(i)$. On original source space choose measurable densities $0\le\omega_l$, $\sum_l\omega_l\le1$. These may depend on the actual finite family and on scale or position. Let $a_i=m_{R_i}\omega_{l(i)}>0$ and define
\[
 d\mu_l=\omega_l\,dy,\qquad p_i'=1-(1-p)/a_i.
\]
For group $l$ price rows in the *same* $R_i$ relative to $\mu_l$ with threshold $p_i'$; a nonpositive target can be dropped. If $P_l$ is that minimum local cost relative to $\mu_l$, then
\[
                           P_p(\mathcal R)\ge\sum_l P_l .                \tag{1}
\]
More strongly, arbitrarily chosen nonnegative local multipliers give a **single original-family** feasible dual vector whose objective is at least their summed local objectives. Indeed put $S_l=\sum_{i:l(i)=l}\lambda_i1_{R_i}$ and $S=\sum_l S_l$. Pointwise $\sum_l\omega_l\min(S_l,1)\le\min(S,1)$, hence
\[
 (S-1)_+\le\sum_l\omega_l(S_l-1)_++\sum_l(1-\omega_l)S_l.
\]
The last term integrates to $\sum_i\lambda_i(1-a_i)|R_i|$. Subtracting in the r096 dual leaves the exact local rows $[p-(1-a_i)]|R_i|=p_i'\mu_l(R_i)$. This proves the certificate claim and (1). One may alternatively put any original feasible $f$ in each local program and sum $\int f\,d\mu_l\le\int f$. Each family is finite, so the cell-atom LP relative to $\mu_l$ gives local strong duality. This proof charges the physical cost once even if the group assignments overlap on source space. It only needs the summed ownership density bound.

When $a_i\ge1-\gamma$ the universal safe local threshold is $p'=1-(1-p)/(1-\gamma)=p-\gamma(1-p)/(1-\gamma)$. This threshold loss cannot be removed from a general local-row implication: fill the whole discarded part of one rectangle and put exactly the still-needed mass on its retained part; then the retained mean is exactly $p_i'$. An arbitrarily small additional occupied mass gives the same test with a strict original inequality. When $\gamma=O(\tau)$, the one-row threshold change is only $O(L/N)$ at $p_N$.

**Low-collision version.** If $p_i'\ge p'>0$, choose in each group any probability mixture $q_{li}$ of its retained normalized box densities relative to $\mu_l$, $W_l=\sum_i q_{li}1_{R_i}/\mu_l(R_i)$. Then $P_l\ge (p')^2/(\int W_l^2 d\mu_l)$ by testing the local constraints and Cauchy-Schwarz. Even arbitrarily concentrated discarded portions cannot spoil the sum in (1); taking a global $L^2$ norm *before* removing these portions has no such protection. The Gram entry is $\mu_l(R_i\cap R_j)/(\mu_l(R_i)\mu_l(R_j))$, with actual intersections. For disjoint panels $U_l=\mathbb T^n\times P_l$ take $\omega_l=1_{U_l}$. Suppose physical covered fraction is $\theta$, $a_i\ge1-\gamma$, and let $\theta_l$ denote the covered fraction in $U_l$ of its own clipped rows. At a fixed horizontal point, each group interval has at least $1-\gamma$ of its length in the union of the clipped group intervals, so
\[
             \sum_l |U_l|\theta_l\ge\theta/H_1(1-\gamma).
\]
If a non-strict endpoint occurs, first reduce the level then take the limit in $H_1(t)=2/t-1$. After rescaling $P_l$ to a circle, clipped intervals are genuine periodic arcs (two edge pieces become one arc), so the local programs can also use true rescaled geometric certificates. For example, if each local cover has a physical receiver-causal assignment with coverage $\theta_l$ in the sense of r096 §2, the local Gram bound gives $P_l\ge |U_l| (p')^2\theta_l^2/3^N$. Jensen across panels yields $P_p(\mathcal R)\ge (p')^2[\theta/H_1(1-\gamma)]^2/3^N$ for a complete partition of unit volume. No price is paid for the number of externally conditioned panels. Near-$Q_N$ coverage would require expensive *within-panel* incidence for this class. Unrestricted covers are not known to admit such a partition with certified local prices.

## 2. Arbitrarily many scale-separated panels paid by the genuine row halos

Partition a circle into finitely many disjoint intervals $P_l$ of lengths $w_l$, and on $\mathbb T^n\times P_l$ take the actual indicator $1_{D_l}(K_l y)$, with periodic finite step sets $D_l\subset\mathbb T^n$, $d_l=|D_l|$, and integer $K_l>0$. Empty rows may be ignored. The exact original volume is $d=\sum_l w_l d_l$. Set $h_l(t)=|\{M_{S,n}1_{D_l}>t\}|$ (volume, not normalized). There are $0<b<T<\infty$ depending on the finite collection of base grids and $N$, and an absolute $c>0$, such that whenever the intervals $[b/K_l,T/K_l]$ are disjoint, for any $0<\delta\le L/N$ and large $N$,
\[
 \sum_l w_l h_l(p)\ \le\ |\{M_{S,N}1_E>p\}|\ \le\
 H_1(1-2(1-p)-2\delta)\sum_l w_l h_l(p-\delta)
 +\delta^{-1} Q_N e^{-cN}\,d.                                         \tag{2}
\]
Implicit adjustment of $c$ absorbs fixed multiplicative constants. This includes arbitrarily large $m$ and arbitrarily large additional gaps in the frequencies. In particular one can take $\delta=e^{-cN/2}$, still with an exponentially negligible second price relative to $Q_Nd$. The horizontal cost is $1+O_B(\sqrt{L/N})$. The input hypotheses do not include a lower density bound: tiny components are paid according to their actual volume. The conditions are on the *entire* receiver-dependent family of true rectangles, not on pointwise counterfactual costs.

For the proof fix $k=\lfloor n/2\rfloor$, $r=\lfloor k/8\rfloor$. Each base $D_l$ has a bad set $B_l$ on which a rectangle of $D_l$-mean greater than $1/2$ can have its first-$k$ median length outside $[b,T]$, with
\[
 |B_l|\le C(8/3)^r[H_{n-r}(1/3)+1]d_l .                              \tag{3}
\]
Here is the r098 random-coordinate argument with the new, fixed slack. For a fixed $r$-subset $J$, at short sides discard grid boundaries in a total neighborhood of measure $\le d_l$. Choose $b$ accordingly, taking the minimum across the finite collections of faces and positive $d_l$. Off this neighborhood all short sampled coordinates are in their receiver's respective step cells. Integrate the true $(n-r)$-dimensional strict indicator halo of the frozen slices to obtain $H_{n-r}(1/3)d_l$ plus the discarded $d_l$. For long sides each coordinate of length $>T$ in the periodic lift gives a density at most $1+1/T$ relative to the uniform circle. Choose $(1+1/T)^r<3/2$. Long success implies true $1/3$ success for the fully integrated $J$ marginal. The weak indicator bound applies to that $[0,1]$ marginal by r092 zero-cost purification, giving the same price. For any offending median, a uniformly selected $r$-subset lands wholly in the short or long coordinates with probability $\ge(3/8)^r$. Average the indicated measurable bad sets over $J$ before Markov; do not union all subsets. This gives (3). At $p_0=1/3$, r091 (6) gives $\log H_{n-r}(1/3)\le\mu(n-r)+O(\sqrt N+L)$, since the inverse Gaussian is bounded. As $\mu>\log(8/3)$ and $r\sim N/16$, (3) is at most $Q_Ne^{-c_0N}d_l$ for a fixed $c_0>0$ and large $N$. No assumption on the size of $d_l$ has entered.

Let $W(y,z)=1_{B_l}(K_l y)$ on panel $P_l$. Its integral is $\le Q_Ne^{-c_0N}d$. Discard receivers with $M_z W(y,\cdot)(z)>2\delta$ at physical cost $\le H_1(2\delta)Q_Ne^{-c_0N}d$. The median of the actual first $k$ side lengths of a remaining receiver's witness can match at most one index $l$, as the $K_l$ rescale the same rectangle exactly. For any other panel its horizontal rectangle average is $\le1/2$ unless that panel's $W(y,\cdot)=1$. Write $a=m_{I_z}1_{P_l}$ and $v$ for the actual horizontal $D_l$ average in the matched panel. With no match take $a=0$. Then the whole witness mean is at most
\[
             1/2+a(v-1/2)_++\tfrac12 m_{I_z}W .
\]
A strict success forces $a>1-2(1-p)-2\delta$ **and** $v>p-\delta$. The same receiver thus has $z$ in the specified genuine one-dimensional indicator halo of a single matched panel, and $K_l y$ in the true $p-\delta$ halo of its own base. Their product volume is at most $H_1(1-2(1-p)-2\delta)w_l h_l(p-\delta)$. Summing over all matches and adding the discarded physical volume proves the upper. The lower uses honest boxes with their $z$ side inside the corresponding panel. This is a true two-sided transfer for each source, keeping track of row prices, horizontal inflation and threshold slack separately.

**A sharp-cost class containing the range of scale codings of r098.** Each $D_l$ may, for instance, be obtained from an arbitrary bounded unconditional log-convex downward-closed *physical* packet of r094 by translation, small isotropic scaling into a unit-cell chart, and periodization. Take each bounding box coordinate width at most $1/4$. Any interval meeting two successive copies of a coordinate support in positive measure has occupancy by their periodic projection at most $2s/(1+s)$ for coordinate width $s\le1/4$: it covers at least one entire gap $1-s$ and at most two width-$s$ pieces per gap. This is $<1/2$, excluding such witnesses at $p-\delta$ for large $N$. Thus its strict periodic halo costs no more than the corresponding normalized physical upper. If the packet is not already step, (3) remains valid at slightly modified absolute thresholds by approximating it from above by a finite grid *superset* of volume at most $2d_l$; a bounded log-convex packet has boundary of measure zero in the coordinate charts and any coordinate faces are null. Use this superset only to price the short/long exceptional sets at the fixed threshold $1/2$. The matched-panel $p-\delta$ halo in (2) still uses the actual original packet and original volume. r094 (17) and r092 (8) give $h_l(p-\delta)\le d_l\exp(q_{0,N}+c_Q L+o_B(L))$ uniformly, because $n=N-1$ and the threshold shift is at most $L/N$. This also covers per-layer different log faces, thicknesses and arbitrarily different densities; arbitrary translations and separate nonsingular axis rescalings preserving the chart-width condition are harmless. Formula (2) yields the sharp static upper for this entire panel class while permitting the unbounded scale range that produces the information obstruction. The matching lower for the class follows by taking one r093 cross (scaled to the chart and repeated) independent of $z$. The r098 base made of many finite translated packet unions is also directly covered by (2), with its *actual* $h_l$ values, independent of whether its geometry is log-convex.

## 3. The optimal Gram certificate fails on a genuine macroscopic near-benchmark cover

Fix $B$ and write $p=1-d_N/N-BL/N$, $L=\log N$. There are finite families $\mathcal R_N$ of positive-volume rational periodic rectangles covering more than $1-1/N$ of the unit torus, and single binary physical sources $E_N$ with every $R\in\mathcal R_N$ strictly successful, such that
\[
 \log (1/|E_N|)=\log Q_N+o_B(L),\qquad
 Q_N^{-1}N^{-O_B(1)}\le P_p(\mathcal R_N)\le Q_N^{-1}N^{o_B(1)}.       \tag{3}
\]
For a mixture $W_q=\sum_{R\in\mathcal R_N}q_R1_R/|R|$, its single-power Hölder certificate and Rényi divergence are
\[
 C_j(q)=\frac{p^{j/(j-1)}}{(\int W_q^j)^{1/(j-1)}}
       =p^{j/(j-1)}e^{-D_j(W_q)},\qquad j>1 .                           \tag{4}
\]
All these are legitimate lower certificates for $P_p$: a common feasible $0\le f\le1$ obeys $p\le\int fW_q\le(\int f)^{1-1/j}(\int W_q^j)^{1/j}$. Nevertheless the *same* family satisfies
\[
 \sup_{q,\,1<j<\infty} C_j(q)\le \exp[-\mu N+O_B(L)]
      = Q_N^{-1}\exp[-d_N+O_B(L)].                                      \tag{5}
\]
The supremum allows both the mixing probabilities and the order to depend arbitrarily on $N$ and on the full family. In particular optimizing a Gram mixture or selecting among arbitrarily many isolated power certificates cannot give a general sharp upper via r096 (8). This is a superpolynomial deficit from the true LP value in (3), not evidence that the LP value itself is smaller than the sought bound. For each *fixed* $j>1$, letting $v=1-1/j$, the stronger statement is
\[
 \log\sup_q C_j(q)\le -[\Lambda_\mu(v)/v]N+o_{j,B}(N),
 \quad \Lambda_c(v)=-\log(1-v)+\log\frac{1+e^{-(1-v)c}}{1+e^{-c}} .     \tag{6}
\]
Strict Jensen gives $\Lambda_\mu(v)/v>\mathbb E Z_\mu=\mu$. At $j=2$ this is a fixed linear-in-dimension loss even for the **best** quadratic Gram certificate.

### One family, one source, true strict margins

Use the r092 §3 shifted iid reference $c_N=\mu-\sigma\sqrt{L/N}/2$, $r_N=e^{-c_N}$, with iid sum $S$ of $N$ increments $Z_c=Y+cB$ (unit exponential and independent Bernoulli of probability $(1+e^c)^{-1}$). Set $\delta=L/N$ and let $s=s_{p+2\delta}(S)-1$. The probability level $p+2\delta$ lies in $(0,1)$ for large $N$; r092 (7)–(8) at this reference give $s=q_{0,N}+c_QL+o_B(L)$ and $\Pr(S>s)>p+2\delta$. Let $T=N^4$. On a small cube $[-a,a]^N$ within one torus chart write $U(y)=\sum_i -\log(|y_i|/a)$. Use the packet $D=\{U>T\}$ and certified receivers $D^- =\{U>T-s\}$, both restricted to this cube, with coordinate hyperplanes ignored. At any such receiver $x$, the true box having coordinate endpoints $x_i$ and $-r_N x_i$ stays in the chart and samples $U(y)=U(x)+S$ in distribution. Thus it has $D$-mean strictly exceeding $p+2\delta$. The actual volume ratio $h=|D^-|/|D|$ satisfies
\[
 h=G_N(T-s)/G_N(T)=\exp[s+O(Ns/T)],\qquad G_N(t)=e^{-t}\sum_{k=0}^{N-1}t^k/k!\quad(t\ge0).
\]
Indeed the polynomial ratio lies between $(1-s/T)^{N-1}$ and one. In particular $h$ differs from $e^s$ by a factor $1+o(1)$.

Take $f=|D^-|\le1/4$ by reducing $a$ if needed. For $M=\lceil5L/f\rceil$ independent uniform torus translates, let $E_N$ be the union of the translated packets. Its raw source count integral is $\lambda=M|D|=(Mf)/h$ with $5L\le Mf\le5L+1$. Consequently $\lambda=o(1)$, $s+\log\lambda=O(\log L)$ and $\lambda=Q_N^{-1}N^{o_B(1)}$. The expected uncovered receiver volume outside all translated $D^-$ is $\le N^{-5}$. The expected loss in volume from source raw count to the union is at most the ordered pair overlap integral in expectation, $\le\lambda^2$. By Markov one realization has covered volume $>1-1/(2N)$ and $\lambda/2\le |E_N|\le\lambda$. Fix it once, before choosing a power or a mixture.

At each of its certified receivers retain the single-packet witness just described. Enlarge it outward to a rational open torus box by a volume factor at most $1+N^{-2}$, with the receiver in its interior. All enlargements are still strictly $p$-successful for the common original $E_N$. The candidates with rational endpoints are countable; a finite subfamily has union volume $>1-1/N$. Associate with each retained box one actual packet and its pre-enlargement witness. This constructs $\mathcal R_N$. For any LP-feasible fractional source, its strict halo at $p-\delta$ contains this full union. Periodic purification at this lower threshold gives $P_p(\mathcal R_N)\ge (1-1/N)/H_N(p-\delta)$. The r091 physical indicator upper (6) is $\log H_N(p-\delta)\le q_{0,N}+C_H L+O_B(1)$ by the critical Gaussian inverse in r092 (3); this proves the lower side of (3). The true binary $E_N$ gives its upper side.

### All orders at once: a physical entropy test

For each translated chart let $U_l$ be its own translated log-product coordinate, extended by $-\infty$ outside that chart. Define the nonnegative physical function
\[
 V(y)=\max\{0,\ \max_{l\le M}(U_l(y)-T)-\log\lambda\}.
\]
It is allowed to take infinite values on volume-null coordinate faces. Lebesgue probability on the whole torus obeys, for $u\ge0$,
\[
 |\{V>u\}|\le \exp[-(1-N/T)u].                                       \tag{7}
\]
Indeed put $x=u+\log\lambda$, so $T+x>0$ for large $N$. The union bound and the exact Gamma volume ratio give $|\{V>u\}|\le\lambda G_N(T+x)/G_N(T)$. At $x\le0$ the Gamma polynomial ratio is $\le1$; for $x>0$ it is at most $(1+x/T)^{N-1}\le\exp(Nx/T)$. Now $x\le u$ proves (7). For $t=1-1/N$, tail integration gives $\int e^{tV}\le1+t/(1-N/T-t)\le3N$ for large $N$.

In contrast, on each true pre-enlargement box, its own chart gives $V(y)\ge U_l(y)-T-\log\lambda$. The expectation on that box is at least $N\mathbb E Z_{c_N}-s-\log\lambda=\mu N-O_B(L)$: the mean function has a smooth maximum at $c=\mu$ so $N\mathbb EZ_{c_N}=\mu N+O(L)$. On every enlarged box positivity retains at least $(1+N^{-2})^{-1}$ times this expectation. Thus $\int W_q V\ge\mu N-O_B(L)$ for *every* original-box mixture. Gibbs' entropy inequality gives
\[
 D(W_q\|1)=\int W_q\log W_q\ge t\int W_q V-\log\int e^{tV}\ge\mu N-O_B(L).
\]
Truncating $V$ and passing to the limit is also valid; the finite mixture is bounded for each finite family. Finally $D_j(W_q)=\frac1{j-1}\log\int W_q^j\ge\int W_q\log W_q$ by Jensen under the measure $W_q\,dy$. As $p^{j/(j-1)}\le1$, (5) follows uniformly, including powers approaching one or infinity. The test $V$ is used to analyze failure of the power certificates; it is not substituted for the physical indicator input.

### A linear exponential loss for each fixed power

For fixed $j>1$ set $v=1-1/j$. The real cumulant generator of $Z_{c_N}$ is $\Lambda_{c_N}$ displayed in (6). Write $M_v=N\Lambda'_{c_N}(v)$ and $a_N=M_v-s-2N^{2/3}$. Since $\Lambda'_\mu(v)>\mu$ and $c_N\to\mu$, $a_N$ is positive of order $N$. In each packet keep its core $U_l>T+a_N$ and let $D'$ be the union of these cores, for the *same* previously fixed rectangles. The Gamma volume bound gives $|D'|\le\lambda\exp(-a_N+O_j(Na_N/T))$. Under the iid tilt by $v$, the sum $S$ has mean $M_v$ and variance $O_j(N)$ uniformly. Chebyshev gives $\Pr_v(|S-M_v|\le N^{2/3})=1-o(1)$. Returning to the original probability,
\[
 \Pr(S>s+a_N)\ge(1-o(1))\exp[N\Lambda_{c_N}(v)-vM_v-vN^{2/3}].
\]
Every original witness gives at least this fraction to its actual core by $U_l(x)>T-s$, and rational enlargement loses at most a factor $(1+N^{-2})^{-1}$. Denote this uniform lower by $u_N$. It also bounds $\int_{D'}W_q$ below for every mixture. Hölder in this small physical set gives $C_j(q)\le |D'|(p/u_N)^{1/v}$. Substitution yields
\[
 \log C_j(q)\le \log\lambda+s-(N/v)\Lambda_{c_N}(v)+O_j(N^{2/3})+o(1),
\]
proving (6). The core can depend on $j$, but the rectangle family, input, density, and macroscopic coverage do not. For $j=2$ the minimum Gram energy $\min_q\int W_q^2$ is therefore at least $\exp(2\Lambda_\mu(1/2)N-o(N))$. The primal cost is only bounded within the band in (3); nothing here proves a new lower or upper coefficient for $C_N$ or unrestricted $H_N$.

## 4. A near-benchmark counterexample to the exact one-dimensional multiplier at the same threshold

For physical sources $E\subset\mathbb R^N$, write $E_z$ for the $(N-1)$-dimensional slice at one coordinate $z$. It is not generally valid, even at scores $Q_NN^{o(1)}$, to upper bound the full halo volume by $H_1(p)$ times the integral of the **actual** sliced halo volumes at the unchanged threshold. There is a sequence of bounded unconditional indicators for which
\[
 \frac{|\{M_{S,N}1_E>p\}|}{\int|\{M_{S,N-1}1_{E_z}>p\}|\,dz}
 \ \ge\ 1+(\sqrt{2\pi}-o(1))\sqrt{L/N},\qquad
 H_1(p)=1+(2\sqrt\mu+o(1))\sqrt{L/N}.                                \tag{8}
\]
Here $\sqrt{2\pi}>2\sqrt\mu$. Both numerator and denominator, divided by the same $|E|$, have logarithm $\log Q_N+o(L)$. This is a strict failure of the proposed exact horizontal factor, not a fixed-power separation or an obstruction to an $N^{o(1)}$ bound.

Take $n=N-2$, suffix $D=E_{T',\mathbf1}$ the true equal-weight cross of r093 in $n$ dimensions, and prefix $A=E_{T,\mathbf1}$ the same construction in two dimensions. The common source is $A\times D$, the sliced coordinate $z$ being one of the two prefix coordinates. Write $b(y)=M_{S,n}1_D(y)$, and $h_{2,T}(t)=|\{M_{S,2}1_A>t\}|/|A|$. Factorization and the fact that every occupied prefix slice $A_z$ is an interval give the **exact** identities
\[
 |\{M_{S,N}1_E>p\}|=|A|\int_{b>p}h_{2,T}(p/b)\,dy,
 \qquad \int|\{M_{S,N-1}1_{E_z}>p\}|dz
       =|A|\int_{b>p}(2b/p-1)\,dy.                                     \tag{9}
\]
For every fixed dimension $j$ and fixed $0<t<1$, the complete static-cross indicator ratio converges as its thickness tends to infinity to $\exp(s_{j,t})$, where $s_{j,t}$ is the supreme independent deterministic-vector hit quantile. This follows even without a special high-threshold assumption: the interior Gamma ratio $G_j(T-s_{j,t})/G_j(T)$ lower and the containing-box upper $e^{jb_0}G_j(T-s_{j,t}+jb_0)/G_j(T)$, $b_0=\log(2/t-1)$, have the same limit by the finite Gamma polynomial. These are complete true-box estimates of r093. The quantile is continuous in $t\in(0,1)$: the independent vector uses a compact parameter space in the original ratios, converges in law under parameter limits (including zero ratios), and each sum cdf is continuous strictly increasing on $(0,\infty)$. Thus at fixed $N$ the distribution of $b$ under the uniform receiver on $\{b>p\}$ has a limit as $T'\to\infty$, with tail $\exp(s_{n,t}-s_{n,p})$ for $p<t<1$. It may have mass at $1$, which is harmless. Fatou as $T\to\infty$ in the first integral, for each fixed $N,T'$, permits the lower response $\exp(s_{2,p/b})$; the second integral needs no such approximation.

For the two-coordinate calibration, each increment density on positive arguments is at least $e^{-t}/2$. The cdf on $Z_1+Z_2\le s$ is consequently at least $[1-(1+s)e^{-s}]/4$. Conversely choose the two identical parameters $c=s$ in the unshifted-plus-shifted exponential formula: throughout this simplex except its null faces the cdf is $[1-(1+s)e^{-s}]/(1+e^{-s})^2$. Thus, uniformly for small $s>0$, the infimum cdf is $s^2/8+O(s^3)$, and
\[
              e^{s_{2,1-u}}=1+\sqrt{8u}+O(u)\qquad(u\downarrow0).       \tag{10}
\]
For $u=1-p/b\in(0,1-p]$, subtract $2b/p-1=1+O(u)$. Dividing the two integrals of (9) therefore yields, in the thickness-lower-limit sense, at least $1+(\sqrt8-o(1))\mathbb E\sqrt{1-p/b}$; the expectation uses the suffix receiver distribution conditioned on $b>p$, and the $O(u)$ cost is at most $O(\sqrt{1-p})$ times the main $\sqrt u$ integrand.

For completeness, on the fixed additive density scale $\xi=L/N$, every fixed $x\ge0$ satisfies
\[
                s_{n,p+x\xi}-s_{n,p}=-x+o(1).                          \tag{11}
\]
Both endpoint optimizing vectors lie in the r092 (7)–(8) near-optimal window: $n=N-2$, and each threshold differs from its $n$-root reference by $O_{B,x}(L/n)$. r092 localization gives mean $\mu n+O_{B,x}(L)$ and variance $\mu n+O_{B,x}(\sqrt{nL})$ for each such independent vector. For all ranks in any fixed $O_{B,x}(L)$ window around $q_{0,n}$ the standardized point is $-z_n+O_{B,x}(L/\sqrt n)$, and the relative Gaussian density change is $o(1)$. The uniform absolute unstandardized density error $O(1/n)$ from r022/r093 (regular unit exponential components) is $o(L/n)$; the Gaussian density at the reference point is $(1+o(1))L/N$ by Mills. Consequently the quantile difference for each endpoint's *fixed* vector is $-x+o(1)$ by integration of the density (the second quantile for that vector lies within $O_x(1)$ by the same local bound). Taking the optimizing vector at either endpoint supplies the two opposite bounds in (11), with no smoothness assumption on an optimizing-vector selection.

It follows under the limiting suffix receiver law above that $(b-p)/\xi$ tends in distribution to a unit exponential. In particular truncation and the tail integral imply $\liminf \mathbb E\sqrt{1-p/b}/\sqrt\xi\ge\int_0^\infty \sqrt x e^{-x}dx=\sqrt\pi/2$, since $1-p/b\ge b-p$. This only uses bounded intervals of $x$ before taking the truncation limit. Combining with (10) proves (8). The estimates can be achieved by finite thicknesses for each $N$: first choose $T'$ to approximate the bounded suffix test, then $T$ to approximate the lower in (9), with an arbitrarily small error compared with $\sqrt\xi$. The suffix halo may also be required to reach at least half its limiting ratio at $p$. Its full upper is $\exp(s_{n,p})$ eventually by r093, and the two-dimensional prefix full ratio is bounded by $\exp(\max\{2\log(2/p-1),s_{2,p}\})=1+O(\sqrt{1-p})$ uniformly in its thickness. Thus the product halo upper from its two factor halos, and the suffix lower alone on the prefix interior, place both normalized scores at $\log Q_N+o(L)$. Small scalings into one torus chart of coordinate width $\le1/4$ also give the same example periodically: near-one witnesses cannot span successive positive copies in any coordinate, as in §2.

### Exponentially small panel threshold and multiplier slack

In (2) one can reduce the fixed unmatched mean ceiling $1/2$ to $u_N=e^{-aN}$ for a sufficiently small fixed $a>0$, while keeping an exponential margin in the source-priced exceptional volume. In the short/long proof of (3) test the reduced dimension at $u_N/2$, taking $(1+1/T)^r<2$. The r001 physical strong $L^t$ bound at a *fixed* large $t$ and periodic transference give $H_{n-r}(u_N/2)\le (2/u_N)^t\exp((n-r)b_t)$, $b_t\to\mu$ as $t\to\infty$. Choose $t$ sufficiently large, then $a$ sufficiently small: as $r\sim N/16$ and $\mu-\log(8/3)>0$,
\[
 (8/3)^r [1+H_{n-r}(u_N/2)]\le Q_N e^{-c_1 N}
\]
with $c_1>2a$. The same self-contained local exception may be handled for nonstep packets using outer finite steps of volume $\le2d_l$ as before. Set $\delta=e^{-aN}$, discard $M_zW>\delta/(1-u_N)$, and use the nonmatching mean bound $u_N$ in the witness calculation. It forces the matched actual base mean $>p-\delta$ and the matched panel fraction $>1-[(1-p)+\delta]/(1-u_N)\ge p-O(e^{-aN})$. Consequently the upper becomes
\[
 |\{M_{S,N}1_E>p\}|\le [H_1(p)+O(e^{-aN})]\sum_l w_l h_l(p-e^{-aN})
                        +Q_N e^{-aN}\sum_l w_l d_l.                    \tag{12}
\]
The choices of $b,T$ may be adjusted for the new ceilings. No uniform same-threshold continuity in $t$ is claimed for arbitrary individual cumulative halo volumes $h_l(t)$. For the log-convex price class of §2, (12) yields the static sharp upper even in the presence of arbitrarily separated scale labels.

**Local scalar checks.** Single-core Python evaluated $\Lambda_\mu(v)/v-\mu$ at $j=1.25,1.5,2,3$ as $0.144895549,0.266473523,0.464033520,0.754865754$. For (8), $\sqrt{2\pi}=2.506628275$ and $2\sqrt\mu=2.261384127$. At $s=0.1,0.01,0.001$ the explicit upper simplex cdf in (10) divided by $s^2$ was $0.1289503074,0.1254145750,0.1250416458$, approaching $1/8$. These scalar checks are not measurements of unrestricted $H_N$ or of primal optimality.
