# Directional receiver deletion and the source cost of retaining benchmark witnesses

For $N\ge1,A>1$ let $K=K_{H,N}^{[A]}=Q_{\Gamma,N}^{[A]}$ be r028/r083's exact one-corner norm at cap $A$. Write $J=J_N$ and $\Phi=\Phi_N$. Let $Q_N$ denote the established scalar benchmark, the supremum over all deterministic ratio vectors in $[0,1]^N$ and all finite caps (r028); it is distinct from $K$, and $Q_N\ge K$ for any fixed cap by the $r=0$ subfamily. Thresholds are strict; all unqualified volumes are ordinary physical Lebesgue volumes. The general $C_N$ coefficient is still open.

## 1. Exact strong norm when all side lengths are free but their directions are fixed

For $\sigma\in\{\pm1\}^N$, $M^\sigma F(x)$ is the supremum over side lengths $l_i>0$ of the mean on the rectangle with endpoints $x_i$ and $x_i+\sigma_i l_i$. Different receivers can choose entirely different lengths. The genuine noncentered operator is $M_S=\max_\sigma M^\sigma$: cutting any rectangle through its receiver along all coordinate hyperplanes expresses its mean as a convex combination of the means of its positive-volume corner pieces. Conversely every corner rectangle is allowed.

**Theorem.** For every fixed $\sigma$, the optimal cap-$A$ constant in $|\{M^\sigma F>1\}|\le K J(F)$ is exactly $K$.

For a single coordinate and fixed sign put $P_i^\sigma$ for its one-sided maximal operator, and $Hg(t)=t^{-1}\int_0^t g(u)du$. For a bounded compactly supported one-dimensional fiber $f$, at $s>0$ the rising-sun identity is $s|E_s|=\int_{E_s} f$, where $E_s=\{P^\sigma f>s\}$. For backward averages, let $P(t)=\int_0^t(f-s)$ and $m(t)=\inf_{v\le t}P(v)$. Then $E_s=\{P>m\}$, and on each component $(a,b)$ the running minimum stays fixed and $P(a)=P(b)$. Finite endpoints follow from $P'=-s$ off a compact interval and $P(t)\to+\infty$ at $-\infty$, $P(t)\to-\infty$ at $+\infty$. For fixed $s>0$ the open set $E_s$ lies in a bounded interval (the compact source interval enlarged by $\|f\|_1/s$ in the backward direction of propagation). Its disjoint open components are at most countable, and $f-s$ is integrable on that bounded interval, so summing the component equalities is legitimate. Forward averages follow by reflection. This is precisely the one-sided component of r010 §1. The additional positive-tail term in r011 is for the union of *both* directions within a coordinate.

Integrating the fiber identity in the other coordinates on $\mathbb R^N$ gives $sD_s=\int_{E_s}f\le\int_0^{D_s}f^*(v)dv$, where $D_s=|E_s|$ and the star is the single global rearrangement. If $s'>s$ and $0<D_{s'}<\infty$, then $Hf^*(D_{s'})\ge s'>s$. Hence $D_{s'}\le|\{Hf^*>s\}|$, since $Hf^*$ decreases. Letting $s'\downarrow s$ gives $(P_i^\sigma f)^*\le Hf^*$ a.e. The zero-volume case is immediate.

For $0\le F\le A1_{[-b,b]^N}$ sequentially use each coordinate once. Intermediate functions are bounded by $A$ times the product of the maximal averages of $1_{[-b,b]}$ in the used coordinates and the original indicators in the unused coordinates. The used-coordinate factors lie in $[0,1]$ and tend to zero at infinity. Every positive level set is therefore bounded, and every next unused-coordinate fiber has compact support. These facts justify the preceding fiber integrals, Fubini, and finite rearrangement integrals even when an intermediate function is not globally $L^1$. The genuine rectangular supremum is bounded above by the composition of the $P_i^{\sigma_i}$ by successively pulling suprema into integrals of nonnegative functions. Equality with the composition is not needed. Its global rearrangement is bounded by $H^NF^*$, by induction and positivity. Rational side lengths suffice for measurability; positive-length box integrals depend continuously on their endpoints for locally integrable input.

At each $u>0$ with $H^NF^*(u)>1$, set $h(s)=F^*(ue^{-s})$. Independent uniform factors in the Hardy averages give $H^NF^*(u)=\mathbb E_{\Gamma(N,1)}h$, and $\int_0^\infty e^{-s}\Phi(h(s))ds=u^{-1}\int_0^u\Phi(F^*(v))dv\le J(F)/u$. The scalar optimum gives $u\le KJ(F)$. Thus the strictly successful volume for the original $M^\sigma$ has that upper bound. Box truncation $F1_{[-b,b]^N}\uparrow F$ extends it to all finite modular cap-$A$ inputs by monotone convergence of every average, supremum, strict level set and modular.

For the converse use the all-backward direction. The r028 §4 cap-preserving radial lower construction at $r=0$ is even in every coordinate. Its fixed-anchor witness box $[0,x]$ stays inside the orthant of its receiver. Restrict the source to the positive orthant and use the positive-orthant subset of the reference receiver set. The certified receiver volume and modular each divide by $2^N$. There each witness is all-backward. Sending r028's spatial parameter to infinity recovers each strictly admissible scalar ratio and hence $K$. Reflect coordinates to obtain any fixed sign. Alternatively the signed fixed-corner operator splits into the separate physical orthants with disjoint sources, so one-orthant near-optimizers exist. The original ingredients and finite-dimensional Gamma sharpness already appear in r010/r011/r028/r083. The direction theorem is a corollary with the physical corner allowed to move for every receiver; no external priority claim is made.

## 2. A physical residual one can delete at the target price

Let $m_F(x)=\#\{\sigma:M^\sigma F(x)>1\}$. Then
\[
 \int m_F\,dx\le2^N KJ(F),\quad |\{m_F\ge q\}|\le 2^NKJ(F)/q. \tag{1}
\]
For any fixed family $U$ of directions, its union of success sets has measure at most $|U|KJ(F)$. Thus $q=\max(1,\lceil 2^NK/Q_N\rceil)$ permits direct payment of the full physical high-multiplicity residual at $Q_N J(F)$. A global list of at most $Q_N/K$ directions likewise has physical residual at most $Q_N J(F)$. The same bounds hold integrating $0\le\rho\le1$ over these sets.

More flexibly, at every source point let $\Lambda(y)$ be a measurable list of available directions (each of the finitely many membership indicators measurable). Define $F_\sigma=F1_{\{\sigma\in\Lambda(y)\}}$. Then
\[
 \left|\bigcup_{\sigma}\{M^\sigma F_\sigma>1\}\right|
 \le K\sum_\sigma J(F_\sigma)=K\int \#\Lambda(y)\Phi(F(y))dy . \tag{2}
\]
This charges the *same* source point once per allowed relative direction, where the allocation is source-dependent, while rectangle lengths and the resulting corners may vary freely with the receiver. All these sets are directly paid residuals, not countable fixed-physical-corner source columns.

## 3. The near-lower-bound physical receiver body has very few successful directions at each point

Let $T_N=\sqrt{\log N}$, $\eta=N^{-2}$. Use the r069 shifted deterministic reference $r_N=e^{-c_N}$, $c_N=\mu-\frac12\sqrt\mu\sqrt{\log N/N}$; the fixed cap-$3/2$ scalar calibration gives bounded nondecreasing $h_N:[0,\infty)\to[0,2]$ with
\[
 {\mathbb E}h_N(S_{N,r_N})>1+\eta,\qquad
 \log j_N=-\log Q_N+o(\log N),\quad j_N=\int_0^\infty e^{-s}\Phi_N(h_N(s))ds. \tag{3}
\]
Indeed choose a strictly admissible profile arbitrarily close to the fixed-cap benchmark (or normalize the exact mean-one optimizer); multiply by $1+2\eta$. The cap stays below two and the modular distortion is at most $D_N(1+2\eta)=e^{O(1/N)}$. All fixed caps $>1$ have the r069 shifted-reference logarithmic calibration; equality of this cap-specific norm with unrestricted $Q_N$ at finite $N$ is not assumed. Extend $h_N$ by zero on the negative axis. With $N$ fixed first, take the physical cap-two source
\[
 F_{N,L}(y)=h_N(-\sum_i\log|y_i|-L)\,1_{[-1,1]^N}(y),\quad
 E_{N,L}=\{x\in[-1,1]^N:0\le t=-\textstyle\sum_i\log|x_i|-L\le T_N\}.
\]
Values on the coordinate planes may be assigned arbitrarily. On all of $E=E_{N,L}$ the rectangle $R_x=\prod_i[{-r_Nx_i},x_i]$ (endpoints unordered) has average at least the mean in (3), since $t\ge0$ and $h_N$ is nondecreasing. Hence $\rho=1$. The shifted Gamma volume formula from r028 follows by product log coordinates (convolution of $N$ unit exponentials). For $V_L=2^N e^{-L}L^{N-1}/(N-1)!$ it gives
\[
 |E|/V_L\to1-e^{-T_N},\qquad J_N(F_{N,L})/V_L\to j_N \quad(L\to\infty\text{ at fixed }N).
\]
For the second limit one integrates $e^{-t}(1+t/L)^{N-1}\Phi_N(h_N(t))$ on $[0,\infty)$; boundedness at fixed $N$ gives domination for $L\ge1$. Choose finite $L=L(N)$ sufficiently large for both limits to relative error $1/N$. This gives $|E|/J(F)=Q_N\exp(o(\log N))$. Only the explicitly counted shell is being valued.

A direction pointing away from zero in coordinate $i$ can be replaced, for this coordinatewise absolute-decreasing original source, by evaluation at $x_i$, with no decrease in average. Fix an integer $k<N$ and drop any $k$ outward coordinates of a successful direction. Relax all other directions. For each remaining coordinate the unrestricted optimal search can be restricted to intervals $[-r_i x_i,x_i]$, $0\le r_i\le1$. Indeed first split at the receiver to use a side interval; the outward piece has mean at most the receiver value. The inward piece ending short of zero can be enlarged to zero by adding higher values. If its far endpoint lies past $-x_i$, trimming that outer part discards only values no larger than any retained value. An inward interval to zero has average at least the receiver value. These replacements are simultaneous across the product by coordinatewise monotonicity and positivity. The remaining sample increments are independent $Z_{r_i}=-\log|U_i|$ for $U_i\sim\operatorname{Unif}[-r_i,1]$. The rectangle side ratios $r_i$ are fixed for this receiver before sampling. There is no feedback. Thus success implies existence of a deterministic vector of $n=N-k$ ratios for which $\mathbb E h_N(t+\sum_{i=1}^n Z_{r_i})>1$ (or obtain a strict feasible rectangle by the definition of a strict supremum). Scalar admissibility yields the necessary inequality
\[
 1\le Q_n\int_0^\infty e^{-s}\Phi_n(h_N(t+s))ds\le Q_n e^t j_N. \tag{4}
\]
Here $Q_n$ is the project's supremum over *vectors* of fixed ratios and over caps as in r028; $\Phi_n\le\Phi_N$ is essential. For $k=\lfloor\epsilon\log N\rfloor$, $\epsilon>0$ fixed, the calibrated expansion of $\log Q$ gives $\log Q_{N-k}-\log Q_N=-\mu k+o(\log N)$. To check the $d$ term, extend the Gaussian root to real $v$ by $\bar\Phi(z)=\sqrt{\mu/v}\,z$: implicit differentiation gives $0<z'(v)\le z(v)/(2v)$, hence $d'(v)=O(\sqrt{\log v/v})$. Since $t\le T_N=o(\log N)$, (4) is impossible for all sufficiently large $N$, uniformly on the shell. Each true successful direction therefore differs from the receiver's all-inward sign in at most $k_N=o(\log N)$ coordinates. This bounds every $m_F(x)$ there by $\sum_{j\le k_N}\binom Nj=\exp(o((\log N)^2))$.

Each sign orthant contributes exactly $2^{-N}|E|$ by symmetry. Consequently any fixed list of $s$ directions can capture at most $s2^{-N}\exp(o((\log N)^2))|E|$. For $s\le (Q_N/K)N^{o(1)}=\exp((\mu-1)N+o(N))$ this fraction is exponentially small. The high-multiplicity deletion in (1) also misses this shell, because $2^N K/Q_N=\exp((\log2+1-\mu)N+o(N))$. Thus a physical body of the full established lower-bound order can consist of locally low multiplicity spread across all spatial orientations.

## 4. Exact exponential cost of retaining canonical rectangle witnesses

For an arbitrary family of direction-dependent physical source pieces, integrate their contributions over actual receiver positions before summing the sign patterns. Write
\[
 r_*=e^{-\mu}=\mu-1,\quad p_N=r_N/(1+r_N)\longrightarrow p=r_*/(1+r_*),\quad h_2(p)=-p\log p-(1-p)\log(1-p)=\mu-1+\log\mu. \tag{5}
\]
For each fixed $0<\theta\le1$, define $b_N(\theta)$ as the infimum of $\sum_\sigma J_N(G_\sigma)/J_N(F)$ over arbitrary measurable pieces $0\le G_\sigma\le F$ such that a subset of the specified physical shell $E$ of volume at least $\theta|E|$ has original witness mean $\fint_{R_x}G_{-\operatorname{sign}x}>1$. Here $F,E,R_x$ are exactly those of §3, at finite $L(N)$.

**Theorem.** $\displaystyle \lim_{N\to\infty}N^{-1}\log b_N(\theta)=h_2(p)$ for every $0<\theta\le1$. More quantitatively, at any budget $\sum_\sigma J_N(G_\sigma)/J_N(F)\le e^{(\lambda+o(1))N}$ with fixed $\lambda<h_2(p)$, the shell average of the retained rectangle means and hence the successfully retained physical fraction are $O(e^{-c_\lambda N})$. The calibrated budget $Q_N/K$ has rate $\mu-1$, a strict exponent gap $\log\mu$. Local computation gives $p\simeq0.217811705720$, $h_2(p)\simeq0.524124324658$, $\log\mu\simeq0.245659781897$.

**Lower bound.** Draw $X$ uniformly in physical volume on $E$, independently of iid $U_i\sim\operatorname{Unif}[-r_N,1]$, and set $Y_i=X_iU_i$, $S=-\sum_i\log|U_i|$. Then $Y$ is conditionally uniform on $R_X$. Put $B=\operatorname{sign}U\in\{\pm1\}^N$. The count $K_B$ of negative entries is binomial $(N,p_N)$. For fixed $\delta>0$, there is a set $\mathcal T$ of typical sign patterns with $\Pr(B\notin\mathcal T)\le e^{-c_\delta N}$, and each $b\in\mathcal T$ satisfies $\Pr(B=b)\le e^{-(h_2(p)-\delta)N}$. Indeed $N^{-1}\log\Pr(B=b)=(K_b/N)\log p_N+(1-K_b/N)\log(1-p_N)$; continuity and Hoeffding on a small fixed interval about $p_N$ give both claims. Likewise $\Pr(S>(\mu+\delta)N)\le e^{-c'_\delta N}$, since $Z=-\log|U_i|$ has mean $m(r_N)=1-r_N\log r_N/(1+r_N)\le\mu$ and density bounded by $2e^{-z}1_{z>0}$. In fact Taylor expansion gives $\log\mathbb E e^{bZ}\le b\mu+2b^2/(1-b)^3$ for $0<b<1$: use $e^t\le1+t+t^2e^t/2$ for $t\ge0$, integrate the density bound and use $\log(1+v)\le v$. The same tail estimate holds uniformly for independent nonidentically distributed $r_i\in[0,1]$, centering each by $m(r_i)$. Correlation between $B$ and $S$ is allowed throughout.

Conditional on any fixed $U$ with $S\le(\mu+\delta)N$, change variables $y_i=x_iU_i$ in the shell integral. The determinant inverse is $e^S$, and $-\operatorname{sign}x=-\operatorname{sign}y\,B$ coordinatewise. Dropping the image-shell indicator gives the upper bound $|E|^{-1}e^{(\mu+\delta)N}\int G_{-\operatorname{sign}y\,B}(y)dy$. Integrate on the typical event. For each fixed $y$ off its null coordinate planes the map $b\mapsto-\operatorname{sign}y\,b$ is a bijection on sign labels. Hence
\[
 \begin{split}
 \mathbb E G_{-\operatorname{sign}X}(Y)
 &\le 2\Pr(B\notin\mathcal T)+2\Pr(S>(\mu+\delta)N)\\
 &\quad+\frac{e^{(\mu+\delta)N}}{|E|}\sum_{b\in\mathcal T}\Pr(B=b)\int G_{-\operatorname{sign}y\,b}(y)dy\\
 &\le O(e^{-cN})+\frac{e^{(\mu-h_2(p)+2\delta)N}}{|E|}\sum_\sigma J_N(G_\sigma).
 \end{split} \tag{6}
\]
Only positivity, exact Jacobians and $\int G_\sigma\le J_N(G_\sigma)$ were used. The two exceptional probabilities are multiplied by the uniform cap two. Since $|E|/J_N(F)=e^{\mu N+o(N)}$, choose $2\delta<h_2(p)-\lambda$ to prove the exponentially small fraction claim by Markov. Applying this along a putative violating subsequence gives the lower bound on the infimum even if no minimum is attained.

**Matching construction.** Put $k_N=\lceil Np_N+2\sqrt{N\log N}\rceil$. At each $y$ retain the full physical height $F(y)$ in just the directions within Hamming distance $k_N$ of $-\operatorname{sign}y$, setting other $G_\sigma(y)$ to zero. Exactly $\sum_{j=0}^{k_N}\binom Nj$ directions receive each nonzero source point, so the actual sum of modular costs divided by $J_N(F)$ equals this integer. Its exponential rate is $h_2(p)$ by Stirling and $p<1/2$. For each original $x\in E$, the relative-direction label $-\operatorname{sign}x$ is discarded precisely when $K_B>k_N$ under the canonical rectangle sampling (up to source-zero points); the pointwise lost height is bounded by two. By Hoeffding the loss of average is at most $2N^{-8}$, uniformly over $x$. The original average exceeds $1+N^{-2}$, hence each of those rectangles still succeeds at threshold one for all sufficiently large $N$. This yields the upper bound simultaneously for $\theta=1$ and every smaller $\theta$. The entropy simplification in (5) follows from $\log(1+r_*)=\log\mu$ and $-p\log r_*=p\mu=\mu-1$.

The exact exponent is for retaining the displayed rectangles with one reusable physical source piece per *relative direction*. The physical residual allowance still can pay this entire shell at lower-bound order. Success on other rectangles after direction-dependent truncation remains a separate question. No universal residual or general endpoint upper bound at the target scale is claimed here.
