# The full physical indicator halo on weighted hyperbolic crosses, and finite common-side tests

Let $N$, $L=\log N$, $\tau=d/N$, $q_0=\mu N-d$, $Q_N$, $c_Q$ have the meanings in r092. The physical critical density is $p_N=1-\tau-BL/N$ for a fixed real $B$; in the endpoint upper application take the positive constant selected in r092. Let $s_{p,w}$ be the supremum over **deterministic** vectors $r_i\in[0,1]$ of the upper-hit quantile at probability $p$ of $\sum_i w_i Z_{r_i}$, with the independent increments of r011/r022, and $s_p=s_{p,\mathbf1}$. Normalize $0\le w_i\le1$, $\max_i w_i=1$. Define
\[
 E_{T,w}=\{x\in[-1,1]^N:\sum_i w_i(-\log|x_i|)>T\},\qquad T\ge0 .
\]
Values on coordinate hyperplanes are immaterial. Each set has finite positive measure. Translations, coordinate permutations and nonsingular diagonal axis scalings preserve all physical volume ratios below.

## 1. Exact extremum for the whole physical halo, with freely receiver-dependent sides

For every fixed $B$, all sufficiently large integers $N$ satisfy
\[
 \boxed{\quad\sup_{T\ge0,\,0\le w_i\le1,\,\max w_i=1}
 \frac{|\{M_S1_{E_{T,w}}>p_N\}|}{|E_{T,w}|}=e^{s_{p_N}}.\quad}         \tag{1}
\]
Consequently the logarithm of this restricted *full* physical halo constant is $q_0+c_QL+o(L)$ (r092 (8)). All rectangles containing the receiver are available, all receivers inside and outside the original cube are counted, and their choices of sides may depend arbitrarily on the receiver. The equal-weight source reaches the indicated supremum as $T\to\infty$ for each qualifying fixed $N$. This is about the entire near-$Q_N$ physical mass of a class of actual sources. The $gL$ gain of the r029 adaptive-hit comparison cannot be realized on this class by choosing different rectangle sides for different receivers.

**Physical reduction.** By sign symmetry take $x_i>0$. Every corner interval may be taken to be $[-r_ix_i,x_i]$, $0\le r_i\le1$. Indeed split an arbitrary rectangle through $x$ at the coordinate planes of $x$, retaining a corner piece of at least its average. For any one coordinate the integral of $1_E$ over the remaining coordinate intervals is even and nonincreasing with coordinate magnitude. An interval pointing outward can be replaced by $[0,x_i]$, one starting between $0$ and $x_i$ can be extended to $0$, and one going left beyond $-x_i$ can be trimmed to $-x_i$, without reducing the average. Applying this successively proves the reduction, including when some $w_i=0$. For any such *one fixed rectangle*, a uniform sampled point has independent coordinates, with $-\log(|y_i|/x_i)$ distributed as $Z_{r_i}$. Thus inside the positive unit cube, writing $u_i=-\log x_i$, the entire halo is precisely, up to a null level surface,
\[
 \sum_i w_i u_i > T-s_{p,w} .                                        \tag{2}
\]
This uses a deterministic vector for a particular receiver, with the supremum over *all* such vectors taken at that receiver. The last-block choices cannot change after seeing the sampled first coordinates. All vector laws have continuous positive sum density on $(0,\infty)$ since one weight is one; hence for $s<s_{p,w}$ some vector strictly exceeds the threshold, while $s>s_{p,w}$ excludes every vector. Inside the cube the averaging intervals stay in the cube. Outside, ignoring the source cube still makes (2) a necessary condition. There is also an independent necessary condition: $M_S1_E(x)\le M_S1_{[-1,1]^N}(x)=\prod_i\min(1,2/(1+|x_i|))$. Consequently $|x_i|<2/p-1=e^b$ whenever outside in that coordinate, where $b=\log(2/p-1)>0$.

Let $G_w(t)=\mathbb P(\sum_i w_i Y_i>t)$ for independent unit exponentials $Y_i$, with $G_w(t)=1$ for $t\le0$. Then $|E|=2^N G_w(T)$ and the enlargement $|x_i|\le e^b$, (2), bounds the whole halo by $2^N e^{Nb}G_w(T-s_{p,w}+b\sum w_i)$. Since a weight equals one, conditioning on the other variables gives $G_w(t-a)\le e^a G_w(t)$ for every $a\ge0$ and every real $t$. This follows from the exact exponential tail, including negative arguments. Writing $D=\sum_i(1-w_i)$, the finite-$N$ geometric bound is therefore
\[
 \log(|\{M_S1_E>p\}|/|E|)\le\max\{Nb,\ s_{p,w}+bD\}.              \tag{3}
\]
For equal weights whenever $s_p\ge Nb$ the right side is $s_p$. The interior (2) gives the lower ratio $G_{\bf1}(T-s_p)/G_{\bf1}(T)$, tending to $e^{s_p}$ as $T\to\infty$ by the integer Gamma tail $G_{\bf1}(t)=e^{-t}\sum_{j=0}^{N-1}t^j/j!$ for $t\ge0$. No assumptions on exceptional outside receivers, fixed anchor counting, or angular share of the original halo are needed.

**Paying the box enlargement for unequal weights.** At $p=p_N$, $b=O(\tau)\to0$, and $s_p=\mu N-d+O(L)$ by r092 (8). If $D>N/4$, for any fixed independent vector the weighted mean is at most $3\mu N/4$, the variance is at most $C N$, and Chebyshev bounds its upper-hit quantile by $3\mu N/4+O(\sqrt N)$ since $p_N\to1$. Thus both entries on the right of (3) are below $s_p$ eventually, uniformly in this range.

For $D\le N/4$, consider a fixed vector and let $q(w)$ be its weighted sum quantile at $p_N$. Uniformly in this region, and on every segment from $w$ to $\mathbf1$ obtained by increasing a coordinate,
\[
                     \partial_{w_i}q(w)\ge c>0.                       \tag{4}
\]
Here is the analytic justification. Write $S_w=\sum w_i Z_{r_i}$ with mean $M$ and variance $V\asymp N$. At least $N/2$ of its coefficients are at least $1/2$; deleting one term leaves at least $N/2-1$. Centered moments through any fixed order are uniformly bounded per term, and each such large coefficient has the exponential characteristic factor $(1+w_i^2 t^2)^{-1/2}$ in modulus. The r022 §5 Fourier estimates therefore extend with absolute constants: the absolute one-term cdf Edgeworth error is $O(1/N)$; the absolute unstandardized density error from the Gaussian with matching $M,V$ is $O(1/N)$ for both the sum and each deleted sum. For details, on the standardized frequency interval $|t|\le\epsilon\sqrt N$, the log characteristic function is $-t^2/2+K_3(it)^3/(6V^{3/2})+O(t^4/N)$, with a uniform Gaussian bound on its modulus. Subtracting the Gaussian and its cubic term and integrating the difference divided by $|t|$ gives $O(1/N)$ as in r022 (10); the density inversion after subtracting the Gaussian alone has standardized $O(N^{-1/2})$ error. Outside this interval at least $N/2-1$ factors give a product bounded by $(1+c t^2/N)^{-(N/2-1)/2}$, with exponentially small integrals. Division by $\sqrt V$ gives the asserted density estimate. The same proof works if some of the other coefficients vanish.

The quantile expansion argument of r092 (7) gives $q(w)=M-\sqrt V z+O_B(L)$. At that point, and at $q(w)-w_i u$ with $1\le u\le2$ for the law with coordinate $i$ deleted, the matching unstandardized Gaussian densities are uniformly comparable to $L/N$. Their positions differ from the corresponding $-z$ by $O_B(L/\sqrt N+z/N)$ in standardized units, and $zL/\sqrt N=o(1)$. The absolute $O(1/N)$ error is smaller. Writing $f_w,f_{-i}$ for sum and deleted densities, differentiating the exact convolution identity for the cdf at its quantile gives
\[
 \partial_{w_i}q(w)=\frac{\mathbb E[Z_{r_i} f_{-i}(q(w)-w_i Z_{r_i})]}{f_w(q(w))}.
\]
The continuous uniformly bounded deleted density, integrable first moment, and positive denominator justify differentiation under the integral and the implicit rule; endpoints use one-sided limits. Each $Z_{r_i}$ has probability at least $(e^{-1}-e^{-2})/2$ in $[1,2]$ from its unshifted exponential component. Restricting the numerator to this interval proves (4), with $N$ sufficiently large depending on $B$. Integrating through successive coordinates gives $q(w)\le q(\mathbf1)-cD$. Taking suprema *after* the coordinatewise estimate yields $s_{p,w}+bD\le s_p$ eventually, since $b<c$. This proves (1).

## 2. A physical integral upper which charges conflicts of last-block sides

Split the coordinates into $m+n=N$. Write $E_z=\{v:(z,v)\in E\}$, $U_z(v)=M_{S,n}1_{E_z}(v)$. For any first-block rectangle $A\ni u$, at a receiver $(u,v)$ the honest best with these first sides is
\[
 J(A,v)=\sup_{B\ni v}m_{A\times B}1_E.
\]
Define a sequence of upper bounds which progressively force last-block sides to be shared. For $k\ge1$, take independent $Z_1,\dots,Z_k$ uniform on the very same $A$, and put
\[
 P_k(A,v)=\mathbb E\left[\sup_{B\ni v}\frac1k\sum_{j=1}^k m_B1_{E_{Z_j}}\right],
 \qquad \mathcal M_kE(u,v)=\sup_{A\ni u} P_k(A,v).                  \tag{5}
\]
For every measurable input set (choose a Borel representative),
\[
 M_{S,N}1_E\le\mathcal M_{k+1}E\le\mathcal M_kE\le\mathcal M_1 E
 =M_{S,m}(M_{S,n}1_E).                                                \tag{6}
\]
The first inequality keeps a fixed honest $B$ before taking the expectation. Convexity of the supremum and averaging the $k+1$ possible leave-one-out $k$-term sums give monotonicity. In particular the two-row upper is the following explicit subtraction from the iterated maximum:
\[
 M_{S,N}1_E(u,v)\le\sup_{A\ni u}\left[m_A U_z(v)-\frac12 m_{A^2}
 \{U_z(v)+U_{z'}(v)-M_{S,n}(1_{E_z}+1_{E_{z'}})(v)\}\right].       \tag{7}
\]
The braces are nonnegative and bound from below $(U_z+U_{z'}-1-M_{S,n}1_{E_z\cap E_{z'}})_+$ by the binary intersection identity. Formula (7) permits the final rectangle to depend on the receiver and the sampled *pair*, while requiring the two integrals within a pair to use one last-block rectangle. This brackets the true geometry between a common last choice for all $z$ and a new choice at every $z$. All the inner suprema have countable rational rectangle representatives, up to harmless outward endpoint approximation; thus the $A^k$ integrands are jointly measurable. For $A$ one may also use outward rational approximation since changing the uniform law on $A$ by $o(1)$ in total variation changes its $k$-fold law by $o(1)$ for fixed $k$. No trace of a general higher-dimensional maximal function on a diagonal is assumed bounded here.

**A finite-grid precision bound.** Suppose $E$ is a binary step source on a finite axis-parallel grid with bounded support and the last-block coordinate grids have $m_i$ bounded cells including the bounding-box faces. Put $M=2^n\prod_{i=1}^n(m_i+1)$. At every last-block receiver off the grid faces and for all first-block rectangles,
\[
 D_{\rm Ber}(P_k(A,v)\|J(A,v))\le\frac{\log M}{k},\qquad
 D_{\rm Ber}(\mathcal M_k E(u,v)\|M_{S,N}1_E(u,v))\le\frac{\log M}{k}.              \tag{8}
\]
Here $D_{\rm Ber}(a\|b)=a\log(a/b)+(1-a)\log[(1-a)/(1-b)]$, using the usual endpoint limits. For each list of sampled rows the last-block response is that of a nonnegative step function on the *same* source grid. The r079 finite-corner argument, adjusting one far endpoint at a time, attains its maximum among the at most $M$ boxes with receiver $v$ as near corner and a grid vertex as far corner. The same applies to the exact $A$-averaged profile. Each candidate has independent row responses $X_{j,B}\in[0,1]$ of mean at most $J=J(A,v)$. For every $t>0$, convexity gives $\mathbb E e^{tX_{j,B}}\le 1-J+Je^t$. Jensen, the finite maximum bound, and independence across sampled rows for fixed $B$ give
\[
 tkP_k\le\log\mathbb E\exp(tk\max_B k^{-1}\textstyle\sum_jX_{j,B})
 \le\log M+k\log(1-J+Je^t).
\]
Taking the optimal $t\ge0$ yields (8) because $P_k\ge J$. The endpoints $J=0,1$ follow directly. For the second inequality take a sequence of outer rectangles approaching the $\mathcal M_k$ supremum. If the limit exceeds $M_{S,N}1_E$, replace $J(A,v)$ by the larger $M_{S,N}1_E$ for all sufficiently large terms, which decreases $D_{\rm Ber}$ there, then pass to the limit. Equality of the two outer suprema makes the assertion immediate.

For instance if $0<\delta<p<1$ and $kD_{\rm Ber}(p\|p-\delta)>\log M$, then $\{\mathcal M_k E>p\}\subset\{M_{S,N}1_E>p-\delta\}$ off the last grid faces. At critical $p=1-\tau+O(L/N)$ and $\delta=o(\tau)$ the denominator is $D_{\rm Ber}(p\|p-\delta)\sim\delta^2/(2\tau)$. This gives a definite accuracy budget for replacing fiberwise adaptation by a finite number of shared observations; e.g. $\delta=cL^2/N$ makes the sufficient $k$ of order $N^{3/2}(\log M)/L^{7/2}$. The finite grid size depends on the source and no dimension-uniform estimate of the weak norm of $\mathcal M_k$ is asserted. An actual upper estimate on the receiver integral of (7), or on a finer $\mathcal M_k$, at this density could be used for $H_N$; reducing it back to $\mathcal M_1$ removes precisely its nonnegative penalty.

## 3. Binary rows where the full physical maximum and the relaxed maxima separate

Here is a complete two-dimensional calculation using *all* physical rectangles, in contrast to a construction with a prescribed outer side. In horizontal coordinates let $I_1=[0,1]$, $I_2=[1-a,1]$ with $0<a<1$. In the vertical interval $[0,1]$ alternate the two indicator rows in the two equal halves of periods of length $\varepsilon=1/j$. Off that interval use empty rows. The resulting actual set has $|E_\varepsilon|=(1+a)/2$. Test at $(x,y)=(1+\delta,y)$ with $\delta>0$, $y>1$. Put $U_1=1/(1+\delta)$, $U_2=a/(a+\delta)$, and $d_0=1/y$. For fixed $k$ and $0<\delta<a$, with $\delta<a/(k-1)$ when $k\ge2$, the limits of the *actual operators* are
\[
 \begin{split}
 M_{S,2}1_{E_\varepsilon}(x,y)&\longrightarrow d_0 U_2,\\
 \mathcal M_k E_\varepsilon(x,y)&\longrightarrow d_0 U_2+
       (U_1-U_2)\,(d_0/2)(1-d_0/2)^{k-1}.                          \tag{9}
 \end{split}
\]
The two-coordinate iterated maximum allowing the horizontal side to change after observing the sampled vertical coordinate is the case $k=1$.

For the proof, fix any actual vertical interval $A$ through $y$. If it meets the nonempty rows its length is at least $y-1$. The categorical row probabilities on it differ entrywise from $(d/2,d/2,1-d)$ by $O(\varepsilon/(y-1))$, where $d=|A\cap[0,1]|/|A|\le d_0$. This estimate is uniform in $A$, by cancellation over complete periods. For $k$ independent samples the error in a bounded expectation is $O(k\varepsilon/(y-1))$ by product coupling. The deterministic row average has the same single-row error uniformly over horizontal intervals. At $x>1$ only horizontal far endpoints $0$ or $1-a$ need to be tested on a nonempty mix, by the step argument. If the conditional proportion of long rows among occupied sampled rows is $t$, their shared density equals the larger of $(a+(1-a)t)/(1+\delta)$ and $U_2$. The long side wins precisely when $t>U_2$. In the deterministic half-half mixture the short side wins for $\delta<a$. In $k$ samples containing any short row, $t\le1-1/k<U_2$ under the displayed condition; if there are only long occupied rows the long side wins. For the multinomial counts $N_1,N_2$ the gain over always taking the short side thus equals $(U_1-U_2)\mathbb E[(N_1/k)1_{N_2=0}]=(U_1-U_2)(d/2)(1-d/2)^{k-1}$. The expected response is nondecreasing in $d$ by coupling each extra occupied row in place of a zero row. As $A=[0,y]$ attains $d_0$, taking the full outer supremum proves (9).

For example take $a=1/2$, $p=1-\eta$, $\delta/\eta\in[0.32,0.34]$ and $(y-1)/\eta\in[0.40,0.42]$. Uniform expansion of (9) for $k=1$ gives the true deficit $(2\delta+y-1)+O(\eta^2)\ge1.04\eta+O(\eta^2)$, and the iterated deficit $(1.5\delta+y-1)+O(\eta^2)\le0.93\eta+O(\eta^2)$. Uniformity over this whole patch follows from $y-1\ge0.40\eta$. For small $\eta$ and sufficiently small $\varepsilon=o(\eta^2)$ the relaxed operator succeeds strictly on the entire positive-area patch and the actual full operator is strictly below $p$ there. In particular one may use $\eta=1-p_N$; the two-dimensional construction detects genuine fiber adaptation at the designated density without changing the source between receivers. Its area is $0.0004\eta^2$; this local example alone does not yield a positive fraction of an $N$-dimensional extremal halo or allow a per-step penalty to be multiplied across $N$ coordinates. Formula (1), by contrast, accounts for a whole high-dimensional benchmark halo.

## 4. The grid-entropy term is sharp for physical binary rows

For fixed integers $k,M\ge1$, $0<q<1$, and every tolerance $\zeta>0$, there is a bounded two-dimensional binary step set and a positive-area receiver patch on which, simultaneously,
\[
 |M_{S,2}1_E-q|<\zeta,\qquad
 \left|\mathcal M_k E-\mathbb E\max_{1\le\ell\le M}B_\ell/k\right|<\zeta,
 \quad B_\ell\ \text{independent }\mathrm{Bin}(k,q).                \tag{10}
\]
The last-block horizontal grid uses $M$ bounded cells, regardless of how many first-block rows are needed.

To construct it, let $X=(X_1,\ldots,X_M)$ have independent Bernoulli$(q)$ coordinates. On each fine vertical period in $[0,1]$ allocate its $2^M$ consecutive subintervals the exact probabilities of these bit vectors. Take a period length $1/j$; outside $[0,1]$ use the all-zero vector. Choose $0=T_0<T_1<\cdots<T_M$ with $(T_{\ell-1}+1)/(T_\ell+1)\le\kappa$ for a freely small $\kappa>0$. In row $X$, occupy the horizontal interval $(T_{\ell-1},T_\ell]$ precisely when $X_\ell=1$. Test horizontal receivers $x\in[-1,0]$. For any list of $k$ rows, every horizontal interval through $x$ averages a convex combination of the empirical bit means $m_\ell=k^{-1}\sum_{j=1}^k X_{j,\ell}$ and zero. Conversely, $[x,T_\ell]$ captures at least $(1-\kappa)m_\ell$ because its last shell has that relative length. Thus the supremum is between $(1-\kappa)\max m_\ell$ and $\max m_\ell$. Exactly the same statement holds for the mean bit vector of a whole vertical interval.

Take vertical receivers in $[1+\gamma,1+2\gamma]$. On any vertical interval $A$ through such a receiver, the row-vector law is within $2/(j\gamma)$ in total variation of the mixture of a true independent-coordinate bit vector with weight $d=|A\cap[0,1]|/|A|$ and a zero vector otherwise. This bound depends only on the two incomplete periods; their probability masses bound total variation for *all* bit patterns. The product error on $k$ independent rows is at most $2k/(j\gamma)$. In the ideal mixture, zeroing a subset of rows reduces the empirical maximum. Hence the expectation is bounded above by the full independent bit-vector value in (10). For $A=[0,y]$, $1-d\le2\gamma$; coupling shows its ideal $k$-row maximum differs from the full law by at most $2k\gamma$. Similarly the deterministic mean-bit vector is the constant $dq$ in the ideal mixture, so the true full rectangle supremum is bounded between $(1-\kappa)(q/(1+2\gamma)-2/(j\gamma))$ and $q+2/(j\gamma)$. First choose $\kappa$ and $k\gamma$ sufficiently small relative to $\zeta$, then $j$ large. This proves (10) uniformly on the stated patch. It permits irrational row heights/lengths and is an ordinary finite axis-parallel measurable step set.

For fixed $q$ and $0<r<\log(1/q)$, take $\log M/k\to r$ and $k\to\infty$. Let $q_r\in(q,1)$ solve $D_{\rm Ber}(q_r\|q)=r$. A union bound and the binomial Chernoff estimate exclude any larger limiting maximum. At any $u\in(q,q_r)$ Stirling applied to the single count $\lceil ku\rceil$ gives a binomial probability at least $\exp[-kD_{\rm Ber}(u\|q)+o(k)]$; independence across $M$ columns makes the probability all miss this count tend to zero. Bounded convergence in probability gives $\mathbb E\max B_\ell/k\to q_r$. Taking $\zeta\to0$ in (10) shows the coefficient one in the grid entropy $\log M/k$ of (8) cannot be decreased uniformly for this physical class; the fixed-rational-grid corner overcount in (8) is only $O(1)$ in the logarithm for these one-dimensional last blocks. One may fix any near-one target $p$ first, choose $q<p$ and $r>D_{\rm Ber}(p\|q)$ below $\log(1/q)$, and take $k,M$ sufficiently large to obtain whole patches of strict relaxed success and honest failure. There is no dimension-only finite sample count making the local gap small uniformly over the source grids: for any given $k$, taking $M\gg q^{-k}$ makes a column of $k$ ones overwhelmingly likely. The construction can use very distant horizontal endpoints and very small receiver volume relative to source size, so it gives no lower bound of this type for the sharp global weak constant. A physical *volume* estimate could still avoid this entropy penalty on most of a large halo.

For $q=1/2$ and $(k,M)=(2,2),(3,3),(4,4)$, local Python exact rational evaluation of $k^{-1}\sum_{i=1}^k[1-\mathbb P(B_1<i)^M]$ gave $11/16$, $47/64$, $49063/65536$ respectively.

## 5. Scope of the remaining upper problem

The physical class in (1) includes variable weights, thickness and arbitrary receiver-specific choices of all rectangle sides. Generic indicator sets need not have hyperplane boundaries in logarithmic coordinates. For them, (7) still involves the actual pair of fibers, their common last-block side and the outside receiver integral. Replacing its pair term by separate fiber suprema returns the relaxation which admits spurious successes already in (9). A bound at the scale required for the general endpoint must estimate how much of a near-critical **physical receiver volume** survives those common-side constraints, or obtain an equivalent geometric covering estimate; (1), (7) and (8) leave that part open. Thus no new coefficient for unrestricted $H_N(p_N)$ or $C_N$ is asserted here.

## 5. Reproducible verification and status

The final [experiment93.py](../../code/experiment93.py) gives exact rational certificates for twenty-one sufficient Gamma premises and 140 explicitly specified Gamma-tail polynomial ratios. Independent [verify93.wl](../../code/verify93.wl) uses a separate rational logarithm bound for the premises. All arithmetic, actual exponential ratios and numerical cross-checks are replayed by [check93.py](../../code/check93.py); scope and resource limits are in [certification93.md](certification93.md) and [validation93.json](validation93.json). The finite premises concern the equal-weight geometric calculation, not a certified onset for the all-weight theorem or numerical values of the optimized quantiles. Extra full-CPU probes are numerical diagnostics only, separately recorded in [extra93](extra93/README.md). No general sharp indicator upper follows from the restricted physical extremum or the local grid-entropy obstruction.
