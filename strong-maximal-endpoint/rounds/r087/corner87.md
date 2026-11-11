# Fixed physical corners have superexponentially small receiver efficiency on the near-$Q_N$ physical shell

Date: 2026-11-11. Strict theorem, independently audited.

## 1. Statement and consequences

Use the cap-two $F=F_{N,L(N)}$ and positive-volume shell $E=E_{N,L(N)}$ from r085 §3. Write $J=J_N$. The entire shell has soft demand $\rho=1$ at $\eta=N^{-2}$ and $|E|/J(F)=Q_N\exp(o(\log N))$. For finite-modular $0\le G\le F$ and an arbitrary physical corner $a\in\mathbb R^N$, $H_aG(x)$ is the true average on the box with opposite corners $a,x$, whichever signed orthant contains $x-a$. Coordinate-plane degeneracies are null. Let
\[
 K_E=\sup_{a,\,0\le G\le F,\,J(G)>0}\frac{|E\cap\{H_aG>1\}|}{J(G)}.
\]
Then for all sufficiently large $N$,
\[
 \boxed{\log K_E\le-\tfrac12N\log N+\tfrac12N\log\log N+O(N).} \tag{1}
\]
The bound is uniform over continuous physical corners, arbitrary measurable true sub-sources, and signed orientations. Only its upper direction is asserted. In particular the general, unrestricted one-corner norm $K_H$ from r083 and this restricted efficiency are different quantities.

For a finite or countably weighted true column cover, $B=\sum_iw_i J(G_i)$, $\Gamma=\sum_iw_i1_{\{H_{a_i}G_i>1\}}$, define $D_E=\int_E(1-\Gamma)_+$. Then
\[
 |E|-D_E=\int_E\min(1,\Gamma)\le K_E B. \tag{2}
\]
This is direct physical Lebesgue integration and Tonelli, including for fractional weights. For any sequence $B/J(F)\le e^{O(N)}$, the fraction of this near-$Q_N$ shell that the columns pay is at most $\exp[-\tfrac12N\log N+\tfrac12N\log\log N+O(N)]$. The remainder must be paid as actual physical deficit. Exact pointwise coverage, or even $D_E\le(1-\theta)|E|$ for fixed $\theta>0$, requires
\[
 \log\frac{B}{J(F)}\ge\log Q_N+\tfrac12N\log N-\tfrac12N\log\log N-O(N). \tag{3}
\]
These formulae concern the shell carrying actual lower-bound-order volume, unlike the tiny r084 slab. Paying the whole shell as residual still fits the contemplated $Q_N\exp(o(\log N))J(F)$ allowance, so this does not disprove a universal decomposition with that allowance, and gives no improved general $C_N$ upper bound.

## 2. Lift the existing uniform logarithmic scalar action deficit to actual rectangles

On $E$, $0\le t=-\sum_i\log|x_i|-L(N)\le T_N=\sqrt{\log N}$. The canonical profile $h_N$ is nondecreasing, zero on negative arguments, capped at two, and satisfies $\log j_N=-\log Q_N+o(\log N)$ for $j_N=\int_0^\infty e^{-s}\Phi_N(h_N(s))ds$. Consider any actual rectangular success from $x$, even for $G\le F$. Normalize its sides to independent multipliers $U_i$ uniform between $1$ and $d_i=a_i/x_i$ and put $W_i=(-\log|U_i|)_+$. Removing the cube cutoff and increasing the radial argument give $F(x_iU_i:i\le N)\le h_N(t+\sum W_i)$. Write $r_*=e^{-\mu}=\mu-1$ and $m(r)=1-r\log r/(1+r)$, with $0\log0=0$.

If $d_i=-r_i$ for $0\le r_i\le1$, $W_i$ has precisely the r022 increment law $Z_{r_i}$. In every other case $W_i$ is stochastically dominated by $Z_0\sim\mathrm{Exp}(1)$. For a short inward interval with $0<d_i<1$ the probability $\Pr(W_i>s)=(e^{-s}-d_i)_+/(1-d_i)\le e^{-s}$ for $s>0$. For $d_i>1$ it vanishes; for $d_i=-r$ with $r>1$ it is $2e^{-s}/(1+r)\le e^{-s}$. These also give the stochastic order at $s=0$ by right continuity. Replace every such irregular coordinate by an independent $Z_0$, keeping the already valid coordinates unchanged. Product quantile coupling and monotonicity of $h_N$ show that
\[
 \mathbb E h_N(t+\textstyle\sum_iZ_{\widetilde r_i})>1
\]
for the resulting deterministic proportion vector $\widetilde r\in[0,1]^N$ attached to this receiver and its one actual rectangle. No common physical input is replaced in the later charging estimate. The domination here is solely a necessary action test on an individual successful rectangle.

The shifted scalar profile's modular is at most $e^{t}j_N$. As $t\le T_N=o(\log N)$ and $\log Q_N=\mu N-d_\mu+\alpha\log N+o(\log N)$, it satisfies the *fixed* r022 test $J_{\rm scalar}\le e^{-s_N}$ at $s_N=\mu N-d_\mu+(\alpha-1)\log N$ for all large $N$. r022 report §§2–6, especially (13), proves uniformly for every feasible deterministic vector and nondecreasing profile at each such fixed constant that
\[
 D(\widetilde r)=\sum_i[\mu-m(\widetilde r_i)]\le C\log N. \tag{4}
\]
This application needs the internal deficit bound of r022, not merely its final formula for $\log Q_N$. In particular the cruder $O(\sqrt{N\log N})$ mean estimate from r086 can be sharpened for this shell. Each physically irregular coordinate replaced by zero incurs the fixed loss $\mu-1$ in (4).

Fix $e=1/50$ and put $r_-=r_*-e$, $r_+=r_*+e$. The elementary root bracket $0.27<r_*<0.28$ places $[r_-,r_+]$ inside $[1/4,3/10]$. For instance $t e^{1+t}$ is strictly increasing and its values at the two rational endpoints bracket one by the exponential series. Call physical sides with $d_i=-r_i$, $r_i\in[r_-,r_+]$ regular and all other sides exceptional. The derivatives satisfy $m^{\prime}(r)=(-\log r-1-r)/(1+r)^2$ and $m^{\prime\prime}(r)=(r-r^{-1}+2\log r)/(1+r)^3$. By $\log r\le r-1$, on $[1/4,3/10]$ we have
\[
 -m^{\prime\prime}(r)\ge\frac{r^{-1}-3r+2}{(1+r)^3}\ge\frac{133000}{65910}>2.
\]
Hence $\mu-m(r)\ge(r-r_*)^2$ on the regular interval. The derivative sign shows that every proportion outside this interval has gap at least $e^2$, and the stochastic replacement of an irregular side also loses $\mu-1=r_*>e^2$. Thus every true success, with exceptional set $I$, $J_0=I^c$, satisfies the two simultaneous constraints
\[
 |I|\le D/e^2\le C_e\log N,\qquad \sum_{j\in J_0}(r_j-r_*)^2\le D\le C\log N. \tag{5}
\]
The required budget includes *all* original irregular coordinates even when their true multiplier has an atom after positive-part truncation. Equation (5) is a necessary condition only and requires no regularity of the sub-source $G$.

## 3. Fixed physical corner: outer rectangle charge on the exceptional fibers

Fix $a,G,I$ and $n=N-k$. If any $a_j=0$ with $j\in J_0$, this success class is empty. Otherwise flip each regular coordinate so that $a_j=-A_j<0$; its receiver must then satisfy $x_j=A_j/r_j>0$, with
\[
 \ell_j=A_j/r_+,\quad u_j=A_j/r_-,\quad x_j\in[\ell_j,u_j],\quad O_j=[-A_j,u_j].
\]
The fixed box $O=\prod_{j\in J_0} O_j$ contains every regular-source box. With this harmless unit-Jacobian reflection, the same true source can be used throughout. The ratio $|O|/\prod_{j\in J_0}(A_j+x_j)$ is at most $R_e^n$ with
\[
 R_e=\frac{r_+(1+r_-)}{r_-(1+r_+)},\qquad
 P_e=\frac{u_j-\ell_j}{|O_j|}=\frac{2e}{r_+(1+r_-)}.\tag{6}
\]
Let $H_{a_I}^{(I)}$ average only in exceptional coordinates from the fixed corner $a_I$, and put $T(x_I)=|O|^{-1}\int_O H_{a_I}^{(I)}G(y_{J_0},x_I)dy_{J_0}$. Positivity and true product averaging imply $H_aG(x)\le R_e^n T(x_I)$ for all regular receiver coordinates in their cells (a.e. suffices). The fixed signed Hardy tensor has $L^2$ norm at most $2^k$: split in exceptional variables at $a_I$ into disjoint orthants, translate and reflect each onto the positive orthant, use the half-line Hardy $L^2$ bound two in each coordinate, and sum the squared integrals over orthants. There is no factor $2^k$ from summing the separate orthants. Thus, by Jensen first in the $O$ probability average and then integration in $x_I$,
\[
 \int T(x_I)^2dx_I\le \frac{4^k}{|O|}\int_{O\times\mathbb R^k}G^2\le\frac{2\,4^k}{|O|}J_N(G).\tag{7}
\]
The second inequality uses $0\le G\le2$ and $G\le\Phi_N(G)$. With $k=0$, use the one-point measure space for $x_I$. All statements concern physical Lebesgue measure before and after reflections.

Let $X_{J_0}\subset\prod[\ell_j,u_j]$ be *any* measurable set containing the regular projections of this success class. The implication above, Markov applied to $T^2$ at threshold $R_e^{-2n}$, and Fubini give
\[
 |E\cap\{H_aG>1,\ \mathrm{exceptional\ set}=I\}|\le
 2\,4^k R_e^{2n}\frac{|X_{J_0}|}{|O|}J_N(G).\tag{8}
\]
No estimate on an atomic receiver space is used. If $X_{J_0}$ is merely the whole cell, this gives the proposed $2\,4^k(P_eR_e^2)^n$ formula exactly. The quadratic constraint in (5) gives a sharper choice of $X_{J_0}$ at the single fixed $e=1/50$.

## 4. Euclidean ball volume gives the quantitative estimate

Use as $X_{J_0}$ the image under $x_j=A_j/r_j$ of the intersection of the parameter cube $[r_-,r_+]^n$ with $\sum_j(r_j-r_*)^2\le C_e\log N$ from (5). The change of variables has Jacobian $\prod A_j/r_j^2\le r_-^{-2n}\prod A_j$, while $|O|=(1+1/r_-)^n\prod A_j$. Dividing gives an upper bound $[r_-(1+r_-)]^{-n}$ times the parameter volume. Hence, for a fixed constant $C_0$ independent of $I,a,G,N$,
\[
 \frac{|X_{J_0}|}{|O|}\le C_0^n\Big(\frac{\log N}{n}\Big)^{n/2}.\tag{9}
\]
Indeed the Euclidean $n$-ball of radius $s$ has volume at most $(\sqrt{2\pi e_{\rm Euler}}s/\sqrt n)^n$, by integration of a Gaussian majorant with optimized scale. Combining (8)–(9), summing $\binom Nk$ possible exceptional subsets for $0\le k\le C_e\log N$, and using $n=N-O(\log N)$ proves (1). More explicitly, the combinatorial sum has logarithm $O((\log N)^2)$ and replacing $n\log(n/\log N)$ by $N\log(N/\log N)$ costs the same order. All fixed geometric factors contribute only $O(N)$. No numerical dimension onset is claimed.

Equations (2) and (3) now follow by weighted Lebesgue integration and the shell calibration, independently of overlap of the column successes or overlap of their source pieces.


## 5. A finite conditional bound with explicit constants

If every successful action on the specified shell has scalar deficit at most a given $D>0$, (5) permits $K=\lceil D/e^2\rceil$ exceptions. Let $C_{\rm geom}=R_e^2/[r_-(1+r_-)]$, choose a rational upper bound $\bar C\ge C_{\rm geom}$, set $n_0=N-K$ and $A=121D/7$; here $121/7>2\pi e_{\rm Euler}$. If $0<K<N/2$ and $n_0\ge\bar C^2A$, then (8), the Gaussian ball bound and a binomial count give
\[
 K_E\le 2\exp[N h_2(K/N)]4^K\bar C^{n_0}(A/n_0)^{n_0/2}. \tag{10}
\]
Indeed $n\mapsto n\log\bar C+\frac n2\log(A/n)$ decreases for $n\ge n_0$ under the displayed gate, so every $k\le K$ uses the same upper factor at $n_0$. Also $\sum_{k\le K}\binom Nk\le e^{Nh_2(K/N)}$: use $(q+1-q)^N=1$ at $q=K/N<1/2$ and note that $q^k(1-q)^{N-k}\ge q^K(1-q)^{N-K}$ when $k\le K$. These steps are finite inequalities.

The certificates test (10) at explicitly assumed $D\le\gamma\log N$. The unknown uniform constant and dimension onset of r022 (4) are not instantiated by a numerical choice of $\gamma$. Thus these conditional numbers do not evaluate the actual finite-dimensional $K_E$.

## 6. Verification and scope

Local Python rational arithmetic and independent Wolfram15 verified the geometric identities, curvature floor and directed constants, followed by 80/160-digit evaluation of (10). Of 48 conditional inputs, 35 satisfy its gates; the other 13 are explicitly excluded rather than counted as validated bounds. The all-$N$ result (1) comes from the proof and the stated r022 uniform theorem, not a fitted numerical slope. See [certification87.md](certification87.md), [validation87.json](validation87.json), and [independent audit](audit87.md). The unrestricted fixed-corner constant remains unchanged; it is the efficiency on the specified receiver shell that is small.
