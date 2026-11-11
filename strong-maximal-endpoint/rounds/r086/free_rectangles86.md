# Exact directional source cost for arbitrary physical rectangles on the benchmark shell

Date: 2026-11-11. Strict proof; independently audited. All logarithms are natural.

This note extends r085's canonical-witness entropy theorem. The charge is **one true source piece for each relative direction**, at its actual Orlicz modular. The pieces can vary arbitrarily in space, overlap, and be smaller in height than their common source. After choosing these pieces, the receiver may optimize over all signs, every positive side length, and every rectangle extending from it in that sign. The result still concerns the canonical cap-two source and the specified physical shell, not the general endpoint constant or r083's signed-fixed-physical-corner column LP.

## 1. Scope and result

Use the $F=F_{N,L(N)}$, $E=E_{N,L(N)}$ of r085 §3. Thus $0\le F\le2$, $J_N(F)>0$, $|E|/J_N(F)=Q_N\exp(o(\log N))=\exp(\mu N+o(N))$. On $E$ the coordinates are nonzero a.e. and $0\le t(x)=-\sum\log|x_i|-L(N)\le\sqrt{\log N}$. Globally $F(y)=1_{[-1,1]^N}(y)h_N(-\sum\log|y_i|-L(N))$, with $h_N$ nondecreasing, zero on $(-\infty,0)$, capped at two, and $j_N=\int_0^\infty e^{-s}\Phi_N(h_N(s))ds=\exp(-\mu N+o(N))$. Write $r_*=e^{-\mu}=\mu-1$, $p=r_*/(1+r_*)$, $h_2(p)=\mu-1+\log\mu$; the logarithm base is $e$. The actual r085 witness used for the matching construction has ratios $r_N\to r_*$ and strict margin $N^{-2}$.

For $0<\theta\le1$, let $b_N^{\rm free}(\theta)$ be the infimum of $\sum_\sigma J_N(G_\sigma)/J_N(F)$ over measurable $0\le G_\sigma\le F$ such that $|E\cap\bigcup_\sigma\{M^\sigma G_\sigma>1\}|\ge\theta|E|$. Here $M^\sigma$ has *every* positive side length at every receiver, with the coordinate sign fixed by $\sigma$ as in r085.

**Theorem.** For every fixed $0<\theta\le1$,
\[
 \lim_{N\to\infty}N^{-1}\log b_N^{\rm free}(\theta)=h_2(p). \tag{1}
\]
At any budget $\sum_\sigma J_N(G_\sigma)/J_N(F)\le e^{(\lambda+o(1))N}$ for fixed $\lambda<h_2(p)$ the relative volume of this true success union in $E$ is $O(e^{-c_\lambda N})$. Source-assignment cost $(Q_N/K_{H,N})N^{o(1)}J_N(F)$ has exponent $\mu-1<h_2(p)$ and therefore cannot use these free direction operators to recover a positive fraction of this physical shell.

## 2. Every strictly successful true witness has few exceptional sides

Fix $x\in E$ off the coordinate planes and an arbitrary true corner rectangle at $x$. Uniform sampling has the form $Y_i=x_iU_i$, with independent $U_i$ uniform in intervals with endpoints $1$ and some $d_i\in\mathbb R\setminus\{1\}$. The normalization allows the other endpoint to lie before zero, short of zero, or arbitrarily far beyond either end of the canonical range, whichever sign $x_i$ has. Let $W_i=(-\log|U_i|)_+$ and $m_i=\mathbb EW_i$. The cube cutoff and monotonicity give $F(Y)\le h_N(t(x)+\sum W_i)$. The scalar cost gives the simple physical envelope $h_N(s)\le e^s j_N$ for $s\ge0$, since $\Phi_N(u)\ge u$ and $h_N$ is nondecreasing on $[s,\infty)$.

* If $d_i=-r$, $0\le r\le1$, then $m_i=m(r)=1-r\log r/(1+r)$ (with $0\log0=0$). The unique maximum is $m(r_*)=\mu$: the derivative on $(0,1)$ is $(-\log r-1-r)/(1+r)^2$ and changes sign once.
* If $0\le d_i<1$, $U_i$ is uniform in $[d_i,1]$ and stochastically dominates a uniform on $[0,1]$, so $m_i\le1$. If $d_i>1$, $m_i=0$.
* If $d_i=-r$ with $r>1$, the mean of the positive part is $2/(1+r)\le1$, by integration over $[-1,1]$.

In **all** cases $\Pr(W_i>s)\le2e^{-s}$ for $s>0$. For $0\le d_i<1$ this follows by the preceding stochastic order, for $d_i>1$ directly, and for $d_i<0$ from the length bound $2e^{-s}/(1+|d_i|)$. This is a tail bound; very short non-crossing intervals need not satisfy an analogous density bound. For $0<b<1$, integration by parts using the tail gives $\mathbb E[W_i^2e^{bW_i}]\le4/(1-b)^3$. Taylor's inequality $e^u\le1+u+u^2e^u/2$ for $u\ge0$ gives
\[
 \log\mathbb E e^{bW_i}\le b m_i+2b^2/(1-b)^3. \tag{2}
\]
Independence holds for each actual fixed rectangle; no random feedback is substituted.

For every fixed $v>0$, all sufficiently large $N$ admit **no** strictly successful $F$-rectangle from $E$ with $\sum m_i\le(\mu-v)N$. Indeed choose small fixed $b>0$ in (2); the chance $\sum W_i>(\mu-v/2)N$ is $O(e^{-c_v N})$. On the complementary event the $h_N$ envelope and $t=o(N)$ give $F(Y)\le e^{-vN/2+o(N)}$. On the exceptional event use the cap two. This proves the assertion uniformly in $x$ and every side choice. In particular it applies to every witness for $G_\sigma\le F$.

Fix $0<a<\min(r_*,1-r_*)/2$. Call a side regular if $d_i=-r_i$, $|r_i-r_*|\le a$, and exceptional otherwise. Every exceptional side loses a fixed positive amount from $\mu$ in the mean $m_i$, by the cases above, compactness of $[0,1]$, and uniqueness of the maximum. Thus the number of exceptional sides is bounded by a deterministic $k_N(a)=o(N)$ uniformly on *all* successes from $E$. Indeed a subsequence with a positive fraction would contradict the preceding uniform assertion at a fixed $v>0$. Only after fixing $a$ is $N$ sent to infinity.

There is also a quantitative uniform count that retains the actual scalar cost. Set $q=-\log j_N-T_N-\log4$, $T_N=\sqrt{\log N}$. On $\{\sum W_i\le q\}$ the envelope gives $F(Y)\le1/4$. Thus $1<\mathbb E F(Y)\le 1/4+2\Pr(\sum W_i>q)$ forces $\Pr(\sum W_i>q)>3/8$. Applying (2) and independence yields, for every $0<b<1$,
\[
 \sum_i m_i> q-\frac{2Nb}{(1-b)^3}-\frac{\log(8/3)}{b}.
\]
For $N>1$ put $b=N^{-1/2}$ and set
\[
 B_N=\mu N+\log j_N+T_N+\log4+
 \frac{2\sqrt N}{(1-N^{-1/2})^3}+\sqrt N\log(8/3).
\]
Every such witness satisfies $0\le\mu N-\sum_i m_i<B_N$. The monotonicity of $m$ on either side of its unique maximum gives the explicit positive gap $g_a=\min\{\mu-1,\mu-m(r_*-a),\mu-m(r_*+a)\}$. Hence its exceptional count satisfies $|I|\le (B_N)_+/g_a$. For the stated r085 profiles $\log j_N=-\log Q_N+o(\log N)$, so the existing calibration and $d_\mu=O(\sqrt{N\log N})$ give $B_N=O(\sqrt{N\log N})$. This improves $o(N)$ to $O_a(\sqrt{N\log N})$ and is uniform over all receivers and actual witnesses for that source. The finite-dimensional formula retains $j_N$; the asymptotic profile calibration supplies no certified numerical onset.

## 3. Control all exceptional side lengths without enumerating their scales

For each index set $I$ of size at most $k_N(a)$ and absolute sign label $d\in\{\pm1\}^I$, put $J=I^c$ and $r^+=r_*+a$. Let $P_I$ be the composition, in increasing order of indices, of the noncentered one-dimensional maximal operators in $I$ (identity if empty). Then $P_I G\le2$ when $G\le2$, and
\[
 \sum_\sigma\int (P_IG_\sigma)^2\le 16^{|I|}\sum_\sigma\int G_\sigma^2\le 2\cdot16^{|I|}\sum_\sigma J_N(G_\sigma).\tag{3}
\]
For an elementary bound, the one-sided global rearrangement by the Hardy operator in r085 §1 and Hardy's $L^2$ constant two imply $\|M_i^\pm\|_{2\to2}\le2$ by slicing (or by the same global rearrangement). For completeness, Hardy's bound follows for nonnegative compact $f$ by writing $A(t)=\int_0^t f$ and integrating $\int(A/t)^2=2\int f(A/t)$; boundary terms vanish by $A(t)^2/t\le\int_0^t f^2$ at zero and compact support at infinity. Cauchy–Schwarz gives $\|Hf\|_2\le2\|f\|_2$, and monotone truncation extends this to general nonnegative $L^2$ input. The maximum of the two one-sided coordinate operators has norm at most four. This bound iterates on $L^2$; no global $L^1$ estimate on $P_I G$ is claimed.

For each $x$, define $\sigma(x;I,d)$ to have coordinates $-\operatorname{sign}x_j$ for $j\in J$ and $d_i$ for $i\in I$. Form the enlarged mean
\[
 A_{I,d}(x)=\mathbb E_U (P_IG_{\sigma(x;I,d)})(x_I,(x_j U_j)_{j\in J}),\qquad U_j\stackrel{\mathrm{iid}}\sim\operatorname{Unif}[-r^+,1].\tag{4}
\]
For a successful rectangle with exceptional set $I$ and assigned absolute exceptional signs $d$, average first in its exceptional coordinates and dominate this average by $P_I$, pointwise on the remaining fibers. Regular rectangles are contained in the enlarged intervals in (4); their volume ratios cost at most $e^{ND_a}$, $D_a=\log[(1+r_*+a)/(1+r_*-a)]$. Therefore its success implies $A_{I,d}(x)>e^{-ND_a}$. This is a pointwise bound on the genuine rectangle supremum, with the usual a.e. measurable coordinate maximal representatives. For a fully simultaneous interpretation, choose rational positive lengths for the original strict witnesses; all such averages form a countable family, and arbitrary positive lengths are recovered by continuity of locally integrable box integrals at positive lengths. The exceptional-coordinate domination holds for a.e. $x_I$ and a.e. values of the other integration variables by Fubini; linear changes in nonzero coordinates preserve null sets. The total number of pairs $(I,d)$ is at most $\sum_{i\le k_N(a)}\binom Ni2^i=\exp(o(N))$ for each fixed $a$.

## 4. Typical sign patterns pay for the regular coordinates

There exist $a_0>0,c>0$ depending on the fixed gap $h_2(p)-\lambda$, such that for every *fixed* $a\in(0,a_0]$ and all sufficiently large $N$ (the lower dimension cutoff may depend on $a$),
\[
 |E|^{-1}\int_E A_{I,d}^2\le O(e^{-cN}) \tag{5}
\]
for all the pairs above. The constant $c$ can be chosen independently of the final value of $a$. This order of choices avoids a circular requirement between action tolerance and the Chernoff rate.

Here are details. Choose $a_0<\min(r_*,1-r_*)/2$. Let $n=|J|=N-o(N)$, uniform over the index sets for this fixed $a$. For $S_J=-\sum_{j\in J}\log|U_j|$ use (2) and $m(r^+)\le\mu$ to bound $\Pr(S_J>(\mu+\delta)N)\le O(e^{-c_1N})$ uniformly for all $a\in(0,a_0]$. For the signs $B_j=\operatorname{sign}U_j$, the probability of a negative entry is $p_a=r^+/(1+r^+)\in[p,p_{a_0}]\Subset(0,1/2)$. Hence $h_2(p_a)\ge h_2(p)$. Choose a small fixed frequency tolerance using uniform boundedness of $|\log(p_a/(1-p_a))|$ on this compact range. By Hoeffding there is a typical set $\mathcal T_{n,a}$ with complement probability $O(e^{-c_2N})$ such that each pattern there has probability at most $e^{-(h_2(p)-\delta)n}\le e^{-(h_2(p)-2\delta)N}$ for all sufficiently large $N$. The exponents $c_1,c_2>0$ are uniform in $a$; one may use $n\ge N/2$ for the tail rate.

Take $X$ uniform on $E$ and set $Y=(X_I,X_jU_j:j\in J)$. By Jensen, $\mathbb E A_{I,d}(X)^2$ is bounded above by the expectation of $(P_IG_{\sigma(X;I,d)}(Y))^2$. The two exceptional sampling events in the previous paragraph contribute at most four times their combined probability, by the cap. Otherwise the determinant inverse of $X\mapsto Y$ at fixed $U$ is $e^{S_J}\le e^{(\mu+\delta)N}$. For fixed $Y=y$ and $d$, as $B$ varies the indices in $J$ of $\sigma(X;I,d)$ vary bijectively over all signs there. Integrating in $U$ and dropping the image-shell indicator, as in r085 §4, bounds the remaining expectation by
\[
 \frac{e^{(\mu-h_2(p)+3\delta)N}}{|E|}\sum_\sigma\int(P_IG_\sigma)^2
 \le e^{(\lambda-h_2(p)+3\delta+o(1))N}.\tag{6}
\]
Here (3), $|I|=o(N)$, $|E|/J_N(F)=e^{\mu N+o(N)}$ and the input cost hypothesis are all necessary. Choose $3\delta<h_2(p)-\lambda$ to get (5) with an exponent $c>0$ uniform in small fixed $a$.

Finally fix $a>0$ sufficiently small that $2D_a<c/2$. On each relevant $I,d$, Markov's inequality and (5) give shell fraction $\le O(e^{(2D_a-c)N})$ for $A_{I,d}>e^{-ND_a}$. Summing the $\exp(o(N))$ pairs from §3 proves the claimed exponential suppression of *all* genuine one-sided maximal success sets simultaneously. The matching upper in (1) is immediate from r085's pass-through Hamming masks retaining every original canonical rectangle: allowing alternative rectangles only increases coverage. As in r085 this holds for each fixed $\theta>0$ without any attained infimum.

For an explicit asymptotic choice at the target rate $\lambda=\mu-1$, take sign-frequency tolerance $\epsilon=1/100$, $\delta=1/20$ and moment parameter $b=1/100$. Uniformly for the regular sampling ratios, $|\log(p_a/(1-p_a))|\le\mu$, so each typical pattern has probability at most $e^{-n(h_2(p)-\mu\epsilon)}$. Hoeffding and $n\ge N/2$ give exceptional sign probability at most $2e^{-N/10000}$. The moment bound gives the other exceptional probability at most $e^{-N/4000}$, since $b\delta-2b^2/(1-b)^3>1/4000$. The remaining budget exponent has margin $\log\mu-\delta-\mu\epsilon>9/50$, before the $o(N)$ costs. Thus the squared-mean exponent in (5) may be taken as $c=1/20000$ for all sufficiently large $N$. Directed constants certify $2D_{10^{-6}}<1/40000$; after the $e^{o(N)}$ enumeration the final success fraction is $O(e^{-N/100000})$ at this target budget rate. This only fixes valid asymptotic constants; no numerical dimension onset is supplied.

## 5. Interpretation and computation

This removes the restriction to preserving the old rectangle: changing all side lengths, choosing a new orientation, or using arbitrary exceptional-coordinate scales and translations of the far endpoint cannot recover this shell at the target per-direction source cost. The common fixed endpoint in $M^\sigma$ is the receiver itself. The theorem does **not** prove failure of r083's separate architecture, where a single fixed physical corner and its one source piece can serve different relative directions at different receivers without a separate per-direction charge. The shell's physical score is at the permitted $Q_N\exp(o(\log N))$ residual scale, so directly paying that shell also remains possible. This is not a contradiction to the target general endpoint upper bound.

The numerical inputs needed for the limiting entropy were checked in r085. The analytic argument here uses only explicit moment and cap estimates, Hardy's $L^2$ bound, binary concentration, and the established r085 physical shell. No numerical proof of an asymptotic onset or fixed-dimension sharp size is claimed.

For the final archive, local Python/SciPy and Wolfram independently checked 72 endpoint/moment combinations, including short, long, crossing and outward edges. Four rational action tolerances have directed mean-gap certificates. Details and the distinction between exact constants and numerical quadrature are in [certification86.md](certification86.md); [audit86.md](audit86.md) records independent proof review. The short MATLAB sanity-check session was closed. No finite-dimensional onset, new endpoint bound, or success claim for the separate fixed-physical-corner architecture is inferred from these calculations.
