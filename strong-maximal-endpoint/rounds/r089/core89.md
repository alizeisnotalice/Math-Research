# A sharp fixed-coefficient height threshold for the general physical endpoint norm

Date: 2026-11-11. Strict asymptotic theorem; independently audited.

Put $L=\log N$, $\rho=\sqrt{L/N}$, $\mu=1+e^{-\mu}$, $\sigma=\sqrt\mu$, and $C_N(A)=\sup_{0\le F\le A}|\{M_SF>1\}|/J_N(F)$ with the original $\Phi_N(t)=t(1+\log^+t)^{N-1}$. Write $C_N$ for the unrestricted norm. Use standard normal density $\phi$, right tail $\bar\Phi$, $z=z_N>0$ solving $\bar\Phi(z)=\sigma z/\sqrt N$, $d=\sigma\sqrt N z$, $q_0=\mu N-d$. All results are asymptotic, and errors for fixed $K$ or $P$ may depend on those fixed parameters.

## 1. General reduction by physically clipping the excess

For every fixed $P>0$,
\[
 K\ge\sigma(1+\sqrt{2P+1})\quad\Longrightarrow\quad
 1\le C_N/C_N(1+K\rho)\le1+O_P(N^{-P}).\tag{1}
\]
At the exact boundary $K=2\sigma$ the ratio is $1+O(L^{-1/4})$. No optimality of the polynomial exponent in (1) is claimed. Below we prove that for *every* fixed $0<K<2\sigma$ the logarithm of this ratio grows faster than $\log N$. In particular $2\sqrt\mu$ is the sharp constant for preserving the endpoint through its $\log N$ scale among caps of the form $1+K\rho$, and at that very constant the ratio actually tends to one. This does not identify the still-unknown unconstrained $\log N$ coefficient.

First, merely tracing the physical low source in r073 §1 already gives (1) with the stricter $K>\sigma+\sqrt{2\mu(P+1)}$ by discarding the whole high source at threshold $N^{-(P+1)}$. The following argument uses r073's proved exponentially accurate reduction $C_N\le(1+O(Ne^{-N/100}))C_N(2)$, the genuine strong $L^p$ upper estimate r001, and a different cost for the clipped excess. In particular an estimate uniform in shrinking caps does not follow merely from r073's statement for each fixed cap.

Assume $F\le2$. Put $s=N^{-(P+1)}$, $a=(1-s)(1+K\rho)$, $f=F\wedge a$, $g=F-f$; then $1<a<2$ eventually. The low part at threshold $1-s$ costs at most $C_N(1+K\rho)(1-s)^{-N}J_N(f)$ by Orlicz dilation. Set $\Delta J=J_N(F)-J_N(f)$, so $J_N(f)+\Delta J$ is *exactly* the original modular, without assuming additivity on overlapping parts. On $[a,2]$, $(\log\Phi_N)'=[1+(N-1)/(1+\log t)]/t\ge cN$ for $c=1/[2(1+\log2)]$. Therefore
\[
 \|g\|_p^p\le \frac{(C p/N)^p}{\Phi_N(a)}\,\Delta J\quad(p\ge2).\tag{2}
\]
Indeed at any height $a+u$ the denominator difference is at least $\Phi_N(a)(e^{cNu}-1)$; $\sup_{v>0}v^p/(e^v-1)\le(C'p)^p$ by splitting $v\le1$ and $v>1$. For an explicit finite coefficient one may take $C=1/c=2(1+\log2)$ in (2): for $v\le1$, $v^p/(e^v-1)\le v^{p-1}\le1\le p^p$; for $v>1$, $e^v-1\ge e^v/2$ and $v^pe^{-v}\le p^pe^{-p}$, so $2e^{-p}\le1$ for $p\ge2$. This extra $p/N$ distinguishes throwing away an excess from throwing away the full original height.

The strong estimate supplies $|\{M_Sg>s\}|\le H_N\Delta J$, $H_N=\exp(N b_p)(Cp/(Ns))^p/\Phi_N(a)$, $b_p=\mu+\mu/(2p)+O(p^{-2})$. Take $p=b\sqrt{N/L}$, $b=\sigma/\sqrt{2P+1}$. For fixed $P,K$, $\log\Phi_N(a)=K\sqrt{NL}+O(L+Ns)$. With $W=\sqrt{N/L}$,
\[
 \log H_N\le\mu N+[-K+\mu/(2b)+(P+1/2)b]\sqrt{NL}-(b/2)W\log L+O(W+L).\tag{3}
\]
At equality in (1) the bracket equals $-\sigma$; larger $K$ only helps. The cap-two physical benchmark satisfies $\log C_N(2)\ge q_0-C_{\rm low}L$, and $z\le\sqrt L$ for all large $N$ by Mills, so (3) is strictly smaller than $\log C_N(2)$ eventually. Since $\{M_S F>1\}\subset\{M_S f>1-s\}\cup\{M_S g>s\}$, suprema and the exact two-part expense give
\[
 C_N(2)\le\max\{(1-s)^{-N}C_N(1+K\rho), H_N\}.
\]
The second term lies strictly below the left side, proving (1). At $K=2\sigma$ take $p=\sigma W$, $s=N^{-1}L^{-1/4}$. The bracket is still $-\sigma$, the secondary term in (3) becomes $-(\sigma/4)W\log L$, and dilation costs $1+O(L^{-1/4})$. This proves the assertion at equality. The general cap-two and $L^p$ results are applied before any source-dependent scalar argument.

**Quadratic exponential rate above the critical window.** If $\delta_N\to0$ and $\delta_N/\rho\to\infty$, then for every fixed $0<c<1/(2\mu)$,
\[
 C_N/C_N(1+\delta_N)\le1+O_c(Ne^{-cN\delta_N^2}).\tag{3a}
\]
The same bound holds for every sufficiently small fixed $\delta>0$ (depending on $c$). For clarity, take $s=e^{-cN\delta^2}$, $a=(1-s)(1+\delta)$, $p=b/\delta$, $b=\sigma/\sqrt{2c}$. In the varying case $N\delta^2\gg L$, so $s=o(\delta)$, $Ns=o(1)$, $p\ge2$, $Cp/N\le1$ for all sufficiently large $N$. The denominator gives $\log\Phi_N(a)=N\delta+O(N\delta^2+Ns+\delta)$. Discard the $Cp/N$ improvement in (2). Then
\[
 \log H_N\le\mu N+N\delta[\mu/(2b)+cb-1+O_c(\delta)+o(1)]
 =\mu N+N\delta[\sigma\sqrt{2c}-1+O_c(\delta)+o(1)].
\]
This is below the physical cap-two benchmark because $N\delta\gg\sqrt{NL}$. The fixed-small-$\delta$ version follows by choosing its upper endpoint sufficiently small to absorb the $O_c(\delta)$ remainder, keep $p$ inside the large-$p$ expansion, and arrange $c\delta^2<1/100$ so that r073’s cap-two error is smaller than the new dilation error. For larger fixed caps, monotonicity or r073 supplies an exponential reduction, without a claim that this particular near-one rate constant is optimal for all caps.

## 2. Below $2\sigma$: an upper estimate that holds for every physical source

The following reductions are used explicitly. r011 bounds every physical cap-$A$ score by the scalar feedback optimum with the same cap: at a normalized receiver one obtains nondecreasing $0\le h\le A$, $j(h)=\int e^{-t}\Phi_N(h(t))dt$, and a predictable policy with $\mathbb E_\pi h(S_N)>1$. Conversely each strictly feasible *fixed*-ratio nondecreasing scalar profile yields a true physical lower bound $1/j(h)$ at the same cap by r028. In the general upper direction r021 §6 (15) says that for every such policy and every real $u$,
\[
 F_\pi(u)\ge F_{G'}(u)-2N^{-4},\quad G'\sim \mathcal N(\mu N+C_H L,\mu(N+t_0)),\quad t_0=O(L).\tag{4}
\]
There is no independence assumption on the physical policy. For any bounded nondecreasing test $V$, this gives $\mathbb E_\pi V\le\mathbb E V(G')+2N^{-4}(\sup V-\inf V)$, by Stieltjes integration (or layer cake). The terminal and comparison Gaussian are continuous in law. Write $G\sim\mathcal N(\mu N,\mu N)$.

Consider $A=1+\epsilon$, $\epsilon=K\rho$, and $\sigma\le K<2\sigma$. Define the bounded slope test
\[
 H_s(x)=1_{\{x\ge s\}}+N^{-1}\min(N\epsilon,(x-s)_+),\qquad U(s)=\mathbb E H_s(G).
\]
If $j(h)\le e^{-s}$, monotonicity gives $\Phi_N(h(x))\le e^{x-s}$ at nonnegative $x$. For $x<s$ this bounds $h$ by $e^{x-s}$. For $0\le u=x-s\le N\epsilon$, the inverse modular estimate r022 (3) or r032 §4 yields $(h-1)_+\le u/N+C u^2/N^2$ uniformly, since $N\epsilon=O(\sqrt{NL})$. For larger $u$ the genuine cap takes over. The single upper test $H_s(x)+e^{x-s}1_{\{x<s\}}$ is nondecreasing, continuous at $s$, and bounded by $1+\epsilon$. One must apply (4) to this *whole* test, not to the isolated lower-side term. This proves
\[
 \mathbb E_\pi h(S_N)\le U(s)+O(L^2/N)\tag{5}
\]
whenever $|s-q_0|=O(\sqrt N/L)$. Indeed the uniform inverse error is $O(\epsilon^2)=O(L/N)$. Since $s<\mathbb E G'$, $\mathbb E[e^{G'-s}1_{G'<s}]\le f_{G'}(s)=O(L/N)$. The $G'$ versus $G$ probability shift at $s$ is $O(L^2/N)$: in this window their standardized negative arguments differ by $O(L/\sqrt N)$ and their densities there are $O(L/\sqrt N)$, by the root Mills relation. The ramp term is $1/N$-Lipschitz, so coupling by the same standard normal costs $O(L/N)$. The strict-level exclusion in (5) is uniform over all feasible profiles and policies.

The same Gaussian model has a matching true-source test at this accuracy. On nonnegative arguments use $h=H_s$ itself, with cap exactly $1+\epsilon$. For an iid sum $S$ of the $r_*=e^{-\mu}$ increments, its mean and variance are $\mu N$ and $\mu N$. The r022 absolute one-term Edgeworth expansion (10) bounds its cdf minus that of $G$ by $C(1+y^2)\phi(y)/\sqrt N + C/N$ with $y=(u-\mu N)/(\sigma\sqrt N)$. At $s$ this is $O(L^2/N)$ in the same window. Integrating that bound between $s$ and $s+N\epsilon$ costs $O(1)+O(\epsilon)$ before division by $N$, since $\int(1+y^2)\phi(y)dy=2$. Consequently
\[
 \mathbb E H_s(S)=U(s)+O(L^2/N),\qquad j(H_s)\le (1+N\epsilon)e^{-s}.\tag{6}
\]
For the modular bound use $\Phi_N(1+u/N)\le e^u$, integrate up to $N\epsilon$, then integrate the constant cap exactly. This handles the expensive cap height by delaying it in the physical lower profile. A constant-height indicator would instead pay it on its full plateau.

Let $f_0=\phi(z)/(\sigma\sqrt N)\sim L/N$, $w=(K/\sigma)\sqrt L-z$, $\Psi(w)=\phi(w)-w\bar\Phi(w)$, and define
\[
 D_K=\frac{\sigma\Psi(w)}{\sqrt N\, f_0}.\tag{7}
\]
The exact mean identity gives $U(q_0)=1+\mathbb E(q_0-G)_+/N-\sigma\Psi(w)/\sqrt N$, since $F_G(q_0)=d/N$; $\mathbb E(q_0-G)_+=O(1)$ by Mills. In particular the last term is the lost Gaussian positive tail imposed by the cap. The derivative is $U'(s)=-f_G(s)-[F_G(s+N\epsilon)-F_G(s)]/N$. For $|s-q_0|\le 2D_K$, $f_G(s)=(1+o(1))f_0$: $w\ge0$ eventually, $D_K=O(\sqrt N/L)$, and the logarithm of the density ratio is $O(zD_K/\sqrt N+D_K^2/N)=o(1)$. The second summand in the displayed first-derivative formula has magnitude at most $1/N=o(f_0)$.

If $K=\sigma$ then $w\to0$ and $D_K\sim\sigma\sqrt N/(\sqrt{2\pi}L)$. If $\sigma<K<2\sigma$, put $a=K/\sigma-1\in(0,1)$. Then $w=(a+o(1))\sqrt L$, $\Psi(w)\sim\phi(w)/w^2$, and $D_K=N^{(1-a^2)/2+o(1)}$. In either case $D_K\gg L$ and $f_0D_K\gg L^2/N$. Testing $s=q_0-(1\mp\delta)D_K$ for any fixed $0<\delta<1$ in (5) and (6) gives respectively a strict uniform exclusion and a strictly successful actual cap-respecting lower profile. Its logarithmic modular overhead $O(L)$ is negligible. Thus
\[
 \log C_N(1+K\rho)=q_0-(1+o(1))D_K,\qquad
 \log[C_N/C_N(1+K\rho)]=(1+o(1))D_K\quad(\sigma\le K<2\sigma).\tag{8}
\]
The second assertion uses the established general bracket $\log C_N=q_0+O(L)$. A smaller fixed $K$ is already ruled out by monotonicity and the $K=\sigma$ result; the next section gives its actual larger loss.

## 3. Smaller caps and the two requested subcritical scales

For any of the following $\epsilon$, let $z_\epsilon=\bar\Phi^{-1}(\epsilon)$, $q_\epsilon=\mu N-\sigma\sqrt N z_\epsilon$. The bounded-profile argument gives the general physical upper bound
\[
 \log C_N(1+\epsilon)\le q_\epsilon+O(L)\quad\left(\epsilon=N^{-1},\ (NL)^{-1/2},\ \text{or }K\rho\text{ for fixed }0<K<\sigma\right).\tag{9}
\]
To see this, $j(h)\le e^{-s}$ forces $h\le N^{-2}$ when $S\le s-2L$. Combining the cap with strict success gives $F_\pi(s-2L)<\epsilon/(1+\epsilon-N^{-2})\le\epsilon$. Test $s=q_\epsilon+(C_H+3)L$ in (4): its Gaussian standardized argument is at least $-z_\epsilon+L/(\sigma\sqrt{N+t_0})$ for large $N$, since $z_\epsilon>0$ and the larger variance helps. The normal cdf exceeds $\epsilon$ by at least $c\phi(z_\epsilon)L/\sqrt N\gg N^{-4}$ (all stated $z_\epsilon\asymp\sqrt L$, $\epsilon\ge N^{-1}$ eventually). This contradicts feasibility, proving (9).

The matching physical lower profile in this range can be even simpler. Put $\ell=\max(0,\log(\Phi_N(1+\epsilon)-1))$ and $h=1_{[s,\infty)}+\epsilon 1_{[s+\ell,\infty)}$. It has $j(h)\le2e^{-s}$; at the fixed iid $r_*$ sum, writing its cdf as $F_*$, strict feasibility is precisely $F_*(s)<\epsilon[1-F_*(s+\ell)]$. Also $\ell\le N\epsilon$ since $\Phi_N(1+\epsilon)\le e^{N\epsilon}$. We use r013/r014's fixed-iid relative moderate-deviation formula $F_*(\mu N-\sigma\sqrt N y)=\bar\Phi(y)(1+O(L^{3/2}/\sqrt N))$ for $0\le y\le 2\sqrt L$; this is only used in the true lower construction. For $\epsilon=N^{-1}$ and $(NL)^{-1/2}$ choose $s=q_\epsilon-C_1 L$. Normal Mills bounds show $F_*(s)\le\epsilon(1-cL^{3/2}/\sqrt N)$ by taking a sufficiently large fixed $C_1$: the normal logarithmic decrement from $C_1L$ is comparable to $C_1 z_\epsilon L/\sqrt N$. Also $F_*(s+\ell)\le C\epsilon$, since $z_\epsilon \ell/\sqrt N=O(1)$ for both cases. The chosen decrement dominates this last probability because $\epsilon=o(L^{3/2}/\sqrt N)$. Thus each profile succeeds strictly, and r028 gives the true lower bound at its actual cap. Together with (9),
\[
 \boxed{\quad\log C_N(1+\epsilon)=\mu N-\sigma\sqrt N z_\epsilon+O(L)\qquad [\epsilon=N^{-1}\ \hbox{or}\ (NL)^{-1/2}].\quad} \tag{10}
\]
In particular cap $1+1/N$ has first negative term $-\sqrt{2\mu N\log N}(1+o(1))$. For $\epsilon=(NL)^{-1/2}$ the leading negative term is the same as for the general norm, but the general-to-capped log loss satisfies
\[
 \log[C_N/C_N(1+(NL)^{-1/2})]\sim \sigma\sqrt{N/L}\,\log L. \tag{11}
\]
Indeed $\bar\Phi(z)/\epsilon=\sigma z\sqrt L$, whose logarithm is $\log L+O(1)$, while $z,z_\epsilon\sim\sqrt L$. Integrating the Gaussian hazard rate $\phi/\bar\Phi\sim x$ between these two quantiles proves the comparison, with $O(L)$ negligible here.

For fixed $0<K<\sigma$, $\epsilon=K\rho$, test the same two-level profile at $s=q_\epsilon-\delta\sqrt{N/L}$ for any fixed $\delta>0$. The relative moderate-deviation formula and the Gaussian hazard give $F_*(s)/\epsilon\to e^{-\delta/\sigma}<1$. Also $F_*(s+\ell)=o(1)$ since $s+\ell\le\mu N-(\sigma-K+o(1))\sqrt{NL}$; Chebyshev at the fixed reference suffices. Hence $h$ succeeds for every such fixed $\delta$ eventually. Letting $\delta\downarrow0$ after the limit, (9) yields $\log C_N(1+K\rho)=q_\epsilon+o(\sqrt{N/L})$. Now $\bar\Phi(z)/(K\rho)\to\sigma/K$, giving
\[
 \log[C_N/C_N(1+K\rho)]\sim\sigma\log(\sigma/K)\sqrt{N/L}. \tag{12}
\]
Equation (12) is not uniform as $K\to\sigma$; (8) resolves the boundary case. Equations (1), (8), and (12) concern all physical inputs and their actual Orlicz denominator. They do not reduce the general problem to indicator functions: the lower profiles already show why delayed heights can matter. No assertion about the ratio's exact polynomial rate above $2\sigma$ or about $N$-dependent transition windows in $K$ is made.


## 4. Independent numerical verification

The final [experiment89.py](../../code/experiment89.py) uses local mpmath to evaluate 35 Gaussian ramp roots, 28 small-cap quantiles and 49 explicit clipped-excess fees, at two precision targets with adaptive guard digits. Independent Wolfram15 evaluates the same 112 parameter cases by different root methods and the unexpanded $L^p$ equation. A directed rational bracket for $2\sqrt\mu$ additionally certifies the constant itself, not finite-dimensional physical norms. All tests, resource limits and numeric qualifications are in [certification89.md](certification89.md) and [validation89.json](validation89.json). The analytical physical cap theorem above is not inferred from a Gaussian fit or unverified dimension onset.
