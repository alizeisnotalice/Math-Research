# Physical height and support laws on successful rectangles

Use the genuine noncentred axis-parallel strong maximal operator, strict threshold one and the physical modular $J(F)=\int\Phi_N(F)$, $\Phi_N(t)=t(1+\log^+t)^{N-1}$. Write $L=\log N$, $\rho=\sqrt{L/N}$, $\mu=1+e^{-\mu}$, $\sigma=\sqrt\mu$, $\bar\Phi(z)=\sigma z/\sqrt N$, $d=\sigma\sqrt N z$, $q_0=\mu N-d$. Here $\bar\Phi$ is the standard normal upper tail. The physical benchmark has $\log Q_N=q_0+\alpha L+o(L)$, $\alpha=(9-\mu)/24$. Constants are asymptotic for sufficiently large integer $N$.

## 1. Exact purification and uniform quadratic cost

For fixed $N$ and $1\le A<\infty$, the exact supremum $C_N(A)$ can be taken over $F$ taking values in $D_A=\{0\}\cup[1,A]$. In fact for every original source with $0<J(F)<\infty$ and $0\le F\le A$, leave $F1_{F>1}$ untouched. In each half-open mesh-$2^{-j}$ cube $I$, choose a measurable subset of $I\cap\{F\le1\}$ of volume $\int_I F1_{F\le1}$, which is possible by atomlessness. On the low domain replace $F$ by the indicator of the selected union. The cost is exactly the same, by linearity of $\Phi_N$ there. The integral of the change is zero on every full mesh cube. On a fixed bounded rectangle the error is bounded by the volume of its intersection with the boundary cubes, which tends to zero. A strict witness therefore persists eventually for its receiver. Fatou bounds the original strict halo volume by the liminf of the new ones. The converse inclusion of source classes proves equality. In r079/rational79 §A the strict-halo density transfer was already proved for ordinary compact rational steps. The new feature here is the gap in the height range with *exact* cost and exactly unchanged superunit source. Purification can alter the original subunit spatial support. All statements below about witness holes concern the actual support of the purified input.

The identities $\log\Phi_N(1+x)=\log(1+x)+(N-1)\log(1+\log(1+x))$ and Taylor's theorem on $0\le x\le1/4$ give
\[
 \log\Phi_N(1+x)=Nx-(N-1/2)x^2+O(Nx^3).
\]
There is also a global directed remainder bound, for all $x\ge0$ and $N\ge1$:
\[
 \frac{x^3}{6}\,f_N'''(x)\le R_N(x):=f_N(x)-Nx+(N-1/2)x^2\le\frac{7N-5}{6}x^3,
 \quad f_N(x)=\log\Phi_N(1+x).
\]
Indeed $f_N'''(x)=[2+(N-1)H(\log(1+x))]/(1+x)^3$, where $H(u)=(2u^2+7u+7)/(1+u)^3>0$ and $H'(u)=-(2u^2+10u+14)/(1+u)^4<0$ on $u\ge0$. Thus $f_N'''$ decreases from $7N-5$, and the integral Taylor remainder has the stated bounds. The constant in the preceding $O(Nx^3)$ is therefore absolute. Thus at $A=1+2\sigma\rho$, uniformly over measurable $E$ and $0\le u\le2\sigma$ on $E$, the physical cost of $F=1_E(1+\rho u)$ is $e^{O(L^{3/2}/\sqrt N)}$ times
\[
 \int_E\exp\{\sqrt{NL}u-(1-1/(2N))Lu^2\}.
\]
Finite nonzero cost is equivalent to $0<|E|<\infty$ here. Applied before taking suprema, this and the independently audited r089 cap theorem give a multiplicative $1+O(L^{-1/4}+L^{3/2}/\sqrt N)$ reduction of the unrestricted norm to this marked-set supremum with the original halo. For any witnessing rectangle $R$, $m_R F>1$ iff $|R\setminus E|/|R|<\rho m_R(1_Eu)$, equivalently $\theta_R<\rho\bar u_R/(1+\rho\bar u_R)$ where $\bar u_R$ averages on $E\cap R$. Here and below $m_R$ denotes averaging on the whole rectangle.

## 2. Physical theorem

Let $F$ be *any* such purified input with critical cap, and $P>0$ fixed. There exists a measurable receiver set $B=B(N,P,F)$ with
\[
 |B|\le Q_N N^{-P}J(F).                                             \tag{1}
\]
Every $x\notin B$ and **every** finite positive-volume rectangle $R\ni x$ with $m_RF>1$ satisfy, simultaneously,
\[
 \left|\frac{|R\setminus E|}{|R|}-\frac dN\right|\le C_P\frac{L^2}N,
 \qquad 0<m_RF-1\le C_P\frac{L^2}N.                                \tag{2}
\]
In particular $m_R[(F-1)_+]=d/N+O_P(L^2/N)$. The rescaled excess on the occupied part has conditional distribution obeying
\[
 \sup_{y\in\mathbb R}\left|
 \frac{|\{v\in E\cap R:\sqrt N(F(v)-1)-d/\sqrt N\le y\}|}{|E\cap R|}
 -\Phi(y/\sigma)\right|\le C_P L N^{-1/4}.                         \tag{3}
\]
These claims compare the actual local uniform volume distribution inside a common witnessing rectangle, for arbitrary measurable input geometry. The exceptional set can contain successful receivers and costs physical volume (1). In particular the estimates apply to all successes outside it without any local cost-efficiency hypothesis. A near-$Q_N$ sequence loses only a vanishing fraction of its halo there if $P$ is chosen greater than its fixed logarithmic efficiency loss. Every retained witness carries at least $c_P N^{1/4}/L$ distinct essential heights in a fixed central rescaled interval, in the sense that fewer than this many distinct values cannot support the whole restricted law there. Indeed every jump of the conditional cdf has size at most twice the right side of (3), whereas $[-\sigma,\sigma]$ has fixed positive limiting mass.

### A common-source perturbation test

Fix $s=q_0-DL$ with $D>0$ fixed. For real $t$ define the nondecreasing envelope
\[
 W_s(t)=\begin{cases}0,&t<s,\\
 \min(1+4\sigma\rho,\Phi_N^{-1}(e^{t-s})),&t\ge s.
 \end{cases}
\]
Here the inverse is the ordinary unconditional inverse of the increasing function $\Phi_N$. For increasing (non-strict monotonicity suffices) transformations $T$ on $D_{1+4\sigma\rho}$ with $T(0)=0$ and $0\le T\le3$, the following is a physical bound. Let $G'\sim\mathcal N(\mu N+C_H L,\mu(N+t_0))$, $t_0=O(L)$ as in r021 §6 (15). If $b>0$ and
\[
 \mathbb E T(W_s(G'))+6N^{-4}\le b,
\quad\text{then}\quad |\{M_S(T(F))>b\}|\le e^s J(F).               \tag{4}
\]
It suffices here to take continuous $T$ on the disconnected domain (all tests below have this property). To prove (4) truncate physically in boxes and increase the boxes to cover the full strict halo; $T(0)=0$ ensures the truncations commute with $T$. For a bounded compact truncation use r011 §§1–2's global rearrangement comparison. If its transformed halo has volume greater than $e^sJ(F)$, choose a normalized volume coordinate $a$ strictly between these two numbers where the r011 upper applies. Its nondecreasing profile before transformation is $h(t)=F_{\rm trunc}^*(a e^{-t})$ for $t\ge0$. Rearrangement of the monotone continuous transform gives the transformed profile $T(h)$ a.e., with consistent quantile representatives. The r011 feedback law puts no mass on a terminal Lebesgue null set (the last increment has a conditional density). Thus some predictable policy obeys $\mathbb E_\pi T(h(S_N))>b$, with
\[
 j(h)=\int_0^\infty e^{-t}\Phi_N(h(t))dt\le J(F)/a<e^{-s}.
\]
Monotonicity gives $\Phi_N(h(t))\le e^t j(h)$. The purified height gap forces $h=0$ before $s$, hence $h\le W_s$ on the terminal domain. The r021 all-policy one-sided cdf inequality in its increasing-test form bounds the expectation by $\mathbb ET(W_s(G'))+2N^{-4}\operatorname{osc}(T(W_s))$. This contradicts the strict threshold. No modular of the artificially perturbed source is substituted for $J(F)$; it is the monotone pushforward of the original rearrangement that permits both upward and downward tests.

### The estimates for the envelope

For this $s$ and $G'$ one has
\[
 p_s=\mathbb P(G'<s)=d/N+O_D(L^2/N),\qquad \mathbb E W_s(G')\le1+O_D(L^2/N). \tag{5}
\]
Indeed $z\sim\sqrt L$, $f_{\mathcal N(\mu N,\mu N)}(q_0)\sim L/N$ by Gaussian Mills. The standardized threshold for $G'$ differs from $-z$ by $O_D(L/\sqrt N)$ and the standard density in that short interval is $O(L/\sqrt N)$, proving the first formula. On $t\ge s$, the quadratic inverse estimate from §1 gives $W_s(t)\le1+(t-s)/N+O(\rho^2)$ until the cap is reached, and the same upper line then continues to bound the constant cap. Thus its expectation is at most $1-p_s+(\mathbb EG'-s+\mathbb E(s-G')_+)/N+O(\rho^2)$. The overshoot $\mathbb E(s-G')_+=O_D(1)$ follows by Mills, since $\sqrt N\phi(z)/z^2=O(1)$ and the small threshold shift changes this by a bounded factor. Now $\mathbb EG'-s=d+O_D(L)$, proving (5).

Write $Y(v)=\sqrt N(v-1)-d/\sqrt N$ on occupied heights. The finite subprobability $\nu_s$ defined by $\mathbb E[1_{G'\ge s}1_{Y(W_s(G'))\in\cdot}]$ satisfies
\[
 \sup_y|\nu_s(({-}\infty,y])-\Phi(y/\sigma)|=O_D(L/\sqrt N).       \tag{6}
\]
On $s\le G'\le\mu N+\sigma\sqrt{NL}$ the cap $1+4\sigma\rho$ is inactive for large $N$, by §1 and $d\le\sigma\sqrt{NL}$. The inverse expansion gives $Y(W_s(G'))=(G'-\mu N)/\sqrt N+O_D(L/\sqrt N)$ uniformly there. The excluded Gaussian regions have total probability $O_D(\rho)$ by (5) and the ordinary upper tail bound. The normal variable $(G'-\mu N)/\sqrt N$ has mean $O(L/\sqrt N)$ and variance $\mu(1+O(L/N))$; the uniform cdf displacement is $O(L/\sqrt N)$ by boundedness of $\phi$ and of $|y|\phi(y)$. Its density is uniformly bounded, so the preceding uniform coupling error implies (6).

### Applying the tests to real rectangles

Fix $D=P+5$. First $T(v)=v$ in (4)–(5) excludes any witness of mean above $1+C_D L^2/N$ at cost $e^sJ(F)$. For $\kappa=1/2$ use both transforms $T_\pm(v)=v\pm\kappa$ on occupied heights and $T_\pm(0)=0$. They are nonnegative and monotone on the purified range, and on every successful rectangle their averages exceed $1\pm\kappa(1-\theta_R)$ respectively (the same constant $1$ comes from its original strict mean). On the scalar side their means are at most $1+O_D(L^2/N)\pm\kappa(1-p_s)$. Testing each sign via (4) and increasing $C_D$ excludes either deviation in (2), each at most $e^sJ(F)$ in volume.

For (3) put $h=L N^{-1/4}$, $a=h/(2\sqrt N)$. For each grid point $y$ of mesh $h$, take the continuous descending ramp $\psi_y(w)=1$ for $w\le y$, $0$ for $w\ge y+h$, affine in between. Use
\[
 T_\pm(v)=v\pm a\psi_y(Y(v))\quad(v\ge1),\qquad T_\pm(0)=0.
\]
Both are increasing since the perturbation has Lipschitz constant at most $1/2$ in $v$. They lie in $[0,3]$ for large $N$. Let $b_y=\int\psi_y\,d\nu_s$, $r_y=m_R(1_E\psi_y(Y(F)))$. If $\pm(r_y-b_y)>C_D L^2/(Na)$ on a successful rectangle, its transformed mean exceeds $1\pm a b_y+C_D L^2/N$. Equations (4) and (5) bound the whole corresponding exceptional halo by $e^sJ(F)$ for each sign and grid point; the witnesses need not be selected measurably. All the exceptional sets used here are fixed-threshold ordinary physical halos, so measurable supersets suffice. Formula (6) yields $|b_y-\Phi(y/\sigma)|\le O_D(h+L/\sqrt N)$ by cdf sandwiching across one mesh step. Since $L^2/(Na)=2h$, all nonexceptional witnesses have $r_y=\Phi(y/\sigma)+O_D(h)$ simultaneously on the grid. Take the finite grid to extend beyond $[-5\sigma\sqrt L,5\sigma\sqrt L]$ by two mesh steps. It has fewer than $N$ points for large $N$ and contains the entire possible $Y(F)$ range on $E$. Monotonicity between adjacent points and the two ramp sandwiches bound the *subprobability* cdf of the occupied heights on every witness by the same Gaussian cdf up to $O_D(h)$; outside the grid use the support range and Gaussian tails. Normalizing by $1-\theta_R=1-O_D(\rho)$ using (2) gives (3). The union of the exceptional physical halos has volume at most $(2N+3)e^{q_0-DL}J(F)\le Q_N N^{-P}J(F)$ for large $N$, proving (1).

## 3. Arbitrary original grey-scale sources at the critical cap

For each fixed $P>0$, all sufficiently large $N$ depending only on $P$, and **every** measurable $0\le f\le 1+2\sigma\rho$ with $0<J(f)<\infty$, there is a measurable $B_f$ with $|B_f|\le Q_NN^{-P}J(f)$. Simultaneously for every $x\notin B_f$ and every finite positive-volume rectangle $R\ni x$ on which $m_Rf>1$, define $\delta_R=m_R(1-f)_+$ and $p_R=m_R\min(f,1)=1-\delta_R$. Then
\[
 \delta_R=d/N+O_P(L^2/N),\qquad m_R(f-1)_+=d/N+O_P(L^2/N),\qquad 0<m_R f-1\le C_P L^2/N.                     \tag{7}
\]
Here $\delta_R$ measures missing mass below one. It agrees with the volume fraction of zeros for purified inputs. For general grey levels it includes the fractional deficits from $0<f<1$.

Extend $Y(v)=\sqrt N(v-1)-d/\sqrt N$ to every $v\ge0$ and put $y_0=Y(1)=-d/\sqrt N$, $h=L N^{-1/4}$. Both the original whole-rectangle height law and the conditional law obtained by purifying the subunit values satisfy
\[
 \sup_y|m_R1_{\{Y(f)\le y\}}-\Phi(y/\sigma)|\le C_P h,
 \qquad \sup_y|V_R(y)/p_R-\Phi(y/\sigma)|\le C_P h,                 \tag{8}
\]
where $V_R(y)=m_R[1_{f>1}1_{Y(f)\le y}+f1_{f\le1}1_{y_0\le y}]$ is the unnormalized occupied mass below the cutoff. The constants in (7)–(8) are uniform in the choice of $f,x,R$.

**Proof.** Apply §1 to $f$ and take $F_j$, with $J(F_j)=J(f)$ exactly. Put $B_f=\liminf_j B(N,P,F_j)$. Fatou gives its asserted volume bound. At any $x\notin B_f$ infinitely many indices obey $x\notin B(N,P,F_j)$. Along these indices every fixed original strict witness $R\ni x$ eventually remains successful because $m_R F_j\to m_R f$. The same boundary-cell cancellation as in §1 yields $m_R1_{F_j=0}\to\delta_R$: the integral on the original subunit portion of the low indicator converges to $\int_R f1_{f\le1}$, and the superunit part is unchanged. Taking limits in (2) proves the deficit and mean statements. The exact scalar identity $m_R(f-1)_+=m_Rf-1+\delta_R$ gives the second assertion of (7).

On occupied heights the superunit contribution to the unnormalized cdf stays fixed and the low contribution is exactly one convergent scalar times $1_{y_0\le y}$. Hence these cdfs converge *uniformly in $y$* to $V_R$. Their total masses converge to $p_R>0$, so (3) passes to the second assertion of (8). Since $\delta_R=O_P(\rho)=o(h)$, the same error bounds $V_R$ itself relative to the Gaussian cdf. For $y\ge y_0$ the original cdf on the whole rectangle equals $V_R(y)+\delta_R$ exactly. For $y<y_0$ it is at most $m_R1_{f\le1}=V_R(y_0)+\delta_R=O_P(h)$, while $\Phi(y/\sigma)\le\Phi(-z)=d/N=O(\rho)$. These two ranges prove the first assertion of (8). No witness-dependent exceptional set has been added.

In particular each retained $R$ puts a fixed positive fraction of its volume in $Y(f)\in[-\sigma,\sigma]$. For all sufficiently large $N$ this entire interval corresponds to heights strictly above one ($z>1$). Each atom has mass at most $2C_P h$ by (8) and continuity of the limiting cdf. Thus either the central superunit range on $R$ is not supported on finitely many heights, or it uses at least $c_P N^{1/4}/L$ distinct values of positive measure. Consequently for **any sequence** of original grey-scale sources within the critical cap whose entire superunit ranges are supported on $o(N^{1/4}/L)$ distinct values modulo null sets, the whole physical halo is contained in $B_f$ for every fixed $P$ and sufficiently large $N$. Its score is then $\le Q_N N^{-P}$, regardless of how many subunit levels the sources use.

## 4. Remaining comparison

The known general logarithmic gap $(\mu-1)L/8+o(L)$ persists. The argument proves (1)–(3) by perturbing individual common physical witnesses, and estimates actual receiver volume; it does not assert pointwise near-optimal local cost as excluded in r069. Nor does it charge fixed-corner columns or use a Gaussian as the physical norm. A concrete next physical quantity is the nested family $E_y=\{v\in E:Y(F(v))>y\}$. On every retained rectangle, $m_R1_{E_y}=\bar\Phi(y/\sigma)+O_P(LN^{-1/4})$ for all $y$ while its hole mass obeys the much sharper (2). Separate constant-threshold indicator bounds for those sets discard their common sides. Even an $O_P(LN^{-1/4})$ threshold error at a central Gaussian quantile shifts the scalar log-volume location by order $\sqrt N\,L N^{-1/4}$, too large on the $L$ scale. One specific next test is a joint physical halo estimate using the *same* sides for base deficit and two central layer densities, with their actual nested volumes charged at their modular increments; its gain would have to come from incidence across receivers. Another is to sharpen the transition-scale tail/hole comparison using common-side perturbations localized at the height-one plateau, to compare the signed $L^2/N$-scale cdf displacement against the fixed-action coefficient. Neither estimate follows from (1)–(3) alone.

## 5. Independent local software verification

The final [experiment90.py](../../code/experiment90.py) constructs directed rational logarithmic bounds for the positive Taylor remainder at $x=2\sigma\sqrt{L/N}$, using the preceding third-derivative formula. It computes $\sqrt{\mu L}$ first and only then divides by $\sqrt N$, so the fixed absolute grid of the square-root helper cannot erase a tiny critical height. Its independent mpmath value is the positive integral Taylor remainder. [verify90.wl](../../code/verify90.wl) instead subtracts the exact quadratic from the original logarithmic expression at adaptive high precision. Both use 80/160 digit targets on thirteen scalar dimensions $N=2^k$ from $2^{16}$ through $2^{1024}$, with extra guard digits, and check the Gaussian quantile and grid-error balances. The grid is counted without enumerating its points. Twelve Wolfram symbolic gates check the derivative, positivity and scale identities.

The complete replay, arithmetic interval qualifications and resource reports are in [certification90.md](certification90.md) and [validation90.json](validation90.json). These computations check the price expansion and scalar scales, not the physical theorem's finite starting dimension or its asymptotic constants. Equations (1)–(8) follow from the physical perturbation and purification proofs above.
