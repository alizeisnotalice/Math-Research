# Steiner compression can decrease a physical strong halo; a different support-plane bound for log-convex sets

Use the actual noncentred axis-parallel strong maximal operator. Write $L=\log N$, $\mu=1+e^{-\mu}$, $\sigma=\sqrt\mu$, $\bar\Phi(z)=\sigma z/\sqrt N=\tau$, $d=N\tau$, $q_0=\mu N-d$, $c_Q=(9-\mu)/24$, $p_N=1-\tau-BL/N$ with fixed real $B$. In the r092 endpoint application $B>0$ is fixed. All $N$-dependent conclusions are for sufficiently large $N$.

## 1. A strict physical counterexample to coordinate Steiner compression

In two dimensions let $E=\{(x,y):0<y<1,\ y<x<y+1\}$. Its area is one. Horizontal-fiber Steiner compression gives the unit rectangle $E^*=[-1/2,1/2]\times[0,1]$ up to null sets. For all $0<\eta\le10^{-4}$ the strict halo volumes satisfy
\[
 |\{M_S1_E>1-\eta\}|-|\{M_S1_{E^*}>1-\eta\}|\ \ge\ \tfrac18\sqrt\eta.       \tag{1}
\]
Indeed the original strict halo contains $E$ up to a null set. Just outside its sloping left side take $x=y-\delta$, $1/4<y<3/4$ and $0<\delta<(1/4)\sqrt{2\eta}$. The actual rectangle $[x,x+1/4]\times[y-1/4,y]$ has only one missing triangle, of area $\delta^2/2$. Its opposite sloping boundary is not met: throughout this rectangle the largest horizontal minus vertical coordinate is $1/2-\delta<1$. The whole patch is strictly successful; its area is $\sqrt{2\eta}/8$, by the unit Jacobian of $(\delta,y)\mapsto(x,y)$. The full halo volume of *any* unit rectangle at threshold $1/A$, $A>1$, is exactly $1+4A\log A$. One coordinate has maximal density $1/(1+\operatorname{dist}(x,I))$ outside a unit interval, one inside. The exterior face regions give $4(A-1)$ and the four exterior corners give $4(A\log A-A+1)$, in addition to the rectangle itself. Thus at $A=(1-\eta)^{-1}$ its excess is at most $4\eta/(1-\eta)^2$. Subtracting yields (1) throughout the stated range.

A rational certificate at $\eta=1/10000$ restricts to $\delta<1/320$. Even the boundary value of the witness mean is $12799/12800>9999/10000$. The patch has area $1/640$. Using $\log(1/(1-\eta))\le\eta/(1-\eta)$, the exact rectangle's entire exterior excess is at most $40000/99980001$, leaving a strictly positive certified difference of $74380001/63987200640$. These fractions were independently evaluated by local Python `fractions.Fraction`.

**Higher dimensions and an efficient receiver-volume scale.** The direction of the difference in (1) holds on a whole near-one level interval. For any measurable suffix set $D\subset\mathbb R^m$ with $0<|D|<\infty$, use the physical products $E\times D$ and $E^*\times D$. Their maximal functions are exactly the products of the two factor maxima, by independent choices of rectangular sides and positivity. Write $b(z)=M_{S,m}1_D(z)$. At $p\ge1-10^{-4}$ all slices with $b(z)>p$ use a prefix threshold $p/b(z)\in[p,1)$. On slices with $b(z)\le p$ neither strict halo contributes. Thus the full-dimensional volume difference is at least
\[
 \frac18\int_{b>p}\sqrt{1-p/b(z)}\,dz>0.                              \tag{2}
\]
The positive lower bound follows even on $D$ where $b=1$ a.e. This proof uses all actual suffix slices and all their varying thresholds; a two-dimensional difference at just one threshold would not suffice.

At $m=N-2$, $p=p_N$, $\Delta=L/N$, choose one actual equal-weight r093 cross $D=E_{T,\mathbf1}$ in $m$ dimensions with $T$ sufficiently large that its halo at $p+\Delta$ reaches at least half the limiting r093 ratio. Both $p$ and $p+\Delta$ in dimension $m$ have the critical form with bounded additive $O(L_m/m)$, because $\tau_N-\tau_m=o(L/N)$; monotonicity brackets them using fixed coefficients in r092 (8). Also $q_{0,m}=q_{0,N}-2\mu+o(1)$, from implicit differentiation of the Gaussian root ($z'_N=O(1/(Nz_N))$). Then
\[
 \frac{|\{M1_{E\times D}>p\}|-|\{M1_{E^*\times D}>p\}|}{|D|}
 \ \ge\ Q_N N^{-1/2+o(1)}                                             \tag{3}
\]
where this notation is a logarithmic lower bound. In (2) retain $b>p+\Delta$, so $1-p/b\ge\Delta/(p+\Delta)$; the suffix ratio there is $\exp(q_{0,N}+c_Q L+o(L))$. Both full ratios to common source volume $|D|$ have logarithm $q_{0,N}+c_QL+o(L)$: the two-dimensional interiors give the suffix lower bound and each product halo is contained in the product of its two individual $p$ halos. Their prefix halo areas are uniformly bounded near one, e.g. by domination by a fixed enclosing rectangle. The suffix full upper is r093 at $p$. This exhibits strict Steiner decrease on actual benchmark-order products, and this particular shear product still has coefficient $c_Q$. It rules out the proposed monotone compression of arbitrary inputs to coordinatewise downsets.

## 2. A radial integral upper for log-convex unconditional downsets

Assume $E$ is invariant under sign in each coordinate, is downward closed in the positive orthant, and has finite positive measure. Suppose additionally its coordinatewise logarithmic image $K=\{\xi:(e^{\xi_i})_i\in E\}$ in that orthant is convex. Replace it by its closure since its boundary has measure zero. In the log chart write $\xi=\theta+s\mathbf1$, $\theta\in V=\{\sum\theta_i=0\}$. The roof $u(\theta)=\sup\{s:\theta+s\mathbf1\in K\}$ is finite everywhere. The positive set contains one nondegenerate origin box giving a lower bound on every ray. Inclusion of successively larger boxes under each radial source point and finiteness of volume give the upper bound. Downward closure implies $|u(\theta)-u(\theta')|\le\|\theta-\theta'\|_\infty$ and convexity makes $u$ concave. At almost every $\theta$, the log supporting hyperplane at $\theta+u(\theta)\mathbf1$ has unit-sum normal
\[
 w(\theta)=N^{-1}\mathbf1-\nabla_V u(\theta),\qquad w_i\ge0 .                         \tag{4}
\]
The formula comes from the tangent inequality for the roof and the normalization; nonnegativity follows by subtracting positive coordinate vectors inside $K$.

For any normalized weight vector let $a_p(w)$ be the supremum, over deterministic side vectors, of the upper-hit quantile at probability $p$ of $\sum_i w_i Z_{r_i}$ as in r093. Define the probability $\nu_E(d\theta)=e^{Nu(\theta)}d\theta/\int_V e^{Nu}$. Then the following whole-volume bound holds without a bound on the number, orientations or spacings of faces:
\[
 \frac{|\{M_S1_E>p\}|}{|E|}\le \int_V \exp\{N a_p(w(\theta))\}\,d\nu_E(\theta),
 \qquad \int_V w_i(\theta)\,d\nu_E(\theta)=1/N.                       \tag{5}
\]
Indeed the honest box at a positive receiver may use intervals $[-r_ix_i,x_i]$ by separate even monotonicity, exactly as in r093. Under its uniform distribution the source log vector is $\xi-Z$ with independent increments of the chosen deterministic vector. For $\xi=\theta+s\mathbf1$ to be successful at some $s>u(\theta)$ its supporting halfspace requires $\mathbb P(\sum_i w_i Z_{r_i}\ge s-u(\theta))>p$. Hence the radial halo extends no further than $u(\theta)+a_p(w(\theta))$. This is also an upper for points below the source roof. Volume in this chart is a common constant times $e^{Ns} ds\,d\theta$. Integrating proves the first inequality. Quantile measurability follows also by using rational side vectors and thresholds, via continuity in distribution of the increments and the absence of point masses when at least one coefficient is positive. To check the second identity, $e^{Nu}$ is integrable, differentiable a.e. with integrable weak derivatives since $u$ is Lipschitz. Integration by parts along each tangent direction in $V$ using cutoffs of radius $R$ and gradient $O(1/R)$ gives $\int\nabla_Vu\,e^{Nu}=0$. Formula (4) proves the claim. This source-dependent angular probability and normalization survive arbitrary switching of the supporting face with the receiver; (5) is not an inequality for a feedback-in-the-integral policy.

## 3. Sharp coefficient for multiple bounded-weight log faces

Fix $A<\infty$. Uniformly over $\lambda_i=Nw_i\in[0,A]$, $\sum_i\lambda_i=N$, and over all deterministic side vectors, the scaled quantile obeys
\[
 \sup_{w:\,Nw_i\le A} N a_{p_N}(w)\le q_0+c_Q L+o_{A,B}(L).          \tag{6}
\]
For $A\ge1$ equality in the logarithmic coefficient is reached by equal weights, by r092. We give the details because r093 normalized instead by $\max w_i=1$ and did not give (6) for this constraint.

Write $m_i,v_i,k_i$ for the mean, variance and third cumulant of $Z_{r_i}$, and
\[
 D=\sum_i\lambda_i(\mu-m_i),\quad H=\sum_i(\lambda_i-1)^2,\quad
 M=\mu N-D,\quad V_2=\sum_i\lambda_i^2 v_i,\quad K_3=\sum_i\lambda_i^3 k_i .
\]
At least $c_A N$ coefficients are $\ge1/2$ since they sum to $N$ and are bounded by $A$; thus $V_2\asymp_A N$. r093's weighted Fourier argument (absolute Edgeworth cdf error $O_A(1/N)$, exponential factor in these $c_AN$ coordinates) and r092's quantile inversion give, uniformly,
\[
 s=M-\sqrt{V_2}z+\frac{K_3}{6V_2}(z^2-1)+O_{A,B}(1)                  \tag{7}
\]
for the $p_N$ hit quantile of $\sum\lambda_i Z_{r_i}$. Per-coordinate stability from r022 implies $|v_i-\mu|+|k_i-\mu(3-\mu)|\le C\sqrt{\mu-m_i}$. Cauchy-Schwarz gives $V_2=\mu(N+H)+O_A(\sqrt{ND})$, since $\sum\lambda_i^3\le A^2N$. Consequently
\[
 \sqrt{V_2}-\sqrt{\mu N}\ge c_A H/\sqrt N-C_A\sqrt D .
\]
For any vector with $s\ge q_0-O(L)$, (7) forces $D+c_A zH/\sqrt N\le O_{A,B}(L)+C_Az\sqrt D$. This yields $D=O_{A,B}(L)$ and $H=O_{A,B}(\sqrt{NL})$. The discarded vectors do not affect the desired upper by taking a fixed sufficiently large coefficient in the localization test. Now $V_2-\mu N=O_{A,B}(\sqrt{NL})$, and $\sum\lambda_i^3=N+O_A(H)$ by expanding about one and using $\sum(\lambda_i-1)=0$. Another Cauchy-Schwarz estimate on the weighted third cumulants gives $K_3/(6V_2)=(3-\mu)/6+o_{A,B}(1)$ uniformly.

Put $\ell=d/(\mu N)$. Taylor gives $\sqrt{V_2}z=d+\ell(V_2-\mu N)/2+o(1)$. The r022 coordinate minimization, uniform for $0\le t\le A\ell$, reads $\mu-m_i+t(v_i-\mu)/2\ge-\mu(\mu-1)t^2/8-Ct^3$. Use it at $t=\ell\lambda_i$, multiply by $\lambda_i$, and sum to obtain
\[
 D+\ell(V_2-\mu N)/2\ge (\mu\ell/2)H-\frac{\mu(\mu-1)}8\ell^2[N+O_A(H)]-O_A(N\ell^3)\ge -g z^2+o(1).
\]
The first positive term absorbs the $O_A(\ell^2H)$ loss. Substituting in (7) proves (6), since $(3-\mu)/6+g=c_Q$ and $z^2/L\to1$.

A second allowed class of supporting normals has at most $r$ nonzero coordinates for a fixed $r$. With probability bounded below by a fixed positive constant depending only on $r$, each of these increments is $\le\mu/2$: take all their unshifted exponential components and restrict each exponential. Thus $a_{p_N}(w)\le\mu/2$ eventually on this class. In particular coordinate cutoff faces of bounded boxes are allowed. Combining with (5), **every** finite-volume unconditional log-convex downset for which almost every supporting normal belongs to either the bounded-weight or this sparse class satisfies
\[
 |\{M_S1_E>p_N\}|/|E|\le\exp(q_0+c_Q L+o_{A,B,r}(L)).              \tag{8}
\]
This covers genuine intersections of arbitrarily many weight-qualified log halfspaces and coordinate cutoff planes, including smooth variation of their normals. Families may change with $N$ and the number of facets need not be bounded. The equal-weight r093 cross, with coordinate cutoff faces, proves optimality of the coefficient within the class. As an explicit family of several distinct dense faces, intersect the unit cube with conditions $\prod |x_i|^{w_i^{(j)}}<e^{-T_j}$, where every row of positive weights sums to one with $Nw_i^{(j)}\le A$. The side choices and the active row may vary between receivers. No union bound over facets appears in (8).

## 4. Unconditional convex bodies are exponentially below the benchmark

There is an absolute $p_0<1$ such that for all $p>p_0$, all dimensions, and every unconditional convex body $E$ of positive volume,
\[
             \{M_S1_E>p\}\subset (1+\sqrt2)E\quad\text{modulo null sets}.
                                                                            \tag{9}
\]
Hence the full volume ratio is at most $(1+\sqrt2)^N$, whose log rate $0.8813735870\ldots$ is strictly below $\mu$. This argument permits all orientations of supporting faces of the unconditional body.

For $x$ in the positive orthant let $a=\|x\|_E>1$, the ordinary radial gauge. The radial boundary point $x/a$ has a supporting functional $n$ with nonnegative coefficients and $n\cdot(x/a)=1$ (unconditional monotonicity, or take absolute coefficients in the polar body). Put $b_i=n_i x_i/a$, so $\sum b_i=1$. For each honest averaged rectangle reduced as in r093, its sampled absolute coordinates are $x_i V_i$, with $V_i$ independent and distributed as the absolute value of a uniform point on $[-r_i,1]$. Hitting $E$ requires $\sum_i b_iV_i\le1/a$. The true mean of every $V_i$ is $(1+r_i^2)/(2(1+r_i))\ge\sqrt2-1$ by one-variable differentiation. Its variance is bounded below by an absolute $v_0>0$: its unshifted uniform $[0,1]$ mixture component has weight at least $1/2$ (the other component is uniform on $[0,r_i]$), so $v_0=1/24$ works. All fourth central moments are at most one.

For $Y=\sum_i b_i(V_i-\mathbb EV_i)$, independence implies $\mathbb E Y^4\le C(\mathbb E Y^2)^2$ with an absolute $C$ (expand the fourth moment and use the variance lower bound). Interpolation $\mathbb E|Y|\ge(\mathbb EY^2)^{3/2}/(\mathbb EY^4)^{1/2}$ and Cauchy-Schwarz for $Y_+$ give $\mathbb P(Y>0)\ge1/(4C)$. If $1/a<\sqrt2-1$, at least this much sample probability misses the body for every box; taking $p_0=1-1/(4C)$ proves (9). Unconditional ordinary convexity is used in this proof. Arbitrarily oriented, nonunconditional convex bodies are not covered by this homothety.

## 5. What is left in the geometric upper

Steiner compression generally decreases the halo in (1)–(3), so it cannot justify taking only unconditional sources in an upper proof. The log-convex case has the source-specific whole-volume integral (5); (8) settles all the weight-qualified multiface sets and (9) treats ordinary unconditional convex bodies regardless of which faces are steep. The unrestricted weight concentration issue in (8) is settled by the independent-tail comparison in §6 below. Arbitrary unconditional downsets can also have nonconvex logarithmic images, for which a globally supporting log halfspace need not exist. Finally, arbitrary physical sources are not covered by either of those classes. The actual full $H_N(p_N)$ upper at coefficient $c_Q$ and the corresponding general $C_N$ gap remain unresolved.

## 6. Removing the normal-weight restriction by a tilted comparison

This section resolves the unbounded-weight question left after §3. It uses the special independent exponential-plus-shift laws and an actual lower-tail probability comparison. For every fixed $B$, uniformly over **all** $\lambda_i\ge0$ with $\sum_i\lambda_i=N$ and all deterministic side actions,
\[
 \sup_{\lambda,\mathbf r}s_{p_N}\!\left(\sum_i\lambda_i Z_{r_i}\right)
       =q_0+c_Q L+o_B(L).                                                \tag{10}
\]
Consequently (5) gives the full sharp logarithmic coefficient $c_Q$ for the supremum over **all** unconditional log-convex positive finite-volume downsets. The matching lower source is the r093 equal-weight cross. The earlier bounded-normal result remains a shorter proof for its class. The normalization $\sum\lambda_i=N$ in (10) differs from bounding each coefficient by one, and is essential in applying it in every direction of (5).

### 6.1 Global Laplace curvature

Write the action as $c=-\log r\in[0,\infty]$, $Z_c=Y+c B_c$, $Y$ unit exponential and $B_c$ independent Bernoulli with probability $1/(1+e^c)$ (no shift at $c=\infty$). For $u\ge0$ define
\[
 F(u,c)=\log\mathbb E e^{-u Z_c}=\log(1+e^{-(1+u)c})-\log(1+e^{-c})-\log(1+u),
 \qquad \ell(u)=\inf_c F(u,c).
\]
For $u>0$ the unique minimizer $c(u)$ satisfies $e^{uc}=1+u(1+e^{-c})$: the left side strictly increases and the right side strictly decreases in $c$, with opposite signs at the endpoints in the resulting root equation. Equivalently $F_c$ changes from negative to positive there. Substitution gives $\ell(u)=-u c(u)$. The smooth extension to $u=0$ has $c(0)=\mu$ (divide the root equation by $u$ and use the implicit function theorem). All needed derivatives are bounded on a compact $u$ interval, including zero. At the minimum write $p=(1+e^c)^{-1}$ and $a=1+u-p$. Then $F_{cc}=up$, $F_{cu}=p(ca-1)/(1+u)$, $F_{uu}=(1+c^2pa)/(1+u)^2$. Consequently
\[
 u(1+u)^2\ell''=u-p+2pca-c^2pa(1-p)=u[1+pc(2-y)]-p(1-y)^2,
 \quad y=c(1-p)=\frac{\log(1+v)}v,\quad v=u(1+e^{-c}).
\]
Since $0\le1-y\le\min(v/2,1)$, $u\ge v/2$ and $p\le1/2$, the last subtractand is at most $u/2$. This proves $\ell''(u)\ge 1/[2(1+u)^2]$. At zero $\ell'(0)=-\mu$, $\ell''(0)=\mu$ by the centered variance at the mean-maximizing action (or direct expansion of the root). Local Wolfram exact simplification of the Hessian, Schur complement numerator, and its final decomposition returned three zero residuals.

For any $s>0$, set $u_i=s\lambda_i$. Define
\[
 P=N\log(1+s)-\sum_i\log(1+u_i),\quad A=\sum_i[F(u_i,c_i)-\ell(u_i)],
 \quad \Psi(s)=\sum_iF(s\lambda_i,c_i).
\]
Here $P,A\ge0$ and $\Psi(s)\ge N\ell(s)+P/2+A$, since $\ell(u)+(1/2)\log(1+u)$ is convex and $\sum u_i=Ns$. This Laplace ordering alone does not order the cdfs.

### 6.2 A lower saddlepoint prefactor uniform in arbitrary weights

Under the tilt $e^{-sS-\Psi(s)}$, each weighted term is, independently, distributed as $a_i Z_{d_i}$ where $a_i=\lambda_i/(1+u_i)$, $d_i=(1+u_i)c_i$. This follows by tilting the exponential and Bernoulli terms separately. Its variance is between $a_i^2$ and $C a_i^2$, and its centered absolute moments of fixed order $j$ are $\le C_j a_i^j$, uniformly even for extreme actions. Indeed $d^j/(1+e^d)$ is bounded. Write $q=-\Psi'(s)$ and $V_s=\Psi''(s)$. We need the absolute lower bound
\[
 \mathbb E_s[e^{s(S-q)}1_{S\le q}]\ge c/(1+s\sqrt{V_s}).             \tag{11}
\]
Decompose $S=G+D$ into the independent sum of *all* exponential components and all Bernoulli shifts. Put $v^2=\operatorname{Var}G=\sum a_i^2$, so $v^2\asymp V_s$ and $\operatorname{Var}D\le C_0 v^2$. Cantelli for the upper side and Chebyshev for the lower side show $\mathbb P_s(-C_1v\le D-\mathbb E_sD\le v/4)\ge c_1>0$ for some absolute $C_1$ (first use $\mathbb P(D-ED\le v/4)\ge1/(1+16C_0)$, then take $C_1$ large).

For the pure exponential sum the density at $EG+xv$ is at least $c_C/v$ for $-1/2\le x\le C$, for each fixed $C<\infty$, uniformly over its weights. Here are details avoiding any unsupported lower bound on tails of general log-concave laws. Normalize weights by $v$, sort them decreasingly, and extract a coordinatewise convergent subsequence $b_{j,n}\to b_j$ from any sequence. All standardized centered sums converge along it to
\[
 \sqrt{1-\sum b_j^2}\,G_0+\sum b_j(Y_j-1)
\]
with independent standard normal $G_0$ and unit exponentials $Y_j$. To see this, after removing the largest $m$ coordinates, each remaining weight is $\le1/\sqrt m$; the log characteristic function of this tail differs from its centered Gaussian quadratic by $O_t(1/\sqrt m)$ for fixed frequency $t$, from the exact exponential factor and $\sum b_{j,n}^2=1$. Take $n$ then $m$ to infinity. Every limit puts positive mass on every nonempty open interval in $(-1,\infty)$. With a nonzero Gaussian component this is immediate. Otherwise $\sum b_j\ge(\sum b_j^2)^{1/2}=1$; sufficiently many first terms have positive density on an interval strictly containing any prescribed target subinterval above $-1$, and the independent remaining centered tail has arbitrarily small variance. Thus by compactness and Portmanteau there are uniform positive masses on one fixed interval strictly between $-1$ and $-1/2$ and another strictly above $C$. The pure exponential sum is log-concave: convolution with an exponential of scale $a$ sends $g$ to $a^{-1}e^{-x/a}\int_{-\infty}^x e^{t/a}g(t)dt$, still log-concave for log-concave $g$ (for log-concave $h$, $(\int_{-\infty}^x h)/h(x)=\int_0^\infty h(x-z)/h(x)dz$ is nondecreasing where $h>0$). The masses furnish a density lower bound at some point in each of the two intervals; log-concavity gives it throughout the intervening compact interval. This proves the asserted uniform pure-exponential density bound. On the indicated Bernoulli event it applies to the convolution at every $q-h$, $0\le h\le v/8$. Integrating this lower density against $e^{-sh}$ proves (11). The full Bernoulli mixture was never assumed log-concave.

### 6.3 Preliminary localization from the physical critical rank

Fix $\epsilon>0$ and test $q=q_0+(c_Q+\epsilon)L$. Suppose some deterministic weighted vector satisfies $\Pr(S\le q)\le 1-p_N=\tau+BL/N$. If $q\ge ES$, the uniform fourth-moment argument for centered sums used in §4, applied to $-S+ES$, gives a fixed positive lower bound on that probability (all summands have sd comparable to $\lambda_i$ and bounded fourth moments in these units). Hence $q<ES$ eventually. There is a unique positive saddlepoint $s$ with $-\Psi'(s)=q$: the derivative increases strictly to zero since every positively weighted component has nonzero density arbitrarily close to zero. Let $s_*$ solve $-N\ell'(s_*)=q$. It exists by strict convexity near zero, and $s_*\sim z/(\sigma\sqrt N)$. The iid reference uses the *single* action $c_*=c(s_*)$, with the same saddlepoint $s_*$ because $F_u(s_*,c_*)=\ell'(s_*)$. Write $I_*=N\ell(s_*)+s_*q=-z^2/2+O_B(\rho L)$. This expansion follows from $\ell(t)=-\mu t+\mu t^2/2+O(t^3)$ in the scalar minimization and $q=\mu N-d+O(L)$.

The true saddle value $\Psi(s)+sq$ is at most zero, the value at $s=0$. Convexity and the lower bound on $\ell''$ imply $s=O(\rho)$: for $t\le1$ the baseline $N\ell(t)+tq$ is at least $-O(d)t+cNt^2$; for $t\ge1$ it stays positive by convexity and its secant slope near $1$. On this small interval $N\ell(s)+sq\ge I_*+cN(s-s_*)^2$. To control the prefactor let
\[
 d(u,s)=\frac{u-s}{1+s}-\log\frac{1+u}{1+s},\qquad P=\sum_i d(u_i,s).
\]
For $s\le1$ and $|u-s|\le1$ this is comparable to $(u-s)^2$, while for $|u-s|>1$ it is bounded away from zero, uniformly (from $x-\log(1+x)$ with $x=(u-s)/(1+s)$). Therefore $s^2 V_s\le C\sum\min(u_i^2,1)\le C(Ns^2+P)$. The exact tilt identity and (11) now yield
\[
 I_*+cN(s-s_*)^2+A+P/2-\log[C(1+\sqrt{Ns^2+P})]\le\log(\tau+BL/N).
\]
Mills gives $\log(\tau+BL/N)=-z^2/2-\log z+O(1)$ and $Ns^2=O(L)$, $z\asymp\sqrt L$. It follows that $X:=N(s-s_*)^2$, $P$ and $A$ are bounded by absolute constants for large $N$ (the right after cancellation is at most $C+(1/2)\log(1+P/L)$). Thus $s\asymp s_*\asymp\rho$, every $u_i$ is bounded by a fixed constant because $d(u,s)$ grows at infinity, and $\sum(u_i-s)^2\le CP$. In particular $\sum(\lambda_i-1)^2\le C N/L$. A fixed positive fraction of the $\lambda_i$ are in $[1/2,2]$, providing smoothing even when other weights grow.

### 6.4 Comparing the exact tilted Mills factors without losing a logarithmic rank

For $u$ in the bounded interval just obtained, put $A_i=F(u_i,c_i)-\ell(u_i)$, $d_i'=\min(1,|c_i-c(u_i)|)$, taking $d_i'=1$ at infinity. Uniformly on this compact interval,
\[
 A_i\ge c u_i(d_i')^2,\quad
 |u_i(F_u(u_i,c_i)-\ell'(u_i))|\le A_i+C u_i^2d_i',\quad
 |F_{uu}(u_i,c_i)-F_{uu}(u_i,c(u_i))|\le C d_i'.                    \tag{12}
\]
Indeed $D(u,c)=(F(u,c)-\ell(u))/u$ extends smoothly in $u$ at zero with $D(0,c)=\mu-m(c)$. Derivatives needed here are uniformly bounded for $c\in[0,\infty]$ since $c^j e^{-c}$ is bounded, and converge at $c=\infty$. Each minimum is unique, stays in the interior, and $D_{cc}(u,c(u))=(1+e^{c(u)})^{-1}$ including the limit at zero. Local Taylor and compactness away from the minimizing curve prove the first inequality. Along this curve $D_u=0$ by differentiating its zero minimum, so $|D_u(u,c)|\le C\min(1,|c-c(u)|)$. Differentiate $uD$ to prove the second. The third follows by the same uniform Lipschitz/boundedness estimate. Terms at $u_i=0$ cause no division in the resulting inequalities.

Write $h_N=Ns^3\asymp L^{3/2}/\sqrt N$. The localization implies $\sum_i u_i^3\le C(Ns^3+P^{3/2})$ by $|u_i|^3\le C(s^3+|u_i-s|^3)$. Thus
\[
 \sum u_i^2 d_i'\le C\sqrt{A}\sqrt{h_N+P^{3/2}}
 \le C(P+A+\sqrt{A h_N})=:C R.
\]
Both $g(u)=u\ell'(u)$ and $k(u)=u^2 F_{uu}(u,c(u))$ have bounded second derivatives on this interval, and $k'(u)=O(u)$ near zero. For either function its sum on the $u_i$ differs from $N$ times its value at $s$ by at most $CP$, by Taylor and $\sum(u_i-s)=0$. From (12) it follows that
\[
 |s\Psi'(s)-Ns\ell'(s)|\le CR,\qquad
 |s^2V_s-Nk(s)|\le CR.
\]
Since $\Psi'(s)=-q=N\ell'(s_*)$, strict curvature gives $Ns|s-s_*|\le CR$. Hence with $b=s\sqrt{V_s}$ and $b_*=s_*\sqrt{NF_{uu}(s_*,c_*)}$,
\[
 |b^2-b_*^2|\le CR,
 \qquad b\asymp b_*\asymp\sqrt L.                                   \tag{13}
\]
The positivity and order of $V_s$ also follow directly from the many regular coordinates and $\sum(\lambda_i-1)^2=O(N/L)$.

We need a local Gaussian approximation which retains the exact main Mills factor. For the standardized tilted centered sum let $f$ denote its density. Its uniformly valid error is
\[
 \sup_t |f(t)-\phi(t)|\le C\gamma_3,
 \quad \gamma_3=C\sum_i a_i^3/V_s^{3/2}
 \le C(N^{-1/2}+P^{3/2}/L^{3/2}).                                   \tag{14}
\]
Here $a_i=\lambda_i/(1+s\lambda_i)$ as above. To verify the first bound, telescope the product of characteristic functions, replacing each standardized summand by its equal-variance centered Gaussian. The single-factor error on the entire line is $\le C a_i^3|t|^3/V_s^{3/2}$ by third absolute moments and Taylor with remainder, also for its matching Gaussian. Every other product, in every replacement order, has at least $cN-1$ regular factors, each either Gaussian or with the genuine tilted exponential component of scale bounded above and below. In standardized frequency their moduli are at most $(1+c't^2/N)^{-1/2}$ (a Gaussian obeys the same bound with adjusted constants). Thus the integral of $|t|^3$ against their product is bounded by an absolute constant for large $N$, by direct integration. Fourier inversion proves the first estimate. For the second use $s^3\sum a_i^3\le\sum u_i^3\le C(Ns^3+P^{3/2})$, $s^2V_s\asymp L$. For the iid tilted reference the same estimate is $O(N^{-1/2})$.

Set $M(b)=\int_{-\infty}^0 e^{bt}\phi(t)dt=e^{b^2/2}\bar\Phi(b)$. Integrating (14) costs $O(\gamma_3/b)$; since $M(b)\asymp 1/b$, the relative error is $O(\gamma_3)$. Also $|\partial_b\log M(b)|\le C/b$ for $b\ge1$, directly from the integral with the change of variable $v=-bt$. By (13) this gives $|\log[M(b)/M(b_*)]|\le CR/L$. Denote by $F_*(q)$ the actual iid cdf at this same $q$ for the action $c_*$. Using both exact tilt identities gives
\[
 \begin{split}
 \log\frac{\Pr(S\le q)}{F_*(q)}
 &\ge A+P/2-C\left[\frac{P+A+\sqrt{A h_N}}{L}
                    +N^{-1/2}+\frac{P^{3/2}}{L^{3/2}}\right]
 \ge -C N^{-1/2}.                                                      \tag{15}
 \end{split}
\]
The rate difference in the first inequality is at least $A+P/2$ by §6.1 and minimization of the baseline at $s_*$. For the last step, $P$ is bounded, and $C\sqrt{A h_N}/L\le A/4+C' h_N/L^2$ with $h_N/L^2=O(1/\sqrt{NL})$. In particular the potentially larger universal $1/L$ remainder of the *asymptotic* Mills expansion never enters: the two exact Mills factors were compared before any large-$b$ simplification.

Finally r092 (8) bounds the actual hit quantile of the unweighted iid action $c_*$ by $q_0+c_QL+o_B(L)$, strictly below $q-\epsilon L/2$ for large $N$. Its unstandardized density throughout $[q-\epsilon L/2,q]$ is at least $cL/N$: $c_*=\mu+O(\rho)$ by differentiability at zero, its sum mean is $\mu N+O(L)$, variance $\mu N+O(\sqrt{NL})$, and its absolute unstandardized Gaussian density error is $O(1/N)$ by r093 §1. The standardized positions differ from $-z$ by $O(L/\sqrt N)$, negligible relatively in the Gaussian density. Consequently $F_*(q)\ge\tau+BL/N+c_\epsilon L^2/N$. Multiplying by $\exp(-C/\sqrt N)$ still exceeds $\tau+BL/N$ (the product loss on the displayed lower bound is $O(\sqrt L/N)$). This contradicts the assumed competitive tail, proves the upper of (10) for every $\epsilon>0$, and the equal-weight lower follows from r092. Combining (5) and (10) removes both the bounded and sparse normal restrictions from (8).

The generalized equal-input Laplace ordering in §6.1 was checked symbolically; (11)–(15) establish the additional probabilistic control necessary for the cdf conclusion. The resulting sharp log coefficient is a theorem for the full **log-convex unconditional** class, allowing unrestricted weight concentration and arbitrary receiver-dependent face changes. Formula (1) still prevents a global Steiner upper reduction, and the argument does not apply to a nonconvex log roof or to arbitrary indicators. The unrestricted endpoint remains open.

### 6.5 Actual angular remainder for an arbitrary unconditional downset

Without convexity, the positive orthant of any finite positive-volume measurable unconditional downset still has a finite Lipschitz log roof $u$: downward closure alone proves the two coordinatewise ray bounds and the Lipschitz estimate in §2. Inclusion or exclusion of its graph has no volume effect. At differentiability points $w=N^{-1}\mathbf1-\nabla_Vu$ is nonnegative with sum one by monotonicity (differentiate the nondecreasing function $s-u(\theta)$ in the original log coordinates). For the actual vector of independent log increments let $\bar Z=N^{-1}\sum Z_i$, $Z_V=Z-\bar Z\mathbf1$. On the ray $\theta$ the sampled membership condition is $s\le\bar Z+u(\theta-Z_V)$. Define its remainder relative to the tangent by
\[
 D_\theta(Z)=u(\theta-Z_V)-u(\theta)+\nabla_Vu(\theta)\cdot Z_V .
\]
Suppose a measurable $b_N(\theta)\ge0$ satisfies, for a fixed $A<\infty$ and a.e. $\theta$,
\[
       \sup_{\text{deterministic side vectors}}\Pr(D_\theta(Z)>b_N(\theta))\le A L/N.
\]
Then the full physical estimate is
\[
 \frac{|\{M_S1_E>p_N\}|}{|E|}
 \le \exp(q_0+c_QL+o_{A,B}(L))\int e^{N b_N(\theta)}\,d\nu_E(\theta).              \tag{16}
\]
Indeed each successful genuine box at $s$ has probability strictly larger than $p_N-AL/N$ of $w\cdot Z\ge s-u(\theta)-b_N(\theta)$, by deleting just the indicated remainder event. The weighted sum has no atoms. Apply (10) at the shifted density with coefficient $B+A$, and integrate rays as in (5). This gives an explicit remaining physical integral even for a nonconcave roof. If its angular exponential moment is $N^{o(1)}$, the same sharp upper follows. For a concave roof one takes $b_N=0$; no uniform control of the moment is asserted for arbitrary nonconcave roofs.

### 6.6 Multiplicative calibration against the finite-dimensional fixed-vector indicator price

Let $s_p^{\rm det}=\sup_{\mathbf r}s_{p_N}(\sum_i Z_{r_i})$ with the original equal weights. The constants in (11)–(15) are uniform for $q=q_0+O_B(L)$, with dependence on the fixed bound implicit in the $O(L)$ window. By r092 (8), $s_p^{\rm det}$ lies in such a window. Test at $q=s_p^{\rm det}+J/\sqrt L$, where $J$ is a large fixed constant. For large $N$ all these tests lie within one fixed larger logarithmic window, so its comparison constant can be fixed before $J$. For the matched iid action at this $q$, $F_*(s_p^{\rm det})\ge1-p_N$ by definition. Its density throughout this short interval is $\ge c L/N$ by the same local estimate used in §6.4. Hence $F_*(q)\ge1-p_N+cJ\sqrt L/N$. In (15) the relative loss $O(1/\sqrt N)$ applied to this lower bound costs at most $C'\sqrt L/N$ (for each fixed $J$ eventually), since $1-p_N\asymp\sqrt{L/N}$. Choosing $J$ large rules out the competing tail. Together with the inclusion of equal weights this proves
\[
 0\le\sup_{\lambda,\mathbf r}s_{p_N}(\sum_i\lambda_iZ_{r_i})-s_p^{\rm det}
    \le O_B(1/\sqrt L),\qquad
 1\le \frac{H_N^{\rm logconv,unconditional}(p_N)}{\exp(s_p^{\rm det})}
    \le 1+O_B(1/\sqrt L).                                             \tag{17}
\]
The left physical inequality uses the whole r093 cross halo as thickness increases; the right uses the angular integral (5). This compares actual finite-dimensional physical and fixed-indicator prices within a factor tending to one, without asserting that their finite-dimensional suprema are identical.

There is also a useful localization for any fixed $C>0$. If a choice of $\lambda,\mathbf r$ reaches the rank $q=q_0-CL$ (equivalently $\Pr(S\le q)\le1-p_N$, allowing equality), then
\[
 H=\sum_i(\lambda_i-1)^2=O_{B,C}(\sqrt{NL}),\qquad
 D=\sum_i\lambda_i(\mu-m_i)=O_{B,C}(L).                              \tag{18}
\]
In fact $F_*(q)=\tau+O_{B,C}(L^2/N)$ by its iid Edgeworth cdf expansion: the mean differs from $\mu N$ by $O(L)$, the variance by $O(\sqrt{NL})$, moving its standardized argument from $-z$ by $O(L/\sqrt N)$ with density $O(L/\sqrt N)$; the cubic cdf term there is $O(L^2/N)$. Thus the assumed tail has $\log[\Pr(S\le q)/F_*(q)]\le O_{B,C}(\rho L)$. Retaining the positive $A+P/2$ in the first line of (15), absorbing its smaller multiples and the mixed term as there, now yields $A+P=O_{B,C}(\rho L)$. Since $H\le CP/s^2$, this proves the first assertion. Also $\mu-m(c_i)\le C[(d_i')^2+u_i^2]$ uniformly: near $\mu$ use its nondegenerate maximum, $|c(u)-\mu|\le C u$, and away use boundedness. Equations (12) and $\sum u_i^3\le C(Ns^3+P^{3/2})$ give $sD\le C(A+\rho L+P^{3/2})=O(\rho L)$, proving the second. For a log-convex roof (18) implies that angular directions with $N^2|\nabla_Vu|^2>C_{B,C}\sqrt{NL}$ have radial halo volume at most $e^{q_0-CL}$ times their own angular source volume. Indeed their normals have no deterministic side vector reaching that threshold, and the support-plane inclusion in §2 applies on each such ray. This statement concerns actual receiver volume while allowing the action to depend on direction.
