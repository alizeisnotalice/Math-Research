# Information required by macroscopic physical covers and a dual for finite periodic rectangle families

On the unit torus use the normalized Lebesgue measure and the genuine periodic rectangle averages of r095. The probability that a uniform receiver has a strict witness is its physical halo fraction $\theta$. Write $d=|E|$ in this file for the ordinary indicator source density; denote the r092 Gaussian correction by $d_N=\sigma\sqrt N z_N$. We are concerned with $p_N=1-d_N/N-B\log N/N$, $\log Q_N=\mu N-d_N+c_Q\log N+o(\log N)$. By r095 (6) a near-one macroscopic periodic coverage tests the unrestricted $H_N$ through its logarithmic coefficient. No scalar feedback lower bound is used as a physical norm.

## 1. Linear information about future receiver coordinates is necessary

Let $X=(X_1,\ldots,X_N)$ be a uniform periodic receiver, with independent uniform coordinates, and $A_i$ a finite positive-length periodic interval containing $X_i$. For a true receiver selection the whole list $A$ is chosen from $X$ (possibly with extra randomness) before sampling $Y_i$ independently and uniformly in their respective intervals. Periodic integration means uniform in the real representative and then reduction mod one. On any uncovered receivers use arbitrary valid default intervals. If a source $E$ of density $d$ has strict averages above $p$ for a receiver fraction at least $\theta$, the resulting law $P$ obeys $P(Y\in E)\ge a:=p\theta$.

Call a sequential reference *causal in the receiver coordinates* if it sees only $X_{\le i}$ in choosing side $A_i$. It may use previous sides, previously sampled points $Y_{<i}$, and independent private randomization. This sampled-coordinate feedback is deliberately allowed in the reference. The true operator may see the entire receiver $X$ when choosing every side, but it must choose a common box for its source integral; the two permissions are distinct.

For every such causal reference $Q$, the full output $Y$ has a density $w_Q$ and
\[
 \int w_Q^2\le3^N,\qquad Q(Y\in E)\le3^{N/2}\sqrt d.                \tag{1}
\]
Thus a receiver-causal selection serving fraction $\theta$ at threshold $p$ has $d\ge (p\theta)^2 3^{-N}$. This is an exponential loss relative to density $\exp(-\mu N+o(N))$ at any fixed positive macroscopic coverage and $p\to1$, because $\log3<\mu$. This concerns a single common source throughout the receiver space. There is no restriction on the shapes of $E$ or on the aspect ratios of its boxes.

For an actual physical selection define its information about the future receiver by
\[
 \mathcal I=\sum_{i=1}^N I_P(A_i;X_{>i}\mid X_{\le i},A_{<i}).        \tag{2}
\]
Information and logarithms are in nats. In particular a valid lower bound, whenever $b:=3^{N/2}\sqrt d<a$, is
\[
 \mathcal I\ge D_{\rm Ber}(a\|b)
       \ge a\,[\tfrac12\log(1/d)-\tfrac N2\log3]-\log2 .             \tag{3}
\]
For $\theta_N=1-1/N$, $p=p_N$ and any sequence $\log(1/d)\ge\mu N-d_N-O(\log N)$ that serves this fraction, **every** measurable physical witness selection satisfies
\[
 \liminf\mathcal I/N\ge(\mu-\log3)/2=0.0899261\ldots .               \tag{4}
\]
The inequality applies in each fixed permutation of the coordinate order, since $X$ has product measure and the proof does not select a privileged axis. Infinite information is permitted for continuous side encodings. There is a direct finite-message implication: if one public label $J$ selected from $X$ suffices to compute every physical side from $(J,X_{\le i})$ (also allowing extra conditionally generated randomization through the stated factorization), the Shannon entropy of this label is at least the right-hand side of (3). In particular a library of $M$ receiver-causal policies combined with an arbitrary full-receiver selector needs $\log M\ge((\mu-\log3)/2-o(1))N$. The genuinely unbounded family of arbitrary physical receiver-dependent boxes is not ruled out by this conclusion.

### The one-dimensional periodic calculation

For any one-dimensional source set $D$ of length $u$, the periodic measure of $\{M1_D>t\}$ is at most $\min(1,u(2/t-1))$ for $0<t<1$. On the line the sharp indicator estimate follows from the r011 §1 rising-sun inequality $t|H|\le\int_H1_D+\int(1_D-t)_+\le(2-t)|D|$; periodizing and truncating in large intervals with bounded witnesses gives the torus upper as in r095 §3. Equivalently one may work directly on the circle. The full periodic maximal output from an arc achieves this distribution: at an outside point at distance $s$ from an arc of length $u$, an interval of length $u+s$ sees the entire arc. This lower and the preceding upper give the distribution, without assumptions about longer intervals.

If $X_1$ is uniform and any legal interval is chosen from it, the sampled output $Y_1$ therefore obeys
\[
 P(Y_1\in D)\le\int M1_D\le G(u),\quad
 G(u)=\int_0^1\min(1,u(2/t-1))\,dt
     =u+2u\log\frac{1+u}{2u} .                                        \tag{5}
\]
Randomization independent before observing the fresh receiver coordinate does not affect the bound. This is the sharp envelope for a *single set test*, since an arc and its maximizing intervals attain it. It need not be an equality for all tested sets under one common interval policy.

Let $g=G'$ on $(0,1)$; it is nonnegative decreasing with integral one,
\[
 g(u)=2\log\frac{1+u}{2u}-\frac{1-u}{1+u},\qquad
 g'(u)=-\frac{2}{u(1+u)^2},\qquad \int_0^1g^2=3 .                   \tag{6}
\]
For the last identity integrate by parts using $u g(u)^2\to0$ at zero and $g(1)=0$. The result is $4\int_0^1g(u)/(1+u)^2du=4\int_0^{1/2}[2\log(1/(2v))-1+2v]dv=3$, with $v=u/(1+u)$. Equation (5) implies absolute continuity of the output and for its density $h$ implies $\int(h-t)_+\le\sup_{0\le u\le1}[G(u)-tu]=\int(g-t)_+$ for $t\ge0$. In particular $\int h^2\le\int g^2$ by integrating the positive parts. More generally convex power moments of order $r>1$ are bounded by $J_r=\int g^r$, using the corresponding positive-part integral.

### Tensorization under the causal reference; comparison with the actual selection

Under $Q$, conditional on $(X_{<i},A_{<i},Y_{<i})$, the next $X_i$ is uniform and independent. With this history fixed the next sampled coordinate is exactly a legitimate one-dimensional interval output under fresh uniform input, including any use of already observed sampled coordinates. After mixing over hidden history conditional only on $Y_{<i}$, the same (5) holds because it bounds each set for each fixed history. Standard Borel conditional densities therefore obey the same positive-part bounds almost everywhere; one can first use a countable generating algebra and extend the inequalities by regularity. Factoring the joint density into conditional densities and integrating backwards gives $\int w_Q^r\le J_r^N$ for every $r>1$. Hölder gives (1) at $r=2$. The same argument and the entropy chain rule give $D(Q_Y\|\mathrm{Unif})\le N\int g\log g$.

For the physical $P$, form $Q$ by sampling $X$ with the same uniform product law, but replacing each full-receiver side kernel by the conditional kernel $P(A_i\in\cdot\mid X_{\le i},A_{<i})$. Conditionals on unused histories can be chosen valid (replace any history not assigning total mass to intervals containing $X_i$ by a default containing interval). In both laws the source samples use the same uniform interval kernels. The chain rule for relative entropy gives exactly $D(P\|Q)=\mathcal I$, including the possibility $+\infty$; equivalently use the $(X,A)$ marginal first, since uniform conditional samples given the boxes agree. This $Q$ is a special causal reference of (1). Applying data processing to the event $Y\in E$ and using $P(Y\in E)\ge a>b\ge Q(Y\in E)$ proves (3) by the two monotonicities of binary KL. For a public message one can alternatively sample its original marginal independently of $X$ under the reference, then follow the same label-specific causal kernels (extend invalid histories as above). The path divergence is at most $I_P(J;X)\le H(J)$ if the extra causal kernels agree conditional on the label; data processing gives the same lower bound. For a deterministic finite label and deterministic side maps this requires no additional randomization assumptions. More generally $b$ in (3) may be reduced to $J_r^{N/r}d^{1-1/r}$ whenever this is smaller, by the same proof.

The available physical selector really can use unrestricted information about all receiver coordinates. Thus (3) rules out replacing it by a subexponential public library or by sample feedback with only causal receiver access at macroscopic coverage. It does not claim an information ceiling for the true selection, nor convert the one-dimensional rearrangement in (5) into a bound on the unconstrained endpoint coefficient.

### The optimized information exponent and its scope

Let $\Lambda(r)=\log\int_0^1 g^r$, $r>1$. There is a unique $r_*>1$ with $r_*\Lambda'(r_*)-\Lambda(r_*)=\mu$. Indeed this left side at one is $\int g\log g\le\log3<\mu$ by (6) and Jensen under $g(u)du$. Its derivative is $r\operatorname{Var}_{g^r/J_r}(\log g)>0$. The tilted distributions concentrate on $\{g>M\}$ for every fixed $M$ as $r\to\infty$; binary relative entropy against Lebesgue measure there shows $r\Lambda'-\Lambda\to\infty$. All derivatives and second moments for the powers used here are finite since $g$ grows only logarithmically at zero and vanishes linearly at one. Therefore
\[
 \kappa=\max_{r>1}\frac{(r-1)\mu-\Lambda(r)}r=\mu-\Lambda'(r_*) .                    \tag{4a}
\]
In the assumptions of (4), using $b=J_{r_*}^{N/r_*}d^{1-1/r_*}$ in (3) yields
\[
 \mathcal I\ge \kappa N-\big[\kappa+(r_*-1)/r_*\big]d_N-O_B(\log N),               \tag{4b}
\]
where the last constant may also depend on the fixed logarithmic slack of $\log(1/d)$. We used $p_N(1-1/N)=1-d_N/N+O_B(L/N)$ and $d_N^2/N=O(L)$. A lightweight local midpoint computation in $-\log u$ gave $r_*\approx1.550294$ and $\kappa\approx0.123243$; these floating values are separate from the exact variational constant.

This exponent is optimal for the enlarged reference class retaining just the one-coordinate Lorenz domination, which allows the product density $w=g^{\otimes N}$. To see its limiting calibration write $S=\sum\log g(u_i)$, $v=\Lambda'(r_*)$, and use the actual product tilt $(g^{r_*}/J_{r_*})^{\otimes N}$. On the event $C_N=\{|S-Nv|\le \sqrt N(\log N)^{1/4}\}$ its probability tends to one by Chebyshev. This window is $o(d_N)$, so its uniform volume also satisfies the density qualification in (4) eventually. Changing measure back to uniform gives $\log |C_N|=-N[r_*v-\Lambda(r_*)]+o(N)=-\mu N+o(N)$. Changing instead to $g^{\otimes N}$ gives $\log\int_{C_N}g^{\otimes N}=-N[(r_*-1)v-\Lambda(r_*)]+o(N)=-\kappa N+o(N)$. Hence the linear exponent in the Hölder small-set estimate for this envelope cannot be increased at rate $\mu$. This calibration does not prove joint realization of $g$ by any physical interval policy and does not claim sharp physical lookahead information.

For $r>1$, integration by parts in $\int g^r$ as in (6) also gives $J_r=r\int_0^1[-2\log t-1+t]^{r-1}dt$. In particular $J_3=16$ exactly. A rigorous finite elementary upper at $r=3/2$ comes from Jensen on each subinterval of the equal $m$-partition:
\[
 J_{3/2}\le\frac{3}{2m}\sum_{j=0}^{m-1}\sqrt{m[F((j+1)/m)-F(j/m)]},\qquad
 F(x)=-2x\log x+x+x^2/2,\quad F(0)=0.
\]
The right side with $m=16$ is approximately $1.58101346$. The analytic bound holds before numerical evaluation and provides a finite certificate route for a stronger explicit coefficient than (4).

The formal 56-case rational computation in `information96.py`, recorded in `information/information96.json`, certifies $\kappa>0.12324$; its strongest row is $r=31/20$ on the squared 1024-cell partition, with lower bound approximately $0.12324068221332726$. The certificate uses directed logarithm intervals, a rational lower bound on $\mu$, and verified integer-power upper inequalities for the concave-power Jensen terms. This numerical certification of a strict lower bound is independent of the preceding approximate quadrature for the variational optimum.

### Two coordinates may be revealed together

As a self-contained strengthening, split the receiver into consecutive blocks of one or two coordinates. The whole current receiver block may be known before its arbitrary true block rectangle is chosen; earlier receivers, choices and source samples may also be used. The next entire receiver block is still fresh uniform. Such a selection with macroscopic coverage $\theta_N\to1$ at $p_N$ satisfies
\[
  \log(1/d)\le [2\log2-(1-\log2)/2+o(1)]N < (1.233+o(1))N.        \tag{4c}
\]
This excludes the leading rate $\mu$ even with genuine two-coordinate within-block lookahead.

For a two-coordinate source of volume $u$, the mass sent by one chosen box, integrated over its uniform receiver, is at most $\int_0^1\min(1,uH_2(t))dt$, by level integration and torus transference. r011 gives $H_2(t)\le\exp[s_t(L_1+L_2)]$ where the $L_i$ are independent positive variables with $\Pr(L_i>s)=2/(1+e^s)$ for $s\ge0$. Indeed that is the one-step maximum hit tail of an increment at threshold $s$, so conditional quantile coupling dominates the full two-step hit law by $L_1+L_2$. The indicator-to-rank bound of r011 then gives the asserted inequality in the weak-constant direction. Thus the integrated upper is bounded by $G_2(u)=\Pr(T_2>-\log u)$ for $T_2=L_1+L_2+Z$, with $Z$ an independent unit exponential; condition on the quantile of $L_1+L_2$ to evaluate the integral. Its derivative $g_2(u)=e^s f_{T_2}(s)$, $s=-\log u$, is nonincreasing since $e^s f_{T_2}(s)=\int_0^s e^v f_{L_1+L_2}(v)dv$. As in (5) this bounds the Lorenz curve of every fresh-block output and gives entropy at most $\int g_2\log g_2=\mathbb ET_2-h(T_2)$. All integrals are finite: $g_2(e^{-s})\le 2s^2$ by bounding each $L_i$ density by $2e^{-t}$.

Tail integration gives $\mathbb E L_i=2\log2$. With $U=\tanh(L_i/2)$ uniform on $(0,1)$, the density at $L_i$ is $(1-U^2)/2$, giving $h(L_i)=2-\log2$. Independence and conditioning give $h(T_2)\ge h(T_2\mid L_2,Z)=h(L_1)$. The two-coordinate entropy per coordinate is consequently at most $2\log2-(1-\log2)/2=1.23286795\ldots<\mu$. The one-coordinate rate is at most $\log3$ from (6). The chain rule bounds total output entropy relative to uniform by the sum of these block bounds. Binary data processing on $E$ gives $(p_N\theta_N)\log(1/d)-\log2$ as a lower bound for that entropy, proving (4c). There is no claim for a selector that reads unopened receiver blocks.

## 2. A finite-family nonanticipative dual for the continuous periodic density problem

There is a distinct route which retains the full physical receiver dependence of the chosen boxes. It gives exact finite dual certificates and clarifies the quantifier required of a continuum upper bound. A periodic rectangle below is a product of arcs of length at most one (full circles are allowed). This suffices for every strict periodic witness: split at the receiver into corner subrectangles, retain a strict piece, then in any longer coordinate split its corner interval into complete periods and the remaining interval adjacent to the receiver modulo one. The average is a convex combination of these two averages, so one of them preserves strictness. Iterate in the coordinates. Nonfull arcs can then be enlarged by arbitrarily small amounts to rational-endpoint arcs making the receiver interior, still preserving a strict average. Full-circle coordinates need no change.

For a finite family $\mathcal R=(R_j)_{j\le m}$ of positive-volume rational rectangles set
\[
 P_p(\mathcal R)=\min_{0\le f\le1}\left\{\int f:\int_{R_j}f\ge p|R_j|\ (1\le j\le m)\right\}.
\]
With the continuous torus density $d_N(p,\theta)$ of r095 (6), one has the exact identity
\[
 d_N(p,\theta)=\inf_{\substack{\mathcal R\ \mathrm{finite\ rational}\\ |\cup_jR_j|>\theta}}
 \ \max_{\lambda_j\ge0}\left[p\sum_j\lambda_j|R_j|
                    -\int\left(\sum_j\lambda_j1_{R_j}-1\right)_+\right].             \tag{7}
\]
In particular the outer infimum is over *covers of actual receiver volume*. It does not permit the side on one coordinate of $R_j$ to change when sampling another coordinate of that same rectangle.

To prove (7), the physical infimum is unchanged by replacing halo coverage $\ge\theta$ by $>\theta$. If equality holds, put into the source an arbitrarily small positive-measure subset of its halo complement; differentiation gives strict success a.e. on the added subset. For a physical strict halo exceeding $\theta$, its successful rational rectangles have a finite subfamily whose union already exceeds $\theta$. Its original source is feasible for $P_p$, proving one direction. In the other direction, partition the torus into the finitely many atoms of the rectangle faces. Cell averaging preserves every integral and total price, so $P_p$ is an attained finite-dimensional LP. For an optimizer, replace $f$ by $(1-\epsilon)f+\epsilon$ to get *strict* constraints and at most $\epsilon$ additional price. In each cell take an ordinary indicator subset of the prescribed cell mass. That one source now succeeds on all of each $R_j$, and $\epsilon\downarrow0$ gives the reverse inequality. Taking nonnegative Lagrange multipliers on its rows and minimizing the linear integrand on $[0,1]$ gives the objective in (7). Feasibility with $f=1$ and ordinary finite LP strong duality prove maximization is exact. These steps extend the fixed-input/fixed-group allocation calculations of r077 to the outer periodic receiver-covering problem via r095's density criterion.

In a geometric form let $q\in\Delta_m$ be a probability vector and $W_q=\sum_j q_j1_{R_j}/|R_j|$, a mixture of the **original true box** laws. Define $a_p(w)=\min\{\int f:0\le f\le1,\ \int fw\ge p\}$; this is the amount of highest-density volume needed to collect mass $p$, allowing a fraction on a constant level. Then
\[
 P_p(\mathcal R)=\max_{q\in\Delta_m}a_p(W_q),\qquad
 P_p(\mathcal R)\ge \frac{p^2}{\min_{q\in\Delta_m}q^T K q},\quad
 K_{ij}=\frac{|R_i\cap R_j|}{|R_i||R_j|}.                             \tag{8}
\]
To check the exact first formula, write $\lambda_j=tq_j/|R_j|$ and optimize the common $t\ge0$. The one-row LP dual gives $a_p(W_q)=\max_{t\ge0}[pt-\int(tW_q-1)_+]$. The second claim follows from Cauchy-Schwarz and $f^2\le f$. The matrix $K$ is a positive semidefinite Gram matrix of physical densities. Rational coordinates yield exact rational entries; any nonnegative rational $q$ summing to one gives a rigorous computable lower certificate. The full dual in (7) also provides a linear rational certificate for rational $p$. A minimizer $q_*$ of the Gram energy satisfies $(Kq_*)_j\ge q_*^TKq_*$ for every box, with equality for every positive $q_{*,j}$, by testing movement to individual vertices. Thus a family feasible for a true source of density $d$ forces *every* mixture to have $q^TKq\ge p^2/d$ and even has one certified original-box mixture with average at least $p^2/d$ on every box. The collision in (8) is between two distinct receiver choices, with their genuine normalized intersection in physical space; no diagonal weak estimate in a product of receiver copies is used.

For example on the unit circle take $R_1=[0,2/3]$, $R_2=[1/3,1]$, $p=3/4$. The primal puts the full shared third into the source and adds $1/6$ from each outside third. The dual $\lambda_1=\lambda_2=1$ gives $P_p=2/3$. The symmetric mixture $q=(1/2,1/2)$ has $\int W_q^2=9/8$, giving only $1/2$ from the quadratic relaxation. The original LP distinguishes these prices; strict feasibility is recovered by an arbitrarily small increase as above.

**What a grid certificate can and cannot settle.** On a finite catalog the outer minimization over subsets covering $>\theta$ still has to be handled, for example by enumeration or a certified disjunctive search with per-family certificates (7)–(8). A near-optimal sampled source alone bounds the density from above, hence $H_N$ from below via r095; it does not certify a continuous upper. If actual witnesses are known to have each side at least $\ell_i>0$, rounding both arc endpoints outward to mesh $h_i$ increases volume by at most $D=\prod_i(1+2h_i/\ell_i)$ and retains each receiver and all the former source mass. For their actual strict density $>p$ the enlarged grid boxes are strictly feasible at $p/D$. A lower certificate for **every** qualifying grid cover at the weakened threshold therefore also bounds the physical source cost for this side-bounded subfamily. A fixed unrestricted source may have arbitrarily small useful sides; in the full problem (7) equivalently takes the decreasing limit of the finite dyadic catalog minima as the mesh tends to zero, using strict source perturbation and outward approximation. No lower bound uniform over all catalogs is proved here.

The lookahead result connects directly to these certificates when a cover contains a causal assignment of its actual boxes to $\theta$ of all receivers. Extend the assignment by a valid causal policy on all receivers if such an extension exists (equivalently require a policy whose covered outputs choose members of the family on that fraction). Weight its successful boxes by their normalized frequencies. For a physical assignment take its genuine independent uniform samples conditional on the full selected box. If the causal assignment uses adaptive sampled feedback, generate its box list using an independent auxiliary sampling history, then draw a fresh independent physical sample from every chosen coordinate interval. Equivalently the latter draw can be made online while the auxiliary history controls subsequent choices. The $(X,A)$ law and the successful assignment frequencies are unchanged, the new output is still from a receiver-causal reference, and it is conditionally product uniform given the final box. For this reference the resulting $W_q$ is bounded by $w_Q/\theta$ because it restricts the nonnegative box mixture to the covered assignments; (1) then gives $P_p(\mathcal R)\ge p^2\theta^2/3^N$ from (8). Arbitrary covers need not admit such an assignment. In particular causally recombining coordinates of disjoint valid witnesses generally produces a new rectangle outside the family, and such mixtures cannot be substituted in the exact dual.

**Local numerical check.** Direct midpoint quadrature in $t=-\log u$ over $0<t<40$ with 40000 cells gave $\int g=1.00000002$, $\int g\log g=0.67815971$, and variance $1.03932411$ for $\log g$ under density $g$. A separate 20000-cell computation gave $\int g^2=2.999999999999946$; (6) proves the value exactly. A numerical power check at $r=1.5$ gave $J_{1.5}\approx1.57663$ and a prospective information coefficient $[(r-1)\mu-\log J_r]/r\approx0.12263$. The explicit certified coefficient in (4) uses the exact second moment, not this floating quadrature.

The remaining full upper requires a density certificate for every macroscopic cover in (7), including those selected with extensive lookahead about actual receivers. Neither the Gram lower bound without a universally controlled mixture nor the information necessary condition alone proves that closure.
