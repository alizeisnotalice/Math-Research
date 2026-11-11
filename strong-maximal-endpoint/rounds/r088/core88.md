# Exact logarithmic cost on an affine patch for translated reflection-center columns

Date: 2026-11-11. Strict patch theorem and a general sufficient residual estimate; independently audited.

## 1. Geometry and theorem

Fix $0\le r\le1$, including the shifted-reference ratio $r=r_N$. Denote by $A_{c,r}$ the genuine average on the rectangle whose opposite endpoints are $x$ and $(1+r)c-rx$. At each fixed $c$ this is a translate of the r028 fixed-ratio operator, with exactly the same cap-restricted weak norm. Conversely $\sup_c A_{c,r}=M_S$ since all rectangles having a receiver as a corner are obtained by $c=(a+rx)/(1+r)$, and splitting any rectangle at its receiver gives a convex combination of such corner averages. Degenerate coordinates are irrelevant to the averages. The r075 reanchoring identity changes the source with its receiver when all centers are moved to zero; our columns use actual common-source pieces $G\le F$. The r078/r081/r084/r087 column obstructions fix the physical endpoint $a$ itself. Those conclusions do not directly apply to this translated-operator family. In particular the r085 near-$Q_N$ reference shell is genuinely covered by its one matching center $c=0$ at cost $J_N(F)$.

Here is a different positive-volume receiver patch where a full pointwise budget for the new columns can be computed to within $e^{O(N)}$. The constants are uniform for all $r\in[0,1]$. Throughout $N\ge1024$, $\eta=N^{-2}$, $g(t)=1-(1-t)/N$ on $[0,1]$, zero outside, and
\[
 F(y)=e^{1/2}\prod_{i=1}^N g(y_i),\qquad 0\le F<2.
\]
Use the actual modular $J=\int\Phi_N(F)$. We have $1/2<J\le 2(1+\log2)^{N-1}$ by the cap/unit support and $\int F=e^{1/2}(1-1/(2N))^N$. For $z\in I=[1/4,3/4]$ define
\[
 v(z)=1+N(e^{-z/N}-1),\qquad x(z)=1+\frac{(1-v)^2}{2N g(v)},\qquad X=\{(x(z_i))_{i=1}^N:z\in I^N\}.
\]
Thus $1/4\le v<0.751$, $1<x<1+1/N$. The coordinate map $z\mapsto x$ has strictly positive continuous derivative on $I$. The pullback to $I^N$ of physical a.e. coverage is legitimate. Define $t=\log M_SF$ on $X$. For a target $0\le d\le1$ there, $B_X(d)$ is the infimum of $\sum_j w_jJ_N(G_j)/J$ such that $\sum_j w_j1_{\{A_{c_j,r}G_j>1\}}\ge d$ a.e. on $X$, over finite or countable physical columns, $w_j\ge0$, $0\le G_j\le F$.

**Theorem.** For the r081 flat cutoff $\chi$ and saturation demand $d_{\rm sat}=1_{\{M_SF\ge1+\eta\}}$ on $X$,
\[
 \log B_X(d_{\rm sat})=N\log N+O(N),\qquad
 \log B_X\!\left(\chi((M_SF-1)/\eta)\right)=\tfrac32N\log N+O(N). \tag{1}
\]
The upper bounds construct genuine columns on **the whole parameter cube** $X$ for the same fixed $F$, not only on a selected depth slice. Applied to a globally requested pointwise budget the estimates give lower bounds by restriction; the upper estimates here are patch estimates only. Dropping this whole patch costs $|X|\le N^{-N}\le2N^{-N}J$ in physical volume. Thus (1) is an obstruction to pointwise coverage and a calibration for flat-boundary gridding, not an obstruction to a physical residual allowance.

## 2. Fixed-center pricing for arbitrary true sub-sources

At $x=x(z)>1$, taking a one-dimensional rectangle outward from $x$ sees no source. Any inward interval of positive average is improved by ending at $x$ and clipping its left endpoint to $[0,1]$. For $b\in[0,1]$ its average is exactly
\[
 K_b(x)=g(v)-\frac{(b-v)^2}{2N(x-b)}.
\]
This follows by integrating the affine function from $b$ to 1 and using $\int_v^1g=(x-v)g(v)$. The unique optimum is $v$, so product factorization yields $t=1/2-\sum z_i/N$, and the maximizing reflection center per coordinate is
\[
 C(z)=\frac{v(z)+r x(z)}{1+r}.
\]
This center is quantitative: $|dv/dz|=g(v)$, while $|dx/dv|=(1-v)/(Ng(v))+(1-v)^2/(2N^2g(v)^2)<2/N$. Therefore
\[
 |C'(z)|=\frac{g(v)|1+r\,dx/dv|}{1+r}>1/3.\tag{2}
\]
It is monotone and the $N$-dimensional inverse Jacobian is at most $3^N$. This nonvanishing derivative is needed for moving-endpoint columns; the corresponding assertion about the optimal *physical endpoint* alone would not price them.

For $b=2$ or 3 form a full-dimensional depth slab $P_b\subset I^N$ as follows. Translate $u_i\in[3/8,5/8]$ by $z_i=u_i-2N^{-b}$, and impose $|\sum_i u_i-N/2|\le N^{1-b}/2$. Then
\[
 \tfrac32N^{-b}\le t\le\tfrac52N^{-b},\qquad |P_b|_z\ge4^{-N}N^{-b}.\tag{3}
\]
The sum of independent centered interval uniforms has symmetric unimodal density, by induction since a centered sliding-interval integral preserves evenness and monotonicity on the positive half-line. The mean density on its central interval is at least that on its whole support of length $N/4$, giving even a factor $4$ beyond the claimed bound. On $P_2$ the saturation demand is one. On $P_3$ the flat demand is at least $q_N=e^{-2N/3}/2$ by $\chi(s)\ge e^{-1/s}/2$ for $0<s<1$ and $e^t-1\ge t$.

Fix any center $c$ and actual $G\le F$. At success on $P_b$ every opposite endpoint $a_i=(1+r)c_i-rx_i$ satisfies $a_i<1$; otherwise both endpoints are exterior. Clipping to $a_i^+=\max(0,a_i)$ can only improve the average of $F$. Dividing the clipped average by $M_SF=e^t$ and using $-\log(1-s)\ge s$ on $[0,1)$ gives
\[
 \sum_i(a_i^+-v_i)^2<(2N+2)t\le6N^{1-b}.
\]
Every $v_i\ge1/4$ and $N\ge1024$ imply that all $a_i^+>0$, hence there was no clipping. Moreover $a_i<0.751+\sqrt{6/N}<7/8$, so each actual successful column has $J_N(G)\ge\int G>8^{-N}$. Finally $a_i-v_i=(1+r)(c_i-C(z_i))$. By (2) the $z$-volume that a single true success column captures in the slab is at most $3^NV_N(6N^{1-b})^{N/2}\le(C_1N^{-b/2})^N$, where $V_N$ is the Euclidean unit ball volume and the last step follows by a Gaussian integral. Integrate the almost-everywhere weighted cover on $P_b$ against $dz$. Tonelli and the lower cost per effective column yield $B_X\ge |P_b|_z q/[8^NJ(C_1N^{-b/2})^N]$, using $q=1$ or $q_N$ respectively. Equation (3) proves the two lower directions in (1).

## 3. Matching coverage of all positive depths

For $0<\delta\le1/4$, use a coordinate grid of centers in $[0,2]^N$ of mesh at most $\sqrt{\delta}/64$ (include endpoints). It has at most $(C_2/\sqrt{\delta})^N$ nodes and contains a node $c$ within that distance in each coordinate of $C(z)$. At the corresponding physical receiver the true endpoint $a_i=(1+r)c_i-rx_i$ then obeys $|a_i-v_i|\le\sqrt{\delta}/32$; hence $0<a_i<0.8$ and $x_i-a_i>0.2$. Since $g(v_i)>0.99$, the exact relative losses $l_i=(a_i-v_i)^2/[2N(x_i-a_i)g(v_i)]$ satisfy $\sum l_i<\delta/300$. By $\log(1-l_i)\ge-2l_i$ for $l_i\le1/2$, the full-source column $G=F$ satisfies
\[
 \log A_{c,r}F(x)\ge t-\delta/150>0\qquad(t\ge\delta/2).
\]
Each grid therefore covers all receivers whose margin exceeds the indicated threshold. In particular one grid with $\delta=\log(1+\eta)$ and unit weights covers the entire saturation target at total cost $\le (C_3N)^N J$. Boundary equality is included by the strict slack in the estimate.

For the full soft target set $\delta_j=2^{-j}/4$, $j\ge0$. On the layer $\delta_j/2<t\le\delta_j$ weight the corresponding grid uniformly by $q_j=\chi((e^{\delta_j}-1)/\eta)$; other nonnegative columns do not interfere. This covers the entire target, by monotonicity of $\chi$ and $-1/4\le t\le1/4$. Its normalized total cost is at most
\[
 C_2^N\sum_{j\ge0}\delta_j^{-N/2}q_j. \tag{4}
\]
For $\delta_j\ge\eta/4$ there are $O(1+\log N)$ terms, each $\le(4/\eta)^{N/2}$. For $\delta_j<\eta/4$, $s=(e^{\delta_j}-1)/\eta\le2\delta_j/\eta<1/2$, and the definition gives $q_j\le e^2e^{-\eta/(2\delta_j)}$. Writing $v_j=\eta/\delta_j$, the remaining sum is at most
\[
 e^2\eta^{-N/2}\sup_{v>0}[v^{N/2}e^{-v/4}]\sum_{v_j>4}e^{-v_j/4}
 \le C_4\eta^{-N/2}(2N/e)^{N/2}.
\]
The residual sum is bounded by a universal constant since $v_j$ doubles. With $\eta=N^{-2}$ this closes (1). The matching construction uses no source mass reconstruction, atomic-receiver weak pricing, or independent relabeling of source height.

For finite certificates one may use the sharper actual-modular bound $J(F)\le\int F+\int F^N\le3(21/20)^N$, proved exactly as in r084 (our height $e^{1/2}$ is smaller than its certified height). Keeping $\Phi_N$ throughout, the preceding constants yield the explicit conservative bounds
\[
\frac{(5N/5544)^N}{3N^2}\le B_X(d_{\rm sat})\le(184N)^N,
\qquad
\frac{(5N^{3/2}/11088)^N}{6N^3}\le B_X(\chi)\le(132N^{3/2})^N.\tag{4a}
\]
For the upper estimates take at most $(130/\sqrt\delta)^N$ nodes per grid. Since $\log(1+\eta)\ge\eta/2$ and $130\sqrt2<184$, the saturation constant follows. The coarse part of (4), divided by $N^{3N/2}$, is at most $(1+2\log_2N)(2/\sqrt N)^N\le1$ for $N\ge1024$. The fine part is at most $e^2(2/e)^{N/2}<8(3/4)^8<1$, using $8/3<e<11/4$ and $N\ge16$. Their sum is at most two, and $2\le(132/130)^N$ follows from Bernoulli for $N\ge65$. For the lower estimates use $|P_b|\ge4^{-N}N^{-b}$, inverse center Jacobian $3^N$, $\sqrt{12\pi e}<11$, the actual-modular upper bound above, and $q_N>2^{-(N+1)}$. Sharper finite lower values keep the exact factorial in the unit-ball formula. These bounds concern the patch budget, not $C_N$.

## 4. A general sufficient estimate for a physical residual

The following finite construction works for every nonnegative finite-modular input and does not use an unknown endpoint constant. Fix packets $0\le G_j\le F$, probability measures $\pi_j$ on physical center space, and integers $k_j\ge1$, $1\le j\le m$. At a receiver $x$ one may give a measurable choice of a genuine positive-side rectangle $R_j(x)=[x,a_j(x)]$ with $\fint_{R_j}G_j>1+\delta_j(x)$, $\delta_j(x)>0$; on any part where no witness is selected set $p_j(x)=0$. For a selected witness put $L_i=|a_{ji}-x_i|$, $d_i=\operatorname{sign}(a_{ji}-x_i)$, $c^*=(a_j+rx)/(1+r)$, and let $p_j(x)=\pi_j(C_j(x))$ for the set
\[
 C_j(x)=\left\{c^*+\frac{(d_i L_i u_i)_i}{1+r}:u_i\ge0,\ \prod_i(1+u_i)\le1+\delta_j(x)\right\}.
\]
Every such center retains the whole original rectangle inside its true averaging box, whose volume is at most $(1+\delta_j(x))|R_j(x)|$. Nonnegativity and the strict input margin therefore imply success. Its exact physical center volume is
\[
 |C_j(x)|=\frac{|R_j(x)|}{(1+r)^N}\int_0^{\log(1+\delta_j(x))}\frac{e^s s^{N-1}}{(N-1)!}\,ds.\tag{5}
\]
This follows from $z_i=\log(1+u_i)$ and the volume of simplex sections. For Borel witness data the center-incidence sets are jointly Borel by their inequalities, hence their $\pi_j$ probabilities are measurable. Rational far endpoints can be selected measurably when an unassigned strict witness is needed, by endpoint continuity and countable first selection. Using a particular full-volume witness with a weaker certified margin is also allowed.

For any integrable physical target $0\le\rho\le1$ there is a deterministic realization of $k_j$ centers for every $G_j$, all with unit weights, whose budget and deficit satisfy
\[
 B=\sum_{j=1}^m k_jJ_N(G_j),\qquad D\le\int\rho(x)\prod_{j=1}^m(1-p_j(x))^{k_j}\,dx.\tag{6}
\]
To see this draw the centers independently from their indicated proposals. At each $x$ any selected center in its certified set gives coverage at least one, so the residual is bounded by $\rho$ times the event that all certificates were missed. Independence and Fubini give (6) in expectation, and some realization achieves at most that bound. For example, $G_1=F$, $k_1=1$, $\pi_1$ a point mass at zero retains the whole canonical near-$Q_N$ reference shell using its actual original witnesses. More generally one may pay $B$ through the overlap multiplicities of physical packets $G_j=F1_{P_j}$. Success of an exact majorization outside these certified sets is not required. Formula (5) uses outward containment alone; it need not give the actual size of the success opportunity, particularly for a smooth optimized witness as in the affine example.

A general theorem would still need packet and proposal choices, for every source, simultaneous at budget $N^{o(1)}J_N(F)$ and physical residual $Q_NN^{o(1)}J_N(F)$. The affine-source computation and (6) do not supply those choices.

Local SymPy returned exactly zero for the residual in the displayed $dx/dv$ identity and verified the two-dimensional integral in (5) equals $(1+\delta)\log(1+\delta)-\delta$ before the physical prefactor. Local SciPy/NumPy evaluated the model tail weight $\delta_j^{-N/2}e^{-\eta/(2\delta_j)}$ on $j=0,\ldots,109$ at $N=1024,4096,16384,65536$. The maximizing indices were $28,34,40,46$, respectively, corresponding exactly to $\delta_j=\eta/N$ at these dimensions. These are checks of the analytic tail calculation, not numerical values of the covering infimum.

## 5. Reproducible verification

Python directed rational arithmetic encloses the full infinite dyadic construction cost: finite terms use stable log-sum-exp and the remainder has an explicit geometric bound. Small exponentials are bounded, not silently lost. Outward-center volume (5) is certified by a positive 81-term series with a rigorous remainder, with large rational logarithms compressed by a 192-bit integer envelope. Independent Wolfram checks use numerical logarithms from the start, bounded negligible terms and per-task memory guards. See [certification88.md](certification88.md), [validation88.json](validation88.json), and [audit88.md](audit88.md). These calculations do not choose packets/proposals for all inputs and do not close the general sharp upper bound.
