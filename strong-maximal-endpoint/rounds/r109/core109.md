# Recovering the actual price of a macroscopic deep-core cover after truncating its load

Let $L=\log N$, $p=1-d_N/N-BL/N$, and use the original unit torus volume and actual rectangle dual of r096. Here $p$ is the unmodified strict success threshold; $q_0=\mu N-d_N$, and $\log Q_N=q_0+c_Q L+o(L)$. The new issue is whether the r102 loss of *every* uncut power certificate persists after deleting the common high-density region and pricing the remaining load. It does not, for the explicit covers below.

## 1. Two deterministic certificate tests

For a finite true family and a probability mixture of its normalized physical box densities write $W=\sum q_R1_R/|R|$. Every feasible $0\le f\le1$ has $\int Wf\ge p$; the actual indicator infimum with a strict margin agrees with the closed finite LP value $P_p$ by r096. For any source-space set $A$ one has
\[
 P_p\ \ge\ \frac{(p-\int_A W)_+^2}{\int_{A^c} W^2}.                              \tag{1}
\]
Use the estimate only when its denominator is positive; a zero denominator with positive numerator is incompatible with feasibility. Indeed the discarded region can supply at most $\int_AW$, and Cauchy-Schwarz on its physical complement uses $f^2\le f$. A nearly constant discarded mass on every box is sufficient, but controlling its mass under the **chosen actual mixture** is all that (1) requires. A deep region that absorbs mass $\ge p$ makes this test vacuous.

Here is a version that can be checked before superposing groups. Choose actual local mixtures $W_l$, group weights $\pi_l\ge0$ summing to one, and local discard sets $A_l$; let $B_l=W_l1_{A_l^c}$. In (1) take $A=\bigcup A_l$, $W=\sum\pi_l W_l$. A sufficient lower bound on the numerator before squaring and an upper on the denominator are
\[
 p-\sum_l\pi_l\int_{A_l}W_l-\sum_{l\ne j}\pi_l\int_{A_j}W_l,
 \qquad \sum_l\pi_l^2\int B_l^2+\sum_{l\ne j}\pi_l\pi_j\int B_l B_j .                \tag{2}
\]
The inequalities are by positivity on the *original* source space. Thus cross hits on other removed sets and pairwise overlaps of surviving loads must be counted. One cannot simply renormalize each row and ignore these two costs.

For completeness, if $S_l$ are nonnegative raw local loads (nonnegative linear combinations of actual box indicators) put $b_l=\min(1,S_l)$. The nonlinear dual applied once to $S=\sum S_l$ gives
\[
 P_p\ge\sum_l\left[p\!\int S_l-\int(S_l-1)_+\right]-\sum_{l<j}\int b_l b_j .        \tag{3}
\]
Indeed $\min(1,\sum S_l)\ge 1-\prod(1-b_l)\ge\sum b_l-\sum_{l<j}b_l b_j$. This charges only the clipped-load pair collision; an arbitrarily high raw spike need not be paid twice. The following construction computes all terms needed for (1) and (3).

## 2. The canonical physical family and its original price

For each sufficiently large $N$ there exists one **finite, rational, genuine rectangle family** $\mathcal R_N$ with covered volume $>1-1/N$, admitting all the uncut-power obstructions in r102, and satisfying
\[
 P_p(\mathcal R_N)=(1+o(1))(5L)e^{-u},\qquad
 u=s_p(S),\quad S=\sum_{i=1}^N Z_{c_N},\quad c_N=\mu-\tfrac{\sqrt\mu}{2}\sqrt{L/N}.    \tag{4}
\]
Here $u$ is the hit quantile of this single shifted iid reference, $\Pr(S>u)=p$. r092 gives $u=q_0+c_Q L+o_B(L)$. On this very same family, one discarded deep set and the sole remaining second moment in (1) already give
\[
 \text{certificate from (1)}\ \ge\ (4e^{-2}-o(1))(5L)e^{-u}.                          \tag{5}
\]
Thus the clipped quadratic certificate suffices through the log scale, whereas the best original Gram lower certificate loses a fixed positive multiple of $N$ in logarithm on the same examples. The full nonlinear certificate recovers the factor $1+o(1)$ in (4). This is a choice of the r102 construction retaining a finite discretization of the actual surface witnesses. It is not asserted for arbitrary subfamilies deleted after the construction.

### Actual receiver mixture on one packet

Use the r102 packet in a small cube $[-a,a]^N$ inside one torus chart, with $U(y)=\sum_i-\log(|y_i|/a)$ and $D_t=\{U>t\}$ restricted to the cube. Put $T=N^4$, $s=s_{p+2L/N}(S)-1$, $C=T-s$, and $H=\{U>C\}$ in the cube. Then $s=q_0+c_Q L+o_B(L)$, $u-s=O_B(1)$ by the shifted iid local density $f_S(u+O_B(L))\sim L/N$, and every box at a receiver with $U(x)=C$, opposite endpoint $-r_Nx$, $r_N=e^{-c_N}$, has true $D_T$-mean $\Pr(S>s)>p+2L/N$. Such boxes cover the entire interior of $H$ up to null faces: when $U(y)>C$ and every $0<|y_i|<a$, take $x$ with the same signs and log coordinates proportional to those of $y$, scaled down to sum $C$. Then $y$ is strictly inside the anchored box at $x$. The volume ratio $|H|/|D_T|=e^{s+O(Ns/T)}$ follows from the integer Gamma polynomial.

Put a concrete probability $q$ on the **genuine boxes**, taking independent uniform signs of $x$ and log magnitudes uniform on the simplex $t_i\ge0$, $\sum t_i=C$. Let $W$ be its mixture density. On a uniformly sampled point in its selected box the log vector is $v_i=t_i+Z_i$, with independent increments. The signs of the sampled point are uniform independent of $v$, since the receiver signs were independent uniform signs and the relative sampled signs do not depend on them. Thus for $V=\sum v_i$ the radial variable has precisely the law $C+S$ under the mixture. Its angular density admits the following pointwise majorant with respect to genuine physical volume:
\[
        W(y)\le K(V):=\frac{(N-1)!}{2^N a^N C^{N-1}} e^V f_S(V-C)\quad(v_i\ge0).      \tag{6}
\]
Outside the cube the density is zero; take $f_S(t)=0$ for $t<0$. To prove (6), the simplex law in the coordinates $t_1,\dots,t_{N-1}$ has density $(N-1)!/C^{N-1}$. In the convolution density at $v$ substitute $z_i=v_i-t_i$ in those coordinates, $z_N=V-C-\sum_{i<N}z_i$. Dropping only the constraints $t_i\ge0$ leaves the usual convolution giving $f_S(V-C)$. Divide by $2^N a^N e^{-V}$ for the signs and Jacobian. No angular prior is replaced by a Gaussian.

The majorant $K$ is nondecreasing in $V$: $e^v f_S(v)$ is proportional to $\sum_{j=0}^N\binom Nj(v-jc_N)_+^{N-1}$, by the shifted exponential/Gamma decomposition. Define
\[
 V_*=C+u,\quad A_0=D_{V_*},\quad A_+=D_{V_*+h},\quad f_*=f_S(u)\sim L/N,
 \quad K_*=K(V_*),                                                      \tag{7}
\]
where $h>0$ is fixed independently of dimension. $A_0$ is the right marginal price location; for the original r102 threshold slack it is shifted inward from $D_T$ by $u-s=O_B(1)$ in radial rank. The Gamma polynomial gives $K_*|A_0|=(1+O_B(N^2/T)) f_*$, and $|H|/|A_0|=e^u[1+O_B(N^2/T)]$. For the first identity compare $\sum_{j<N}V_*^j/j!$ to $V_*^{N-1}/(N-1)!$ with relative error $O(N/T)$, then $V_*^{N-1}/C^{N-1}=1+O(N^2/T)$. These estimates are uniform in the small chart scale $a$.

### One clipping followed by a finite moment

The mass of the discarded physical $A_+$ is *exactly the same on every surface box*, equal to
\[
 \alpha=\Pr(S>u+h)=p-(h+o(1)) f_* .
\]
This uses the uniform reference local density $f_S(u+O_B(L))/f_*=1+o(1)$ in the bounded shift. The remaining second moment is
\[
 \int_{A_+^c} W^2\le (e^h+o(1))\,\frac{f_*^2}{|A_0|}.                              \tag{8}
\]
Indeed squaring (6) and integrating over log-simplex sections yields the upper
\[
 \frac{(N-1)!e^{C}}{2^N a^N C^{N-1}}
 \int_0^{u+h}e^v f_S(v)^2(1+v/C)^{N-1}\,dv.
\]
The polynomial factor is $1+O_B(N^2/T)$ in this range. On $[u-2L,u+h]$ the density ratio to $f_*$ tends uniformly to one: the iid action is $\mu+O(\rho)$, standardized arguments differ by $O(L/\sqrt N)$, and the unstandardized density error from the matching Gaussian is $O(1/N)=o(f_*)$ by r093. On $[0,u-2L]$ the uniform density upper $O(1/\sqrt N)$ makes the remaining integral $O(e^{u-2L}/N)=o(e^u f_*^2)$. Multiplication by the prefactor and the identity following (7) prove (8). Therefore (1) already prices the single packet by at least $(h^2e^{-h}-o(1))|A_0|$, maximized at $h=2$. The logarithmic price recovery uses a true *truncated* moment; a full moment computed first and then multiplied by a clipping fraction does not give (8).

### Single-packet nonlinear calibration

For the same mixture use raw multiplier density $S_0=W/K_*$. By monotonicity of the majorant,
\[
 \int(W-K_*)_+\le \int_{A_0}K-K_*|A_0|,
 \qquad \int_{A_0}K=\int_u^\infty f_S(v)(1+v/C)^{N-1}dv=p+O_B(N^2/T).
\]
The weight error is at most $\mathbb E[e^{(N-1)S/C}-1]=O_B(N^2/T)$ by the uniform iid mgf near zero, so this estimate includes the entire far tail. It follows directly that $p\int S_0-\int(S_0-1)_+\ge |A_0|(1-O_B((N^2/T)/f_*))=(1-o(1))|A_0|$. This exhibits actual multipliers rather than inferring their existence from a separate capacity theorem.

### Finite boxes, real coverage, and simultaneous translations

The following limiting/discretization details are important because a surface of receivers alone has zero physical volume. For fixed $N,T,a$, the compact simplex parameter and finite signs map continuously into $L^2$ of normalized true box indicators on each sign component: side lengths are bounded away from zero for this finite problem, and small endpoint perturbations have small symmetric difference. Discretize the probability mixture to a finite convex combination $W^{(0)}$ of its true surface boxes within $N^{-4}$ in $L^2$ (hence also $L^1$ on the torus). Additionally include finitely many *zero-weight* genuine surface boxes whose interiors cover at least $(1-N^{-2})|H|$, possible by the interior covering above and inner regularity. The finite common constraints include both lists. This gives a fixed covered set $H'$ of volume $(1-o(1))|H|$ inside $H$; for definiteness choose an inner compact portion of that volume and its finite cover. Every retained surface row has the same stated original means on each radial core. The numerical count of boxes is not claimed to be uniform. One may take the chosen chart sufficiently small so $|H'|\le1/4$.

Write $M=\lceil5L/|H'|\rceil$. Independently translate the entire finite list $M$ times on the torus. Total raw $A_0$ price is
\[
             D=M |A_0|=(1+o(1))(5L)e^{-u},\qquad M/K_*=(1+o(1))D/f_* .                \tag{9}
\]
Every receiver in any translated $H'$ has a true physical witness. Expected uncovered fraction is at most $e^{-M|H'|}\le N^{-5}$, so union volume exceeds $1-1/N$ with probability tending to one. The union of translated $D_T$ is a feasible near-benchmark binary source as in r102; if its density lower calibration is needed, total raw volume is $M|D_T|=Q_N^{-1}N^{o(1)}$ and the expected raw-to-union loss is at most its square, so actual volume is asymptotic to the raw volume with probability tending to one.

For the sharper primal upper in (4), set $u_-=s_{p+3N^{-2}}(S)$ and use the common true union $\bigcup_l D_{C+u_-}^{(l)}$. Each untranslated surface box hits its own such core with probability $p+3N^{-2}$. We have $u-u_-=O(1/(NL))$ by the local lower density, hence its total raw volume is $(1+o(1))D$. The same source is feasible with strict margin even after rational enlargements of volume factor $\le1+N^{-3}$. Thus $P_p\le(1+o(1))D$.

To obtain the clipped quadratic lower in this *same* random experiment use uniform weights $1/M$ across translated copies of $W^{(0)}$. Its mass in its own $A_+$ is at most $\alpha+N^{-4}$; the exact unperturbed rows in fact give equality before rational approximation. For distinct independent translations, $\mathbb E\int_{A_j^+}W_l^{(0)}=|A_+|$. Therefore the expected cross contribution in the first line of (2) is $\le M|A_+|=O(D)=o(f_*^2)$. Markov makes it $o(f_*)$ with probability tending to one. For $B_l=W_l^{(0)}1_{(A_l^+)^c}$, independent translations have $\mathbb E\int B_lB_j=(\int B_1)^2\le1$. The total weighted cross second moment in (2) has expectation at most one, whereas the own contribution is bounded by $(e^h+o(1))f_*^2/D$ by (8) and the chosen $L^2$ approximation. Here $f_*^2/D\to\infty$ exponentially, so the cross contribution is $o(f_*^2/D)$ with probability tending to one. This yields (5) by (2) at $h=2$.

For the exact nonlinear lower take all raw loads $S_l=W_l^{(0)}/K_*$ in (3). The local objectives differ from the computed surface value by at most $N^{-4}/K_*=o(|A_0|)$ because the positive part is Lipschitz and the mixture has total mass one. The expected clipped cross interference is at most $\tfrac12(M/K_*)^2=o(D)$ using $D/f_*^2\to0$. Markov gives total interference $o(D)$ with probability tending to one, hence $P_p\ge(1-o(1))D$. Together with the preceding upper this proves (4). All coverage, feasibility, and collision estimates hold simultaneously for one fixed realization, by the union bound on their vanishing failure probabilities.

**Ensuring the same global clipping mass on *every* retained box.** This strengthening can be imposed in the same realization. Use $h=2$, $A=\bigcup A_l^+$ and $\varepsilon=N^{-3}$. Within the original finite list of one packet assign its certified compact $H'$ points measurably to one containing box, by first selection in the fixed list. After translation, discard any box whose normalized mass on the union of the *other* copies' $A_j^+$ exceeds $\varepsilon$. Conditionally on the translation of the box's own packet, this event has probability at most $M|A_+|/\varepsilon=O(D/\varepsilon)$ by independence, Fubini and Markov, uniformly over its box index. Hence the expected lost normalized mixture weight (uniform across packets) is $O(D/\varepsilon)$; the expected sum of lost assigned receiver volumes is $O(M|H'|D/\varepsilon)=O(LD/\varepsilon)$. The former is $o(f_*)$ and the latter $\le N^{-2}$ with probability tending to one, by Markov and the exponential smallness of $D$. Before removal the uncovered volume is $\le N^{-2}$ with probability tending to one as well. Thus the retained finite family still covers $>1-1/N$. Deleting these raw weights costs at most $p(M/K_*)\,o(f_*)=o(D)$ in (3) or the single global dual objective, since the subtracted positive part only decreases. For (1) normalize the remaining mixture to probability; its clipped second moment increases by at most the square of the $1+o(f_*)$ renormalization factor by pointwise positivity, and every remaining box now has global discarded mass between $\alpha$ and $\alpha+\varepsilon$. Later finite outward rational enlargement changes these masses by at most $O(\varepsilon)$. Thus uniformly on *all* final rows $\alpha_R=m_R1_A=\alpha+O(N^{-3})$. If one wishes to normalize each surviving conditional row, its lower threshold is exactly $1-(1-p)/(1-\alpha_R)=(p-\alpha_R)/(1-\alpha_R)$, equal to $(2+o(1))f_*/(1-p)$ here. The proof above also records the collision cost on that remaining part. This deletion does not affect the feasible upper source or the uncut r102 obstructions, which apply to every retained box individually.

Finally enlarge each of these finitely many translated boxes outward to a rational periodic box, preserving all receivers and keeping the enlargement factor $\le1+N^{-3}$. They may be chosen so close that the total convex mixture changes by $\le N^{-4}$ in $L^2$ and $L^1$, by finite-dimensional endpoint continuity. The original raw dual loads correspond to genuine nonnegative multipliers on the *new* boxes by multiplying each normalized box coefficient by its reciprocal volume. The nonlinear value changes by at most $(M/K_*)N^{-4}=o(D)$; in (1) the discarded-set mass changes by $o(f_*)$ and the surviving $L^2$ norm changes by $o(f_*/\sqrt D)$. Strict feasibility was already checked for the slight upper enlargement. If repeated boxes occur, merge their nonnegative multipliers. All members of this family still satisfy the r102 uncut entropy and fixed-order lower-mass tests, since each originated from a true chart/surface witness and was enlarged by at most $1+N^{-3}$. In particular its raw quadratic failure and its clipped quadratic recovery coexist on the same finite rational cover.

### Finite multilevel moment certification of the matching nonlinear lower

Even the $1-o(1)$ lower in (4) can be certified by finitely many *truncated first and second moments* per packet, with a growing number of layers. Here are the details, to contrast this statement with the one-cut constant $4/e^2$ in (5). Work initially with the exact surface mixture, write $S_0=W/K_*$, $v=V-C$, $H=\log L$. Declare the whole outside of the chart part of the low region (there $S_0=0$). Use three regions:

- On $v<u-H$ use $(S_0-1)_+\le S_0^2/4$. The proof of (8) at upper integration endpoint $u-H$ gives $\int_{v<u-H}S_0^2=O(e^{-H}|A_0|)$. In detail, the part below $u-2L$ divided by $e^u f_*^2$ is $O(1/(NL^2))=o(e^{-H})$ and on $[u-2L,u-H]$ the density ratio is $1+o(1)$.
- On $v>u+H$ use $(S_0-1)_+\le S_0$. Compared with integrating $K/K_*-1$ there, this loses at most its geometric volume $(1+o(1))e^{-H}|A_0|$ by the same Gamma polynomial ratio. This uses the actual first moment of the deep part; its uncut second moment need never be computed.
- Divide each side of $[u-H,u+H]$ at $u$ into equal bins of widths $\Delta\asymp L^{-1/2}$. Let $C_j$ denote the physical shell of a bin, $m_j=|C_j|^{-1}\int_{C_j}S_0$, $v_j=|C_j|^{-1}\int_{C_j}S_0^2-m_j^2$. The moment-only check
\[
 \int_{C_j}(S_0-1)_+\le |C_j|[(m_j-1)_++\sqrt{v_j}]
\]
holds by Lipschitz continuity and Cauchy-Schwarz. The true mass of a bin is $\int_{C_j}W=\int_{\rm bin}f_S(v)dv$. The majorant's first mass is the same integral weighted by $(1+v/C)^{N-1}=1+O(N^2/T)$. Furthermore $\max_{C_j}K/\min_{C_j}K=1+O_B(\Delta+1/L)$ on these bins: the Gaussian density with the fixed action's mean and variance changes by $1+O_B(\rho\Delta)$ within a bin, the absolute true-density error $O(1/N)$ has relative size $O(1/L)$, and the exponential factor changes by $e^{\Delta}$. It follows that $\int_{C_j}S_0^2\le |C_j|m_j^2[1+O_B(\Delta+1/L+N^2/T)]$, since $S_0^2\le(\sup_{C_j}K/K_*)S_0$. The summed standard-deviation errors are $O_B(L^{-1/4})\sum_{\rm bins}|C_j|m_j=O_B(HL^{-1/4})|A_0|$, using $\int_{u-H}^{u+H} f_S=O(H f_*)$. In bins below $u$ the mean-positive-part term is zero; in bins above $u$ it is at most the integral of $K/K_*-1$, since $K/K_*\ge1$ there and dominates $S_0$.

Altogether the computable upper on $\int(S_0-1)_+$ exceeds the majorant bound $\int_{A_0}(K/K_*-1)$ by at most $O_B(e^{-H}+H L^{-1/4})|A_0|$. This is $o(|A_0|)$. Thus $O(\sqrt L\log L)$ physical bins, their first two moments, one outside second moment and one inside first moment certify the same $(1-o(1))|A_0|$ local dual value. The majorant above is how these moments are bounded for this example, not an assumption that the true output is exactly radial. At each fixed $N$ the finite-bin bound is continuous in $L^2$ of the mixture (including bins of zero variance). Hence finite surface discretization and rational endpoint approximations in the proof may be taken sufficiently fine to keep its change $o(|A_0|)$; no uniform bound on box count is asserted. For distinct packet copies the interference in (3) is also bounded by their raw **cross** moment $\int S_lS_j$, whose expectation under independent translations is at most $1/K_*^2$. Their sum is $o(D)$ with probability tending to one even with this looser bound; the destructive large raw **self** second moments are never substituted. The previous small pruning loss is paid after obtaining this unpruned load objective. Thus a finite, location-resolved moment certificate recovers (4) as well as the one-cut bound (5).

## 3. A concrete conditional-price interface for unrestricted covers

The weights of several different local geometries need not be equal. Suppose one can select actual groups of boxes from a given cover, a normalized local mixture $W_l$ in each group, a local discard set $A_l$, a residual mass $g_l>0$, and a local target price $d_l>0$ with
\[
 p-\int_{A_l}W_l\ge g_l,\qquad \int B_l^2\le g_l^2/d_l,\qquad B_l=W_l1_{A_l^c}.
\]
Write $D=\sum_l d_l$, $J=\sum_l d_l/g_l$, and choose $\pi_l=d_l/(g_l J)$. In (2) the uncorrupted residual is at least $D/J$, and the own clipped second moment is at most $D/J^2$. In particular if the measured cross terms satisfy
\[
 \sum_{l\ne j}\pi_l\int_{A_j}W_l\le(1-\alpha)D/J,
 \qquad \sum_{l\ne j}\pi_l\pi_j\int B_l B_j\le (C-1)D/J^2,
 \quad 0<\alpha\le1,\ C\ge1 ,
\]
then the original common-source LP costs at least $\alpha^2 D/C$. This is a deterministic finite-moment test with different source locations and different conditional residual masses. It does not require bounding the Gram energy before clipping. To prove the sharp general upper by this route one would need such groups for **every** finite covering family at the macroscopic coverage of r095/r096, with $\alpha^2 D/C\ge Q_N^{-1}N^{-o(1)}$ (at fixed positive covered fraction). The displayed cross bounds are sufficient for this particular route; other proofs are not required to realize this decomposition. The exact necessary-and-sufficient dual condition remains the top-$p$-mass bound of r096 (7)–(8) for some original-box mixture in every such family. Independent translates and the radial formula (6) verify the needed costs for the constructed near-benchmark example. Neither property is available just from arbitrary rectangle coverage or from small raw source density.

**Local calculation.** For $B=0$ and $N=16,64,256$, local single-core Python used the actual finite binomial mixture of integer Gamma laws with the shifted $c_N$. At $h=2$ the rank $u$ was $16.85396317,72.40657342,303.64101337$ respectively. In the fixed-$N$, infinite-thickness scalar expression $[F_S(u+2)-F_S(u)]^2/\int_0^{u+2}e^{v-u}f_S(v)^2dv$, deterministic numerical quadrature gave $0.54568612,0.53618476,0.53805475$, compared with $4e^{-2}=0.54134113\ldots$. The integration used 400-panel Simpson quadrature above $\max(0,u-50)$; it is only a scale check. The actual physical bounds use the analytic estimates, not numerical quadrature or an assumed uniform angular law for the output.
