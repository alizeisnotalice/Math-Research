# Near-optimal indicators can lose a prescribed polynomial under one Steiner compression

Use the notation of r092–r094: $L=\log N$, $\tau=d/N=\sigma z_N/\sqrt N$, $p=p_N=1-\tau-BL/N$ for fixed $B$, $q_0=\mu N-d$, and $\log Q_N=q_0+c_Q L+o(L)$. All scores below are genuine strict strong indicator halo volumes divided by original indicator volume. The construction preserves this denominator exactly under one-coordinate fiber compression.

## 1. Full-volume compression counterexample at near-$Q_N$ efficiency

Let $\delta=L/N$ and $n=N-1$. For large $N$ there is $h_N$ with $\log h_N=q_0+c_Q L+o_B(L)$ and the following property. For every integer $1\le m\le c\delta h_N$, with an absolute small $c>0$, there is a bounded measurable indicator source $E\subset\mathbb R\times\mathbb R^n$ such that its first-coordinate Steiner compression $E^*$ satisfies
\[
 \frac{|\{M_S1_E>p\}|}{|E|}\asymp h_N,\qquad
 \frac{|\{M_S1_{E^*}>p\}|}{|E^*|}\asymp h_N/m,
 \qquad |E^*|=|E| .                                                     \tag{1}
\]
All implicit comparison factors are absolute for sufficiently large $N$. Each nonempty first-coordinate fiber is a single interval of length one both before and after compression. Taking $m=\lfloor N^a\rfloor$ gives a log-score loss $aL+O(1)$ for any fixed $a>0$. Taking $m=\lfloor c\delta h_N\rfloor$ leaves the compressed score comparable just to $N/L$, while the original full score remains $Q_N N^{o(1)}$. Thus even restricting to physical inputs which attain the known benchmark through its logarithmic coefficient does not permit a universal $N^{o(1)}$ charge for a single coordinate compression. This does not assert that the uncompressed score exceeds the benchmark by a power.

### Choosing a genuine physical packet

Take a bounded r093 equal-weight hyperbolic cross $D$ in dimension $n$. Choose thickness so that its complete physical indicator halo $H$ at $p+2\delta$ has $h=|H|/|D|$ between half the limiting r093 price and that price. The density $p+2\delta$ in dimension $n$ is of the critical form with a bounded additive $O(L_n/n)$ perturbation, bracketable by fixed r092 coefficients. Also $q_{0,n}=q_{0,N}-\mu+o(1)$, as in r094 §1. Thus $\log h=q_0+c_QL+o_B(L)$; take $h_N=h$. All $y\in H$ admit a real rectangular witness of mean strictly larger than $p+2\delta$. There is a finite bound $T_0$ on each side length of such witnesses: in any coordinate they cannot exceed the corresponding width of a box containing $D$ divided by $p+2\delta$. Hence $H$ is bounded. Rectangles can be chosen measurably by rational outward perturbation and first selection, increasing $T_0$ by a fixed factor when needed.

### Packing colors without losing their strict witnesses

Choose a torus with side $S$ and volume $V=S^n$, sufficiently large to contain translated copies of $D,H$ within one chart and so $|H|/V\le1/4$. In each of $m$ colors put $M=\lceil V/|H|\rceil$ independent uniform translates of $D$ on the torus; all translates are mutually independent. Periodize to $\mathbb R^n$. Write $D_j$ for the union in color $j$, $A_j=D_j\setminus\bigcup_{i\ne j}D_i$. These $A_j$ are measurable and pairwise disjoint. Put $u=M|D|/V$, $\lambda=mu$; then $1/h\le u\le5/(4h)$, $\lambda\le(5/4)c\delta$.

For any torus point $y$, the probability it belongs to some translated $H$ of color $j$ is $1-(1-|H|/V)^M\ge1-e^{-1}$. Condition on that color and pick the first such translate and its bounded witness, lifting periodically. At any spatial point the probability of occupancy by another color is at most $\lambda$. Fubini on the chosen witness bounds the expected removed fraction by $\lambda$; this conditioning uses only the selected color, not the independent other colors. Markov shows the loss fraction exceeds $\delta$ with probability at most $\lambda/\delta$. Therefore, for sufficiently small absolute $c$, the expected sum of the torus volumes of color-success sets $G_j$, each with an actual witness of all sides $\le T_0$ and $A_j$-mean $>p+\delta$, is at least $0.6mV$. They may equivalently be defined by the successful bounded translated certificates, or by a countable family of short rational witnesses with $T_0$ enlarged. Since the sum is at most $mV$, it exceeds $0.3mV$ with probability at least $0.3/0.7$.

The integral of the raw count of all copies is exactly $\lambda V$. The loss from this count down to the surviving disjoint indicator total $\sum1_{A_j}$ is at most the sum of ordered pairwise intersections over distinct copies: if $l\ge2$ copies cover a point its loss is at most $l\le l(l-1)$. The expected integral of these intersections is at most $\lambda^2 V$ by independence. Markov bounds the probability of $\sum|A_j|<3\lambda V/4$ by $4\lambda=o(1)$, also when $m$ varies in the stated range. We can consequently fix a single realization with $\sum|G_j|\ge0.3mV$ and $3\lambda V/4\le\sum|A_j|\le\lambda V$.

### Preventing escape of the complete halo after finite truncation

Let $\eta=1-p$, choose a multiple $P$ of $S$ with $P\ge2nT_0/\eta$, and remove the same periodic set from every color:
\[
 U_\phi=\{y:\{\phi+\textstyle\sum_{i=1}^n y_i/P\}\in[0,2\eta]\},\qquad 0\le\phi<1.
\]
For a uniformly random phase the expected fraction of removed source mass equals $2\eta$. If the phase of a receiver avoids this arc expanded on both sides by $nT_0/P$, its complete chosen short witness avoids the removed set. On a $P$-period cell, the expected loss from the sum of certified receiver volumes is thus at most $3\eta mP^n$. Markov and $\eta\to0$ allow one fixed phase to retain summed receiver volume $\ge0.2mP^n$ and total source volume at least $\lambda P^n/2$ on that cell, simultaneously; the upper source bound $\lambda P^n$ is automatic. This diagonal removal costs only one small fraction, despite controlling long intervals in *every* suffix coordinate.

Take a cube $[0,R]^n$ with $R$ an arbitrarily large integer multiple of $P$. Let $A'_j=(A_j\setminus U_\phi)\cap[0,R]^n$. Their total volume lies in $[\lambda R^n/2,\lambda R^n]$ exactly by periodicity of the surviving source. Choose $R$ so large that discarding the $T_0$-thick boundary layer removes at most $0.05mR^n$ from the summed certified receiver volume. Then a sum of at least $0.15mR^n$ of transverse receiver volumes remains, with honest short witnesses of $A'_j$-mean $>p$. All of this concerns actual finite inputs.

For *any* rectangle, in any one suffix coordinate, an interval of length greater than $T_1=4P/\eta$ has removed fraction at least $2\eta-2P/T_1>\eta$ along every transverse line; at most two endpoint period fragments cause the stated error. Therefore a positive-volume rectangle witnessing any subset of $\mathbb R\times([0,R]^n\setminus U_\phi)$ at threshold $p$ has every suffix side $\le T_1$. All its receivers must lie over $[-T_1,R+T_1]^n$. Increasing $R$ so that $(1+2T_1/R)^n\le2$ bounds this entire transverse receiver volume by $2R^n$. The argument works for arbitrarily eccentric boxes, including receivers outside the original source cube. It avoids bounding the full halo by a saturated torus computation alone.

### Original and compressed volumes

Take length-one disjoint horizontal bands $I_j=[3j,3j+1]$, $E=\bigcup_{j=1}^m I_j\times A'_j$. Each nonempty fiber belongs to exactly one band, so $E^*=I_0\times\bigcup A'_j$ with $I_0=[-1/2,1/2]$ and both source volumes equal $\sum|A'_j|$. For each certified suffix receiver in color $j$ every point of $I_j$ succeeds by taking that interval and its true suffix witness. Thus $|\{M_S1_E>p\}|\ge0.15mR^n$. The original horizontal support is within an interval of length at most $3m$; its one-dimensional indicator halo has length exactly that length times $2/p-1$. The complete transverse bound above gives $|\{M_S1_E>p\}|\le6m(2/p-1)R^n$. After compression the analogous upper is $2(2/p-1)R^n$. The compressed halo also contains $I_0$ over the union of the certified suffix receivers, which has volume at least $0.15R^n$ by the summed lower bound divided by $m$. These are bounds on the full Lebesgue volumes at the same strict threshold. Dividing by the common price in $[\lambda R^n/2,\lambda R^n]$ proves (1).

**Completing all coordinate compressions.** With the saturated choice $m=\lfloor c\delta h_N\rfloor$, the once-compressed true source occupies a fraction $\asymp\delta$ of a box of side lengths $1,R,\ldots,R$. Any further coordinate Steiner compression preserves its volume and keeps it inside a translate of a box of exactly those same side lengths: individual slice lengths are at most the containing box side, and the other coordinate ranges are unchanged. This statement persists through any number of further compressions. For any axis-parallel box $B_0$ with side widths $\ell_i$, the exact product of its one-dimensional maximal indicators implies that $\{M_S1_{B_0}>p\}$ is within $\sum_i\operatorname{dist}(x_i,I_i)/\ell_i<1/p-1$. Splitting by the $j$ outside coordinates and integrating simplex sections shows
\[
 \frac{|\{M_S1_{B_0}>p\}|}{|B_0|}\le\sum_{j=0}^N {N\choose j}\frac{[2(1/p-1)]^j}{j!}
 \le \exp(2\sqrt{2N(1/p-1)}).
\]
The last estimate uses $\sum x^j/(j!)^2\le(\sum (\sqrt x)^j/j!)^2$. Thus the score after *all* coordinate compressions, in any order following the first, is bounded by $C\delta^{-1}\exp(O(N^{1/4}L^{1/4}))=\exp(o(N))$. The original score has logarithm $\mu N-d+c_Q L+o(L)$. For any prescribed order relabel the encoded coordinate to be its first. The full-compression conclusion uses actual volume density in the bounding box, not persistence of a periodic guard after later symmetrizations; it does not assert that every possible ordering on one fixed labeled example has this behavior.

## 2. True full-volume limits on multicenter block steering

The counterexample can also be written as a single arbitrarily sheared fiber product: at suffix point $y\in A=\bigcup A'_j$, put the column $I_0+b(y)$ where $b(y)$ selects the center of its color. Its price is exactly $|A|$. Thus geometry which merely assumes that fibers are translates, without controlling where they go, does not restore the comparison to the full halo of the *same centered source*. The following three price statements specify what translated columns can and cannot realize.

**Separated one-coordinate centers.** For arbitrary pairwise disjoint measurable suffix sets $A_j$ of finite total volume, arbitrary length-one horizontal intervals $I_j$ with mutual gaps at least two, and $p>1/2$, one has the exact full formula
\[
 |\{M_S1_{\bigcup I_j\times A_j}>p\}|
 =\sum_j\int_{\{b_j>p\}}(2b_j(y)/p-1)\,dy,\qquad b_j=M_{S,n}1_{A_j}.   \tag{2}
\]
Indeed each row has length at most one. A true witness's horizontal side has length strictly below $1/p<2$, so it cannot see two different bands in positive measure. The original halo is exactly the union of the individual product halos. Those halos are horizontally disjoint because the length-one interval's $p$ halo extends only $(1/p-1)<1$ on either side. Factorization gives the formula, with all actual suffix rectangles still free. Consequently the true global score is between the $|A_j|$-weighted average of the individual transverse scores at $p$ and $(2/p-1)$ times this average. This applies to §1 itself. The full $m$-fold difference there measures lost receiver volume after aligning channels. It cannot be interpreted as an $m$-fold improvement over the component weak scores before alignment.

**Unrestricted center locations for a box of fixed side lengths.** In $N=n+k$ dimensions let
\[
 F(x,y)=1_D(y)\,1_{b(y)+[-1/2,1/2]^k}(x)
\]
with arbitrary measurable $b:D\to\mathbb R^k$, and $0<|D|<\infty$. The honest source volume is $|D|$. Positive coordinate scalings give the identical assertion for any fixed axis box. Without a bound on the number, speed of variation, or separation of the centers,
\[
 \frac{|\{M_{S,N}F>p\}|}{|D|}\le(1+2/p)^k H_n(p).                    \tag{3}
\]
In fact every side length of the first block in any strict witness is $<1/p$: integrating along that coordinate, even with all other samples fixed, sees occupied length at most one. At receiver $(x,y)$, define $D_x=\{v\in D:|b_i(v)-x_i|\le1/2+1/p\ (1\le i\le k)\}$. Only suffix rows in this set can supply positive measure to the witness. Hence its *same* suffix rectangle through $y$ has $D_x$-density strictly greater than $p$. Integrate $|\{M_{S,n}1_{D_x}>p\}|\le H_n(p)|D_x|$ over $x$; Fubini gives $\int|D_x|dx=(1+2/p)^k|D|$. This argument permits fibers to encode arbitrary suffix samples by moving the centers, while retaining the one physical preceding box at each chosen receiver. All receiver points are included; no one-witness simulation has been substituted for the global halo.

For $p=p_N$, $1\le k\le N/2$, the r091 one-indicator upper (6) in dimension $n=N-k$ and $1-p_N+2n^{-4}=\tau_N+O_B(L/N)$ give
\[
 \log H_n(p_N)\le \mu(N-k)-\sigma\sqrt{N-k}\,z_N+C_H L+O_B(1),
\]
where any change $C_H\log n-C_H L$ is bounded in this range. The inverse Gaussian changes from $-z_N$ by $O_B(1/\sqrt N)$, by the density $\phi(z_N)\asymp L/\sqrt N$ and uniformity on such a short standardized interval. The r021 variance shift contributes $O(L^{3/2}/\sqrt N)$. Thus every source of this shape satisfies
\[
 \log\frac{|\{M_{S,N}F>p_N\}|}{|D|Q_N}
 \le-(\mu-\log3)k+(C_H-c_Q)L+O_B(k\sqrt{L/N})+o(L).                \tag{4}
\]
For $k/\log N\to\infty$ the relative loss is at least $\exp[(\mu-\log3-o(1))k]$ in inverse score. In particular *no* number of center channels restores near-$Q_N$ performance when more than logarithmically many, and at most half, of the coordinates have such translated fixed-width fiber shape. With $k=tL+O(1)$, (4) gives a deficit of any specified power by taking fixed $t$ sufficiently large compared to $C_H-c_Q$. This conclusion does not hold by the same proof for arbitrary varying fiber widths.

**Variable hyperbolic widths, with bounded relative center displacement.** Let the centered fiber on a suffix set $D$ be the box $\prod_{i\le k}[-a_i(y),a_i(y)]$, $0<a_i(y)<\infty$, with finite total volume, and let $E_b$ be its translate by a measurable vector $b(y)$. Suppose $|b_i(y)-b_i^0|\le C_i a_i(y)$ on every occupied row up to null sets for fixed anchors $b_i^0$ and constants $C_i\ge0$. Then $E_b$ is contained in the translate by $b^0$ of the global coordinate dilation of the centered source $E_0$ by factors $1+C_i$ in the first $k$ coordinates. The true complete halo and actual-source price consequently obey
\[
 \frac{|\{M_S1_{E_b}>p\}|}{|E_b|}\le
 \prod_{i\le k}(1+C_i)\,\frac{|\{M_S1_{E_0}>p\}|}{|E_0|}.           \tag{5}
\]
The inclusion uses the real original intervals row by row, the maximal operator commutes with this *one common global dilation*, and $|E_b|=|E_0|$ by Fubini. For any centered hyperbolic cross or general unconditional log-convex $E_0$ of r094, and $\sum_i\log(1+C_i)=o(L)$, this whole shear class therefore retains the sharp static upper coefficient $c_Q$. One arbitrary sample-dependent shift per row bounded by a fixed multiple of that row's width, including all intervals containing a common point, cannot realize the r029 log gain. The union of these shear classes allowing $k=1$ has the matching lower from the centered r093 cross. More generally for prescribed $k=o(L)$ one may take a unit $k$-box times a centered cross in $N-k$ suffix dimensions: its prefix interior contributes the full suffix halo, whose logarithm relative to source volume is $q_0+c_Q L+o(L)$ by r092–093, since replacing $N$ by $N-k$ costs $O(k)+o(L)$ and changes $\tau$ by $o(L/N)$. For larger prescribed $k$ only the displayed upper is claimed. Unbounded relative displacement is genuine: (1) gives complete counterexamples even at one fixed width if centers can spread among too many distinct locations.

**Elementary checks.** Local Python `fractions.Fraction` enumerated all $2^6$ overlap patterns for three colors with two independent copies per color and single-point coverage probability $1/31$. Raw mass was $6/31$, surviving mass $148230000/887503681$, and lost mass $23544906/887503681\le30/961$, the expected ordered pair count. At $\eta=1/24$ the worst certified removed fraction on a length $4P/\eta$ coordinate interval is $2\eta-2/(4/\eta)=1/16>1/24$. The existence proof uses the volume and independence estimates in §1, not numerical evidence about a sampled large-dimensional halo.

The near-optimal one-coordinate Steiner reduction fails in the strong form (1). The translated-box estimates control source families with the stated fixed widths or bounded relative shifts. They do not construct a full physical input scoring $Q_N N^{c}$ for fixed $c>0$, and do not supply the missing general upper for completely unconstrained indicators.

## 3. An all-source periodic density criterion at macroscopic receiver coverage

On the unit $N$-torus let $d_N(p,\theta)$ be the infimum density of a measurable indicator whose **periodic strong rectangle halo** at strict threshold $p$ covers at least a fraction $\theta$, $0<\theta<1$. Each finite axis-parallel rectangle is averaged in the periodic extension. Then, with the unrestricted physical $H_N(p)$,
\[
       \frac{\theta}{H_N(p)}\le d_N(p,\theta)
       \le \frac{-\log(1-\theta)}{H_N(p)}.                              \tag{6}
\]
For the first inequality, first restrict periodic witnesses to sides $\le K$. Truncate the actual periodic source to $[0,R]^N$ for integer $R$, losing only a boundary layer of thickness $K$ among the periodic good receivers. Original source mass is precisely $R^N$ times the torus density. Apply the true physical indicator constant and let $R\to\infty$, then $K\to\infty$; strict finite witnesses exhaust the original periodic halo. This proves the lower density bound.

For the other direction take any bounded positive-volume physical source $D$; such sources suffice by monotone box truncation and strictness. Its physical strict halo is bounded, since each successful rectangle side is at most the corresponding width of a bounding box of $D$ divided by $p$. Dilate $D$ isotropically so that it and its strict halo fit one torus chart, with halo fraction $f\downarrow0$ and unchanged ratio $h$. Take $M=\lceil\log(1-\theta)/\log(1-f)\rceil$ independent uniform translates of this source, using their union as one actual periodic indicator. Its density is at most $Mf/h$ deterministically. Each receiver has probability at least $1-(1-f)^M\ge\theta$ of having a true witness from one of the translated copies. Some realization has covered fraction at least $\theta$, by integration on the torus. Let $f\downarrow0$ and $h\uparrow H_N(p)$ to conclude (6).

For instance with $\theta_N=1-1/N$, the negative logarithm of the density and $\log H_N(p_N)$ differ by at most $O(\log\log N)$, so they have exactly the same coefficient at the $\log N$ scale. A lower bound on the source density of a single almost-everywhere-covering periodic indicator is therefore a sufficient actual geometric upper module, without isolated receiver realizations. No such sharpened lower density estimate is supplied here.
