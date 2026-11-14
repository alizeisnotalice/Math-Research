# Symmetric count-profile sources: exact rectangle formulas

This note gives a non-enumerative formula for the complete uncentered strong-maximal halo of a two-level, unconditional source. It also gives exact source budgets and a finite-dimensional integration scheme. It does **not** identify a near-extremizer or assert a new asymptotic lower bound.

## 1. Source and exact budget

Fix an outer radius (R>1). For (y\in[-R,R]^n), let

\[
N(y)=\#\{i:|y_i|\le1\},\qquad f(y)=w_{N(y)},
\]

and set (f=0) outside the cube. Assume (0\le w_0\le w_1\le\cdots\le w_n<\infty). The indicator source in the question is (w_j=\mathbf 1_{j\ge k}). Each coordinate is inner on a set of length (2), outer on a set of length (2(R-1)), and supported on length (2R). Thus, exactly,

\[
\int f=(2R)^n\mathbb E[w_B],\qquad B\sim\mathrm{Bin}(n,1/R).
\]

In particular,

\[
|E_{R,k}|=(2R)^n\Pr\{\mathrm{Bin}(n,1/R)\ge k\}
=\sum_{j=k}^n\binom nj2^j[2(R-1)]^{n-j}.
\]

This is the actual Lebesgue source cost, not a Boolean-grid normalization. For (k/n\to\alpha>1/R), the binomial tail has exponential rate \(-D(\alpha\Vert1/R)\), so \(n^{-1}\log|E_{R,k}|\to\log(2R)-D(\alpha\Vert1/R)\).

If the comparison denominator is the strong-maximal JMZ Orlicz cost, use

\[
J_n(f)=\int f(1+\log^+ f)^{n-1}
=(2R)^n\mathbb E\!\left[w_B(1+\log^+ w_B)^{n-1}\right],\quad B\sim\mathrm{Bin}(n,1/R).
\]

For the indicator (w_j=\mathbf1_{j\ge k}), (J_n(f)=|E_{R,k}|). For soft profiles bounded by (1), (J_n(f)=\int f=(2R)^n\mathbb E w_B\); for profiles taking values above (1), retain the displayed Orlicz expectation and do not replace it by \(\mathbb E w_B\).

For a general axis-parallel rectangle (Q=\prod_i I_i), put

\[
a_i=|I_i\cap[-1,1]|,\qquad b_i=|I_i\cap([-R,-1]\cup[1,R])|.
\]

The exact average is

\[
\fint_Q f=\sum_{j=0}^n w_j[z^j]\prod_{i=1}^n\left(\frac{b_i}{|I_i|}+\frac{a_i}{|I_i|}z\right).
\]

The indicator case is the Poisson-binomial upper tail \(\sum_{j\ge k}[z^j]\prod_i(q_i+p_i z)\), where (p_i=a_i/|I_i|), (q_i=b_i/|I_i|). The polynomial coefficients are computed by the standard recurrence (c_j\leftarrow q_i c_j+p_i c_{j-1}), in (O(n^2)) arithmetic (or (O(nk)) for one threshold). No (2^n) patterns are listed.

## 2. Full endpoint reduction

The formula below is for the **full** strong maximal operator, not just a preselected witness family. By sign symmetry, take a receiving point (x\in[0,\infty)^n). Fix all rectangle sides but one. If the other-coordinate count is (J), the one-dimensional score for an interval (I\ni x_i) is

\[
\frac{A|I\cap[-1,1]|+B|I\cap([-R,-1]\cup[1,R])|}{|I|},
\quad A=\mathbb E w_{J+1},\quad B=\mathbb E w_J.
\]

Monotonicity of (w) gives (A\ge B\ge0). Hence each side can be optimized by maximizing a linear functional in its inner and outer lengths.

The following canonical intervals suffice, coordinate by coordinate:

| Coordinate of (x_i\) | Maximizing interval(s) | Inner/outer probabilities or support factor |
|---|---|---|
| (0\le x_i\le1) | any interval containing (x_i) wholly in ([-1,1]) | inner probability (1) |
| (1<x_i\le R) | ([-1,x_i]) | (p_i=2/(x_i+1),\ q_i=(x_i-1)/(x_i+1),\ p_i+q_i=1) |
| (x_i>R) | one of ([-1,x_i]), ([-R,x_i]) | types A/B below |

Here is the endpoint check for the last two rows. For (1<x\le R), the interval can be shortened to end at (x); with full support, the average is (B+(A-B)p), where (p) is its inner fraction. Among intervals containing (x), (p\le2/(x+1)), attained by ([-1,x]). An interval leaking outside the support has no larger score. For (x>R), an interval can be shortened to end at (x). As its left endpoint (a) ranges through ([-R,-1]), its score is

\[
\frac{2A+(R-2-a)B}{x-a},
\]

whose derivative has constant sign; its maximum is at (a=-1) or (a=-R). On ([-1,1]) the score is \(((1-a)A+(R-1)B)/(x-a)\), again with constant-sign derivative, so only (a=-1,1) need be checked. On ([1,R]), the score (B(R-a)/(x-a)) is maximized at (a=1). But the (a=-R) value dominates the (a=1) value for (x>R): using (A\ge B), their difference is at least

\[
\frac{B(R+1)(x-R)}{(x+R)(x-1)}\ge0.
\]

Endpoints (a<-R) add only zero-valued length and are dominated by (a=-R). Thus no other interval endpoint improves the score. Applying this conditional one-coordinate argument successively proves the reduction for all rectangles. The proof also applies to every nonnegative nondecreasing count profile (w), not only an indicator threshold.

For (x_i>R), define the two choices

\[
\begin{array}{c|cc}
&s(x)=p(x)+q(x)&\theta(x)=p(x)/s(x)\\ \hline
A:\ [-1,x]&\dfrac{R+1}{x+1}&\dfrac2{R+1}\\[4pt]
B:\ [-R,x]&\dfrac{2R}{x+R}&\dfrac1R.
\end{array}
\]

These are the fraction of the interval inside the source support and the conditional inner probability given support. Notice type B has larger support fraction for (x>R); their ratio is

\[
\rho(x)=\frac{s_B(x)}{s_A(x)}=\frac{2R(x+1)}{(R+1)(x+R)},
\]

which is strictly increasing in (x).

## 3. Symmetry-compressed exact full halo

Partition the positive coordinates into

\[
L=\{x_i\le1\},\quad M=\{1<x_i\le R\},\quad O=\{x_i>R\},
\]

with sizes \(\ell,m,o\), and sort the outer coordinates as (R<x_{(1)}\le\cdots\le x_{(o)}). Let (Y_M) be a sum of independent Bernoulli variables with parameters (2/(x_i+1)), (i\in M). For each (r=0,\ldots,o), set

\[
S_r(x)=\left(\prod_{i\in O}s_A(x_i)\right)
\left(\prod_{j=o-r+1}^{o}\rho(x_{(j)})\right)
\mathbb E\,w_{\ell+Y_M+\mathrm{Bin}(o-r,2/(R+1))+\mathrm{Bin}(r,1/R)}.
\]

The two binomials in this expression are independent of each other and of (Y_M). For fixed (r), the count-profile expectation depends only on (r), not on which outer coordinates are type B. The support product is maximized by assigning type B to the (r) largest (x_i), because \(\rho\) increases. Therefore

\[
M_S f(x)=\max_{0\le r\le o} S_r(x).
\]

For (w_j=\mathbf1_{j\ge k}), replace the expectation by the corresponding Poisson-binomial tail. This is an exact complete-halo formula. For each point, all (o+1) candidates can be computed with coefficient recurrences; no endpoint-choice enumeration (2^o) is needed.

For the indicator source, a successful point at threshold (q>0) must satisfy (x_i<2R/q-R) on every outer coordinate: the average is at most that coordinate's support factor, and the largest available factor is (2R/(x_i+R)). For (q=1/2), the full halo lies in ([-3R,3R]^n). For soft profiles normalized by (w_n\le1), the same bound holds. In general replace (q) in this support argument by (q/w_n).

Consequently, the complete half-halo volume is the finite integral

\[
2^n\sum_{\ell+m+o=n}\frac{n!}{\ell!m!}
\int_{[1,R]^m}\int_{R<x_1<\cdots<x_o<3R}
\mathbf1\{\max_r S_r(x)>1/2\}\,dx_1\cdots dx_o\,dx_M.
\]

The factor (2^n) accounts for signs. The multinomial choice of (L,M,O) combined with the (o!\) outer ordering gives (n!/(\ell!m!)). Empty groups use unit integrals. This reduces a full halo computation to binomial/multinomial sums, Poisson-binomial recurrences, and a bounded integral; it does not convert the problem to a Boolean-cube halo.

On each fixed (L/M/O) stratum and outer ordering chamber, (F(x)=\max_rS_r(x)) is coordinatewise nonincreasing: middle-coordinate success probabilities decrease, and both outer support factors decrease. Thus, for a rectangular integration cell with lower and upper corners (x^-,x^+), (F(x^+)\le F(x)\le F(x^-)). If (F(x^+)>1/2), the whole cell is in the halo, so it contributes to an **inner (lower) volume bound**. If (F(x^-)\le1/2), the whole cell is outside the strict halo, so it can be discarded from the **outer (upper) volume bound**. Cells with (F(x^-) >1/2\ge F(x^+)) are undecided and contribute their whole volume only to the outer bound. Interval arithmetic for the rational score factors and directed enclosures for the finite probability recurrence yield certified finite-dimensional bounds. Adaptive subdivision should focus on undecided cells.

## 4. Low-rate and plausible parameter choices

If (R=2), the simple containment ( [-1,1]^n\subset E_{2,k}\subset[-2,2]^n\) and cube scaling give

\[
\frac{|\{M_S\mathbf1_{E_{2,k}}>1/2\}|}{|E_{2,k}|}
\le 2^n G_{\rm cube}(n),\qquad
G_{\rm cube}(n)\le2e^{2\sqrt{2n\log2}}.
\]

Thus its exponent is at most \(\log2<\mu\): the (R=2) count-threshold family cannot test the target rate, regardless of (k). More generally the exact budget and cube containment give

\[
H(E_{R,k})\le\frac{G_{\rm cube}(n)}{\Pr\{\mathrm{Bin}(n,1/R)\ge k\}}.
\]

For (k/n\to\alpha>1/R), this only permits an exponent as large as (D(\alpha\Vert1/R)). Thus choices with (D(\alpha\Vert1/R)<\mu) are already low-rate by this coarse bound. For (R>e^\mu\), values just above the solution of (D(\alpha\Vert1/R)=\mu) are the most plausible count thresholds to screen; this is only a necessary screening heuristic, not evidence of a large halo. MATLAB numerically solved that scalar equation (diagnostic only):

| (R) | (1/R) | root (\alpha_*\) with (D(\alpha_*\Vert1/R)=\mu) |
|---:|---:|---:|
| 4 | 0.25 | 0.982410043628 |
| 8 | 0.125 | 0.825938268772 |
| 16 | 0.0625 | 0.679803933918 |

For (k=n), the source is exactly the inner cube and its actual halo ratio is subexponential; the cube-containment bound above is then very wasteful. The endpoint formula still computes it exactly. The soft profile (w_j=t^{n-j}), (0\le t\le1), is a particularly cheap factorized control case: (f=\prod_i(\mathbf1_{[-1,1]}+t\mathbf1_{[-R,R]\setminus[-1,1]})), so (M_S f=\prod_i M_1g(x_i)). It should serve as a calibration/low-rate family, not as a claimed median-halo extremizer. General monotone (w_j) remains computable by the recurrence and complete formula above.
