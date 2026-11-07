# r061 — a quantified product-family obstruction for the weighted mesh

Let `N` be the total dimension and put `m=N-2`.  This note tests only the
following *product construction*: the r060 two-dimensional weighted mesh,
at its distinguished receiving point and with its actual common-rectangle
value, is multiplied by the r024 hyperbolic step core in the remaining `m`
coordinates.  It is **not** an upper bound for arbitrary inputs or for the
full product halo.

The result is a useful calibration of the log-scale: a mesh pointwise
adaptive advantage of size `R=(log N)^2/N` cannot be recovered at that
receiver without moving the near-critical hyperbolic core by
`sqrt(mu) log N+o(log N)`.  The real `Phi_N` denominator of the product is
only changed by a constant factor, so it cannot cancel this log loss.

Here `mu` is the positive solution of

\[
 \mu=1+e^{-\mu},\qquad B=e^\mu,
\]

and `Phi_N(t)=t(1+log^+t)^{N-1}`.

## 1. Exact common-action deficit from the mesh

Use r060 with

\[
 R_N={\kappa(\log N)^2\over N},\quad \kappa>0\text{ fixed},\qquad
 M_N=\max\{3,\lfloor N/(\log N)^5\rfloor\}.
\]

Its outer crossing scale is

\[
 \eta_N={2(M_N-1)R_N^2\over1-(2M_N-1)R_N}=O((N\log N)^{-1}). \tag{1}
\]

After the r060 strict scale
\(b_N=(1+\eta_N)(1+N^{-2})\), the *actual common rectangle* at the
distinguished receiver has value

\[
 \gamma_N=b_NS^*_N=(1+N^{-2})(1-R_N)
 =1-R_N+o(R_N). \tag{2}
\]

The larger value `b_N A_N^*=1+N^{-2}` is an adaptive/composed value.  It is
not a strong-maximal rectangle and must not be inserted into a product
strong-maximal witness.  Indeed, for a product input the true factorization
is

\[
 M_{m+2}(g\otimes h)(x,y)=M_mg(x)M_2h(y), \tag{3}
\]

so the second factor at this receiver is exactly the common-action value in
(2), not the fibrewise reoptimized value.

## 2. Exact mixed-dimensional modular accounting

Let `h_N=b_N f_N^{mesh}` on the normalized mesh support (area one), and let
the core terminal be `A 1_[q,infinity)` with `A>1`.  Its product denominator,
relative to a bare-core denominator evaluated with the **same** `Phi_N`, is exactly

\[
 D_N(A):=\int \frac{\Phi_N(Ah_N(y))}{\Phi_N(A)}\,dy. \tag{4}
\]

This is the actual `Phi_N` expression; it is neither a `Phi_m` cost nor a
product of lower-dimensional modulars.

Write `c_N=b_N(1+1/N)`.  The r060 cell bookkeeping gives
\(\int f_N^{mesh}=1-R_N\), and `f_N^{mesh}` is equal to one away from its
low and recovery cells.  Their total complementary (unchanged) area is

\[
 U_N=1-L_N-H_N,
\]
\[
 L_N={p_1\over Q^{M-1}}+2R_N(1-Q^{-(M-1)}),\qquad
 H_N={p_1(Q-1)\over Q^{M-1}}+2R_N(1-Q/Q^{M-1}). \tag{5}
\]

Consequently `U_N=1-O(R_N+N^{-1})`.  Cellwise monotonicity gives the
two-sided, completely explicit estimate

\[
 U_N\,{\Phi_N(Ab_N)\over\Phi_N(A)}
 \le D_N(A)\le
 b_N(1-R_N)\left({1+\log(Ac_N)\over1+\log A}\right)^{N-1}. \tag{6}
\]

For the scale below, (1) implies `N log b_N=O(1/log N)`. A sharper
cellwise bound, retaining recovery-cell mass instead of charging it to
all cells, is
\[
 U_N\frac{\Phi_N(Ab_N)}{\Phi_N(A)}
 \le D_N(A)\le\frac{\Phi_N(Ab_N)}{\Phi_N(A)}[1+(e-1)H_N]. \tag{6a}
\]
Indeed every low or unchanged value is at most `Ab_N`, and
`Phi_N(Ab_N(1+1/N))/Phi_N(Ab_N)<=e`. Also
`0<=log[Phi_N(Ab_N)/Phi_N(A)]<=N log b_N=o(1)`.
Since `U_N->1` and `H_N->0`, this proves
\[
 \boxed{D_N(A_N)\to1,\qquad\log D_N(A_N)=o(1).} \tag{7}
\]
The ordinary `m=N-2` dimensional core uses `Phi_m(A_N)`. Its conversion
factor is exactly `Phi_N(A_N)/Phi_m(A_N)=(1+log A_N)^2=1+o(1)`.
Thus this dimensional bookkeeping supplies no logarithmic compensation.
The unchanged-cell lower bound is essential: a small upper bound on the
cost alone would not exclude a compensating decrease of the denominator.

## 3. The quantile price is exactly logarithmic

Take the r024 core amplitude

\[
 A_N=1+\sqrt{\log N/N}.
\]

For the fixed r024 action, let

\[
 S_m=\sum_{i=1}^m\{Y_i+\mu J_i\},
\quad Y_i\sim\operatorname{Exp}(1),\quad
 J_i\sim\operatorname{Bernoulli}((1+e^\mu)^{-1}),
\]

with independent summands, and write `F_m` for its CDF.  Define the two
critical quantiles by

\[
 F_m(q_N)=1-A_N^{-1},\qquad
 F_m(\widetilde q_N)=1-(A_N\gamma_N)^{-1}. \tag{8}
\]

The first is the zero-margin core threshold; the second is the threshold
needed after the *actual* mesh factor (2).  For example choose a positive strict margin `epsilon_N=N^(-2)=o(R_N)`.
It changes the required CDF mass by `O(epsilon_N)`, hence the core
threshold by `O(epsilon_N N/log N)=o(1)` using the estimates below.
After fixing `N` and this decreased threshold, take the core depth large
enough that the average error is less than `epsilon_N/2` and its relative
modular error is less than `epsilon_N`. This order preserves strictness
and avoids exchanging the dimension and depth limits.

Use the self-contained uniform density/CDF estimates in
[local_md61.md](local_md61.md), including its quantile localization. At lower-tail mass
`1-A_N^{-1}~sqrt(log N/N)`, they give
\[
 f_m(q_N)=\frac{\log N}{\sqrt\mu\,N}(1+o(1)). \tag{9}
\]
Here the standardized lower quantile has magnitude
`sqrt(log N)(1+o(1))`. To avoid a circular use of density uniformity on an
unknown quantile interval, first fix any `C>kappa sqrt(mu)`. The same local
estimate on the **fixed** interval `[q_N-C log N,q_N]` gives
`f_m(q_N+w)/f_m(q_N)=1+o(1)` uniformly there: its standardized translation
is `O(log N/sqrt N)`, and the resulting log-density change is
`O((log N)^(3/2)/sqrt N)=o(1)`. Consequently
\[
 F_m(q_N)-F_m(q_N-C\log N)
 =\frac{C(\log N)^2}{\sqrt\mu\,N}(1+o(1)).
\]
The required CDF mass change is
`A_N^{-1}(gamma_N^{-1}-1)=kappa(log N)^2/N(1+o(1))`, which is smaller.
It is also `o(1-A_N^{-1})`, so the target remains positive. Monotonicity
therefore brackets `qtilde_N` in that fixed interval. Only now apply the
mean-value theorem and the uniform density estimate to get
\[
 q_N-\widetilde q_N
 =\frac{A_N^{-1}(\gamma_N^{-1}-1)}{f_m(\xi_N)}
 =\kappa\sqrt\mu\,\log N\,(1+o(1)). \tag{10}
\]

The physical hyperbolic truncation from r024 converts these step terminals
into bounded compactly supported inputs after taking the depth parameter
sufficiently large for each fixed `N`; it preserves a chosen strict margin.
At the distinguished gadget receiver, compare the **core sectional scores**,
using the same `Phi_N` normalization. Their difference is

\[
 [q_N-\log\Phi_N(A_N)]-
 [\widetilde q_N-\log\Phi_N(A_N)-\log D_N(A_N)]
 =\kappa\sqrt\mu\log N+o(\log N).\tag{11}
\]

The gadget coordinate has been fixed here, so this is a core section,
not a positive-volume `N`-dimensional halo witness. No gadget receiving-box
volume has been inserted. The mixed denominator tends to its bare-core
value, while the **common-action** defect spends a definite logarithmic
coefficient in that section.  It cannot be used at this receiver to convert its adaptive
pointwise advantage into a gain in the unresolved log-`N` coefficient.

This is a product-family and distinguished-receiver result only.  It does
not rule out a non-product source, a non-product reception union, or a
different mesh region whose full strong halo has been separately analyzed.

### Efficiency limitation of this particular core

This step core is efficient only at the leading exponential-rate scale.
In fact (8) of local_md61 gives
`q_N=N mu-sqrt(mu N log N)+o(sqrt(N log N))`, whereas
`log Phi_N(A_N)=sqrt(N log N)(1+o(1))`. Its bare sectional score is thus
\[
 N\mu-(\sqrt\mu+1)\sqrt{N\log N}+o(\sqrt{N\log N}),
\]
which already loses an additional `sqrt(N log N)` compared with `Q_N`.
The quantified logarithmic product price above is a calibration of this
specific step core. It does not exclude tensor constructions using the
actually near-optimal soft core, which need a separate mean-response and
mixed-modular calculation.

## 4. Local Wolfram receipt

`efficiency61.wl` evaluates (4) cell by cell and evaluates (8) through the
exact finite identity

\[
 F_m(s)=\sum_{j=0}^m\Pr\{\operatorname{Bin}(m,p)=j\}
 F_{\Gamma(m,1)}(s-\mu j).
\]

The saved JSON records the completed finite rows.  It is floating numerical
corroboration of (6), (9), and (10), not an interval certificate or a proof
of a statement about general `C_N`.  The receipt uses `kappa=1/4` and
`N=128,256,512,1024`; all two quantile residuals are zero to the displayed
working precision, and the actual ratios lie between the two bounds in (6).
