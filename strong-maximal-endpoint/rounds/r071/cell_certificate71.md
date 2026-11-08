# Certified Cartesian-cell lower bounds for the local-cost kernel

## Conclusion

The proposed cell coefficient is a valid **uniform continuous lower bound**
after two restrictions:

1. each source cell must lie strictly above its receiver cell in every
   coordinate, `k_i>=j_i+1`; and
2. when `c=d\Delta`, the uniformly available marked-label factor is
   `2^{#{i:k_i-j_i>=d+1}}`, not the collocation factor with `>=d`.

With piecewise-constant, label-invariant source heights (or with `G_k` a
pointwise lower bound over every source label), a finite LP satisfying these
inequalities with interval/rational certificates gives a genuine continuum
weak-type lower bound. The existing floating r070 MILP does not do so: it
uses collocation coefficients, includes unsafe coordinate differences, and
uses the unsafe `>=d` label test. Its values remain exploratory only.

## Cartesian form and the Jacobian

Let `v` be the receiver logarithmic coordinate and `w=v+z` the source
coordinate. Write `S=\sum_i z_i`. The coarea factor in the conditional
kernel cancels when `(s,Pz)` is returned to Cartesian `z` coordinates:

\[
 dz={ds\,dH(Pz)\over\sqrt N},\qquad
 e^{-s}{\mathbf1_{D_{s,J}}(Pz)\over\sqrt N H_{N,c}(s)}
 ds\,dH(Pz)={e^{-S}\over H_{N,c}(S)}\,dz.                \tag{1}
\]

Thus, after summing the increment labels that are allowed by `z`, the
joint-Cartesian density is

\[
 {e^{-S}\over H_{N,c}(S)}
 \prod_{i=1}^N\bigl(1+\mathbf1_{\{z_i\ge c\}}\bigr),
 \qquad z_i\ge0.                                         \tag{2}
\]

For a fixed receiver label, (2) is exact on label-invariant sources. More
generally it is a lower kernel multiplied by a common lower height over all
source labels; it is not a pointwise lower bound in terms of an arbitrary
source's label-average. This qualification is essential for the factor in
(2).

The physical cell mass in log coordinates is

\[
 \omega_k=\int_{C_k}e^{-\sum_iw_i}\,dw
 =e^{-\Delta|k|}(1-e^{-\Delta})^N,                        \tag{3}
\]

where `C_k=\prod_i[k_i\Delta,(k_i+1)\Delta)`. If a source has constant
height `g_k` on `C_k` and on every label, its physical input mass (with Haar
probability in the label variable) is `\sum_k\omega_kg_k`.

## Uniform lower cell matrix

Take a receiver cell `C_j` and a source cell `C_k`, and suppose

\[
 k_i\ge j_i+1\quad(1\le i\le N).                         \tag{4}
\]

Then `w_i\ge v_i` for all `w\in C_k`, `v\in C_j`. The function

\[
 B(s)={e^{-s}\over H_{N,c}(s)}
\]

is decreasing, because `H_{N,c}` is increasing. Hence its value at
`v=j\Delta` is a lower bound throughout `C_j`. Put
`M=|k-j|=\sum_i(k_i-j_i)`. Integrating the source cell at that receiver
corner gives

\[
 q_M=\Delta^N\int_0^N
 {e^{-\Delta(M+x)}\,f_{\rm IH,N}(x)
  \over H_{N,c}(\Delta(M+x))}\,dx,                        \tag{5}
\]

where `f_{\rm IH,N}` is the density of the sum of `N` independent
`\operatorname{Unif}[0,1]` variables. This formula is simply the exact
change of variables `w_i=(k_i+a_i)\Delta`, followed by
`x=\sum_i a_i`; no numerical quadrature assertion is involved.

Now assume `c=d\Delta` with integer `d`. For all points of both cells,

\[
 z_i=w_i-v_i\ge(k_i-j_i-1)\Delta.
\]

Therefore the marked choice `J_i=1` is available **throughout** the two
cells only if

\[
 k_i-j_i\ge d+1.                                         \tag{6}
\]

The unmarked choice is always available. The certified retained multiplier is

\[
 m_{jk}=2^{\#\{i:k_i-j_i\ge d+1\}}.                      \tag{7}
\]

Consequently

\[
 L_{jk}=\begin{cases}
 m_{jk}q_M,&\text{if (4) holds},\\
 0,&\text{otherwise}
 \end{cases}                                               \tag{8}
\]

has the pointwise consequence

\[
 (\mathcal J_cG)(v,\sigma)\ge\sum_k L_{jk}g_k
 \quad(v\in C_j),                                         \tag{9}
\]

provided `G(w,\tau)=g_k` on `C_k` for every label `\tau`. The same
statement holds if `g_k` is a pointwise lower bound for `G` on that full
joint cell. For `c=\infty` there are no marked choices and `m_{jk}=1`.

Why the r070 collocation factor is unsafe: it used `k_i-j_i>=d`, but when
equality holds the lower corner separation is only `(d-1)\Delta<c`; the
marked label fails on a positive part of the product cell. It also allowed
coordinate differences zero, for which a source point can lie below a
receiver point. Neither issue affects a numerical heuristic, but both rule
out a uniform lower certificate.

## From a certified finite matrix to a continuum weak lower bound

Let finite, disjoint source cells and receiver cells be chosen, and use
half-open cells so their boundaries have zero physical measure. A certificate
requires rational interval data with directions chosen correctly:

\[
 \underline L_{jk}\le L_{jk},\qquad
 \underline\omega_j\le\omega_j\le\overline\omega_j.
\]

Choose rational nonnegative heights `g_k`, a rational `\lambda>0`, and a
strict rational margin `\eta>0`. If for every selected receiver cell

\[
 \sum_k\underline L_{jk}g_k\ge(1+\eta)\lambda,             \tag{10}
\]

then (9) puts every selected cell (apart from null boundaries) inside the
strict set `{\mathcal J_cG>\lambda}`. Therefore

\[
 \|\mathcal J_c\|_{1\to1,\infty}
 \ge {\lambda\sum_{j\in E}\underline\omega_j
       \over\sum_k\overline\omega_k g_k}.                 \tag{11}
\]

For a finite action with irrational `c`, alignment `c=d\Delta` can be kept
symbolically, or certified by common rational interval enclosures that
preserve the exact comparison in (6). The integral in (5), the exponential
cell weights in (3), and the LP feasibility/objective must be enclosed with
these directed inequalities; an arbitrary floating MILP solution is not an
interval certificate. Equation (10), rather than a mere equality to
`\lambda`, is needed because the usual weak norm uses a strict superlevel
set.

## Check of the fixed-dimensional tail and critical-column notes

The statements in `tail71.md` have no evident label-normalization defect.
With **counted** labels,

\[
 H_{N,c}(T+w)={2^N\over\Gamma(N)}T^{N-1}(1+O(T^{-1}))
\]

for fixed `N,c,w`: every one of the `2^N` label shifts contributes the same
leading polynomial term. In the outer cone all increment labels are allowed,
and for a fixed receiver label the map from an increment label to a source
label is a bijection. Hence the field uses exactly the counted source mass
`m=\sum_\tau\int e^{-\sum y}G_\tau(y)dy`, with no additional `2^N`.
The receiver shell has the compensating counted-label factor `2^N`; these
are precisely the factors that cancel in the tail identity
`t|\{j_cG>t\}|\to m`. With Haar-normalized labels, both counting factors
are removed and the same identity results.

Likewise, for the critical-column quantity

\[
 D_N(s)=\int_0^s e^{-r}{H_{N,c}(s-r)\over H_{N,c}(s)}dr,
\]

the fixed-`N` tail gives limit one, whereas the critical window has
`\partial_s\log H=1+o(1)`. For fixed `r` this makes the ratio tend to
`e^{-r}`, and monotonicity supplies domination by `e^{-r}`. Thus the
critical limit is `\int_0^\infty e^{-2r}dr=1/2`. This is a valid
noncommuting-limits observation, provided the stated compact-window
logarithmic-derivative estimate is available. It is not a
dimension-uniform deep-tail assertion.

