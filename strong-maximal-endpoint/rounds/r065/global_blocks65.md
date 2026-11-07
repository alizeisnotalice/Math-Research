# Global localization for fixed equal blocks

## The fixed-block globalization statement

This note globalizes the critical-strip common-action theorem to arbitrary
receiver locations for a fixed number \(r\) of equal blocks. It is a
bounded-overlap localization argument, not a claim about within-block
angular actions or arbitrary physical rectangles.

Let \(N=rb\), with fixed \(r\), and let \(H:\mathbb R^r\to[0,A]\) have
finite global modular budget
\[
 J(H)=\int_{\mathbb R^r}e^{-\sum_iw_i}\Phi_N(H(w))\,dw<\infty. \tag{1}
\]
At every receiver \(u\), allow an arbitrary equal-block action vector
\(\mathbf c(u)\in I_K^r\). We count only receivers for which the
corresponding block average is strictly larger than one. Using the
fixed-block vector theorem in universal_action65.md and the fixed-anchor
transfer of r028, the conclusion is
\[
 \bigl|\{\hbox{such receivers}\}\bigr|_{e^{-\sum u_i}du}
 \le(1+o(1))Q_NJ(H).                                      \tag{2}
\]
The cap \(A\) and \(r\) are fixed. The result permits a receiver-dependent
action vector, but not a receiver-dependent within-block/angular kernel.

## 0. Fixed-action block transfer

For later use, fix one equal-block action vector \(\mathbf c\) and define
\[
 \mathcal R_{\mathbf c}K(u)=E_{\mathbf c}K(u+T).
\]
The exact block-model transfer lemma is
\[
 \nu\{u:\mathcal R_{\mathbf c}K(u)>1\}\le Q_NJ(K)          \tag{3}
\]
for every bounded \(K:\mathbb R^r\to[0,A]\) of finite global budget.  This
is not a general physical-rectangle assertion.

To prove it first take \(K\) compactly supported and a compact receiver
set \(E\) with a uniform strict margin.  In the physical deep realization,
put
\[
 f_D(x)=K(v(x)-D\alpha)\mathbf1_{[-1,1]^N}(x).
\]
At the physical point \(v=D\alpha+u\), the fixed-anchor rectangle average
is exactly \(\mathcal R_{\mathbf c}K(u)\).  The coarea limits of
r026/report.md give
\[
 K_D^{-1}\int\Phi_N(f_D)\to J(K),\qquad
 K_D^{-1}|E_D|\to\nu(E).
\]
The fixed-global-vector full-input weak theorem of r028/report.md bounds
the physical ratio by \(Q_N\), proving (3) for this compact strict case.

For a general bounded \(K\), take
\(K_M=K\mathbf1_{[-M,M]^r}\).  Translation continuity of the block density
makes \(\mathcal R_{\mathbf c}K_M\) continuous; these functions increase to
\(\mathcal R_{\mathbf c}K\).  Dini convergence on a compact strict subset
of the limiting halo, followed by exhaustion and monotone convergence,
proves (3).  Thus (3) is an unconditional lemma of the fixed block model.

## 1. Angular localization

Use angular and radial coordinates
\[
 x_i(u)=u_i-u_r\quad(1\le i<r),\qquad y(u)=\sum_{i=1}^ru_i. \tag{4}
\]
For \(\mathbf c\in I_K^r\), the block increment vector \(T\) satisfies
\[
 E\,x_i(T)=O(L),\qquad \operatorname{Var}x_i(T)=O(N),       \tag{5}
\]
uniformly in \(\mathbf c\). The first assertion follows from
\(m'(\mu)=0\): every block mean differs from \(b\mu\) by \(O(L)\).
Uniform exponential moments of the Gamma--Bernoulli increments then give,
for any prescribed \(K_0\), a fixed \(C\) such that
\[
 \sup_{\mathbf c\in I_K^r}
 \Pr\{\|x(T)\|_\infty>R\}\le\eta,\qquad
 R=C\sqrt{NL},\quad \eta=N^{-K_0}.                         \tag{6}
\]
Choose \(K_0>3\).

Set \(W=N\). For a shift \(\tau\in[0,W)^{r-1}\), let
\(\{C_k^\tau:k\in\mathbb Z^{r-1}\}\) be the half-open \(W\)-cube partition
of angular space, and define
\[
 H_k(w)=H(w)\,
 \mathbf1_{\{x(w)\in C_k^\tau+[-R,R]^{r-1}\}}.             \tag{7}
\]
If \(x(u)\in C_k^\tau\), then (6) gives
\[
 E_{\mathbf c}H_k(u+T)\ge E_{\mathbf c}H(u+T)-A\eta.       \tag{8}
\]
Put \(\beta=(1-2A\eta)^{-1}\). A strictly successful receiver for \(H\)
is strictly successful for \(\beta H_k\), because
\(\beta(1-A\eta)>1\). Also, for \(a\ge1\) and \(h\ge0\),
\[
 \Phi_N(ah)\le a^N\Phi_N(h).                               \tag{9}
\]
For \(h\le1\), this follows from \(1+\log a\le a\); for \(h\ge1\),
divide by \(\Phi_N(h)\) and use the same inequality. Thus
\(\beta^N=1+o(1)\), and \(\beta H_k\) has a fixed cap
\(A_*=\max\{2A,2\}\).

## 2. The random-shift budget identity

For a fixed source point \(w\), the expected number of enlarged grid cells
in (7) which contain \(x(w)\) equals
\[
 (1+2R/W)^{r-1}=1+o(1).                                   \tag{10}
\]
Tonelli's theorem, applied to the nonnegative modular integrand, therefore
gives the exact expectation identity
\[
 E_\tau\sum_kJ(H_k)
 =(1+2R/W)^{r-1}J(H).                                     \tag{11}
\]
There is a shift \(\tau\) for which the right-hand side bounds the sum.
Fix it. By (9),
\[
 \sum_kJ(\beta H_k)\le(1+o(1))J(H).                        \tag{12}
\]
This handles the countably infinite collection of cells; cells of zero
budget have no successful receiver after the cut and are omitted.

## 3. Cellwise normalization and radial truncation

Normalize the global budget once, by a diagonal translation, to
\[
 J(H)=e^{-q_0},\qquad q_0=N\mu-d_N.                        \tag{13}
\]
This does not alter the target-volume/budget ratio. For each nonzero cell
source \(K_k=\beta H_k\), choose a cellwise diagonal translation
\(a_k\mathbf1\) so that
\[
 \widetilde K_k(w)=K_k(w+a_k\mathbf1),\qquad
 J(\widetilde K_k)=e^{-q_0}.                              \tag{14}
\]
Diagonal translation leaves \(x\) unchanged. If
\(\widetilde u=u-a_k\mathbf1\), then success and weighted volume transform
by
\[
 E_{\mathbf c}\widetilde K_k(\widetilde u+T)
 =E_{\mathbf c}K_k(u+T),\qquad
 d\nu(\widetilde u)=e^{ra_k}\,d\nu(u),                    \tag{15}
\]
where \(d\nu(u)=e^{-y(u)}du\).

For every successful receiver of \(\widetilde K_k\), the r026 point-budget
theorem gives
\[
 J_+\!\left(\widetilde K_k(\widetilde u+\cdot)\right)^{-1}
 \le C_r\sqrt{\frac{N}{b^r}}\,Q_N.                        \tag{16}
\]
The left budget is at most \(e^{y(\widetilde u)}e^{-q_0}\).
Since \(r\) is fixed and \(\log Q_N=q_0+O(L)\), (16) implies
\[
 y(\widetilde u)\ge-C_1L.                                 \tag{17}
\]

The Jacobian of the coordinate change (4) is \(r^{-1}\), so inside one angular cube the
weighted volume above a radial level is exactly
\[
 \int_{\{x\in C_k^\tau,\ y>T\}}e^{-y}\,du
 =\frac{W^{r-1}}r e^{-T}.                                 \tag{18}
\]
Take \(T=C_2L\). The transformed receiver set for cell \(k\) has no mass
below \(-C_1L\) by (17). Its portion with \(y>T\) is bounded by (18).
Undoing (15) and using (14) gives the budget-weighted form
\[
 \nu\{u:x(u)\in C_k^\tau,\ y(u-a_k\mathbf1)>T\}
 \le \frac{W^{r-1}}r e^{-T+q_0}\,J(K_k).                  \tag{19}
\]
Summing (19) and using (12)--(13), the raw weighted volume of all these
far parts is at most
\[
 (1+o(1))\,\frac{W^{r-1}}r e^{-T}.                        \tag{20}
\]
Choosing \(C_2\) sufficiently large makes (20) negligible compared with
\(Q_Ne^{-q_0}\), which is polynomial in \(N\). Notice that this summation
uses \(J(K_k)\): summing the cellwise bound (18) alone over infinitely many
cells would be invalid.

## 4. The near strip and conclusion

For the remaining receivers,
\[
 -C_1L\le y(\widetilde u)\le C_2L.                         \tag{21}
\]
Each normalized cell source has actual critical budget, fixed cap, and a
successful receiver in the fixed block-action box. The vector
universal-action theorem applies to all its receivers in (21): one
cell-independent action \((c_N,\ldots,c_N)\), followed by the common
\(a_N=O(L^4/N)\) repair, receives them all. The fixed-action block transfer
lemma (3) then controls their union in cell \(k\) by
\[
 \nu(\text{near receivers in cell }k)
 \le(1+o(1))Q_NJ(K_k).                                    \tag{22}
\]
The \(o(1)\) is uniform in \(k\), because normalization has removed the
radial location and all caps and block counts are fixed.

Sum (22), use (12), and add (20). This proves (2).

The argument is robust to infinitely many nonempty cells because all
interchanges are Tonelli interchanges of nonnegative quantities. Its
essential restrictions are also clear: fixed \(r\), equal macroscopic
blocks, a fixed height cap, and actions confined to \(I_K^r\). It neither
proves an upper bound for general physical rectangles nor bypasses the
action-dependent angular transport identified in r063--r064.
