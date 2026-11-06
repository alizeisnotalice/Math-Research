# r026: block-microcoordinate deep-layer reduction

Date: 2026-11-06. This concerns actual anchored rectangles and block-sum
inputs. It is a scoped reduction and convex program, not a claim that
correlated profiles improve the general endpoint constant.

Partition the coordinates into fixed blocks \(B_j\), \(|B_j|=b_j\),
\(\sum_jb_j=n\), and put
\[
 V_j(x)=\sum_{i\in B_j}-\log|x_i|,\qquad
 \alpha_j>0,\quad\sum_j\alpha_j=1.
\]
The signed-cube change of variables is exactly
\[
 dx=2^n e^{-\sum_jv_j}\prod_{j=1}^m\frac{v_j^{b_j-1}}{\Gamma(b_j)}\,dv_j.
 \tag{1}
\]
For \(v_j=L\alpha_j+w_j\), define
\[
 K_L=2^ne^{-L}L^{n-m}\prod_j\frac{\alpha_j^{b_j-1}}{\Gamma(b_j)}.
 \tag{2}
\]

## Exact compact-window reduction

Let \(H\) be bounded, nonnegative and compactly supported on
\(\mathbb R^m\), and use the genuine input
\[
 f_L(x)=H(V_1(x)-L\alpha_1,\ldots,V_m(x)-L\alpha_m).
\]
Then (1) and dominated convergence give
\[
 \int\Phi_n(f_L)dx=K_L\{J(H)+o(1)\},\qquad
 J(H)=\int_{\mathbb R^m}e^{-\sum_jw_j}\Phi_n(H(w))\,dw.\tag{3}
\]
The remaining factor is \(\prod_j(1+w_j/(L\alpha_j))^{b_j-1}\to1\).
Wolfram 15 independently returned this power-ratio limit as one; it also
returned \(\Phi_n'(a)=(1+\log a)^{n-2}(n+\log a)\) for \(a>1\), with right
derivative \(n\) at \(a=1\).

For the anchored rectangle increment vector, write
\(T_j=\sum_{i\in B_j}Z_{r_i}\). If a compact receiver set
\(E\subset\mathbb R^m\) satisfies
\[
 \inf_{u\in E}\mathbb EH(u+T)>1+\eta,\tag{4}
\]
then the receiver block-sum set \(V_j=L\alpha_j+u_j\), \(u\in E\), is an
actual strict witness for every sufficiently large \(L\). In fact the
rectangle average is exactly \(\mathbb EH(u+T)\), not an asymptotic proxy.
Its measure is
\[
 K_L\left\{\int_Ee^{-\sum_ju_j}\,du+o(1)\right\}.
\]
Thus the genuine fixed-\(n\) deep-layer lower score is
\[
 \liminf_{L\to\infty}\mathcal L(f_L)
 \ge\frac{\int_Ee^{-\sum u}\,du}{J(H)}.\tag{5}
\]
Signs and block-angular variables were integrated in (1), so (5) is not a
composition surrogate.

## Two-block orthant receiver

For \(m=2\), replace compact \(E\) by an orthant. Assume \(H\) is bounded,
coordinatewise nondecreasing, has a common coordinate lower cutoff (zero if
some \(w_j<-A\)), and has finite \(J(H)\). If for \(q\in\mathbb R^2\)
\[
 \mathbb EH(T-q)>1+\eta,\tag{6}
\]
then every receiver with \(u_j\ge-q_j\) is a strict actual witness by
coordinate monotonicity. The receiver volume divided by \(K_L\) tends to
\(e^{q_1+q_2}\); boundedness, the lower cutoff, and exponential weights
justify dominated convergence in both (3) and this unbounded receiver tail.
Hence
\[
 \liminf_{L\to\infty}\mathcal L(f_L)\ge e^{q_1+q_2}/J(H).\tag{7}
\]
The same assertion holds for fixed \(m\), but the two-block form makes the
first genuine correlated candidate explicit.

## Fixed-receiver convex dual: why the optimizer is correlated

Fix the block increment law with density \(f_T\) on \(\mathbb R_+^m\), and
let \(d\nu(w)=e^{-\sum_jw_j}dw\). Its likelihood ratio relative to \(\nu\)
is
\[
 R(t)=e^{\sum_jt_j}f_T(t)=\prod_{j=1}^mR_j(t_j),\qquad
 R_j(t)=e^tf_{T_j}(t).\tag{8}
\]
For one receiver at the origin, the exact convex program is
\[
 \inf_{H\ge0}\left\{\int\Phi_n(H)d\nu:\int R H\,d\nu\ge1\right\}.
 \tag{9}
\]
On finite boxes with \(R\) truncated it has Fenchel dual
\[
 \sup_{\lambda\ge0}\left\{\lambda-\int\Phi_n^*(\lambda R)d\nu\right\}.
 \tag{10}
\]
Monotone exhaustion gives the same variational candidate whenever the primal
value is finite. Its pointwise KKT form is
\[
 H=0\ (\lambda R<1),\qquad H=1\ (1\le\lambda R\le n),\qquad
 \Phi_n'(H)=\lambda R\ (\lambda R>n).\tag{11}
\]
Thus for two or more nondegenerate blocks the candidate depends on the
product likelihood ratio in (8), not only on \(t_1+\cdots+t_m\). It is a
real correlated block profile. Equations (9)--(11) do not show that its
score (5) or (7) exceeds scalar \(Q_n\): the receiver set, full Orlicz cost,
and block saddle prefactors still need comparison. They do reduce the next
search to a precise finite-dimensional integral problem rather than a formal
iterated-maximal construction.

## A rigorous localized equal-block obstruction

Here is a direct explanation for the observed half-log loss, with deliberately
limited scope.  Let \(m\) be fixed, \(b=n/m\) an integer, and use the same
parameter \(c\) in every coordinate, with \(c\) in a fixed small compact
neighborhood of \(\mu\).  Let \(T_1,\ldots,T_m\) be the independent block
sums and \(S=\sum_jT_j\).  There are constants \(0<a<\mu\), \(C_m,c_m>0\)
such that, for all large \(n\) and all \(s\in I_n=[an,n\mu]\), the
conditional density of \((T_1,\ldots,T_{m-1})\) given \(S=s\) obeys
\[
 \|k_s\|_\infty\le C_m n^{-(m-1)/2}.\tag{12}
\]

Proof: choose the common negative tilt \(\lambda_s\) for which one tilted
increment has mean \(s/n\).  The tilt range is compact.  Conditioning on
the sum is invariant under this common exponential tilt, hence
\[
 k_s(t_1,\ldots,t_{m-1})=
 \frac{g_{b,\lambda_s}(t_1)\cdots
 g_{b,\lambda_s}(t_{m-1})g_{b,\lambda_s}(s-\sum_{j<m}t_j)}
 {g_{n,\lambda_s}(s)}.\tag{13}
\]
The Gamma factor in every tilted characteristic function gives uniformly
\(\|g_{b,\lambda}\|_\infty\le C/\sqrt b\).  The same Fourier local-limit
argument at its tilted mean gives \(g_{n,\lambda_s}(s)\ge c/\sqrt n\).
Substitution in (13) proves (12).  This is the same elementary
Gamma--Bernoulli Fourier input used in r022; no moderate-deviation estimate
is required because the tilt recenters at \(s\).

Let \(H\ge0\) be supported in the total-sum strip \(\sum_jt_j\in I_n\),
and define \(h(s)=\mathbb E[H(T)\mid S=s]\).  Conditional Jensen and (12)
give the pointwise coarea estimate
\[
 \int_{\sum t=s}\Phi_n(H(t))\,d t_1\cdots dt_{m-1}
 \ge c_m n^{(m-1)/2}\Phi_n(h(s)).\tag{14}
\]
Consequently
\[
 J_m(H)\ge c_m n^{(m-1)/2}
 \int_{I_n}e^{-s}\Phi_n(h(s))\,ds.\tag{15}
\]
Both profiles have the identical single-receiver constraint
\(\mathbb E H(T)=\mathbb Eh(S)\).  If \(D_{m,I}\) and \(D_{1,I}\) denote
the infima of the two respective costs under this constraint, then
\[
 D_{m,I}\ge c_m n^{(m-1)/2}D_{1,I},\qquad
 -\log D_{m,I}\le-\log D_{1,I}-\frac{m-1}{2}\log n+O_m(1).\tag{16}
\]
This rigorously supports the \(m=2\) numerical \(\frac12\log n+O(1)\)
loss for profiles whose useful expectation is in the left-central strip.
It does not exclude a profile extracting a material contribution outside
\(I_n\), nor does it compare the resulting fixed-receiver score with the
full scalar weak constant \(Q_n\).  Those are precisely the remaining gaps.

## Height-capped full corner barrier, uniform in all deterministic ratios

The strip restriction can be removed for a bounded-height single-corner
problem, without a common block parameter or balanced blocks. Assume fixed
\(m\), arbitrary block sizes \(b_j\ge1\), arbitrary deterministic
\(r_i\in[0,1]\), and \(0\le H\le A_n\), where
\(\log A_n=o(n)\). If the corner is translated to the origin and
\(\mathbb EH(T)>1\), then
\[
 \frac1{J_m(H)}\le C_m\sqrt{\frac n{\prod_{j=1}^m b_j}}\,Q_n.\tag{17}
\]
A corner at \(-q\) has score
\(e^{\sum q_j}/J_m(H)=1/J_m(H(\cdot-q))\), so (17) applies to it as well.

Write \(c_i=-\log r_i\), with \(r_i=0\) treated as a degenerate zero
Bernoulli shift. Under common tilt \(e^{\eta Z_i}\), \(\gamma=1-\eta>0\),
\[
 Z_i\overset{\eta}=\operatorname{Exp}(\gamma)+
 c_i\operatorname{Bern}\bigl((1+e^{\gamma c_i})^{-1}\bigr),
 \qquad \frac1\gamma\le\mathbb EZ_i^\eta\le\frac\mu\gamma.\tag{18}
\]
The upper bound is \(1+x/(1+e^x)\le\mu\); Wolfram 15 maximized it at
\(x=\mu\), returning \(\mu=1.278464542761074\ldots\). For every
\(s/n\in[1/2,3]\), choose the common tilt whose total mean is \(s\). Then
\(\gamma\in[1/3,2\mu]\), uniformly in the ratio vector.

All tilted moments are uniform because \(c^ke^{-\gamma c}\) is bounded,
the exponential variance is uniformly positive, and every characteristic
function has the Gamma damping factor \(|\gamma/(\gamma-it)|\). The
triangular-array Fourier local limit theorem at its tilted mean therefore
gives \(g_{n,\eta}(s)\ge c/\sqrt n\). For each block no Fourier argument
is even needed: it is a \(\Gamma(b_j,\gamma)\) variable convolved with a
nonnegative Bernoulli-shift mixture, hence
\[
 \|g_{b_j,\eta}\|_\infty\le\|g_{\Gamma(b_j,\gamma)}\|_\infty
 \le C/\sqrt {b_j},\tag{19a}
\]
including \(b_j=1\). Common tilting leaves the conditional law unchanged,
so its conditional \((m-1)\)-density is bounded by
\[
 \|k_s\|_\infty\le C_m\sqrt{\frac n{\prod_jb_j}}.
 \tag{19b}
\]
This is uniform on \([n/2,3n]\) and in every deterministic ratio vector.

The original tails satisfy
\[
 \mathbb P\{S\notin[n/2,3n]\}\le e^{-\kappa n}.\tag{19}
\]
The lower tail follows from \(S\ge\Gamma(n,1)\); a uniform small positive
exponential moment gives the upper tail. Set \(I=[n/2,3n]\) and
\(h(s)=\mathbb E[H(T)\mid S=s]\). By (19),
\[
 \mathbb E[h(S)1_I(S)]\ge1-\varepsilon_n,
 \qquad\varepsilon_n=A_ne^{-\kappa n}=e^{-\kappa n+o(n)}.
\]
Conditional Jensen now gives
\[
 J_m(H)\ge c_m\sqrt{\frac{\prod_jb_j}{n}}
 \int_Ie^{-s}\Phi_n(h(s))\,ds.\tag{19c}
\]
With \(a_n=(1-2\varepsilon_n)^{-1}\), the profile
\(\widetilde h=a_nh1_I\) has strict scalar average exceeding one and
\[
 \Phi_n(a_nt)\le a_n(1+\log a_n)^{n-1}\Phi_n(t)=(1+o(1))\Phi_n(t).\tag{20}
\]
For every vector, \(S=\Gamma(n,1)+B\) gives
\[
 e^sf_S(s)=\mathbb E[e^B(s-B)_+^{n-1}/\Gamma(n)],
\]
which is nondecreasing. Weighted increasing rearrangement turns
\(\widetilde h\) into a nondecreasing scalar profile with no larger cost and
no smaller average. Its cost is at least \(Q_n^{-1}\). Combining this fact,
(19c), and (20) proves (17).

For balanced \(b_j=n/m\), (17) is the genuine
\((m-1)\log n/2+O_m(1)\) barrier observed numerically. If one block remains
small, (17) gives only a constant loss; this precisely quantifies an
unbalanced-block escape that is not ruled out. The statement does not cover
unbounded profiles, growing \(m\), multiple receiver components, or the
full nonanchored strong maximal operator.
