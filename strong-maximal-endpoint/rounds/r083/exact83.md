# Exact capped norm for one signed physical corner

## Scope and governing reduction

For the signed fixed-corner tensor Hardy average \(H_a\), r028 (the
\(\mathbf r=0\) case) already proves the physical equality
\[
 \mathfrak C_{N,H_a}^{[A]}=Q_{0,N}^{[A]}.                \tag{1}
\]
The left side permits arbitrary physical inputs \(0\le F\le A\); the right
side is the scalar Gamma problem below. This note does not assert a new
rearrangement theorem. It solves that scalar variational problem exactly and
then invokes (1).

Let \(S\sim\operatorname{Gamma}(N,1)\), write
\(\overline Q_N(s)=\mathbb P(S>s)\), and minimize
\[
 J(h)=\int_0^\infty e^{-s}\Phi_N(h(s))\,ds
\tag{2}
\]
over \(0\le h\le A\), nondecreasing, subject to
\[
 \mathbb E h(S)=1.                                      \tag{3}
\]
If \(J_*\) denotes the minimum, then
\[
 \boxed{\mathfrak C_{N,H_a}^{[A]}=J_*^{-1}.}             \tag{4}
\]

## 1. KKT profile and unique mean root

For \(b>0\), set \(\lambda_b=\Gamma(N)/b^{N-1}\) and let \(h_b(s)\)
minimize pointwise
\[
 \Phi_N(t)-\left({s\over b}\right)^{N-1}t,
 \qquad0\le t\le A.                                    \tag{5}
\]
Since \(\partial\Phi_N(1)=[1,N]\), this profile is
\[
 h_b(s)=
 \begin{cases}
 0,&s<b,\\
 1,&b\le s\le bN^{1/(N-1)},\\
 [\Phi_N']^{-1}((s/b)^{N-1}),&
 bN^{1/(N-1)}<s<b[\Phi_N'(A)]^{1/(N-1)},\\
 A,&s\ge b[\Phi_N'(A)]^{1/(N-1)}.
 \end{cases}                                             \tag{6}
\]
Endpoint conventions are immaterial. In particular \(h_b\) is
nondecreasing.

Let \(m(b)=\mathbb E h_b(S)\). It is continuous and strictly decreasing,
with limits \(A\) and zero as \(b\downarrow0\) and \(b\uparrow\infty\).
Thus there is a unique \(b_*=b_{N,A}\) such that
\[
 m(b_*)=1.                                                \tag{7}
\]
For every feasible \(h\), integrate the pointwise Lagrange inequality from
(5), with physical weight \(e^{-s}ds\):
\[
 J(h)-\lambda_b(\mathbb Eh(S)-1)
 \ge J(h_b)-\lambda_b(\mathbb Eh_b(S)-1).               \tag{8}
\]
At \(b=b_*\), both means equal one, so \(J(h)\ge J(h_{b_*})\).
Thus \(h_{b_*}\) is the global minimizer even before imposing monotonicity,
and strong duality follows directly from this supporting inequality.

## 2. Layer-cake formulas

For \(t>1\), define the activation location
\[
 s_b(t)=b[\Phi_N'(t)]^{1/(N-1)}.                         \tag{9}
\]
Layer-cake under the Gamma law and under the physical exponential measure
give the exact identities
\[
 m(b)=\overline Q_N(b)+\int_1^A\overline Q_N(s_b(t))\,dt,\tag{10}
\]
\[
 J(b)=e^{-b}+\int_1^A\Phi_N'(t)e^{-s_b(t)}\,dt.         \tag{11}
\]
Indeed, \(\{h_b>t\}=\{S>s_b(t)\}\) for \(t>1\), whereas every
\(0<t<1\) is activated at \(s=b\). Formula (11) uses
\(\Phi_N(h)=\int_0^h\Phi_N'(t)dt\).

The true fixed-corner logarithmic charge can equivalently be written
\[
 \log\mathfrak C_{N,H_a}^{[A]}
 =-\log J(b_*)=b_*-\log\bigl(e^{b_*}J(b_*)\bigr).        \tag{12}
\]

## 3. Strict halos and thresholds

The scalar minimizer has mean exactly one, while the endpoint score uses a
strict halo. For every \(\sigma\in(1,A)\), let \(b(\sigma)\) be the unique
root \(m(b)=\sigma\). The profile \(h_{b(\sigma)}\) is strictly feasible.
As \(\sigma\downarrow1\), continuity of (10)--(11) gives
\(J(b(\sigma))\to J(b_*)\). The physical lower construction in r028,
applied to these strictly feasible profiles and then passed to the limit,
proves the reverse inequality in (4) without replacing a weak threshold by
a strict one.

Envelope differentiation of the pointwise Lagrangian gives
\[
 {J'(b)\over m'(b)}=\lambda_b.                           \tag{13}
\]
Consequently, for the score at mean \(\sigma\),
\[
 {d\over d\sigma}\log[J(b(\sigma))^{-1}]
 =-{\lambda_{b(\sigma)}\over J(b(\sigma))}.            \tag{14}
\]
One-sided versions cover activation endpoints.

## 4. Relation to r082

The capped Gamma-Fenchel argument of r082 supplies a valid fixed-corner
upper coefficient \(e^{q_{N,A}}\). By the already exact physical-scalar
reduction (1) and the minimization above,
\[
 J(b_*)^{-1}\le e^{q_{N,A}}.
\]
The KKT profile, not the r082 radial Fenchel profile, is the actual optimum
for this fixed physical corner. These identities do not provide a packing
theorem for variable corners, so they do not by themselves change a general
strong-maximal endpoint bound.

## 5. A one-sided finite-dimensional certificate

The layer-cake formulas give a particularly simple rigorous upper certificate
for the actual strict fixed-corner norm.  Suppose an interval computation (or
an exact estimate) supplies \(b_{\rm hi}>0\) with
\[
 m(b_{\rm hi})<1.                                       \tag{15}
\]
Since \(m\) is strictly decreasing, the mean-one root satisfies
\(b_*<b_{\rm hi}\).  From (13), \(J'(b)=\lambda_bm'(b)<0\), hence
\[
 J(b_*)>J(b_{\rm hi})\ge e^{-b_{\rm hi}},
\]
where the last inequality is immediate from the positive integral term in
(11).  Therefore
\[
 \boxed{\log\mathfrak C_{N,H_a}^{[A]}=-\log J(b_*)<b_{\rm hi}.}\tag{16}
\]
This requires no separate certification of the scaled integral
\(e^{b_*}J(b_*)\).  The strict-halo issue is already resolved in Section 3:
r028 identifies the physical strict norm with the \(\sigma\downarrow1\)
limit \(J(b_*)^{-1}\), so (16) applies to the genuine cap-\(A\) corner
constant.  In particular, a certified cap-three, dimension-64 value
\(b_{\rm hi}\) may replace a coarser r082 fixed-corner coefficient wherever
only this one-sided corner upper bound is needed. It does not alter the cap
or threshold hypotheses, and it does not supply variable-corner sharing.
