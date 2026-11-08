# Capped fixed-corner tensor-Hardy bound

## Result

Let \(N\ge2\), \(A>1\), and let \(S\sim\operatorname{Gamma}(N,1)\), with
\[
 d\nu_N(s)=e^{-s}{s^{N-1}\over(N-1)!}\,ds.
\]
Define
\[
 t_q(s)=\min\{A,\Phi_N^{-1}(e^{s-q})\}.
\]
The map \(q\mapsto\mathbb E t_q(S)\) is continuous and nonincreasing, strictly decreasing in the regime
\(q+\log\Phi_N(A)>0\) that contains every root of level below A, and has limits \(A,0\) at \(-\infty,+\infty\). Thus there is a unique
\(q=q_{N,A}\) such that
\[
 \mathbb E t_q(S)=1.                                   \tag{1}
\]
For every cap-\(A\) input and every fixed physical corner,
\[
 \boxed{|\{H_aF>1\}|\le e^{q_{N,A}}J_N(F).}            \tag{2}
\]
All signed receiver orientations cost no extra modular factor. The cap in
the definition of \(t_q\) is essential.

## Supporting functions

Except at the null singleton where \(t_q=1\), set
\[
 a_q(s)=
 \begin{cases}
 1,&t_q(s)<1,\\
 1/\Phi_N'(t_q(s)),&1<t_q(s)<A,\\
 0,&t_q(s)=A,
 \end{cases}
 \qquad
 b_q(s)=\sup_{0\le t\le A}\{t-a_q(s)\Phi_N(t)\}.
\tag{3}
\]
Then, for \(0\le t\le A\),
\[
 t\le a_q(s)\Phi_N(t)+b_q(s).                           \tag{4}
\]
On the lower branch \((a_q,b_q)=(1,0)\); on the middle branch the convex
supporting point is \(t_q\); on the capped branch \((a_q,b_q)=(0,A)\).
Consequently
\[
 t_q=a_q\Phi_N(t_q)+b_q,\qquad
 1-\mathbb E b_q=\mathbb E[a_q\Phi_N(t_q)]>0.           \tag{5}
\]
Moreover
\[
 M_q:=\mathbb E[e^Sa_q(S)]
 =e^q\mathbb E[a_q\Phi_N(t_q)],                         \tag{6}
\]
because \(\Phi_N(t_q)=e^{S-q}\) off the capped region, and \(a_q=0\) on
the capped region. In particular,
\[
 {M_q\over1-\mathbb E b_q}=e^q.                         \tag{7}
\]
The cap makes \(M_q\) finite; this is precisely where the uncapped version
can fail.

## Translation proof

First take the positive orientation with the corner translated to zero. If
\(T_i\) are independent \(\operatorname{Exp}(1)\) variables, then
\[
 H_aF(x)=\mathbb E F(x_1e^{-T_1},\ldots,x_Ne^{-T_N}).   \tag{8}
\]
For \(\Omega=\{H_aF>1\}\), (4) gives pointwise on \(\Omega\)
\[
 1-\mathbb E b_q
 <\mathbb E[a_q(S)\Phi_N(F(xe^{-T}))].                  \tag{9}
\]
Integrate, put \(y_i=x_ie^{-T_i}\), and use \(dx=e^Sdy\). Enlarging the
image of \(\Omega\) to all space gives, by Tonelli,
\[
 (1-\mathbb E b_q)|\Omega|\le M_qJ_N(F).                \tag{10}
\]
Equations (7) and (10) prove (2). For all signed orientations, split \(F\)
into the disjoint source orthants relative to the corner. A rectangle in one
receiver orthant meets only its matching source orthant, and
\[
 \sum_\eta J_N(F{\bf1}_{O_\eta})=J_N(F).
\]
This proves the asserted orientation-free cost accounting.

## Thresholds and scope

For \(0<\sigma<A\), choose \(q_{N,A}(\sigma)\) from
\(\mathbb E t_q(S)=\sigma\). Replacing one by \(\sigma\) in (5), (9), and
(10) gives
\[
 |\{H_aF>\sigma\}|\le e^{q_{N,A}(\sigma)}J_N(F).       \tag{11}
\]
Thus \(\sigma=1+o(\log N/N)\) is allowed for fixed \(A\ge2\) at large
\(N\), but the cap remains present. Estimating the resulting shift of \(q\)
requires a separate mean-slope and capped-tail analysis; the coefficient and constant-order expansion are proved separately in `asymptotic82.md`.

There is, however, an exact derivative formula. At values of \(q\) where
\(\mathbb E t_q=\sigma\), differentiation under the gamma integral gives
\[
 {dq_{N,A}\over d\sigma}=-{1\over\mathcal D(q)},        \tag{12}
\]
where
\[
 \mathcal D(q)=
 \mathbb E\!\left[t_q(S){\bf1}_{\{S<q\}}\right]
 +\mathbb E\!\left[
 {t_q(S)(1+\log t_q(S))\over N+\log t_q(S)}
 {\bf1}_{\{1<t_q(S)<A\}}\right].                      \tag{13}
\]
Indeed, \(\partial_qt_q=-t_q\) below one,
\(\partial_qt_q=-\Phi_N(t_q)/\Phi_N'(t_q)\) in the middle branch, and it
vanishes in the capped branch. The middle quotient is exactly the second
integrand in (13). Dominated differentiation is valid because
\(|\partial_qt_q|\le A\).

Thus no source rescaling, and hence no Orlicz dilation loss, is needed to
move the success threshold. The critical-window analysis in
`asymptotic82.md` proves
\(\mathcal D(q_{N,A}(1))\sim(r_N+1)/N\) uniformly in the required local
neighborhood. Consequently, for \(0\le\tau\le N^{-2}\),
\[
 q_{N,A}(e^{-\tau})-q_{N,A}(1)
 ={N\tau\over r_N+1}(1+o(1)).
\]

## Optimality within this proof class

Consider any measurable scalar \(a(s)\ge0,b(s)\) satisfying (4), with
\(D=1-\mathbb Eb>0\) and \(M=\mathbb E[e^Sa(S)]<\infty\). The same proof
has coefficient \(M/D\). Evaluating its inequality at \(t_q\) gives
\[
 D\le\mathbb E[a\Phi_N(t_q)].                           \tag{14}
\]
Also \(e^S\ge e^q\Phi_N(t_q)\): this is equality off the capped region,
and follows from \(e^{S-q}\ge\Phi_N(A)\) on it. Hence
\[
 M\ge e^q\mathbb E[a\Phi_N(t_q)]\ge e^qD.               \tag{15}
\]
So \(M/D\ge e^q\), with equality for (3). This is exact optimality only
among scalar radial Fenchel/translation kernels. It does not furnish a
packing over variable corners or a new general strong-maximal bound.

## Numerical check and asymptotic caution

A standard-library Python quadrature/bisection check gave
\(q_{2,2}\approx1.852480\) and verified (7) to relative error about
\(10^{-15}\) at \((N,A)=(2,2),(8,2)\). This is a sanity check, not proof.
The observed high-dimensional cap-two savings are promising, but proving a
constant-order expansion of \(q_{N,2}\) requires an Edgeworth remainder
smaller than the mean-root slope. The exact formula is proved here; `asymptotic82.md` supplies the separate Gamma-density remainder and constant-order expansion.
