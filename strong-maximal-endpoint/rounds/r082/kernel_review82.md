# Radial Fenchel kernel charge at one fixed physical corner

## Verdict

The proposed charge is correct.  For \(N\geq2\), a fixed cap \(A>1\), and
\(
 \Phi=\Phi_N(t)=t(1+\log^+t)^{N-1},
\)
it gives the fixed-corner weak estimate
\[
 \bigl|\{x:H_aF(x)>1\}\bigr|\leq e^{q_{N,A}}J_N(F),
 \qquad 0\leq F\leq A. \tag{1}
\]
It works simultaneously in all signed receiver orientations around the
*same* physical corner, without a factor \(2^N\).  Moreover \(e^{q_{N,A}}\)
is the exact minimum among proofs which use (i) a nonnegative coefficient
depending only on the total logarithmic radius, (ii) its pointwise capped
Fenchel majorant, and (iii) Markov after the one Fubini calculation below.

This is a fixed-corner theorem only.  It neither identifies the actual best
fixed-corner weak norm nor charges a moving-corner family or the strong
maximal operator.

## 1. Calibration and its elementary boundary checks

Let \(S\sim\Gamma(N,1)\), with density
\[
 \gamma_N(s)={e^{-s}s^{N-1}\over\Gamma(N)}\mathbf1_{s>0}.
\]
For \(q\in\mathbb R\), put
\[
 t_q(s)=\min\{A,\Phi^{-1}(e^{s-q})\}. \tag{2}
\]
The map \(q\mapsto\mathbb E t_q(S)\) is continuous and nonincreasing, with
limits \(A\) and \(0\) at the two ends of the real line.  It is strictly
decreasing once \(q+\log\Phi(A)>0\); this is the only regime in which its
value can be below \(A\).  Thus there is a unique \(q=q_{N,A}\) such that
\[
 \mathbb E t_q(S)=1. \tag{3}
\]

Except at the two one-point transition levels (which have Gamma measure
zero), set
\[
 a(s)=
 \begin{cases}
 1,&t_q(s)<1,\\
 1/\Phi'(t_q(s)),&1<t_q(s)<A,\\
 0,&t_q(s)=A,
 \end{cases}
 \qquad
 b(s)=\sup_{0\leq t\leq A}\{t-a(s)\Phi(t)\}. \tag{4}
\]
At \(t_q=1\), one may use \(a=1,b=0\), or any compatible one-sided
subgradient choice; at \(t_q=A\), the choice \(a=0,b=A\) is valid.  Neither
choice affects an integral.

The cases in (4) have the required equality
\[
 t_q(s)=a(s)\Phi(t_q(s))+b(s). \tag{5}
\]
Indeed, if \(t_q<1\), then \(a=1\) and
\(t-\Phi(t)\leq0\) on \([0,A]\), so \(b=0\).  If \(1<t_q<A\), the
derivative of \(t-a\Phi(t)\) vanishes at \(t_q\), and strict convexity of
\(\Phi\) on \((1,A)\) makes this its maximizer.  If \(t_q=A\), the
choice is simply \(0\cdot\Phi(t)+A\geq t\).  In particular, for every
\(0\leq t\leq A\),
\[
 t\leq a(s)\Phi(t)+b(s). \tag{6}
\]

Put \(P=\mathbb E[a(S)\Phi(t_q(S))]\).  On the non-capped set,
\(\Phi(t_q)=e^{S-q}\), while \(a=0\) on the capped set.  Equations
(3)--(5) therefore give the two exact identities
\[
 M:=\mathbb E[e^Sa(S)]=e^qP,
 \qquad 1-\mathbb Eb(S)=P. \tag{7}
\]
There is no tail-integrability defect: \(a=0\) once
\(S\geq q+\log\Phi(A)\).  The threshold is positive (otherwise
\(t_q=A\) almost surely, contradicting (3)), hence \(M<\infty\), and of
course \(0\leq b\leq A\).  This is precisely where the fixed cap removes
the uncapped inverse-Gamma tail problem.

## 2. Fixed-corner proof, including signs

First take the positive orthant and corner \(0\).  If
\(U_1,\ldots,U_N\) are independent \({\rm Exp}(1)\) variables and
\(S=\sum U_i\), then
\[
 HF(x)=\mathbb E\,F(x_1e^{-U_1},\ldots,x_Ne^{-U_N}). \tag{8}
\]
Define the auxiliary positive operator
\[
 \mathcal JF(x)=\mathbb E\,[a(S)\Phi(F(xe^{-U}))]. \tag{9}
\]
Applying (6) inside (8) yields
\[
 HF(x)\leq\mathcal JF(x)+\mathbb Eb(S).
\]
Thus \(HF(x)>1\) implies \(\mathcal JF(x)>P\).  Tonelli and the change
of variables \(y_i=x_ie^{-U_i}\) give
\[
 \int_{\mathbb R_+^N}\mathcal JF(x)\,dx
 =\mathbb E[e^Sa(S)]\int_{\mathbb R_+^N}\Phi(F(y))\,dy
 =M J_N(F). \tag{10}
\]
Markov and (7) now prove (1).

For a general fixed physical corner \(a_0\in\mathbb R^N\), split both
source and receiver into the \(2^N\) signed orthants relative to \(a_0\),
and reflect each part to (8).  A box with \(a_0\) as a vertex and receiver
in one such orthant uses only the matching source orthant.  Those source
orthants are disjoint up to null coordinate hyperplanes, so summing (10)
adds their modulars to \(J_N(F)\), not to \(2^N J_N(F)\).  Hence (1) holds
for the union of all orientations at this one corner.  This assertion does
not cover boxes with different opposite vertices.

## 3. Exact optimality in the radial Fenchel--Markov class

Consider any measurable radial multiplier \(\alpha(S)\geq0\), with
\[
 b_\alpha(s)=\sup_{0\leq t\leq A}\{t-\alpha(s)\Phi(t)\},\quad
 B=\mathbb Eb_\alpha<1,\quad
 M_\alpha=\mathbb E[e^S\alpha(S)]<\infty. \tag{11}
\]
The identical argument would charge the fixed-corner halo by
\(M_\alpha/(1-B)\).  Apply its pointwise majorant at \(t_q(S)\):
\[
 1=\mathbb Et_q\leq\mathbb E[\alpha\Phi(t_q)]+B.
\]
Since \(\Phi(t_q)\leq e^{S-q}\), with equality before the cap,
\[
 1-B\leq e^{-q}M_\alpha,
 \qquad {M_\alpha\over1-B}\geq e^q. \tag{12}
\]
The multiplier (4) has equality in both steps of (12), by (7).  Therefore
it is optimal in this precisely stated class.  The assertion does not
exclude a nonradial multiplier, a non-Fenchel proof, or an argument using
additional geometry/overlap from producing a smaller true weak constant.

## 4. It strictly improves the r075 capped Hardy charge

The comparison can be made analytically, not only numerically.  Set
\[
 \zeta_A={\log(1+\log A)\over\log A},\qquad
 p=1+\zeta_A(N-1),\qquad D=(p')^N. \tag{13}
\]
Concavity of \(\log(1+u)-\zeta_Au\) on \([0,\log A]\) gives
\(t^p\leq\Phi(t)\) on \([0,A]\).  Thus the r075 cap-\(A\) version of
the tensor-Hardy calculation has charge \(D^p=(p')^{Np}\) (for \(A=2\),
this is exactly its displayed \(K_N\)).

That same number is available within the class (11).  Take
\[
 \alpha(s)=c e^{-s/p'},\qquad c=p^{-1}D^{p-1}. \tag{14}
\]
Since \(\mathbb E e^{S/p}=D\),
\[
 M_\alpha=cD=p^{-1}D^p.
\]
The elementary conjugate bound
\[
 \sup_{t\geq0}(t-\alpha t^p)
 =(p-1)p^{-p'}\alpha^{-1/(p-1)}
\]
shows \(B\leq(p-1)/p\), hence
\[
 {M_\alpha\over1-B}\leq D^p. \tag{15}
\]
The inequality is strict: the unconstrained maximizer in the displayed
conjugate runs through a nonempty interval inside \((1,A)\) as \(s\)
ranges over an interval of positive measure (it starts at \(D^{-1}<1\)).
On that interval \(\Phi(t)>t^p\), so the capped Fenchel intercept is
strictly smaller than the power intercept.  Combining (12) and (15),
\[
 e^{q_{N,A}}<(p')^{Np}. \tag{16}
\]
For \(A=2\), (16) is the promised strict improvement over the r075
constant \(\exp\{N+(2\zeta_2)^{-1}+O(N^{-1})\}\).  It does not assert
that either side is the sharp fixed-corner operator norm.

## Scope audit

The argument needs a fixed finite cap, a fixed corner, nonnegative input,
and the physical Lebesgue measure used in (10).  It does not apply directly
to a conditional action kernel, a receiver-dependent corner, a global
moving-anchor selector, or an uncapped input.  In particular, it is a
stronger and optimally calibrated replacement for the *r075 fixed-corner
Hardy charge*, not a resolution of the source-budget/anchor-packing gap.
