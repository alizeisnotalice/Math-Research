# Audit of the canonical-collapse bounded source

## Verdict

**PASS in the stated critical-window scope.**  The proposed product source
is bounded by two, compact, has strict averages above one throughout a
positive cube of receivers, and has a unique maximizing action in the
specified critical window.  The associated reanchors send that entire cube
to one point.  This is a genuine physical realization of the reanchor
collapse, but it is deliberately not asserted to satisfy the local
near-optimal cost qualification, and it says nothing about optimization
over actions outside the critical window.

Put \(r_r=e^{-c_r}\), \(\ell_0=1+r_r\), and \(l_0=-r_r\).  Assume
\(0<\rho<1/10\), set
\[
 H=2\rho<\ell_0/4,\quad
 \delta={p_{c_r}\rho\over8},\quad
 A=1+{1\over20N},\quad b={1\over20N},\quad
 d={b\ell_0H\over\ell_0-2H}.
\tag{1}
\]
Define \(g\) to be zero outside the enclosing interval [l_0-H,1+H] and, within that interval,
\[
 g(t)=
 \begin{cases}
 A+b(t-l_0),&l_0-H\le t\le l_0+H,\\
 A-d,&l_0+H<t<1-H,\\
 A+b\ell_0+b(t-1),&1-H\le t\le1+H.
 \end{cases}
\tag{2}
\]
The jumps at the interfaces are intentional and do not affect the exact
integrals.  Let \(F_N(y)=\prod_{i=1}^Ng(y_i)\).

## 1. Exact one-dimensional average

Fix \(x\in[1-\delta,1+\delta]\) and
\(c\in[c_r-\rho,c_r+\rho]\).  Write
\[
 \ell=x(1+e^{-c}),\qquad a=x-\ell=-e^{-c}x,
 \qquad u=a-l_0,\quad v=x-1.
\tag{3}
\]
The endpoints lie in their advertised bands.  Indeed \(|v|\le\delta<H\),
and, using \(\delta\le\rho/8\) and \(e^\rho-1\le(10/9)\rho\) for
\(\rho<1/10\),
\[
 |e^{-c}x-r_r|
 \le r_r\max\{1-(1-\delta)e^{-\rho},(1+\delta)e^\rho-1\}
 <2\rho=H.
\tag{4}
\]
Thus the integration interval contains the central constant segment and
has only its two endpoint pieces in the linear bands.  Direct integration
gives
\[
\begin{aligned}
 \int_a^xg(t)dt
 &=A\ell+{b\over2}(v^2-u^2)+b\ell_0(H+v)-d(\ell_0-2H)\\
 &=A\ell+b\ell_0v+{b\over2}(v^2-u^2),
\end{aligned}
\tag{5}
\]
because the choice of \(d\) cancels \(b\ell_0H\).  Since
\(\ell-\ell_0=v-u\), division by \(\ell\) yields the claimed exact
identity
\[
 {1\over\ell}\int_{x-\ell}^{x}g(t)dt
 = A+b(x-1)-{b(\ell-\ell_0)^2\over2\ell}.
\tag{6}
\]
In particular, for fixed \(x\), the unique maximizer among actions whose
endpoints remain in the two bands is \(\ell=\ell_0\).  The uniqueness is
strict because \((\ell-\ell_0)^2/\ell\ge0\), with equality only at
\(\ell_0\).

## 2. The admissible unique critical-window policy

The required action is exactly
\[
 c_*(x)=-\log\!\left({\ell_0\over x}-1\right),
\qquad x(1+e^{-c_*(x)})=\ell_0.
\tag{7}
\]
It lies in the critical window.  Writing \(h=x-1\), one has
\[
 c_*(x)-c_r=\log x-\log(1-h/r_r).
\tag{8}
\]
Here \(\delta/r_r=\rho/[8(1+r_r)]<\rho/8\).  The elementary bounds
\(|\log x|\le2\delta\) and
\(|\log(1-h/r_r)|\le2\delta/r_r\) therefore give
\[
 |c_*(x)-c_r|
 \le2\delta+2\delta/r_r
 ={\rho\over4}\left(p_{c_r}+{1\over1+r_r}\right)<\rho.
\tag{9}
\]
Thus (7) is admissible, and (6) proves unique maximization on the full
critical action interval, not merely at a formal stationary point.

At that action the one-dimensional average is
\[
 m_*(x)=A+b(x-1)\ge A-b\delta>1,
\tag{10}
\]
since \(\delta<1\).  For every
\(x=(x_1,\ldots,x_N)\in[1-\delta,1+\delta]^N\), take the independent
coordinate actions \(c_i=c_*(x_i)\).  Product Fubini and (10) give
\[
 {1\over\prod_i x_i(1+e^{-c_i})}
 \int_{\prod_i[-e^{-c_i}x_i,x_i]}F_N(y)dy
 =\prod_{i=1}^Nm_*(x_i)>1.
\tag{11}
\]

## 3. Boundedness, compactness, and exact collapse

The minimum of \(g\) on its nonzero part is positive: from
\(H<\ell_0/4\), \(d<2bH\), and the possible endpoint corrections are at
most \(bH\).  Its maximum obeys
\[
 \|g\|_\infty\le A+b(\ell_0+H)
 \le1+{3.2\over20N},
\tag{12}
\]
because \(\ell_0<2\) and \(H<0.2\).  Hence
\[
 0\le F_N\le(1+3.2/(20N))^N\le e^{0.16}<2.
\tag{13}
\]
The support is a bounded closed interval in each coordinate, so its product
support is compact.

Finally set
\[
 \lambda(c)={1+e^{-c}\over1+e^{-c_r}}.
\tag{14}
\]
At the maximizing policy, (7) gives
\[
 \lambda(c_*(x))={\ell_0/x\over\ell_0}={1\over x}.
\tag{15}
\]
Coordinatewise, the canonical reanchor is therefore
\[
 T_i(x)=\lambda(c_i)x_i=1,
\tag{16}
\]
on the entire receiver cube.  This is the exact point-collapse envisioned
in the reanchoring discussion.

## Scope

The source in fact fails the wide-window local-cost qualification
**exponentially**, not merely by lack of a proof.  Write \(B=1/20=.05\),
so \(A=1+B/N\) and \(b=B/N\).  On the middle segment,
\[
 Nd={B\ell_0H\over\ell_0-2H}<2BH\le.02,
\tag{17}
\]
while the left-band minimum has increment \(B-BH\ge.04\).  The right
band is larger.  Hence every point of the full source interval
\([l_0-H,1+H]\) obeys
\[
 g(t)\ge1+{.03\over N}.                                  \tag{18}
\]
By (4), every coordinate interval
\([-e^{-c}x,x]\), for the receiver cube and critical actions under
consideration, is contained in this full source interval.  Thus every
conditional angular fibre of the actual rectangle lies in the positive
part of the source, and
\[
 F_N\ge(1+.03/N)^N\ge1.03,
 \qquad K_c^s\Phi_N(F_N)\ge\Phi_N(1.03)\quad(s>0).
\tag{19}
\]
For a wide window whose lower endpoint is
\(N\mu-7\sqrt{N\log N}\), this implies
\[
 j_W^c\Phi_N(F_N)
 \ge \Phi_N(1.03)\int_W e^{-s}\,ds.
\tag{20}
\]
With \(\nu=\log(1+\log1.03)>0\),
\(\log\Phi_N(1.03)=\nu N+O(1)\).  Comparing (20) with the critical
reference scale
\(J_{\rm ref}=\exp\{-N\mu+O(\sqrt{N\log N})\}\) gives
\[
 \log{j_W^c\Phi_N(F_N)\over J_{\rm ref}}
 \ge\nu N+O(\sqrt{N\log N}).                             \tag{21}
\]
It therefore cannot satisfy
\(j_W^c\Phi_N(F_N)\le J_{\rm ref}N^\varepsilon\) for any fixed
\(\varepsilon\).  The r072 cost screen removes this example completely
from the efficiency-qualified receiver class.

The proof only establishes the unique
maximizer inside \([c_r-\rho,c_r+\rho]\); a larger action range may access
the source differently.  Accordingly it is a canonical-source realization
of the Jacobian obstruction, not a full strong-maximal optimization or a
counterexample to any qualified global upper bound.
