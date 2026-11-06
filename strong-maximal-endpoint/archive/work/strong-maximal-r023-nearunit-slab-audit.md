# r023 audit: near-unit cube/slab certificate

Date: 2026-11-06.  This note audits only the stated cube construction.  It is
not an upper bound for the strong maximal constant and does not compare to the
larger hyperbolic-profile mechanism.

Let \(Q=[0,1]^n\), \(E=[0,r]^n\), and
\[
 f_\varepsilon=1_Q+\varepsilon 1_{Q\setminus E},\qquad
 \kappa(r)=r+(1-r)\log(1-r),\qquad
 A_{n,r}=1-r\kappa(r)^{n-1}.
\]
For a point just outside the \(i\)-th upper face, write \(x_i=1+u\) and
write its other coordinates as \(y_j\in[0,1]\).  The rectangle
\[
 [0,1+u]\times\prod_{j\ne i}[y_j,1]
\]
has average
\[
 \frac{1+\varepsilon G_i(y)}{1+u},\qquad
 G_i(y)=1-r\prod_{j\ne i}\frac{(r-y_j)_+}{1-y_j}.
\]
The lower-face rectangle \([-u,1]\times\prod_{j\ne i}[y_j,1]\) gives the
same expression.  Thus \(0<u<\varepsilon G_i(y)\) is already a *strict*
level-one witness.  Moreover every point of \(Q\) is a strict witness, by
averaging over \(Q\).  The single-exterior-coordinate face sets are disjoint,
and
\[
 \int_0^1\frac{(r-y)_+}{1-y}\,dy=\kappa(r).
\]
Consequently the exact witness certificate and exact Orlicz denominator are
\[
 |\{M_Sf_\varepsilon>1\}|\ge 1+2n\varepsilon A_{n,r},\tag{A}
\]
\[
 D_{n,r}(\varepsilon)=1+(1-r^n)[\Phi_n(1+\varepsilon)-1].\tag{B}
\]
If a literal positive margin is wanted, use
\(0<u<(1-\eta)\varepsilon G_i(y)\); no trim is logically needed.

## Fixed dimension expansion

Since \(\Phi_n'(1+)=n\), the normalized lower certificate is
\[
 \frac{1+2n\varepsilon A_{n,r}}{D_{n,r}(\varepsilon)}
 =1+n\{1+r^n-2r\kappa(r)^{n-1}\}\varepsilon+O_n(\varepsilon^2).\tag{C}
\]
Thus the coefficient in the normalized ratio contains the \(+1\) displayed
in (C).  The expression \(n[r^n-2r\kappa^{n-1}]\) is (C)'s coefficient
minus \(n\), and hence needs an explicitly stated baseline subtraction.
For \(n=2,r=1/4\),
\[
 2\{r^2-2r\kappa(r)\}=\frac{6\log(4/3)-1}{8}
 =0.0907615543\ldots.
\]
Here `log 4/3` must mean \(\log(4/3)\), not \((\log4)/3\).
Wolfram 15 symbolically expanded (A)/(B) and returned exactly (C).

## The \(\varepsilon=a/n\) near-unit regime

Set \(r_n=1-c_n/(n\log n)\), with \(c_n\to\infty\) and
\(c_n/\log n\to0\).  Then
\[
 r_n^n\to1,\qquad (n-1)\log\kappa(r_n)=-c_n(1+o(1)),\qquad A_{n,r_n}\to1.
\]
The middle assertion is logarithmic equivalence; the stronger multiplicative
claim \(r_n\kappa(r_n)^{n-1}\sim e^{-c_n}\) requires, for example,
\(c_n\log\log n/\log n\to0\).  For each fixed \(a>0\),
\[
 \Phi_n(1+a/n)\to e^a,qquad D_{n,r_n}(a/n)\to1,qquad
 \frac{1+2aA_{n,r_n}}{D_{n,r_n}(a/n)}\to1+2a.\tag{D}
\]
In particular \(a=1\) yields the strict slab certificate \(3+o(1)\).
The constant-amplitude rectangular level-one witness has volume
\(\Phi_n(1+a/n)\to e^a\), so this slab geometry exceeds it at \(a=1\).

## Larger amplitudes: a sharp scope-limited obstruction

This is only an upper bound on the certificate (A)/(B), not on all witnesses.
Put \(\delta=1-r\), \(t=n\varepsilon\), \(P=\Phi_n(1+\varepsilon)\),
and \(p=1-r^n\).  Elementary Bernoulli and exponential estimates give
\[
 A_{n,r}\le n\delta(1-\log\delta),\qquad
 p\ge\frac{n\delta}{1+n\delta}.\tag{E}
\]
If \(\delta\ge1/(2n)\), then \(p\ge1-e^{-1/2}\) and
\[
 \frac{1+2tA_{n,r}}{1+p(P-1)}
 \le\frac{1+2t}{(1-e^{-1/2})P}.
\]
If \(\delta<1/(2n)\) and the normalized certificate is greater than one,
then \(2tA_{n,r}>p(P-1)\).  By (E),
\[
 \delta<\exp\!\left(1-\frac{P-1}{4t}\right),
\]
and monotonicity of \(\delta(1-\log\delta)\) yields the quantitative bound
\[
 0<\frac{1+2tA_{n,r}}{1+p(P-1)}-1
 \le\frac{n(P-1)}2\exp\!\left(1-\frac{P-1}{4t}\right).\tag{F}
\]
Thus a positive gain is not literally impossible: an extremely tiny shell can
give one, but (F) makes it super-small whenever \(P/t\gg\log n\).
For both \(\varepsilon=1/\sqrt{n\log n}\) and
\(\varepsilon=\sqrt{\log n/n}\), \(\log P\asymp n\varepsilon=t\), so
\(P/t\gg\log n\) and the single-coordinate slab effect is super-small.
This leaves open higher-codimension exterior sets, other inputs, other levels,
and every question about \(C_n\).

## All-codimension near-unit halo at \(\varepsilon=1/n\)

The stronger extension is correct.  Let a witness point have precisely the
coordinates in \(K\), \(|K|=k\), outside \(Q\), with outside extensions
\(u_i>0\), and use the full input interval in these coordinates.  In the
remaining coordinates use \([y_j,1]\).  The average is
\[
 \frac{1+\varepsilon G_k(y)}{\prod_{i\in K}(1+u_i)},\qquad
 G_k(y)=1-r^k\prod_{j\notin K}\frac{(r-y_j)_+}{1-y_j}.
\]
For a fixed \(y\), its strict exterior volume in that orthant is
\[
 H_k(\varepsilon G_k(y))
 =\frac1{(k-1)!}\int_0^{\log(1+\varepsilon G_k(y))}e^s s^{k-1}\,ds.
\]
The cases of different \(K\) and different signs partition the exterior by
its set of outside coordinates, hence do not overlap.  There are
\(2^k\binom nk\) choices.

Take \(r_n=1-1/(n\sqrt{\log n})\) and \(\varepsilon=1/n\).  For fixed
\(k\),
\[
 \mathbb E(1-G_k)=r_n^k\kappa(r_n)^{n-k}\to0,
 \qquad n^kH_k(G_k/n)\to G_k^k/k!
\]
in expectation (bounded convergence is available).  Thus the (k)-face
contribution tends to \(2^k/(k!)^2\).  Fatou gives the lower limit
\[
 \liminf |\{M_Sf_{1/n}>1\}|\ge
 \sum_{k\ge0}\frac{2^k}{(k!)^2}=I_0(2\sqrt2).
\]
For the matching upper limit, use
\(f_{1/n}\le(1+1/n)1_Q\).  The entire level-one halo of the latter is
exactly
\[
 V_n=\sum_{k=0}^n2^k\binom nkH_k(1/n).
\]
Since \(H_k(1/n)\le e^{1/n}n^{-k}/k!\), its summands are dominated by
\(e,2^k/(k!)^2)\), while each fixed summand tends to \(2^k/(k!)^2\).
Dominated convergence gives \(V_n\to I_0(2\sqrt2)\).  The input denominator
tends to one because \(1-r_n^n\to0\) and
\(\Phi_n(1+1/n)\to e\).  Therefore this cube-shell construction has the
sharp asymptotic normalized halo
\[
 \frac{|\{M_Sf_{1/n}>1\}|}{\int\Phi_n(f_{1/n})}
 \longrightarrow I_0(2\sqrt2)=4.2523508795026238\ldots.
\]
The strict inequalities in the lower witness can again be retained by a
factor \(1-o(1)\) in \(\varepsilon G_k\).  This result is still confined to
the near-unit cube-shell family.

## r024: full capped-cube class and constant rectangles

Let \(\mathcal C_n\) be the class of functions supported in \(Q\) satisfying
\(1\le f\le1+1/n\) on \(Q\).  Their Orlicz mass obeys \(J(f)\ge1\), while
positivity and domination give
\[
 \{M_Sf>1\}\subseteq\{M_S[(1+1/n)1_Q]>1\}.
\]
The right-hand halo has volume \(V_n\) above, hence tends to
\(S:=I_0(2\sqrt2)\).  The shell construction has both halo tending to \(S\)
and \(J\to1\).  Therefore, for the normalized level-one quotient,
\[
 \sup_{f\in\mathcal C_n}\frac{|\{M_Sf>1\}|}{J(f)}\longrightarrow S.
\]
This is a statement about the supremum and its limit; it does not say every
member has halo tending to \(S\).

For a constant-amplitude cube/rectangle \(A1_Q\) with
\(A=1+a/n\), \(0\le a\le1\), the normalized quotient has limit
\[
 F(a)=e^{-a}\sum_{k\ge0}\frac{(2a)^k}{(k!)^2}.
\]
It is increasing on \([0,1]\).  If \(S(a)=\sum_{k\ge0}(2a)^k/(k!)^2\),
then
\[
 S'(a)-S(a)=1-\sum_{k\ge2}\frac{k-1}{k+1}\frac{(2a)^k}{(k!)^2}>0.
\]
For example, at \(a\le1\), the \(k=2\) term is at most \(1/3\), and the
remaining unweighted tail is at most
\((2/9)/(1-1/8)=16/63\); their sum is less than one.  Hence its best capped
constant-rectangle limit is
\[
 F(1)=e^{-1}I_0(2\sqrt2)=1.5643524652\ldots.
\]
Wolfram 15 independently evaluated \(S(1)=4.2523508795\ldots\),
\(F(1)=1.5643524652\ldots\), and \(F'(1)=0.1976672593\ldots\).

Finally, if a cube-supported function satisfies \(0\le f\le1+\varepsilon\),
then its level-one halo is bounded by the constant-cube halo.  With
\(\ell=\log(1+\varepsilon)\),
\[
 |\{M_Sf>1\}|\le\sum_{k=0}^n2^k\binom nkH_k(\varepsilon)
 \le2\sum_{k\ge0}\frac{(2n\ell)^k}{(k!)^2}.
\]
The logarithm of the latter is \(O(\sqrt{n\ell})=o(n)\) for
\(\varepsilon=1/\sqrt{n\log n}\) and
\(\varepsilon=\sqrt{\log n/n}\).  This is only a cube-support/near-flat
amplitude obstruction; it says nothing analogous for arbitrary supports.
