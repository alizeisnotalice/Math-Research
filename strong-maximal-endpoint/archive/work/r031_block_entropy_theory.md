# r031: convolution regularizes the moving-jump entropy cusp

Date: 2026-11-06.  Let \(f_{b,c}\) denote the \(b\)-fold density of
\(Z_c\).  On \(t>0\),
\[
 f_{b,c}(t)=\frac{e^{-t}}{(1+e^{-c})^b\Gamma(b)}
 \sum_{k=0}^b{b\choose k}(t-kc)_+^{b-1}. \tag{1}
\]
For a.e. \(t\),
\[
 s_{b,c}(t)=\partial_c\log f_{b,c}(t)
 =bp_c\left(1-\frac{f_{b-1,c}(t-c)}{f_{b,c}(t)}\right),
 \qquad p_c={e^{-c}\over1+e^{-c}}. \tag{2}
\]

## Fixed blocks: a rigorous regularity result

For each fixed \(b\ge2\), and uniformly on a compact
\(J\Subset(0,\infty)\), \(c\mapsto\sqrt{f_{b,c}}\) is differentiable in
\(L^2(dt)\), with derivative
\[
 \frac12s_{b,c}\sqrt{f_{b,c}}. \tag{3}
\]
In particular the Fisher information
\[
 I_b(c)=\int s_{b,c}^2f_{b,c}<\infty \tag{4}
\]
and the block KL expansion is quadratic:
\[
 D(f_{b,c+h}\Vert f_{b,c})=\frac12I_b(c)h^2+o(h^2). \tag{5}
\]

Proof mechanism.  Away from the finite moving knots \(t=kc\), (1) is
classically differentiable.  At zero, the \(k=0\) term is comparable to
\(t^{b-1}\), with a smooth-in-\(c\) coefficient.  At infinity all terms
are \(e^{-t}\) times polynomials uniformly over \(J\).  The remaining
\(O(|h|)\)-wide knot strips are harmless for \(b\ge2\): for \(b=2\) the
new positive-part contribution changes by \(O(|h|)\), so the squared
error of its difference quotient is \(O(|h|)\); for \(b>2\) it is smaller.
The complementary region is handled by dominated convergence.  This is
quadratic-mean differentiability; it fails at \(b=1\), where the moved
jump has order-one discrepancy on its knot strip.

There are useful effective global bounds.  For \(b=2\), direct use of
(1) gives
\[
\frac{f_{1,c}(t-c)}{f_{2,c}(t)}\le {1+e^c\over c},\qquad
\|s_{2,c}\|_\infty\le\max(2p_c,2/c-2p_c). \tag{6}
\]
For \(b=3\), the polynomial forms give
\[
\frac{f_{2,c}(t-c)}{f_{3,c}(t)}\le {8(1+e^c)\over c},\qquad
\|s_{3,c}\|_\infty\le3p_c+24/c. \tag{7}
\]
These are uniform on \(J\).  Decompose any \(b\) into independent blocks
of size 2, with at most one block of size 3.  If \(\xi_j\) are their
scores, convolution differentiation gives
\[
s_{b,c}(S_b)=E_c[\sum_j\xi_j\mid S_b]. \tag{8}
\]
Since block scores are centered,
\[
\|s_{b,c}\|_\infty\le C_Jb,\qquad
I_b(c)\le E_c(\sum_j\xi_j)^2\le C_Jb. \tag{9}
\]
The first estimate, though crude, supplies a genuinely global score tail
bound; no central-limit assertion is being used for it.

## Stationary large-block Fisher limit

At \(c=\mu\), with \(p=p_\mu=(\mu-1)/\mu\),
\[
m'(\mu)=0,\qquad v(\mu)=\mu,\qquad v'(\mu)=\mu-1=\mu p. \tag{10}
\]
The claim is
\[
\boxed{I_b(\mu)\longrightarrow p^2/2.} \tag{11}
\]
The normal-family value is \((v'/v)^2/2=p^2/2\), but a plain CLT is
insufficient: the ratio in (2) has to be known to order \(1/b\).

Here is a direct common-tilt proof of the required central expansion.  Put
\[
t=b\mu+\sqrt{b\mu}\,z,\qquad m_\eta=t/b,\qquad
d=m_\eta-\mu,
\]
and tilt one increment by \(e^{\eta Z}/M(\eta)\).  Let
\(g_k^\eta\) be the centered tilted \(k\)-sum density.  Exact untilting
gives
\[
\frac{f_{b-1,\mu}(t-\mu)}{f_{b,\mu}(t)}
=\frac{e^{\eta\mu}}{M(\eta)}
 \frac{g_{b-1}^\eta(d)}{g_b^\eta(0)}. \tag{12}
\]
For \(|z|\le\sqrt{K\log b}\), the characteristic function has a uniform
two-term Fourier expansion
\[
g_k^\eta(0)=(2\pi kv_\eta)^{-1/2}
\{1+A_\eta/k+O_K(k^{-3/2})\}, \tag{13}
\]
where \(A_\eta=\kappa_{4,\eta}/(8v_\eta^2)
-5\kappa_{3,\eta}^2/(24v_\eta^3)\).  The same Fourier calculation gives
\[
\frac{g_{b-1}^\eta(d)}{g_{b-1}^\eta(0)}
=1+O_K(|d|/b+d^2/b). \tag{14}
\]
The \(A_\eta\) terms cancel between \(b-1\) and \(b\), while
\[
\frac{e^{\eta\mu}}{M(\eta)}
=1-\frac{\mu\eta^2}{2}+O(\eta^3),\qquad
\eta=\frac{z}{\sqrt{b\mu}}+O(z^2/b). \tag{15}
\]
Thus
\[
\frac{f_{b-1,\mu}(t-\mu)}{f_{b,\mu}(t)}
=1+\frac{1-z^2}{2b}
+O_K((1+|z|^d)b^{-3/2}), \tag{16}
\]
and (2) yields
\[
s_{b,\mu}(t)=\frac p2(z^2-1)
+O_K((1+|z|^d)b^{-1/2}). \tag{17}
\]
This is the needed \(1/b\) density-ratio precision.  The only analytic
input beyond r029's absolute local limit estimate is (13) through the
first Edgeworth correction.  It follows by expanding the tilted
log-characteristic function to fourth order on the Fourier scale
\(u=w/\sqrt k\); the global bound
\(|\varphi_\eta(u)|\le(1+u^2/(1-\eta)^2)^{-1/2}\) controls the complement.
It should be recorded explicitly in a formal proof rather than replaced
by “CLT”.

There is a uniform-integrability argument.  Exponential moments give
\[
P_\mu(|S_b-b\mu|>\sqrt{K b\log b})\le Cb^{-M}
\]
for \(K\) sufficiently large and \(M>3\).  From (9), the score-square
tail is at most \(C b^2\cdot Cb^{-M}=o(1)\).  On the moderate band, (17)
and the tilted local limit expansion converge in \(L^2\) to
\(\frac p2(Z^2-1)\), \(Z\sim N(0,1)\).  Hence
\[
I_b(\mu)\to\frac{p^2}{4}E(Z^2-1)^2=\frac{p^2}{2}. \tag{18}
\]
This tail argument is the required replacement for an unjustified
central-CLT integration.

## Exact moment lower bound: the liminf needs no Edgeworth theorem

There is a short independent proof of the sharp lower limit.  Let
\(X=S_b-bm\), \(Y=X^2-bv\), where \(m,v,\kappa_3,\kappa_4\) are the
one-step mean, variance, third and fourth cumulants.  The score covariance
identities are
\[
 {\rm Cov}(s_{b,c},X)=bm'(c),\qquad
 {\rm Cov}(s_{b,c},Y)=bv'(c). \tag{19}
\]
At \(c=\mu\), \(m'=0\), so Cauchy--Schwarz with
\(Y-(\kappa_3/v)X\) gives
\[
 I_b(\mu)\ge
 {v'(\mu)^2\over
 2v(\mu)^2+\{\kappa_4-\kappa_3^2/v\}/b}. \tag{20}
\]
Indeed,
\[
{\rm Var}\{Y-(\kappa_3/v)X\}
=b\{2bv^2+\kappa_4-\kappa_3^2/v\}.
\]

Wolfram 15.0 symbolically differentiated the one-step log mgf
\[
-\log(1-x)+\log(1-p_c+p_ce^{cx})
\]
and, under \(e^\mu=1/(\mu-1)\), returned
\[
\kappa_1=\mu,\quad \kappa_2=\mu,\quad
\kappa_3=\mu(3-\mu),\quad
\kappa_4=\mu(3-\mu)(4-\mu),\quad v'(\mu)=\mu-1.
\]
Consequently
\[
\kappa_4-\kappa_3^2/v=\mu(3-\mu),\qquad
 I_b(\mu)\ge {p_\mu^2\over2+(3-\mu)/(b\mu)}. \tag{21}
\]
This already proves
\[
\liminf_{b\to\infty}I_b(\mu)\ge p_\mu^2/2
\]
with no local limit or score-tail argument.  The Fourier expansion
(12)--(17), together with the block/Chernoff tail estimate, is needed only
for the matching \(\limsup\).

## Actual numerical diagnostic (not proof)

I independently ran deterministic SciPy piecewise quadrature of (1)--(2)
at \(\mu=1+W(e^{-1})=1.278464542761074\).  It gives
\(p_\mu^2/2=0.0237209695742844\) and

| \(b\) | 2 | 3 | 4 | 8 | 16 | 32 |
|---:|---:|---:|---:|---:|---:|---:|
| \(I_b(\mu)\) | .10963008 | .04790451 | .03832874 | .03063139 | .02723543 | .02549538 |

These are consistent with the limit.  The main MATLAB computation through
\(b=2048\) gives the stronger numerical diagnostic
\(I_{2048}(\mu)=.0237489670\), again above the target.  Neither numerical
calculation substitutes for the Fourier-plus-tail proof above.
