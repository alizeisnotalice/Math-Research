# A superexponential obstruction to full pointwise soft column covers

Date: 2026-11-11. This is a theorem about the actual column architecture from r081–r083, not an endpoint lower-bound improvement. All logarithms are natural and all endpoint/cover statements use physical Lebesgue measure.

Let $J_N(F)=\int\Phi_N(F)$, $\Phi_N(t)=t(1+\log^+t)^{N-1}$, and let $H_aG(x)$ be the average over the axis-parallel rectangle with opposite physical corners $a,x$, allowing all signed orientations. Set $\eta=N^{-2}$ and
\[
\rho_F(x)=\chi((M_SF(x)-1)/\eta),\qquad
\chi(t)=\frac{e^{-1/t}}{e^{-1/t}+e^{-1/(1-t)}}\quad(0<t<1),
\]
with $\chi=0$ on $t\le0$ and $\chi=1$ on $t\ge1$. For each $0\le F\le A$ of finite positive modular, minimize $\sum_jw_jJ_N(G_j)/J_N(F)$ over finite or countable families with constant $w_j\ge0$, arbitrary continuous corners $a_j$, $0\le G_j\le F$, and $\sum_jw_j1_{\{H_{a_j}G_j>1\}}\ge\rho_F$ a.e. The supremum over $F$ is $\mathcal W_{N,\eta}(A)$. An empty feasible family has infinite infimum.

## Theorem

For every integer $N\ge1024$,
\[
\boxed{\quad\mathcal W_{N,N^{-2}}(2)\ge
 \frac{(5N^{3/2}/1848)^N}{18N^3}.\quad} \tag{1}
\]
Even the relaxed problem that demands coverage only on $\{M_SF>1+N^{-2}\}$ satisfies
\[
\boxed{\quad\mathcal W^{\rm sat}_{N,N^{-2}}(2)\ge
 \frac{(5N/924)^N}{9N^2}.\quad} \tag{2}
\]
Thus $\log\mathcal W(2)\ge\frac32N\log N-O(N)$ and $\log\mathcal W^{\rm sat}(2)\ge N\log N-O(N)$. Both statements also hold at cap three. They exclude every $\exp(O(N))$ full-cover budget, in particular the proposed matching upper target (11) of r082. They do not show infinite budget in any fixed dimension and do not change the known bounds for $C_N$.

## 1. A source with continuously moving optimal corners

Treat $p=2,3$ together. On $[0,1]$ let $g(t)=1-1/N+t/N$, zero elsewhere, and put
\[
a_p=\frac12+\frac2{N^p},\qquad B_p=e^{a_p},\qquad F_p(y)=B_p\prod_{i=1}^N g(y_i).
\tag{3}
\]
In particular $F_p<33/20<2$. For $z\in[1/4,3/4]$ set
\[
v(z)=1+N(e^{-z/N}-1),\qquad
x(z)=1+\frac{(1-v(z))^2}{2Ng(v(z))}.
\tag{4}
\]
The inequalities $1-u\le e^{-u}\le1-u+u^2/2$ give
$1/4\le v<3/4+9/(32N)<19/25$. We have $g(v)=e^{-z/N}$ and $1<x<1+1/N$. Also
\[
x'(v)=-\frac{1-v}{Ng(v)}-\frac{(1-v)^2}{2N^2g(v)^2}<0,
\quad v'(z)=-g(v)<0.
\]
Consequently $T(z_1,\ldots,z_N)=(x(z_i))$ is a $C^1$ diffeomorphism on the interior parameter cube, with nonzero Jacobian there. Physical a.e. coverage pulls back to $z$-a.e. coverage. The small physical Jacobian does not reduce a pointwise demand under this change of test measure.

For $b\in[0,1)$, direct integration gives the exact identity
\[
K_b(x):=\frac{\int_b^1g(t)dt}{x-b}
=g(v)-\frac{(b-v)^2}{2N(x-b)}. \tag{5}
\]
Any interval containing $x>1$ with positive source mass can be shortened at its right endpoint to $x$, and at a negative left endpoint to zero. Its average is at most $g(v)$, with equality at $[v,x]$. Product factorization therefore yields the actual, unrestricted strong maximum
\[
M_SF_p(Tz)=\exp(a_p-\textstyle\sum_i z_i/N). \tag{6}
\]
There is no finite-anchor restriction in this formula.

## 2. Positive-measure slabs and their cutoff demands

Define
\[
P_p=\left\{z\in[1/4,3/4]^N:
 |\textstyle\sum_i z_i-N/2|\le\frac1{2N^{p-1}}\right\}.
\]
Then
\[
|P_p|_z\ge2^{-N}N^{-p},\qquad
\frac3{2N^p}\le\log M_SF_p(Tz)\le\frac5{2N^p}.
\tag{7}
\]
For the volume assertion, the sum of centered uniforms on $[-1/4,1/4]$ has an even density nonincreasing on the positive half-line. Convolution with a centered uniform preserves this property: the derivative of its interval integral is the right density value minus the left one, nonpositive a.e. The central interval of length $N^{-(p-1)}$ has average density at least the average $2/N$ on the whole support of length $N/2$. This in fact proves the stronger probability lower bound $2/N^p$; (7) deliberately discards its factor two.

For $p=2$, $M_SF_p>1+\eta$, so $\rho_{F_p}=1$ on $T(P_2)$. For $p=3$, put $t=N^2(M_SF_p-1)\ge3/(2N)$. If $0<t<1$, the denominator defining $\chi$ is at most two and $\chi(t)\ge\frac12e^{-1/t}\ge\frac12e^{-2N/3}$. If $t\ge1$, $\chi(t)=1$ gives the same inequality. Hence the respective demands are
\[
w_2=1,\qquad w_3=\tfrac12e^{-2N/3}>2^{-(N+1)}.
\tag{8}
\]
The last inequality follows from $e^{2/3}<2$.

## 3. Complete pricing of all true corners and measurable source pieces

Fix any $a\in\mathbb R^N$ and measurable $0\le G\le F_p$. If $H_aG(Tz)>1$ at some $z\in P_p$, then every $a_i<1$; otherwise the rectangle misses the source in that coordinate. Let $b_i=\max(0,a_i)$. Clipping only removes zero-source portions. By (5)–(6),
\[
1<H_aG(Tz)\le H_bF_p(Tz)
 =M_SF_p(Tz)\prod_i(1-\ell_i),\quad
\ell_i=\frac{(b_i-v_i)^2}{2Ng(v_i)(x_i-b_i)}.
\]
All $0\le\ell_i<1$. Since $-\log(1-\ell_i)\ge\ell_i$,
\[
\sum_i(b_i-v_i)^2\le(2N+2)\sum_i\ell_i
 <\frac{5(1+1/N)}{N^{p-1}}<\frac6{N^{p-1}}.
\tag{9}
\]
Write $R_p=\sqrt6\,N^{-(p-1)/2}<1/10$. In particular $b_i<19/25+1/10<7/8$. Every side of the actual rectangle $[a,Tz]$ is greater than $1/8$, even when $a_i<0$. Strict success consequently forces
\[
J_N(G)\ge\int G>8^{-N}. \tag{10}
\]
The entire success section of this one column on $P_p$ lies in a $v$-ball of radius $R_p$. The Jacobian of $z\mapsto v$ has absolute determinant $e^{-\sum z_i/N}\ge e^{-3/4}$. Thus its $z$-measure is at most
\[
U_{N,p}=e^{3/4}V_N R_p^N,
\qquad V_N=\frac{\pi^{N/2}}{\Gamma(1+N/2)}. \tag{11}
\]
This bound prices every continuous physical corner and every true sub-source. It does not use sampled receivers or a grid of corners.

For any countable feasible cover of $\rho_{F_p}$, integrate its pulled-back inequality on $P_p$. Tonelli applies because the weights and indicator functions are nonnegative. Discard columns with zero-measure success section there; every retained column satisfies (10). If $\mathcal B=\sum_jw_jJ_N(G_j)$, then
\[
w_p2^{-N}N^{-p}\le w_p|P_p|_z
 \le\sum_jw_j|\{z\in P_p:H_{a_j}G_j(Tz)>1\}|_z
 \le U_{N,p}8^N\mathcal B.
\tag{12}
\]
Infinite total cost is harmless. The auxiliary $z$ measure pushes forward to an absolutely continuous physical measure on $T(P_p)$; no weak estimate for arbitrary receiver measures has been substituted for a Lebesgue theorem.

## 4. The actual Orlicz denominator and rational constants

For every $t\ge0$, $\Phi_N(t)\le t+t^N$, because $1+\log t\le t$ when $t\ge1$. Hence
\[
J_N(F_p)\le B_p(1-1/(2N))^N+(B_pI_N)^N,
\quad I_N=\frac N{N+1}\left[1-(1-1/N)^{N+1}\right].
\tag{13}
\]
This is a bound on the actual Orlicz integral, not its exact value. Since $\log(1-u)\ge-u/(1-u)$ for $0<u<1$,
\[
(1-1/N)^{N+1}\ge e^{-1-2/(N-1)}\ge e^{-1-2/1023}>73/200.
\]
Together with $B_p<33/20$ this gives $I_N<127/200$ and
$B_pI_N<4191/4000<21/20$. Therefore
\[
\tfrac12<J_N(F_p)\le3(21/20)^N. \tag{14}
\]
The lower bound follows from $J_N(F_p)\ge\int F_p$ and Bernoulli's inequality $(1-1/(2N))^N\ge1/2$.

The Gaussian integral bound $V_NR^N\le(2\pi e/N)^{N/2}R^N$ follows by integrating $e^{tR^2-t|u|^2}$ and minimizing in $t>0$. Since $e^{3/4}<3$ and $\sqrt{12\pi e}<11$,
\[
U_{N,p}\le3(11/N^{p/2})^N. \tag{15}
\]
The constants are elementary: $e<11/4$, and
$22/7-\pi=\int_0^1t^4(1-t)^4/(1+t^2)dt>0$, so $12\pi e<121$. Dividing (12) by (14) and using (8), (15) proves (1), (2) exactly.

A stronger finite-dimensional certified expression keeps the ball volume:
\[
\mathcal W\ge L_{N,p}:=
\frac{w_p2^{-N}N^{-p}}
 {e^{3/4}V_N(6/N^{p-1})^{N/2}8^N\,3(21/20)^N}.
\tag{16}
\]
In the finite certificates, replace $\pi$ by its rational upper bound $22/7$ and evaluate $(N/2)!$ exactly for the tested even dimensions. Formula (13) additionally gives sharper floating diagnostics, separately labelled as numerical.

## 5. Small physical residual and remote gluing

The obstructing physical receiver region is contained in $(1,1+1/N)^N$. Thus
\[
|T(P_p)|\le N^{-N}\le2N^{-N}J_N(F_p). \tag{17}
\]
It can be discarded by a separately controlled physical-volume residual. No estimate here bounds the residual of a general input.

The same local obstruction can be attached to cap-two almost extremizers. Take compact nonzero $H\le2$ of hard score $s$; spatially dilate $F_p$ to a packet $F'$ with $J_N(F')=J_N(H)/N$. Translate this packet far from $H$ in every coordinate. Any rectangle containing a copied receiver and meeting $\operatorname{supp}H$ then has volume greater than $\|H+F'\|_1$, so it cannot succeed for any $G\le H+F'$. At copied slab receivers every successful column can consequently be restricted to $\operatorname{supp}F'$ without losing **success on that copied slab** or increasing its modular. Applying the translated/dilated local counting argument (12), rather than asserting any new full-space cover of $F'$, gives
\[
\frac{\mathcal B}{J_N(H+F')}\ge\frac{L_{N,p}}{N+1}.
\tag{18}
\]
The hard score of $H+F'$ is at least $s/(1+1/N)$. Thus global near-efficiency does not rescue full pointwise covers. Compact near-extremizers and the cap-two passage are available from r079 and r073.

## Verification and scope

The proof above holds for every integer $N\ge1024$. Python rational arithmetic and Wolfram 15 independently verify its constants and identities; exact factorial/log interval checks cover $N=2^{10},\ldots,2^{21}$, and exact Irwin–Hall checks cover 20 small-dimensional slabs. These are verification of the formulas, not a replacement of the all-$N$ argument. See [certification84.md](certification84.md).

The obstruction arose in an explicitly authorized native Astra core task; an independent agent and the parent checked the proof. No public priority claim is made. The general $C_N$ logarithmic-coefficient gap remains $(\mu-1)/8$. The remaining plausible column route must allow independently controlled physical residuals or change the charging mechanism.
