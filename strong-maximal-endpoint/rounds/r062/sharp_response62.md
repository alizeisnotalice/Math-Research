# A sharp common-location response for fixed balanced blocks

## Scope and statement

This note proves a local, fixed-block statement. It concerns a common
parameter `c`, a fixed number `r` of equal blocks, and one monotone capped
terminal. It does **not** give an upper bound for the strong maximal
constant, for a varying number of blocks, or for a general family of
rectangles.

Let
\[
p_c=(1+e^c)^{-1},\qquad Z_c=E+cJ,
\]
where \(E\sim\operatorname{Exp}(1)\),
\(J\sim\operatorname{Bernoulli}(p_c)\), independently. Put
\[
m(c)=EZ_c=1+cp_c,\qquad \mu=1+e^{-\mu},
\qquad \Phi_N(h)=h(1+\log^+h)^{N-1}.
\]
Fix integers \(r\ge1\) and \(N=rb\), a cap \(A<\infty\), and constants
\(C_0,K>0\). Constants implicit below may depend on this fixed data, but
not on \(H,c,N\).

Let \(T_1,\ldots,T_r\) be independent block sums, each with the density
\(f_{b,c}\) of a sum of \(b\) independent copies of \(Z_c\). Write
\(S=\sum_jT_j\), so \(S\) has density \(f_{N,c}\). Consider every
coordinatewise nondecreasing measurable \(H:\mathbb R_+^r\to[0,A]\) with
\[
\mathbb EH(T)\ge1,\qquad
J(H):=\int_{\mathbb R_+^r}e^{-\sum_jt_j}\Phi_N(H(t))\,dt\le e^{-s_N},
\tag{1}
\]
where
\[
s_N=N\mu-d_N+O_{C_0}(\log N).\tag{2}
\]
Here \(d_N\) is defined canonically by
\[
 d_N=\sqrt{\mu N}\,z_N,
 \qquad \overline\Phi(z_N)=\frac{d_N}{N}
 =\sqrt{\frac\mu N}\,z_N,                                \tag{3}
\]
where \(\overline\Phi\) is the standard Gaussian upper tail.  Thus
\(z_N^2=\log N-2\log\log N+O(1)\).  The uniform
moderate-deviation consequences of `r061/local_md61` are
\[
d_N\asymp\sqrt{N\log N},\quad
F_{N,c}(q_N)=\frac{d_N}{N}(1+o(1)),\quad
f_{N,c}(q_N)=\frac{\log N}{N}(1+o(1)),\tag{4}
\]
uniformly for
\[
|c-\mu|\le C_0\sqrt{\frac{\log N}{N}},\qquad
q_N=s_N-K\log N.\tag{5}
\]

Define the one-sided common-location response
\[
G_d(H,c):=\left.\frac{d}{da}\right|_{a=0+}
\mathbb EH(T+a\mathbf1).
\]
Then, uniformly in all the preceding data,
\[
\boxed{G_d(H,c)=\bigl(r+o(1)\bigr)\frac{\log N}{N}.}\tag{6}
\]
The proof uses neither a KKT equation nor a smoothness assumption on \(H\).

## 1. Cancellation in the common parameter

Since \(p_c'=-p_c(1-p_c)\),
\[
m'(c)=p_c\{1-c(1-p_c)\}.
\]
The defining equation for \(\mu\) gives
\[
p_\mu=\frac{\mu-1}{\mu},\qquad1-p_\mu=\frac1\mu,
\]
and hence
\[
m(\mu)=\mu,\qquad m'(\mu)=0,\qquad m''(\mu)=-p_\mu.\tag{6}
\]
Taylor's formula gives, throughout (5),
\[
m(c)-\mu=-\frac{p_\mu}{2}(c-\mu)^2+O(|c-\mu|^3),
\qquad N\{m(c)-\mu\}=O(\log N).\tag{7}
\]
Also \(\operatorname{Var}Z_c=\mu+O(\rho_N)\), where
\(\rho_N=\sqrt{\log N/N}\). The mean displacement in (7) changes the
moderate-deviation exponent by \(O(d_N\log N/N)=o(1)\), and the variance
displacement changes it by \(O((d_N^2/N)\rho_N)=o(1)\). Therefore the
uniform tilted-Fourier proof in `r061/local_md61` supplies (4) uniformly
in `c`.

## 2. Envelope and mass rigidity

For every \(t\in\mathbb R_+^r\), coordinate monotonicity of \(H\) and
monotonicity of \(\Phi_N\) imply
\[
J(H)\ge\int_{u_j\ge t_j}e^{-\sum_ju_j}\Phi_N(H(u))\,du
\ge e^{-S(t)}\Phi_N(H(t)).
\]
Thus
\[
\Phi_N(H(t))\le e^{S(t)-s_N}.\tag{8}
\]
Set \(P=(H-1)_+\) and \(D=(1-H)_+\). On \(\{S\le q_N\}\), (8) gives
\(\Phi_N(H)\le N^{-K}\), hence
\[
H\le N^{-K}\quad\hbox{on }\{S\le q_N\}.\tag{9}
\]
Indeed \(\Phi_N(h)=h\) on \([0,1]\) and \(\Phi_N(h)\ge1\) above one.

For \(y\ge0\),
\[
\log\{1+\log(1+y)\}\ge\frac{y}{1+y}.\tag{10}
\]
Indeed, with \(v=\log(1+y)\), the difference is
\(\log(1+v)-1+e^{-v}\); it vanishes at zero and has derivative
\((1+v)^{-1}-e^{-v}\ge0\).
Consequently, if \(H=1+y>1\) and
\(\delta=(S-s_N)_+<N-1\), then (8) yields
\[
\frac{(N-1)y}{1+y}\le\log\Phi_N(1+y)\le\delta,
\qquad y\le\frac{\delta}{N-1-\delta}.\tag{11}
\]

Choose a sufficiently large fixed box-width \(B\), independently of the
cap, and put
\[
\mathcal C_B=\left\{\max_j|T_j-bm(c)|\le B\sqrt{b\log N}\right\}.\tag{12}
\]
Uniform Bernstein bounds for \(Z_c\) give, for every prescribed `M`, a
choice of `B` such that
\[
\mathbb P(\mathcal C_B^c)\le N^{-M}.\tag{13}
\]
On \(\mathcal C_B\), (11) bounds \(P\) by \(C\rho_N\). Off this box use
the cap. Since
\(S-s_N=d_N+(S-Nm(c))+O(\log N)\), (11), (13), and elementary first and
second moment estimates give
\[
\mathbb EP\le(1+o(1))\frac{d_N}{N},\qquad
\mathbb EP^2\le(1+o(1))\frac{d_N^2}{N^2}=O(\rho_N^2).\tag{14}
\]
For clarity, \(\mathbb E(S-s_N)_+=d_N(1+o(1))\) and
\(\mathbb E(S-s_N)_+^2=d_N^2(1+o(1))\); divide (11) by \(N\) and charge
the complement of (12) by (13).

The target condition gives \(\mathbb ED\le\mathbb EP\), while (9) and
(4) give
\[
\mathbb ED\ge(1-N^{-K})F_{N,c}(q_N)
=(1-o(1))\frac{d_N}{N}.\tag{15}
\]
Combining (14)--(15), and using (9) on the lower core once more, yields
\[
\left\|D-\mathbf1_{\{S\le q_N\}}\right\|_{L^1}=o(\rho_N).\tag{16}
\]
No stronger unweighted error estimate is needed.

## 3. Blockwise score control

Let
\[
\ell_{b,c}(x)=\frac{f_{b,c}'(x)}{f_{b,c}(x)},\qquad
\ell(T)=\sum_{j=1}^r\ell_{b,c}(T_j).\tag{17}
\]
Then uniformly in the stated parameter window,
\[
\sup_{t\in\mathcal C_B}|\ell(t)|\le C\rho_N,
\qquad \mathbb E\ell(T)^2\le\frac r{b-2}=O(N^{-1}).\tag{18}
\]

To prove the first assertion, let \(K_c(\theta)=\log\mathbb Ee^{\theta Z_c}\).
For a one-block coordinate \(x\) in (12), choose the unique small tilt
\(\theta\) such that \(K_c'(\theta)=x/b\). Then
\[
|\theta|=O\!\left(\sqrt{\frac{\log N}{b}}\right)=O(\rho_N),\tag{19}
\]
and the exact tilted-density identity is
\[
f_{b,c}(x)=e^{-\theta x+bK_c(\theta)}f_{b,c,\theta}(x),
\qquad
\ell_{b,c}(x)=-\theta+\frac{f_{b,c,\theta}'(x)}{f_{b,c,\theta}(x)}.\tag{20}
\]
At the tilted mean, `r061/local_md61` gives the uniform characteristic
function bounds
\[
|\chi_{c,\theta}(u)|\le\left(1+\frac{u^2}{(1-\theta)^2}\right)^{-1/2},
\qquad |\chi_{c,\theta}(u)|\le e^{-c_1u^2}\ (|u|\le\varepsilon).\tag{21}
\]
Fourier inversion therefore gives, without any Edgeworth-parity argument,
\[
\begin{aligned}
|f_{b,c,\theta}'(bK_c'(\theta))|
&\le\frac1{2\pi}\int_{\mathbb R}|u||\chi_{c,\theta}(u)|^bdu\\
&\le C\int_0^\varepsilon ue^{-c_1bu^2}\,du+O(\varrho^b)
\le\frac Cb.\tag{22}
\end{aligned}
\]
The first bound in (21), retaining a fixed number of powers for
integrability, controls the outer integral. The same Fourier split gives
the uniform local lower bound
\[
f_{b,c,\theta}(bK_c'(\theta))\ge c_1b^{-1/2}.\tag{23}
\]
Equations (19)--(23) prove the pointwise part of (18), block by block.

For the second assertion in (18), condition on the Bernoulli count in a
block. The location score of the resulting Gamma\((b,1)\) component has
squared mean \(1/(b-2)\). The mixture score is its conditional expectation,
so conditional Jensen gives the same upper bound for one block. Independence
and zero mean scores finish the proof of (18).

### Gamma--Binomial regularity

The differentiations used below are legitimate for an arbitrary bounded
measurable terminal, not only for a KKT response.  Indeed, with
\[
 g_b(x)=\frac{x^{b-1}e^{-x}}{\Gamma(b)}\mathbf1_{\{x\ge0\}},
 \qquad w_{b,k}(c)=\binom bk p_c^k(1-p_c)^{b-k},
\]
the exact density is
\[
f_{b,c}(t)=\sum_{k=0}^b w_{b,k}(c)g_b(t-ck).                \tag{R}
\]
For \(b\ge4\), \(g_b,g_b',g_b''\in L^1(\mathbb R)\), with the relevant
axis traces equal to zero.  Termwise differentiation of the finite sum
shows that \(c\mapsto f_{b,c}\) is \(C^2\) as an \(L^1\)-valued map on
every compact parameter interval, while \(f_{b,c}\in W^{1,1}\) spatially
and \(f_{b,c}(0)=0\).  The product density inherits these properties.
Thus translation and parameter differentiation can be passed through
\(\int H\,d\mathbb P_c\) for bounded measurable \(H\), with no
regularity of \(H\) and no optimization hypothesis.

## 4. Differentiation and the sharp response

The \(W^{1,1}\) part of the preceding regularity statement permits
translation differentiation against bounded measurable \(H\). After substituting
\(u=t+a\mathbf1\), the axis terms vanish and
\[
G_d(H,c)=-\mathbb E[H(T)\ell(T)]
=\mathbb E[D\ell]-\mathbb E[P\ell].\tag{24}
\]
By (14) and (18),
\[
|\mathbb E[P\ell]|\le\|P\|_2\|\ell\|_2
=O\!\left(\frac{\sqrt{\log N}}{N}\right)
=o\!\left(\frac{\log N}{N}\right).\tag{25}
\]

Put \(I=\mathbf1_{\{S\le q_N\}}\). On \(\mathcal C_B\), (16) and
(18) give
\[
|\mathbb E[(D-I)\ell;\mathcal C_B]|
\le C\rho_N\|D-I\|_1=o(\rho_N^2).
\]
Off this box, \(|D-I|\le1\), so (13) and (18) yield
\[
|\mathbb E[(D-I)\ell;\mathcal C_B^c]|
\le\|\ell\|_2\mathbb P(\mathcal C_B^c)^{1/2}
=O(N^{-(M+1)/2})=o\!\left(\frac{\log N}{N}\right).\tag{26}
\]
Thus
\[
\mathbb E[D\ell]=\mathbb E[I\ell]+o\!\left(\frac{\log N}{N}\right).\tag{27}
\]

For every coordinate \(j\), integration by parts in that axis gives the
exact identity
\[
\begin{aligned}
\int_{\{S\le q_N\}}\ell_{b,c}(t_j)\prod_{i=1}^rf_{b,c}(t_i)\,dt
&=\int_{\{\sum_{i\ne j}t_i\le q_N\}}
f_{b,c}\!\left(q_N-\sum_{i\ne j}t_i\right)
\prod_{i\ne j}f_{b,c}(t_i)\,dt_{-j}\\
&=f_{N,c}(q_N).\tag{28}
\end{aligned}
\]
The lower endpoint is zero because \(f_{b,c}(0)=0\). Summing (28), using
(4), and applying (24)--(27), proves (6).

## 5. Optional two-block common-parameter consequence

This is only a consequence in the two-block residual framework of `r054`.
Let \(r=2\), \(b=N/2\).  With
\(R_{b,c}(t)=e^tf_{b,c}(t)\),
\[
 s=\partial_c\log R_{b,c},\qquad k=\partial_t\log R_{b,c},
 \qquad e=\frac{R_{cc}}R-bm''(c)(1-k),
\]
the **common** second derivative uses the plus residual
\[
 Z_+=e(T_1)+e(T_2)+2s(T_1)s(T_2).                           \tag{29}
\]
This is not the contrast residual with a minus sign.  The ordinary
one-block estimates needed here follow directly from the kernel calculation
in `r054/full_curvature54.md`, Section 6: on its central band,
\(|e|=O((\log N)^2)\), while its explicit conditional Gamma--Binomial
formula for \(R_{cc}/R\), together with inverse-Gamma moments, gives
\(\|e\|_4=O(b^2)\).  The complement has arbitrarily high polynomial
probability, so Hölder yields \(\|e\|_2=O((\log N)^2)\).  Likewise the
Fisher estimate there gives \(\|s\|_2=O(\log N)\).  Independence now gives
\[
 \|Z_+\|_2\le2\|e\|_2+2\|s(T_1)s(T_2)\|_2
 =2\|e\|_2+2\|s\|_2^2=O((\log N)^2),                     \tag{30}
\]
independently of `H`.
The cap and (14)--(15) give
\[
\|H-1\|_2=O(\rho_N^{1/2}),\qquad
|\mathbb E[(H-1)Z_+]|=O((\log N)^2\rho_N^{1/2})=o(1).\tag{31}
\]
The exact common-mode kernel identity isolated in
`r054/full_curvature54.md`, Section 6, and used on the common window in
`r057/common_window57.md`, Section 3, is
\[
\frac{d^2}{dc^2}\mathbb E_cH(T)
=b\,m''(c)G_d+\mathbb E[(H-1)Z_+],\tag{32}
\]
(6)--(7) and (30)--(31) give
\[
\boxed{\frac{d^2}{dc^2}\mathbb E_cH(T)
=-p_\mu\log N\,(1+o(1)).}\tag{33}
\]
Both the identity and the norm estimate are Gamma--Binomial kernel facts:
they are independent of `H` and do not require a KKT characterization.
This still makes no assertion about the unresolved general strong-maximal
upper bound.

## 6. Propagation across the full common-parameter window

This final section upgrades the preceding conclusion from a terminal which
has target value at one common parameter to the entire fixed common window.
It is still a statement about the fixed terminal `H` and the fixed balanced
block model only.

Fix `K_1` and write
\[
 I_N=[\mu-K_1\rho_N,\mu+K_1\rho_N],
 \qquad \mathscr G(c)=\mathbb E_cH(T),                     \tag{33}
\]
where the subscript means that every block is formed from \(Z_c\). Assume
that (1)--(2) hold for the same terminal `H`, and that
\[
 \mathscr G(c_*)\ge1\qquad\hbox{for some }c_*\in I_N.    \tag{34}
\]
The budget in (1) is independent of `c`.

### Sharpened block scores

The estimate in (22) can be improved at the tilted mean.  Make the Fourier
change of variables \(u=v/(\sigma_{c,\theta}\sqrt b)\), and compare the
centered characteristic-function power with its Gaussian approximation.
The uniform third-moment Taylor bound from `r061/local_md61` gives, on the
small-frequency range,
\[
 \left|\chi_{c,\theta}\!\left(\frac{v}{\sigma_{c,\theta}\sqrt b}\right)^b
       -e^{-v^2/2}\right|
 \le Cb^{-1/2}|v|^3e^{-c_2v^2}.                              \tag{35}
\]
The outer and middle ranges are exponentially small by the same split as
in (21)--(22).  In the inversion formula for the derivative, the Gaussian
term has integrand \(-iv e^{-v^2/2}\), which is odd and integrates to zero.
Thus (35), now multiplied by \(|v|\), gives
\[
 |f_{b,c,\theta}'(bK_c'(\theta))|=O(b^{-3/2}),
 \qquad
 \ell_{b,c}(x)=-\theta+O(b^{-1})                            \tag{36}
\]
uniformly for the one-block interval in (12).  This is the only use of
Gaussian parity in the argument.

Let \(s_{b,c}=\partial_c\log f_{b,c}\) and
\(k_{b,c}=\partial_t\log(e^tf_{b,c})=1+\ell_{b,c}\).  The exact Euler
identity for the Gamma--Bernoulli mixture is
\[
 t k_{b,c}(t)+c s_{b,c}(t)=b m(c)-1.                         \tag{37}
\]
Put \(t=bm(c)+\Delta\).  From tilt inversion and (36),
\[
\theta=\frac{\Delta}{b\nu(c)}+O\!\left(\frac{\Delta^2}{b^2}\right),
\qquad \nu(c)=\operatorname{Var}Z_c=m(c)-cm'(c),
\]
and (37) becomes
\[
 c s_{b,c}(t)=-\Delta+t\theta+O(1)
 =\left(\frac{m(c)}{\nu(c)}-1\right)\Delta
   +O\!\left(\frac{\Delta^2}{b}+1\right).                 \tag{38}
\]
At \(\mu\), \(m(\mu)=\nu(\mu)=\mu\); throughout `I_N`,
\(m(c)-\nu(c)=O(\rho_N)\).  On the block box (12),
\(|\Delta|=O(\sqrt{b\log N})\), so (38) proves
\[
 \sup_{t\in\mathcal C_B}|s_{b,c}(t)|=O(\log N),
 \qquad
 \sup_{t\in\mathcal C_B}\left|\sum_{j=1}^rs_{b,c}(t_j)\right|
 =O(\log N).                                                 \tag{39}
\]

### Near-level version of the response lemma

All of Sections 2--4 remain valid if the target condition is weakened to
\[
 \mathbb EH(T)\ge1-\eta_N,\qquad \eta_N=o(\rho_N).          \tag{40}
\]
Indeed, the only changed mass inequality is
\[
 \mathbb ED\le\mathbb EP+\eta_N.
\]
Together with (14)--(15), this gives (16) unchanged:
\[
 \|D-\mathbf1_{\{S\le q_N\}}\|_1=o(\rho_N),
 \qquad \|H-1\|_2=O(\rho_N^{1/2}).                         \tag{41}
\]
The same score decomposition therefore proves (6) under (40), uniformly.
The cap is used in the second assertion in (41).

On a near-level set satisfying (40), differentiate (33) in `c`. Since the
parameter score has mean zero,
\[
 \mathscr G'(c)=\mathbb E_c\left[(H-1)\sum_{j=1}^rs_{b,c}(T_j)\right].
\]
On \(\mathcal C_B\), use (39) and
\(\mathbb E|H-1|=O(\rho_N)\), which follows from (14), (15), and (40).
For an explicit tail estimate, write the Bernoulli count in a block as
\(K\sim\operatorname{Bin}(b,p_c)\), independently of
\(G\sim\operatorname{Gamma}(b,1)\).  Differentiating the Gamma--Binomial
mixture and conditioning on its sum gives the exact score representation
\[
 s_{b,c}(T)=\mathbb E\left[bp_c-(b-1)\frac K G\ \middle|\ T\right].
\]
Conditional Jensen and the inverse-Gamma moments give, for \(b\ge4\),
\[
 \mathbb Es_{b,c}(T)^2
 \le\frac{b^2p_c^2+b(b-1)p_c(1-p_c)}{b-2}=O(b).              \tag{42}
\]
Thus \(\|\sum_js_{b,c}(T_j)\|_2=O(\sqrt N)\).  The cap, (13), and
Cauchy--Schwarz show that the contribution of \(\mathcal C_B^c\) is at
most \(O(\sqrt N\,N^{-M/2})\).  Taking `M` sufficiently large gives
\[
 \sup_{c:\,\mathscr G(c)\ge1-o(\rho_N)}|\mathscr G'(c)|
\le C\rho_N\log N.                                        \tag{43}
\]
All constants in (41)--(43) depend only on the fixed cap, window, and
number of blocks, not on the source terminal.

### Bootstrap to the full interval

Set \(\eta_0=N^{-1/2}=o(\rho_N)\), and let `U` be the connected component
containing \(c_*\) of
\[
 \{c\in I_N:\mathscr G(c)\ge1-\eta_0\}.
\]
The map \(c\mapsto\mathscr G(c)\) is continuous by the total-variation
continuity of the block densities. On `U`, (43) and the length of `I_N`
give
\[
 \mathscr G(c)\ge\mathscr G(c_*)-2K_1C\rho_N^2\log N
 \ge1-C_1\frac{(\log N)^2}{N}.                              \tag{44}
\]
The last error is `o(eta_0)`. Hence no endpoint of `U` in the interior of
`I_N` can lie on the level \(1-\eta_0\); by connectedness, \(U=I_N\).
In particular (44) holds throughout `I_N`. Since
\((\log N)^2/N=o(\rho_N)\), the near-level response statement applies at
every common parameter in the full window. Therefore
\[
 G_d(H,c)=\bigl(r+o(1)\bigr)\frac{\log N}{N}
 \qquad(c\in I_N),                                          \tag{45}
\]
uniformly for every fixed capped monotone terminal satisfying (1), (2),
and (34).

For two equal blocks, combine (45) with (29)--(31). The exact identity
now holds at every `c` in `I_N`, and \(m''(c)=-p_\mu+O(\rho_N)\), so
\[
 \boxed{\mathscr G''(c)=-p_\mu\log N\,(1+o(1))
 \quad\hbox{uniformly for }c\in I_N.}                      \tag{46}
\]
This is the full common-window strong concavity needed for a sharing
argument in this fixed model. It is not a conclusion about a general
physical parameter map or the unresolved strong-maximal upper bound.
