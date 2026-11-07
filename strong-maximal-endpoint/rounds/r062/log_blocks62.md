# Logarithmically many balanced blocks: a common-parameter concavity lemma

## Scope and notation

This is a fixed-terminal result for a growing, but logarithmic, number of
equal blocks. It does not prove a strong-maximal upper bound or the missing
physical-parameter bridge.

Let
\[
p_c=(1+e^c)^{-1},\quad Z_c=E+cJ,\quad m(c)=1+cp_c,
\quad \Phi_N(h)=h(1+\log^+h)^{N-1},
\]
where \(E\sim\operatorname{Exp}(1)\) and
\(J\sim\operatorname{Bernoulli}(p_c)\) are independent.  Let
\(\mu=1+e^{-\mu}\), \(L=\log N\), and \(\rho=\sqrt{L/N}\).  Define the
canonical critical deficit by
\[
d_N=\sqrt{\mu N}\,z_N,\qquad
\overline\Phi(z_N)=d_N/N=\sqrt{\mu/N}\,z_N.                \tag{1}
\]
Thus \(z_N^2=L-2\log L+O(1)\), \(d_N\asymp\sqrt{NL}\), and the
critical lower-tail density is \((1+o(1))L/N\).

Fix \(C,C_0,K,A<\infty\), and choose a fixed sufficiently large \(K_{\rm cut}>10\), let \(N=rb\) be integral with
\(1\le r\le CL\), and put
\[
\mathcal I_N=\{c=(c_1,\ldots,c_r): |c_i-\mu|\le K\sqrt{L/b}
=K\sqrt r\,\rho\}.                                       \tag{2}
\]
For `c` in this box, let \(T_i\) be independent sums of \(b\) copies of
\(Z_{c_i}\).  Let one fixed terminal \(H:\mathbb R_+^r\to[0,A]\) be
coordinatewise nondecreasing and satisfy
\[
J(H):=\int e^{-\sum_it_i}\Phi_N(H(t))\,dt\le e^{-s_N},
\qquad |s_N-(N\mu-d_N)|\le C_0rL.                          \tag{3}
\]
Write \(\mathcal G(c)=\mathbb E_cH(T)\).  Assume
\[
\mathcal G(c_*)\ge1\quad\hbox{for some }c_*\in\mathcal I_N.\tag{4}
\]
We prove
\[
\boxed{\nabla_c^2\mathcal G(c)
\preceq-(p_\mu-o(1))\frac LrI_r\quad(c\in\mathcal I_N).}\tag{5}
\]

### Finite-mixture regularity

All terminal differentiations below are ordinary density differentiations,
not KKT regularity.  Let
\[
 g_b(x)=\frac{x^{b-1}e^{-x}}{\Gamma(b)}\mathbf1_{\{x\ge0\}},
 \qquad w_{b,k}(c)=\binom bk p_c^k(1-p_c)^{b-k}.
\]
The exact one-block density is the finite shifted-Gamma mixture
\[
 f_{b,c}(t)=\sum_{k=0}^b w_{b,k}(c)g_b(t-ck).               \tag{R1}
\]
For \(b\ge4\), \(g_b,g_b',g_b''\in L^1(\mathbb R)\) and their endpoint
traces needed in integration by parts vanish.  Differentiating the finite
sum shows that \(c\mapsto f_{b,c}\) is \(C^2\) into \(L^1\) on compact
`c` intervals; spatially it lies in \(W^{1,1}\) with \(f_{b,c}(0)=0\).
The product density has the same coordinatewise properties.  Hence both
location and parameter derivatives of \(\mathbb E_cH(T)\) are legitimate
for every bounded measurable `H`.

## 1. Triangular-array moderate deviations

Put \(S=\sum_iT_i\) and \(q=s_N-K_{\rm cut}L\). Uniformly in (2),
\[
\mathbb P_c(S\le q)=\frac{d_N}{N}
\left[1+O\!\left(L^{-1}+\frac{rL^{3/2}}{\sqrt N}\right)\right],
\qquad f_{S,c}(q)=\frac LN(1+o(1)).                         \tag{6}
\]

For completeness, this is a genuinely uniform triangular-array version of
`r061/local_md61`.  The exact cancellation \(m'(\mu)=0\) gives
\[
\sum_{i=1}^rb\{m(c_i)-\mu\}=O(rL),\qquad
N^{-1}\sum_{i=1}^rb\operatorname{Var}(Z_{c_i})
=\mu+O(\sqrt r\,\rho).                                    \tag{7}
\]
The tilted characteristic functions have bounded third cumulants uniformly.
Moreover their exponential factors give, for a common small tilt,
\[
\prod_{i=1}^r|\chi_{c_i,\theta}(u)|^b
\le\left(1+\frac{u^2}{(1-\theta)^2}\right)^{-N/2}.         \tag{8}
\]
Consequently the three-frequency proof in `r061/local_md61` has the same
local error \(O((1+z^3)/\sqrt N)\) as in the i.i.d. case, as follows.
Set
\[
 K_N(\theta)=\frac1N\sum_{i=1}^rbK_{c_i}(\theta),
 \qquad v_{N,\theta}=K_N''(\theta).
\]
For a common tilt in a fixed compact subinterval of \(( -\infty,1)\),
write \(\chi_{c_i,\theta}\) for the corresponding centered characteristic
function.  Summing the uniform third-cumulant expansion gives, uniformly
for \(|u|\le N^{1/8}\),
\[
 \left|\prod_{i=1}^r
 \chi_{c_i,\theta}\!\left(\frac{u}{\sqrt{Nv_{N,\theta}}}\right)^b
 -e^{-u^2/2}\right|
 \le\frac{C|u|^3}{\sqrt N}e^{-c u^2}.                       \tag{R3}
\]
On \(N^{1/8}<|u|\le\varepsilon\sqrt N\), the small-frequency Gaussian
bound gives an exponentially small integral; beyond that range (8), with a
fixed number of powers retained for integrability, does the same.  Fourier
inversion of (R3) therefore yields, uniformly for \(z\asymp\sqrt L\),
\[
 f_{S,c}(Nm_N-\sigma_N\sqrt N z)
 =\frac{e^{-z^2/2}}{\sigma_N\sqrt{2\pi N}}
 \left\{1+O\!\left(\frac{1+z^3}{\sqrt N}\right)\right\}, \tag{R4}
\]
where \(m_N=K_N'(0)\), \(\sigma_N^2=K_N''(0)\).

The same Fourier inversion, with the phase of modulus one retained, gives
uniformly in the common compact tilt window
\[
\sup_y\left|\sqrt{Nv_{N,\theta}}\,
f_{S,c,\theta}(NK_N'(\theta)+\sqrt{Nv_{N,\theta}}y)
-\varphi(y)\right|\le C N^{-1/2}.
\]
This phase-uniform local CLT justifies integration over the shifted tilted
density in (R5), rather than only evaluation at its mean.

For the integrated tail, choose the common saddle \(\theta<0\) from
\(K_N'(\theta)=x/N\), and put \(I_N=x\theta/N-K_N(\theta)\).  The exact
common-tilt Laplace identity is
\[
 F_{S,c}(x)=e^{-NI_N}\int_0^\infty e^{\theta v}
 f_{S,c,\theta}(x-v)\,dv.                                  \tag{R5}
\]
Insert (R4) at the tilted mean and evaluate the Gaussian Laplace integral
exactly as in `r061/local_md61`, Section 2.  This gives the uniform Mills
formula and its relative error \(O(L^{-1}+(1+z^3)/\sqrt N)\).
Finally, the mean shift from (7), including (3), changes the critical
exponent by \(O(rL^{3/2}/\sqrt N)\), while the variance shift is smaller.
Applying (R4)--(R5) at the canonical level proves (6).

## 2. Envelope and coordinate location derivatives

The upper-orthant argument yields
\[
\Phi_N(H(t))\le e^{S(t)-s_N}.                               \tag{9}
\]
Set \(P=(H-1)_+\), \(D=(1-H)_+\), and
\(I=\mathbf1_{\{S\le q\}}\).  On the core, \(H\le N^{-K_{\rm cut}}\).
For \(H=1+y>1\),
\[
y\le\frac{(S-s_N)_+}{N-1-(S-s_N)_+},                       \tag{10}
\]
because \(\log[1+\log(1+y)]\ge y/(1+y)\).

Choose a fixed `B` sufficiently large and define the per-block box
\[
\mathcal C_B=\{\max_i|T_i-bm(c_i)|\le B\sqrt{bL}\}.        \tag{11}
\]
Bernstein plus a union bound gives \(\mathbb P(\mathcal C_B^c)\le N^{-M}\)
for any prescribed `M`.  The boxwise maximum in (10) is
\(O(\sqrt{rNL})\).  Its first two moment bounds and (6) give
\[
\mathbb EP\le\frac{d_N}{N}(1+\varepsilon_N),\qquad
\mathbb EP^2\le\left(\frac{d_N}{N}\right)^2(1+\varepsilon_N),
\tag{12}
\]
where
\[
\varepsilon_N=C\left\{L^{-1}+\frac{rL^{3/2}}{\sqrt N}
+\sqrt r\,\rho\right\}=o(1).                             \tag{13}
\]
The cap charges the complement of the box.

If \(\mathcal G(c)\ge1-\eta\), \(\eta=o(\rho)\), then
\(\mathbb ED\le\mathbb EP+\eta\).  With the lower core mass from (6),
\[
\|D-I\|_1\le C\rho\varepsilon_N+O(\eta),
\qquad \|H-1\|_2=O(\sqrt\rho).                             \tag{14}
\]

Let \(\ell_i=f_{b,c_i}'/f_{b,c_i}\).  Tilted Fourier inversion gives
\[
|\ell_i|\le C\sqrt r\,\rho\quad\hbox{on }\mathcal C_B,
\qquad \mathbb E\ell_i(T_i)^2\le(b-2)^{-1}.                \tag{15}
\]
The latter follows also by conditional Jensen from the Gamma score.  Since
`P` is coordinatewise nondecreasing, Stieltjes integration by parts gives
\[
\mathbb E[P\ell_i]=-\mathbb E[\partial_iP]\le0.            \tag{16}
\]
The individual location derivative is therefore bounded below by its
negative-part contribution:
\[
G_{d_i}:=\left.\partial_a\right|_{0+}\mathbb E_cH(T+ae_i)
=\mathbb E[D\ell_i]-\mathbb E[P\ell_i].                    \tag{17}
\]
Integration by parts in the `i`th coordinate gives
\(\mathbb E[I\ell_i]=f_{S,c}(q)\).  Combining (14)--(15) on the box and
Fisher Cauchy--Schwarz off it gives
\[
G_{d_i}\ge\frac LN(1-o(1)),                                 \tag{18}
\]
because its relative error is at most
\[
C\left\{\frac{\sqrt r}{L}+
\frac{r^{3/2}L^{3/2}}{\sqrt N}+r\rho+
\frac{\sqrt r\,\eta}{\rho}\right\}.                     \tag{19}
\]
For \(r\le CL\), the first three terms vanish.  This is deliberately
only a lower location-response estimate.

## 3. Propagation from one successful vector

Let \(s_{b,c}=\partial_c\log f_{b,c}\).  The parity refinement of the
tilted Fourier calculation gives, at the tilted mean,
\[
\ell_i=-\theta_i+O(b^{-1}),\qquad
\theta_i=\frac{T_i-bm(c_i)}{b\nu(c_i)}
+O\!\left(\frac{(T_i-bm(c_i))^2}{b^2}\right).              \tag{20}
\]
The Gaussian derivative vanishes by oddness; the weighted Fourier error is
\(O(b^{-1/2})\), so \(f'_{b,c,\theta}(bK_c'(\theta))=O(b^{-3/2})\),
which proves (20).  The exact Euler identity
\[
T_i(1+\ell_i)+c_is_{b,c_i}(T_i)=bm(c_i)-1                    \tag{21}
\]
then gives \(|s_{b,c_i}|=O(L)\) on \(\mathcal C_B\): indeed
\(m(c_i)-\nu(c_i)=O(\sqrt r\rho)\), and the linear and quadratic
terms are respectively \(O(L)\).

For an explicit tail bound, if \(K_i\sim\operatorname{Bin}(b,p_{c_i})\)
and \(Y_i\sim\operatorname{Gamma}(b,1)\) are independent, then
\[
s_{b,c_i}(T_i)=\mathbb E\left[bp_{c_i}-(b-1)\frac{K_i}{Y_i}
\ \middle|\ T_i\right],
\]
and conditional Jensen gives
\[
\mathbb Es_{b,c_i}(T_i)^2
\le\frac{b^2p_{c_i}^2+b(b-1)p_{c_i}(1-p_{c_i})}{b-2}=O(b).   \tag{22}
\]
Thus on the near-level region of (14),
\[
|\partial_{c_i}\mathcal G(c)|
=|\mathbb E_c[(H-1)s_{b,c_i}(T_i)]|\le C\rho L.             \tag{23}
\]
Use the central \(O(L)\) score and \(\mathbb E|H-1|=O(\rho)\) on the
box, and (22), the cap, and Cauchy--Schwarz off it.

Let \(\eta_0=N^{-1/2}=o(\rho)\), and let `U` be the component containing
\(c_*\) of \(\{c\in\mathcal I_N:\mathcal G(c)\ge1-\eta_0\}\).  A straight
path within the box has \(\ell^1\) length at most \(2Kr^{3/2}\rho\).
For any proposed exit point, integrate (23) only up to the *first* exit
along the segment from \(c_*\); that initial segment remains in `U`.
Thus
\[
\mathcal G(c)\ge1-O\!\left(\frac{r^{3/2}L^2}{N}\right).   \tag{24}
\]
The error is `o(eta_0)` for \(r\le CL\), contradicting the defining exit
level \(1-\eta_0\).  The first-exit argument gives \(U=\mathcal I_N\),
so (24) holds on the full box.  It is `o(rho)`,
and the final term in (19) becomes
\[
\frac{\sqrt r\,\eta}{\rho}
=O\!\left(\frac{r^2L^{3/2}}{\sqrt N}\right)=o(1).           \tag{25}
\]
Therefore (18) holds uniformly throughout \(\mathcal I_N\).

## 4. Hessian closure

For one block define the density residual explicitly by
\[
 e_i=\frac{\partial_{c_ic_i}^2f_{b,c_i}}{f_{b,c_i}}
 +bm''(c_i)\frac{\partial_tf_{b,c_i}}{f_{b,c_i}}.
\]
The completed kernel estimates of `r054/full_curvature54.md`, Section 6,
in the natural window \(|c_i-\mu|\le K\sqrt{L/b}\), give
\[
\mathbb Ee_i(T_i)=0,\qquad
\|e_i(T_i)\|_2\le CL^2,\qquad \|s_{b,c_i}(T_i)\|_2\le CL. \tag{26}
\]
Here is the ordinary density proof of these bounds.  On the central band,
the compact-tilt calculation of `full_curvature54.md`, Section 6, gives
\(|s_{b,c_i}|=O(L)\) and \(|e_i|=O(L^2)\).  For the tail, condition on the
Gamma--Binomial representation \(T_i=Y_i+c_iK_i\), where
\(Y_i\sim\operatorname{Gamma}(b,1)\) and
\(K_i\sim\operatorname{Bin}(b,p_{c_i})\) are independent.  The exact
formula used there is
\[
 \frac{R_{c_ic_i}}R(T_i)=\mathbb E\left[
 b^2p_{c_i}^2-bp_{c_i}(1-p_{c_i})
 -2bp_{c_i}(b-1)\frac{K_i}{Y_i}
 +(b-1)(b-2)\frac{K_i^2}{Y_i^2}\ \middle|\ T_i\right].   \tag{R6}
\]
Inverse-Gamma moments through order eight imply
\(\|R_{c_ic_i}/R\|_4=O(b^2)\); the corresponding mixture formulas for
the spatial and parameter scores give \(\|e_i\|_4=O(b^2)\) and
\(\|s_{b,c_i}\|_4=O(b)\).  Since
\(\mathbb P(\mathcal C_B^c)\le N^{-M}\), Hölder gives
\[
 \|e_i\mathbf1_{\mathcal C_B^c}\|_2\le Cb^2N^{-M/4}=o(L^2),
 \qquad
 \|s_{b,c_i}\mathbf1_{\mathcal C_B^c}\|_2\le CbN^{-M/4}=o(L),
\]
after one fixed sufficiently large choice of `M`.  Together with the
central bounds, this proves (26) uniformly in `i` and in the growing block
count.  Finally \(\mathbb Ee_i=0\) follows from normalization and
\(\mathbb E(1-k_i)=0\), again without an `H` or KKT assumption.

The exact one-block density identities are
\[
\partial_{c_ic_i}^2\mathcal G
=bm''(c_i)G_{d_i}+\mathbb E_c[(H-1)e_i],
\qquad
\partial_{c_ic_j}^2\mathcal G
=\mathbb E_c[(H-1)s_{b,c_i}s_{b,c_j}]\quad(i\ne j).          \tag{27}
\]
By (14), \(\|H-1\|_2=O(\sqrt\rho)\).  Block independence and (26) make
every residual or off-diagonal error in (27) at most \(C\sqrt\rho L^2\):
for a cross term use
\(\|s_{b,c_i}s_{b,c_j}\|_2=\|s_{b,c_i}\|_2\|s_{b,c_j}\|_2\).
Since \(m''(c_i)=-p_\mu+o(1)\), (18), and \(b=N/r\),
\[
\partial_{c_ic_i}^2\mathcal G
\le-(p_\mu-o(1))\frac Lr+C\sqrt\rho L^2.                  \tag{28}
\]
Gershgorin bounds the upper spectral error by \(Cr\sqrt\rho L^2\).
Relative to the main scale \(L/r\), it is
\[
Cr^2\sqrt\rho L=o(1)\qquad(r\le C\log N).                \tag{29}
\]
This proves (5).
