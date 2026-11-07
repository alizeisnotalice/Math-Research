# Exact angular kernel for a deep hyperbolic input

## Result

Within-block angles can be retained with an action-independent raw terminal.
The resulting physical rectangle average has an exact, but action-dependent,
conditional angular kernel. Its action derivatives contain angular transport
and Bernoulli-reweighting terms of indefinite sign. This rules out an
unqualified curvature theorem for all bounded raw terminals. It does not
refute a curvature theorem restricted to critical, near-efficient physical
inputs; that sharper statement remains open. A useful bridge would require
a separate angular rigidity/transport estimate for efficient physical inputs.

This is relevant even if the radial argument is upgraded to remove monotonicity
of its terminal: the obstruction below comes from the retained angular kernel,
not from monotonicity of a radial profile.

## 1. Exact Jacobian and Orlicz cost

Partition the `N` coordinates into blocks \(B_j\), of sizes \(b_j\). On
each open orthant of \([-1,1]^N\), set
\[
a_{ji}=-\log|x_{ji}|,\qquad v_j=\sum_{i\in B_j}a_{ji},
\qquad \omega_{ji}=a_{ji}/v_j\in\Delta_{b_j-1}.
\]
Let \(\sigma_j\) be uniform probability on the simplex; in standard
simplex coordinates, \(d\sigma_j=\Gamma(b_j)d\omega_j\). The full
Lebesgue Jacobian is
\[
 dx=2^N\prod_j e^{-v_j}\frac{v_j^{b_j-1}}{\Gamma(b_j)}
 \,dv_j\,d\sigma_j(\omega_j).                             \tag{1}
\]
This is the exact angular refinement of the radial formula in `r026`,
Section 2.

Fix \(\alpha_j>0\), \(\sum_j\alpha_j=1\), and a bounded raw terminal
\(F(w,\omega)\) that is independent of all actions. At depth `L`, define
the actual input
\[
 f_L(x)=F(v-L\alpha,\omega),\qquad w=v-L\alpha.            \tag{2}
\]
For compact `w` support (or a dominated analogue), its exact full Orlicz
cost is
\[
 \begin{aligned}
 \int_{[-1,1]^N}\Phi_N(f_L)\,dx
 &=K_L\int_{\mathcal D_L} e^{-\sum_jw_j}
 \prod_j\left(1+\frac{w_j}{L\alpha_j}\right)^{b_j-1}\\
 &\quad\times\int\Phi_N(F(w,\omega))\,d\sigma(\omega)\,dw,
 \end{aligned}                                             \tag{3}
\]
where \(\mathcal D_L=\{w_j>-L\alpha_j:1\le j\le m\}\), and
\[
K_L=2^Ne^{-L}L^{N-m}\prod_j\frac{\alpha_j^{b_j-1}}{\Gamma(b_j)}.
\]
Consequently the normalized limit is
\[
\mathcal J(F)=\int e^{-\sum w_j}\int\Phi_N(F(w,\omega))\,d\sigma(\omega)\,dw.
\tag{4}
\]
The matching exact receiver-volume formula is
\[
 |E_L|=K_L\int_Ee^{-\sum u_j}
 \prod_j\left(1+\frac{u_j}{L\alpha_j}\right)^{b_j-1}
 \,du\,d\sigma(\xi).                                     \tag{5}
\]
Hence neither the numerator nor the denominator permits free deletion of
the angles.

## 2. Actual noncentered action and conditional angular kernel

For an anchored interval ratio \(r_{ji}=e^{-c_{ji}}\), the exact log
increment is
\[
 Z_{c_{ji}}=E_{ji}+c_{ji}J_{ji},\qquad
 E_{ji}\sim\operatorname{Exp}(1),\quad
 J_{ji}\sim\operatorname{Bernoulli}(p_{c_{ji}}),
 \quad p_c=(1+e^c)^{-1}.                                   \tag{6}
\]
Equivalently its density is
\[
q_c(z)=\frac{e^{-z}(1+\mathbf1_{\{z>c\}})}{1+e^{-c}},
\qquad z>0.                                                 \tag{7}
\]
This is the logarithm of the absolute value of a uniform point in
\([-e^{-c},1]\), so it is the actual noncentered-rectangle law.

At receiver radial/angle coordinates \(v=L\alpha+u,\xi\), let
\[
 S_j=\sum_{i\in B_j}Z_{c_{ji}},\qquad
 \Theta_{ji}=\frac{(L\alpha_j+u_j)\xi_{ji}+Z_{c_{ji}}}
 {L\alpha_j+u_j+S_j}.                                      \tag{8}
\]
The average of the physical input (2) over the physical anchored rectangle
is exactly
\[
 \boxed{\mathcal A_{L,c}F(u,\xi)=\mathbb E_cF(u+S,\Theta).}\tag{9}
\]

Let \(f_{S,c}\) be the joint density of the block sums. Given one sum,
the exact conditional angular pushforward is
\[
 \begin{aligned}
 \kappa_{j,c}^{V_j,\xi_j}(B\mid s_j)
 &=\frac1{f_{S_j,c}(s_j)}
 \int_{\substack{z_i\ge0\\\sum_{i\in B_j}z_i=s_j}}
 \mathbf1\!\left\{\frac{V_j\xi_j+z}{V_j+s_j}\in B\right\}
 \prod_{i\in B_j}q_{c_{ji}}(z_i)\,dz_{-i},\\
 &\hspace{38mm}V_j=L\alpha_j+u_j .
 \end{aligned}                                             \tag{10}
\]
Thus radial disintegration gives
\[
 \boxed{\mathcal A_{L,c}F(u,\xi)
 =\int f_{S,c}(s)\int F(u+s,\eta)
 \prod_j\kappa_{j,c}^{V_j,\xi_j}(d\eta_j\mid s_j)\,ds.}\tag{11}
\]
Common exponential tilting does not alter (10).  In particular the raw
terminal is action independent, but the conditional kernel is not: it
depends on all \(c_{ji}\) both through the factors \(q_{c_{ji}}\) and the
transport map.

## 3. Transport derivatives

For a smooth raw terminal, fix \(i\in B_j\), write
\(R_j=L\alpha_j+u_j+S_j\), and set
\[
 \mathcal V_{ji}=\partial_{\eta_{ji}}-
 \sum_{k\in B_j}\eta_{jk}\partial_{\eta_{jk}},
 \qquad D^{\rm fr}_{ji}=\partial_{w_j}+R_j^{-1}\mathcal V_{ji}.
\]
Holding the latent Bernoulli variables fixed gives
\[
\partial_{c_{ji}}(u+S,\Theta)
=J_{ji}\left(e_j,\frac{e_i-\Theta_j}{R_j}\right).
\]
Differentiating their weights too yields the exact formula
\[
 \begin{aligned}
 \partial_{c_{ji}}\mathcal A_{L,c}F
 &=\mathbb E_c[J_{ji}D^{\rm fr}_{ji}F(u+S,\Theta)]\\
 &\quad+p_{c_{ji}}'\{\mathbb E_c[F\mid J_{ji}=1]
 -\mathbb E_c[F\mid J_{ji}=0]\}.
 \end{aligned}                                             \tag{12}
\]
The first line contains radial and angular transport, while the second is
Bernoulli reweighting.  To avoid any convention ambiguity about squaring a
vector field with variable coefficients, the same-coordinate pathwise
second derivative is written explicitly as
\[
 J_{ji}\left\{\partial_{w_jw_j}F
 +2R_j^{-1}\mathcal V_{ji}\partial_{w_j}F
 +R_j^{-2}(\mathcal V_{ji}^2-\mathcal V_{ji})F\right\}.   \tag{13}
\]
The \(-\mathcal V_{ji}\) inside the last term includes the derivative of
the angular vector-field coefficients and of \(R_j^{-1}\).  Equivalently,
if \(D_{ji}=\partial_{w_j}+R_j^{-1}\mathcal V_{ji}\) is treated as its
full variable-coefficient vector field, this expression is \(J_{ji}D_{ji}^2F\).
Add the first/second Bernoulli-weight derivatives. Mixed derivatives in a
block contain the analogous two transport fields and derivatives of
\(R_j^{-1}\). These terms have no fixed sign.

## 4. Explicit obstruction

The direct conditional-radialization example in `r062` is already an
angular-kernel obstruction. For two increments, take `V=0`, condition on
`S=s`, and use the action-independent angular terminal
\[
F(\eta)=\mathbf1_{\{\eta_1<a/s\}}.
\]
Assume 0<a<c<s/2 (for instance s=3μ, a=μ/2, and c in a sufficiently small window about μ). The exact calculation is
\[
\mathbb E_c[F(\Theta)\mid S=s]=\frac{a}{2(s-c)},\qquad
\partial_c=\frac{a}{2(s-c)^2}>0,\qquad
\partial_c^2=\frac{a}{(s-c)^3}>0.                          \tag{14}
\]
Thus conditional common-action curvature can be positive for a raw
terminal independent of the action. Smooth cutoffs preserve the strict
sign. At positive deep radius, a fixed smooth angular profile transports by
`O(1/R_j)`, but an angular layer of width `O(1/R_j)` retains the effect.
The Orlicz cost (3) does not penalize that layer.

There is also a simple unconditioned sign obstruction.  For the monotone
raw terminal \(h(y_1,y_2)=\mathbf1_{\{y_1>a\}}\), and a common action
\(c>a\), the full increment law gives
\[
 \mathbb E_c h(Z_{c,1},Z_{c,2})
 =e^{-a}+(1-e^{-a})p_c,
\]
so
\[
 \partial_c^2\mathbb E_c h
 =(1-e^{-a})p_c(1-p_c)(1-2p_c)>0.                            \tag{15}
\]
Thus even coordinate monotonicity does not create an unqualified
negative-curvature principle for the full law.

More generally, any nonzero signed second derivative of kernel (11) can be
paired with a bounded smooth angular test function of either sign; adding a
constant background makes the terminal nonnegative without changing this
curvature conclusion. Hence there is no universal negative Hessian for all
bounded raw terminals.  Neither example is a critical-budget or
near-efficiency counterexample; the critical-efficient physical theorem is
unresolved.

## 5. Conditional route that remains available

On a region \(R_j\ge R\), compare (12) with the frozen-angular radial
response, obtained by replacing \(\Theta\) by \(\xi\) while retaining the
actual radial sums and Bernoulli weights.  Then (12)--(13) give,
schematically,
\[
 |\text{first angular transport error}|\lesssim
 R^{-1}\|\mathcal V_{ji}F\|_\infty,
\]
\[
 |\text{pure second angular error}|\lesssim
 R^{-2}(\|\mathcal V_{ji}^2F\|_\infty+\|\mathcal V_{ji}F\|_\infty).
 \tag{16}
\]
Only the difference between the actual and frozen-angular Bernoulli
reweighting terms has the same `R^{-1}` bound for an angular-Lipschitz
terminal, since the `J=0,1` images in (8) are that far apart.  The entire
Bernoulli reweighting term is part of the radial action response and need
not be small. BV versions use the Stieltjes form of (12).

Thus a physically meaningful bridge could assume or prove small weighted
transport seminorms, for example
\[
\mathfrak T_1(F)=\int e^{-\sum w}\sum_{j,i}|\mathcal V_{ji}F(w,\omega)|
\,dw\,d\sigma(\omega),
\]
and a second-order analogue, relative to the available deep radii and the
radial curvature scale. But a cap, a near-target average, and the complete
Orlicz quantity (4) do not control these seminorms: angular
measure-preserving rearrangements preserve
\(\int\Phi_N(F(w,\omega))d\sigma\) at each `w`, while arbitrarily fine
angular oscillations make transport variation arbitrarily large.

The remaining route toward a physical noncentered-rectangle upper bound is
therefore an angular rigidity theorem for near-efficient actual inputs, or
a quantitative result that large angular transport forces an Orlicz or
receiving-volume loss. The fixed-radial common-curvature lemma alone cannot
provide this bridge.
