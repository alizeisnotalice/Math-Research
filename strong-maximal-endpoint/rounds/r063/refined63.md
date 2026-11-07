# Refined fixed-block response without monotonicity

## Verdict

For a fixed number of equal blocks, the proposed `1/N` refinement is valid
for every bounded measurable terminal satisfying the actual critical budget
and the exact success condition \(EH\ge1\).  It does not need coordinate
monotonicity.  The key is a budget-probability envelope plus an exact
Gaussian tilted CDF integral.  This proof uses two logarithmic buffers for a pointwise off-exceptional-set core bound; alternatively the direct budget integral in nonmono63 (3) yields the needed core expectation with one buffer.

Let \(N=rb\), with fixed `r`, let the blocks be sums of independent
\(Z_c=\operatorname{Exp}(1)+c\operatorname{Bern}((1+e^c)^{-1})\), and
assume \(|c-\mu|=O(\sqrt{\log N/N})\), \(0\le H\le A\),
\[
 J(H)\le e^{-s_N},\quad EH\ge1,\quad
 s_N=N\mu-d_N+O(\log N),                                  \tag{1}
\]
where
\[
 d_N=\sqrt{\mu N}z_N,\qquad
 \overline\Phi(z_N)=\sqrt{\mu/N}\,z_N.                    \tag{2}
\]
Then each coordinate location response satisfies
\[
 G_{d_i}=f_{N,c}(s_N)+\frac1N+o(N^{-1}).                    \tag{3}
\]
For a common action of the fixed equal blocks, the ordinary plus-residual
identity consequently gives
\[
 \frac{d^2}{dc^2}E_cH
 =-p_\mu\left\{\frac{z_N\varphi(z_N)}{\overline\Phi(z_N)}+1\right\}
 +o(1).                                                      \tag{4}

For a near-success level \(EH\ge1-\eta\), the same argument has an
additional location error `O(eta rho)`.  Thus `eta=o((N rho)^{-1})` is
needed for (3); the weaker `eta=o(rho)` only gives the leading response. The already-proved first-exit bootstrap from one successful action gives η=O(L²/N) on the whole fixed common window, which satisfies this stronger condition. Thus (3)–(4) hold uniformly on that entire window.

## 1. Cramér CDF without a Mills error

Put \(L=\log N\).  The common-tilt local CLT in `r061/local_md61` has the
global tilted density bound
\[
 \sup_y\left|\sigma_\theta\sqrt N f_\theta(Nm_\theta+
 \sigma_\theta\sqrt N y)-\varphi(y)\right|\le C N^{-1/2}.\tag{5}
\]
For a lower-tail point `q`, let \(\theta<0\) solve \(m_\theta=q/N\), put
\(\lambda=-\theta\sigma_\theta\sqrt N\), and let \(I\) be the Cramér
rate.  The exact tilted identity is
\[
 F(q)=e^{-NI}\int_0^\infty e^{\theta v}f_\theta(q-v)\,dv. \tag{6}
\]
The Gaussian part of (6) is evaluated exactly as
\(e^{\lambda^2/2}\overline\Phi(\lambda)\), not by replacing it with
`1/lambda`.  With t=v/(σ_θ√N), the local CLT error in the normalized Laplace integral is at most C/(λ√N), whereas the Gaussian integral is comparable with 1/λ; hence its relative error is O(N^−1/2).  Moreover
\[
 \lambda=z_N+O(z_N^2/\sqrt N),\qquad
 NI-\lambda^2/2=O(z_N^3/\sqrt N).
\]
Consequently, uniformly for fixed `r`, \(|c-\mu|=O(\rho)\), and
\(w=O(L)\),
\[
 \begin{aligned}
 F_{N,c}(N\mu-d_N+w)&=\overline\Phi(z_N)
 \{1+O(L^{3/2}/\sqrt N)\},\\
 f_{N,c}(N\mu-d_N+w)&=\frac{\varphi(z_N)}{\sqrt{\mu N}}\{1+O(L^{3/2}/\sqrt N)\}.
 \end{aligned}                                             \tag{7}
\]
The common-parameter mean shift and the slack in (1) only add the same
order for fixed `r`.  Thus
\[
 F_{N,c}(s_N+O(L))=d_N/N+O(L^2/N),
 \qquad f_{N,c}(s_N+O(L))=f_{N,c}(s_N)+o(N^{-1}).             \tag{8}
\]

The same identity proves the required overshoot estimate.  If
`q=s_N-O(L)`, then
\[
 E(q-S)_+=e^{-NI}\int_0^\infty v e^{\theta v}f_\theta(q-v)dv
 \le e^{-NI}\frac{C}{\sqrt N\theta^2}=O(1),               \tag{9}
\]
because \(e^{-NI}=O(L/\sqrt N)\) and \(\theta^2\asymp L/N\).

## 2. Probability envelope and scalar comparison

The joint block density is at most one: each increment density is at most
one, and convolution preserves this bound.  Choose fixed large buffers
`K0,K1`, set
\[
 B=\{H>1,\ \Phi_N(H)>e^{S-s_N+K_0L}\},
 \qquad q=s_N-(K_0+K_1)L.                                  \tag{10}
\]
Since the integrand of `J` on `B` is at least
\(N^{K_0}e^{-s_N}\),
\[
 P_T(B)\le |B|\le N^{-K_0}.                               \tag{11}
\]
On \(B^c\cap\{S\le q\}\), \(H\le N^{-K_1}\), whence
\[
 E[H1_{S\le q}]=O(N^{-K_1})+O_A(N^{-K_0}).                  \tag{12}
\]

Define
\[
 U(S)=\min\{A-1,[\Phi_N^{-1}(e^{S-s_N+K_0L})-1]_+\}.       \tag{13}
\]
With \(P=(H-1)_+\), \(D=(1-H)_+\), one has \(P\le U\) off `B`.  On
the uncapped branch,
\[
 U'(x)=\frac{(1+U)(1+\log(1+U))}{N+\log(1+U)}.             \tag{14}
\]
The cap event is exponentially remote.  On the central box, with δ=(S−s_N+K0 L)_+, inversion gives U=δ/N+O(δ²/N²+δ/N²). The overshoot bound (9) gives Eδ=d_N+O(L), while Eδ²=O(NL); the cap and a polynomial tail handle the complement. Consequently
\[
 EU=d_N/N+O(L/N).                                           \tag{15}
\]
Exact success implies `EP>=ED`; (12) and (8) give
\(ED\ge F(q)-O(N^{-\min(K_0,K_1)})\).  Therefore, using `P<=U` off `B`,
\[
 \|P-U\|_1=O(L^2/N),\qquad
 \|D-1_{S\le q}\|_1=O(L^2/N).                              \tag{16}
\]
The second bound includes the strip between `q` and the envelope cutoff;
its mass is `O(L f_{N,c}(s_N))=O(L^2/N)`.

## 3. The `+1/N` response term

Let \(\ell_i=f_{b,c}'/f_{b,c}\).  On the fixed-block central box
\(|\ell_i|=O(\rho)\), and its complement is negligible by the Fisher
bound and a polynomial tail.  Thus (16) yields
\[
 E[D\ell_i]=E[1_{S\le q}\ell_i]+o(N^{-1})=f_{N,c}(q)+o(N^{-1}).\tag{17}
\]
The scalar envelope is absolutely continuous, so integration by parts
gives
\[
 E[U(S)\ell_i]=-EU'(S)=-N^{-1}+o(N^{-1}).                  \tag{18}
\]
Indeed `U=O(rho)` except on a negligible tail and the omitted lower tail
has probability `O(rho)`.  Replacing `P` by `U` in (18), using (16), gives
\[
 E[P\ell_i]=-N^{-1}+o(N^{-1}).                              \tag{19}
\]
Finally translation differentiation of the Gamma--Binomial product gives
\[
 G_{d_i}=-E[H\ell_i]=E[D\ell_i]-E[P\ell_i].                \tag{20}
\]
Equations (8), (17)--(20) prove (3).

## 4. Common curvature and the explicit normal constant

Summing (3) over the fixed equal blocks and using the ordinary common-mode
plus-residual identity from `r062` gives
\[
 \frac{d^2}{dc^2}E_cH=Nm''(c)f_{N,c}(s_N)+m''(c)+o(1).       \tag{21}
\]
By (2), (7), and \(m''(c)=-p_\mu+o(1)\), this is (4).
The canonical equation and the normal-tail expansion also give
\[
 z_N^2=L-2\log L-\log(2\pi\mu)+o(1),
 \qquad \frac{z_N\varphi(z_N)}{\overline\Phi(z_N)}
 =z_N^2+1+o(1).
\]
Hence the bracket in (4) equals
\[
 L-2\log L-\log(2\pi\mu)+2+o(1).                          \tag{22}
\]
This is a fixed-block auxiliary refinement only; it supplies neither the
action-dependent angular bridge nor a general noncentered-rectangle upper
bound.
