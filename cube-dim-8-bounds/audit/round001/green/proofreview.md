# R001 Green interface: independent analytic audit

2026-10-27. **Verdict: PASS for the analytic Green/source/obstacle interface and its stated sixteen-direction normalization.** This is not acceptance of a new numerical bound on \(C_8\), nor certification of any numerical \(g_0\) or \(I_0\). Source identities and frozen hashes are in `manifest.json`. The parent retains responsibility for the final project proof and adoption.

The author artifact keeps its `upper/round001/green/` source path and is being received with the parent's R002 materials as a backup analytic interface. It is not used to justify the currently accepted numerical \(C_8\) record. The finalized author proof SHA-256 is `67a6314c9d3a0da74c569e0301791e740cb0c03588aececa76a9b0ce5e991cdf`.

I read the new author proof in full and reconstructed its arguments below. Prior low-dimensional PASS labels were not used as mathematical premises. The audit includes an independently written small exact check, actually run with standard-library Fraction arithmetic; it does not execute the author checker or repeat its independent-subset enumeration. All newly authored files are in this audit directory. The previously frozen `audit/round001/upper/` package is unchanged by this work. No source-author files, commits, pushes, subagents, or MATLAB runs are part of this audit.

## 1. Model and generator normalization

Let \(V=\{v_j\}\subset\mathbb R^8\setminus\{0\}\) be finite and spanning, with \(c_j>0\). The directions remain unnormalized. Independent Cauchy coefficients of scales \(c_jt\) give

\[
X_t=\sum_jv_jZ_{j,t},\quad
\psi(\xi)=\sum_jc_j|v_j\cdot\xi|,\quad
P_tf=p_t*f,\quad L=-\psi(D).
\]

For any eight-direction basis \(B\), its partial sum has density

\[
p_t^B(x)=|\det B|^{-1}\prod_{i\in B}
\frac{c_it}{\pi[(c_it)^2+(B^{-1}x)_i^2]}.
\]

Convolution with the remaining probability law proves that the actual density is strictly positive, continuous and vanishes at infinity, and

\[
\|p_t\|_\infty\le Ct^{-8},\quad
C=(\pi^8|\det B|\prod_{i\in B}c_i)^{-1}.
\]

Scaling of every coefficient gives \(p_t(x)=t^{-8}p_1(x/t)\), joint continuity for \(t>0\), an approximate identity, and positive-time \(L^1\) continuity. Symmetry of the independent Cauchy laws gives the self-adjoint convolution semigroup. Spanning implies \(a|\xi|\le\psi(\xi)\le b|\xi|\) for some positive finite constants, by compactness of the unit sphere. Hence the generator's \(L^2\) domain is \(H^1\), its form domain is \(H^{1/2}\), with the \(L^2\) component included in the norms.

The correct test generator is

\[
L\phi(x)=\frac1{2\pi}\sum_jc_j\int_\mathbb R
\frac{\phi(x+sv_j)+\phi(x-sv_j)-2\phi(x)}{s^2}\,ds.
\]

Indeed the integral of \((2\cos(as)-2)/s^2\) is \(-2\pi|a|\); the Fourier symbol is exactly \(-\psi\). The coefficient is linear in \(c_j\), and the vector length remains in the shifts. Near zero, second derivatives control the integrand; at infinity the denominator is integrable. The same bounds in \(L^1\) and the sup norm give \(L\phi\in L^1\cap L^\infty\cap C_0\). The scaling identities used later are therefore valid:

\[
\|L\chi_R\|_\infty=R^{-1}\|L\chi\|_\infty,
\qquad\|L\chi_\varepsilon\|_1=\varepsilon^7\|L\chi\|_1.
\]

Fourier dominated convergence on Schwartz tests gives \((P_t\phi-\phi)/t\to L\phi\) uniformly, also in \(L^2\), and \(\partial_tP_t\phi=P_tL\phi\). No assertion of high-dimensional \(H^1\) continuity is needed or used.

## 2. Exact eight-dimensional singular set

Define

\[
S_3=\bigcup_{\operatorname{rank}(A)\le3}\operatorname{span}(A),
\qquad G(x)=\int_0^\infty p_t(x)\,dt.
\]

This is a finite union of closed proper linear subspaces, so it is closed and Lebesgue null. It includes the origin and all low-rank spans, not just individual direction lines.

**Outside \(S_3\).** For compact \(E\Subset\mathbb R^8\setminus S_3\), let \(d=\operatorname{dist}(E,S_3)>0\), \(D_V=\sum_j|v_j|\), and \(\delta=d/(4D_V)\). If a representation \(x=\sum_js_jv_j\) has its large coefficients \(J=\{j:|s_j|>\delta\}\) of direction rank at most 3, its large part belongs to \(S_3\) and its small part has norm at most \(d/4\). Such a point cannot lie in the \(d/4\)-neighborhood of \(E\). Hence representations landing in that neighborhood have at least four **independent** large directions. Counting four large coefficients without this rank argument would not suffice.

For each independent four-tuple \(I\), extend it to an eight-direction basis \(B_I\). Restrict those four coefficient laws to \(|s|>\delta\) without normalizing their lost mass. Their densities are bounded by \(c_jt/(\pi\delta^2)\); the remaining four basis densities are bounded by \(1/(\pi c_jt)\). Basis change and convolution with the remaining probability law give a density bound independent of time because \(t^4t^{-4}=1\). Sum over the finite possible \(I\). Event domination first gives a measure inequality for every Borel subset of the open neighborhood, hence a density bound there a.e. Continuity of the true \(p_t\) promotes this to every point, separately for each \(t>0\).

Combining this uniform bound for \(0<t\le1\) with \(Ct^{-8}\) for \(t\ge1\) gives a locally uniform integrable time majorant. Consequently \(G\) is finite and locally continuous throughout the complement.

**On \(S_3\).** Write \(x=\sum_{k=1}^ra_kv_{i_k}\), where the directions are independent and \(0\le r\le3\), and extend them to a basis. The remaining random sum \(Y_t\) scales as \(tY_1\); the event \(\|B^{-1}Y_t\|_\infty\le Mt\) has some positive time-independent probability. On that event the first \(r\) basis densities have lower bounds proportional to \(t\), while the \(8-r\) zero-coordinate factors have lower bounds proportional to \(t^{-1}\). Thus

\[
p_t(x)\ge k_x t^{2r-8},\quad0<t\le1.
\]

For \(r=0,1,2,3\), the powers are \(-8,-6,-4,-2\), all nonintegrable at zero. The conclusion is therefore the exact pointwise statement

\[
G(x)=+\infty\quad\Longleftrightarrow\quad x\in S_3.
\]

The argument permits redundant or parallel directions. No enumeration of independent four-tuples is needed to prove the analytic theorem.

## 3. Whole-space integrability and the positive probability source

Split \(G=U_T+W_T\), where \(U_T=\int_0^Tp_tdt\) and \(W_T=\int_T^\infty p_tdt\). Then \(\|U_T\|_1=T\) and \(\|W_T\|_\infty\le C/(7T^7)\). Choosing \(T=(2C/(7\lambda))^{1/7}\) gives

\[
|\{G>\lambda\}|\le 2(2C/7)^{1/7}\lambda^{-8/7}.
\]

For finite-measure \(E\), integrating \(\int_Ep_t\le\min(1,C|E|t^{-8})\) in time proves local integrability. Layer integration above a positive threshold then gives \(H_c=(G-c)_+\in L^1(\mathbb R^8)\). This controls both singular neighborhoods and the far field; homogeneity alone would not justify it.

For a compact smooth test \(\phi\),
\(\int_0^\infty\int p_t|L\phi|\le\int_0^\infty\min(\|L\phi\|_\infty,Ct^{-8}\|L\phi\|_1)dt<\infty\).
Thus semigroup differentiation and the two time endpoints legitimately give \(LG=-\delta_0\). Tonelli and the semigroup law give \(P_tG=\int_t^\infty p_sds\), so \(G=P_tG+U_t\) at its a.e. finite points. Applying Jensen to the positive-part function yields

\[
P_tH_c\ge(P_tG-c)_+\ge H_c-U_t.
\]

Pair against nonnegative compact smooth tests and let \(t\downarrow0\). Uniform convergence of the test generator and the \(L^1\) property of \(H_c\) treat the left side; \(U_t/t\) is an approximate identity on the right. Hence

\[
\mu_c:=LH_c+\delta_0\ge0
\]

as a distribution on **all of \(\mathbb R^8\)**. The positive-distribution theorem makes it a positive Radon measure, including all low-rank singular spaces and level interfaces.

For monotone large cutoffs \(\chi_R\uparrow1\),
\(\int\chi_Rd\mu_c=1+\int H_cL\chi_R\to1\), because \(\|L\chi_R\|_\infty=O(R^{-1})\). Positivity and monotone convergence prove total mass exactly 1, without an a priori finite-mass assumption.

Let \(h=\min(G,c)\). It is bounded, and the previous integrability statements justify \(\mu_c=-Lh\). The author proves no atom at the origin. Its same argument proves the stronger **absence of every atom**: for any \(x_0\), test with \(\chi((x-x_0)/\varepsilon)\), equal to 1 near \(x_0\). Translation invariance and scaling give

\[
0\le\mu_c(\{x_0\})\le
\int\chi((x-x_0)/\varepsilon)d\mu_c
\le c\varepsilon^7\|L\chi\|_1\longrightarrow0.
\]

This is a strict analytic consequence, not a sampled assertion. The result is an atomless positive probability source satisfying \(LH_c=\mu_c-\delta_0\) globally. The comparison theorem below itself requires only a positive finite probability source; it does not rely on atomlessness.

## 4. Obstacle decomposition, compactness, and no loss of mass

For nonnegative bounded compactly supported \(f\), put \(F=\|f\|_\infty\), \(m=\|f\|_1\), and fix \(\kappa>0\). The trivial case \(m=0\) can be omitted. For \(0<t\le1\),

\[
V_t=t\sum_{j\ge0}P_{jt}f,\quad
u_{t,0}=0,\quad u_{t,k+1}=(P_tu_{t,k}+t(f-\kappa))_+.
\]

The positive-time series has a uniform sup-norm tail bounded by \(t\sum_{j\ge1}\min(F,Cm(jt)^{-8})\), at most the corresponding integral. Consequently \(V_t\) is uniformly bounded, and \(V_t-tf\in C_0\). Tonelli proves \(P_tV_t=V_t-tf\); monotone iteration therefore has a fixed point with

\[
0\le u_t\le(V_t-\kappa t)_+.
\]

For each fixed \(t\), this upper bound has essential compact support because \(f\) is compactly supported and the other part of \(V_t\) is in \(C_0\). It follows that \(u_t\in L^1\cap L^2\). Neither \(V_t\in L^1\) nor a common compact support for different \(t\) is required.

Set \(L_t=(P_t-I)/t\) and \(\sigma_t=f+L_tu_t\). On \(\{u_t>0\}\), the fixed-point equality gives \(\sigma_t=\kappa\); on \(\{u_t=0\}\), it gives \(0\le\sigma_t\le\kappa\). Because \(u_t\) is now known to be in \(L^1\), semigroup mass preservation gives \(\int\sigma_t=m\). Pairing with \(u_t\), self-adjoint contractivity gives

\[
\langle u_t,-L_tu_t\rangle+\kappa\|u_t\|_1=\int fu_t.
\]

The resulting uniform \(L^1,L^2,L^\infty\) bounds on \(u_t\), \(L^2\) bound on \(\sigma_t\), and \(L^2\) bound on \(L_tu_t=\sigma_t-f\) are valid. With the unitary Fourier transform, the multiplier is \(-\psi_t=-(1-e^{-t\psi})/t\). For \(tQ\le1\) and \(\psi\ge Q\),

\[
\psi_t\ge(1-e^{-1})Q.
\]

Thus high-frequency \(L^2\) tails tend uniformly to zero along \(t\downarrow0\), at rate \(O(Q^{-2})\) in squared norm. A low-frequency projection lies in a bounded \(H^1\) set, because \(\psi\asymp|\xi|\); local Rellich and a diagonal subsequence give local strong \(L^2\) and a.e. convergence of \(u_t\). Global weak \(L^2\) limits for \(u_t\) and \(\sigma_t\) can be taken on the same subsequence.

Fatou gives \(u\ge0\), \(u\in L^1\cap L^\infty\). Test convergence gives \(Lu=\sigma-f\in L^2\), and the multiplier ellipticity implies \(u\in H^1\). Equivalently, one may first extend the test identity by the standard density of \(C_c^\infty\) in the \(H^1\) generator graph norm, then use the closed self-adjoint operator. No Sobolev embedding into continuous functions is invoked.

The weak limit has \(0\le\sigma\le\kappa\). At a.e. point in \(D=\{u>0\}\), the approximating \(u_t\) eventually becomes positive, so \(\sigma_t=\kappa\) eventually there. Dominated convergence on finite-measure portions of \(D\), combined with weak convergence, identifies \(\sigma=\kappa\) on \(D\).

Positive cutoff tests first give \(\sigma\in L^1\) and \(\int\sigma\le m\). To exclude escaped mass, test \(f=\sigma-Lu\) with \(\chi_R\). Since

\[
|\langle Lu,\chi_R\rangle|=|\langle u,L\chi_R\rangle|
\le R^{-1}\|u\|_1\|L\chi\|_\infty\to0,
\]

one gets \(\int\sigma=m\), hence \(\kappa|D|\le m\). The final decomposition is therefore valid with a **single receiving set** \(D\), independent of the comparison radius:

\[
f=\sigma-Lu,\quad u\ge0,\quad
u\in L^1\cap L^\infty\cap H^1,\quad
0\le\sigma\le\kappa,\quad\sigma=\kappa\text{ on }D,\quad
\int\sigma=\int f.
\]

## 5. Finite-source settlement and every positive scale

Suppose \(K\ge0\), \(K\in L^1\), \(M_K=\int K>0\), and
\(LK=a(\mu-\delta_0)\), where \(a>0\) and \(\mu\) is a positive probability measure. Compact test identities extend to Schwartz tests by smooth cutoffs converging in the Schwartz topology. The second-difference estimate gives uniform convergence of their \(L\)-images, which justifies passage against \(K\in L^1\) and finite \(\mu\). Fourier inversion then gives the continuous function identity

\[
-\psi\widehat K=a(\widehat\mu-1).
\]

Multiplication by the unitary Fourier transform of \(u\in H^1\), and Young bounds for \(L^1*L^2\) and finite-measure convolution, give the legitimate \(L^2\) identity

\[
K*Lu=a(\mu*u-u).
\]

This avoids defining an unsupported convolution of arbitrary distributions. For \(K_r(x)=r^{-8}K(x/r)\) and the probability dilation \(\mu_r\), first-order homogeneity gives

\[
LK_r=\frac ar(\mu_r-\delta_0),\quad
K_r*Lu=\frac ar(\mu_r*u-u).
\]

The factor \(1/r\) is essential and is retained. For the obstacle decomposition of bounded compactly supported \(f\), choose \(\kappa=\lambda/M_K\). For each fixed \(r\), \(\mu_r*u-u\ge0\) a.e. on the **same** \(D^c\); Fubini with \(dx\,d\mu_r\) guarantees that null-set choices of the \(L^2\) representative do not matter even when \(\mu_r\) is singular. Thus

\[
K_r*f\le K_r*\sigma\le\kappa M_K=\lambda
\quad\text{a.e. on }D^c.
\]

Combine exceptional null sets first for positive rational \(r\). The map \(r\mapsto K_r\) is \(L^1\)-continuous; boundedness of \(f\) then gives uniform continuity in the radius of \(K_r*f\). This upgrades the countable conclusion to all positive real radii without reselecting \(D\). It proves

\[
\lambda|\{\sup_{r>0}K_r*f>\lambda\}|\le M_K\|f\|_1.
\]

For general \(f\in L^1\), use \(f_j=\min(|f|,j)1_{[-j,j]^8}\uparrow|f|\). Positive convolution gives monotone convergence for every fixed \(x,r\) with the canonical extended integral. The two suprema commute:

\[
\sup_r\sup_jK_r*f_j=\sup_j\sup_rK_r*f_j.
\]

For each bounded approximant the all-real-radius supremum equals its rational-radius supremum; hence the limit is measurable and equals the all-real-radius supremum for \(|f|\). Strict level sets and input masses converge monotonically. The same weak bound therefore holds for arbitrary \(L^1\), without presuming radius-uniform continuity for the unbounded final input.

If \(0<g_0\le\operatorname{essinf}_{[-1,1]^8}G\) and \(0<c<g_0\), set

\[
K=\frac{(G-c)_+}{256(g_0-c)}.
\]

Then \(K\ge256^{-1}1_{[-1,1]^8}\) a.e., and its source has coefficient \(a=[256(g_0-c)]^{-1}\). This a.e. kernel domination suffices for every translated Lebesgue convolution and every dilation: integrating a null-set change of a kernel against a Lebesgue density changes no convolution value. The full cube maximal operator is therefore bounded by the kernel maximal operator, giving \(C_8\le M_K\). The argument uses one common scale in all eight coordinates; independent rectangle scales are not included.

## 6. Essential face interface, homogeneity, and cost

Scaling gives the actual pointwise identity \(G(rx)=r^{-7}G(x)\), allowing infinite values. Let \(Q=[-1,1]^8\) and

\[
J=\sum_{F\subset\partial Q}\int_FG(y)^{8/7}\,dS(y).
\]

On each open face the ray map \((r,y)\mapsto ry\) has volume Jacobian \(r^7\). Overlaps of face cones have zero volume. Tonelli therefore gives

\[
|\{G>1\}|=J/8,\quad
|\{G>s\}|=(J/8)s^{-8/7},\quad
\int(G-c)_+=\frac78Jc^{-1/7}.
\]

The previously proved weak bound implies \(J<\infty\); trace integrability is not an assumption. Strict positivity of the density on compact \(Q\times[1,2]\) gives a positive lower bound for \(G\) on \(Q\); local integrability makes its essential infimum finite.

The essential cube and essential boundary infima are equal. One direction follows from \(r^{-7}\ge1\) for \(r\le1\). For the other, a positive-surface-measure set of face values \(G(y)\le a<g_{\rm cube}\), followed along an interval of rays with \(r<1\) sufficiently close to 1, would produce a positive-volume subset of \(Q\) with \(G<g_{\rm cube}\). This contradiction is the needed trace argument; lower semicontinuity alone would not justify it. In the finite-frame model, continuity off \(S_3\) and infinite values on \(S_3\) also upgrade an essential cube lower bound to a pointwise cube lower bound.

The kernel mass is

\[
M_K=\frac{7J}{8\cdot256}\frac{c^{-1/7}}{g_0-c}.
\]

The logarithmic derivative is \((8c-g_0)/(7c(g_0-c))\), so the unique optimum is \(c=g_0/8\). At that value,

\[
C_8\le\frac{8^{1/7}J}{256g_0^{8/7}}.
\]

If the actual model has individual coordinate-reflection symmetry and is transitive on all sixteen faces, then for
\(I=\int_{[0,1]^7}G(1,y)^{8/7}dy\) one has \(J=16\cdot128 I=2048I\). A genuine continuous certificate \(I\le I_0\) would consequently yield

\[
C_8\le8^{8/7}I_0/g_0^{8/7},\qquad
C_8^7\le8^8I_0^7/g_0^8.
\]

These are conditional interfaces. This package supplies neither numerical inequality \(G\ge g_0\) nor \(I\le I_0\) for a proposed numerical pair. Arbitrary frames need not have the symmetry needed to replace \(J\) by \(2048I\).

## 7. Short time, fixed-slice Fatou, and long time

On the face \(x_1=1\), put \(J_1=\{j:(v_j)_1\ne0\}\), \(D_1=\sum_{j\in J_1}|(v_j)_1|\), and \(A=\sum_jc_j|(v_j)_1|\). Spanning makes both positive. For \(0<\beta<1\), a representation landing in \(|x_1|>\beta\) must have some normal coefficient larger than \(\beta/D_1\) in absolute value. Complete each such direction to a basis. Its restricted density costs \(O(t)\), while the remaining seven basis factors cost \(O(t^{-7})\). An event union gives an a.e. density bound on the open region. Continuity gives the pointwise bound, and only then may one let \(\beta\uparrow1\) at the fixed face. The result is

\[
\sup_{y\in\mathbb R^7}p_t(1,y)\le Bt^{-6},\quad
B=D_1^2\sum_{j\in J_1}
\frac{c_j}{\pi^8|\det B_j|\prod_{i\in B_j\setminus\{j\}}c_i}.
\]

This avoids conditioning on a zero-probability slice.

The one-dimensional marginal is Cauchy of scale \(At\), but its density identification directly gives only an a.e. identity in \(s\). Fix \(t\), choose good \(s_k\to1\), use pointwise continuity in \((s,y)\), and apply Fatou over \(y\). The valid direction is

\[
\int_{\mathbb R^7}p_t(1,y)dy
\le\frac{At}{\pi(1+A^2t^2)}
\le\min\left(\frac A\pi t,\frac1{\pi At}\right).
\]

No unproved slice equality is needed. Interpolating the short-time sup and mass bounds gives

\[
\|p_t(1,\cdot)\|_{8/7}
\le B^{1/8}(A/\pi)^{7/8}t^{1/8},
\]

because \(-6/8+7/8=1/8\). Minkowski first on finite time/space regions, then monotone convergence, gives the **whole-plane** short-time budget

\[
\left\|\int_0^\varepsilon p_t(1,\cdot)dt\right\|_{8/7,\mathbb R^7}
\le\frac89B^{1/8}(A/\pi)^{7/8}\varepsilon^{9/8}.
\]

The bound includes the low-rank singular neighborhoods and spatial infinity. The set \(S_3\cap\{x_1=1\}\) has zero seven-dimensional measure, so its infinite point values do not contradict this finite norm. Individual free-coordinate reflection symmetry improves the orthant norm by \(2^{-49/8}\); restricting further to the unit positive face preserves that upper bound.

The global long-time pointwise bound is \(R_T\le C/(7T^7)\). It yields an \(L^{8/7}(E)\) bound \(|E|^{7/8}C/(7T^7)\) for finite face domains \(E\). It cannot be directly integrated over all of \(\mathbb R^7\). If a whole-plane long-time estimate is desired, combine instead \(Ct^{-8}\) with the fixed-slice mass \((\pi At)^{-1}\), giving

\[
\|p_t(1,\cdot)\|_{8/7}\le C^{1/8}(\pi A)^{-7/8}t^{-15/8},
\quad
\|R_T(1,\cdot)\|_{8/7,\mathbb R^7}
\le\frac87C^{1/8}(\pi A)^{-7/8}T^{-7/8}.
\]

The distinct \(T^{-7}\) and \(T^{-7/8}\) rates are used on their appropriate domains.

## 8. Sixteen unnormalized directions and the 4+4 model

The model consists of the eight vectors \(W=\{(1,\epsilon_2,\epsilon_3,\epsilon_4)\}\) in each four-dimensional block, with \(c_j=1\). Each direction has length 2. The two blocks share the same semigroup time:

\[
p_t^{(8)}(x',x'')=p_t^{(4)}(x')p_t^{(4)}(x''),\quad
G_8(x',x'')=\int_0^\infty p_t^{(4)}(x')p_t^{(4)}(x'')dt.
\]

The latter is a single integral, not \(G_4(x')G_4(x'')\). Four-dimensional weak constants or previously computed four-dimensional \(g/I\) values therefore cannot simply be squared or reused as eight-dimensional certificates.

A four-dimensional Hadamard basis has determinant magnitude 16; the block basis has determinant magnitude 256. Thus

\[
C=1/(256\pi^8),\quad R_T\le1/(1792\pi^8T^7).
\]

The exact frame Gram matrix is \(8I_8\); therefore \(\sqrt8|\xi|\le\psi(\xi)\le8\sqrt2|\xi|\). These are coarse structural bounds, not numerical \(g/I\) estimates.

Each of the eight normal directions on the first face admits a determinant-256 completion. An explicit construction is \(\operatorname{diag}(w)H_4\): its first column is \(w\), the other columns remain in \(W\), and determinant magnitude stays 16. Add an independent Hadamard basis in the other block. Hence \(A=D_1=8\) and

\[
B=8^2\frac8{256\pi^8}=2/\pi^8,
\quad
\|T_\varepsilon(1,\cdot)\|_{8/7,\mathbb R^7}
\le\frac89\,2^{1/8}8^{7/8}\pi^{-15/8}\varepsilon^{9/8}.
\]

Coordinate flips preserve the unoriented body directions and their Cauchy laws, even if a representative vector changes sign. Within-block permutations, coordinate flips, and block exchange make all faces equivalent and justify \(J=2048I\). The author uses the correct unnormalized convention consistently.

## 9. What the actual finite check establishes

`independent_check.py` was written for this audit and actually exited 0; `actualreceipt.json` records its source hash, interpreter, UTC timing, runtime and status. It verifies the Hadamard/block determinants, Gram matrix, vector lengths, eight explicit normal completions, \(A,D_1,B,C\), source dilation exponent \(-1\), short/long time exponents and prefactors, orthant exponent, face factor 2048, and the raised cost coefficient \(8^8=16777216\). All computations use integers or Fractions; there are no Green nodes, optimizers, subset enumerations, or approximate integration results.

The analytic proof above establishes the singular classification, source positivity and probability mass, atomlessness, obstacle compactness, and full-operator comparison. A finite algebra check cannot prove those infinite-dimensional statements. The source author's reported finite checks are reviewed for consistency and provenance separately in the manifest; their prior labels are not substituted for this analysis.

The finalized author package consists of `proof.md`, `check.py`, `checks.json`, and `sources.json`; the last file is its source manifest. All four are locked by hashes before and after the final audit freeze. The external source hashes in `sources.json` were independently matched to the available files. The saved `checks.json` describes its arithmetic and platform but does not contain command timing or a separate execution receipt. Consequently those saved bytes alone are not an independent historical attestation of the author's run. This audit instead supplies its own actual run receipt, and neither the analysis nor its finite checks require rerunning the author's 1796 independent-four-subset enumeration.

Fixed Lebesgue background is essential throughout: it supplies self-adjoint translation convolution, Fourier diagonalization, mass preservation, null-set invariance under translations/dilations, the measure-convolution Fubini argument, \(K*\sigma\le\kappa\int K\), \(\kappa|D|\le\int\sigma\), and cube-ray Jacobians. This audit makes no extension to variable background measures.

No unresolved substantive gap was found in the stated analytic interface. The remaining numerical work is to construct compatible continuous \(g_0\) and \(I_0\) certificates for the true eight-dimensional model. The Green material currently provides no new accepted numerical value of \(C_8\).

## Reproduction

From the clone root:

```text
python3 cube-dim-8-bounds/audit/round001/green/independent_check.py
```

This reruns only the small exact checks and overwrites their result/receipt in this audit directory. The analytic assessment is contained in this report. Base repository HEAD and all frozen author source hashes are recorded in `manifest.json`; no assertion is made about unrelated missing historical blobs or full-repository integrity.
