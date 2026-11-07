# A terminal-independent common action on the critical fixed-block strip

## Result, and what it does not say

This gives a common-action theorem for a fixed number of equal radial
blocks. It avoids differentiating a moderate-deviation formula in the
parameter: every efficient terminal is compared to one fixed scalar ramp.
It does **not** prove the general strong-maximal upper bound. The fixed
equal-block far-band localization is supplied separately by
global_blocks65.md; action-dependent angular physical kernels remain
outside the present scope.

Let \(N=rb\), with fixed \(r\), \(L=\log N\), and
\[
 Z_c=E+cB_c,\qquad E\sim{\rm Exp}(1),\quad
 B_c\sim{\rm Bernoulli}((1+e^c)^{-1}).
\]
The \(r\) block sums have total \(S\). Put
\[
 q_0=N\mu-d_N,\qquad
 d_N=\sqrt{\mu N}\,z_N,\qquad \overline\Phi(z_N)=d_N/N. \tag{1}
\]
For the local radial theorem, \(H\) is defined on \(\mathbb R_+^r\). Write
\[
 J_+(H)=\int_{\mathbb R_+^r}e^{-\sum_jt_j}\Phi_N(H(t))\,dt.
\]
We consider capped \(0\le H\le A\) with
\[
 J_+(H)\le e^{-q_0+C_*L},\qquad
 E_{c_H}H\ge1\quad\hbox{for some }c_H\in
 I_K=[\mu-K\sqrt{L/N},\mu+K\sqrt{L/N}].                 \tag{2}
\]
The fixed \(C_*\) includes the \(O(L)\) source-budget slack of all critical
translations below. Set \(A_*=\max\{A,2\}\); every imported cap-dependent
estimate is used with this one fixed cap. Constants are uniform in \(H\)
and its translate.

## 1. A fixed scalar reference

Define
\[
 h_0(s)=
 \begin{cases}
 0,&s<q_0,\\
 \min\{2,1+(s-q_0)/N\},&s\ge q_0,
 \end{cases}
 \qquad H_0(t)=h_0\!\left(\sum_jt_j\right).              \tag{3}
\]
The jump at \(q_0\) is intentional. Coarea on the simplex gives
\[
 J_+(H_0)=\frac1{(r-1)!}\int_{q_0}^{\infty}
 e^{-s}s^{r-1}\Phi_N(h_0(s))\,ds
 =e^{-q_0}N^{r-1/2}\Theta_{r,A_*}(1).                   \tag{4}
\]
Indeed, with \(s=q_0+x\sqrt N\),
\[
 \log\Phi_N(1+x/\sqrt N)-x\sqrt N
 =-x^2+O(N^{-1/2}(1+|x|^3)),
\]
while the other ranges are exponentially integrable. The simplex factor
is \(q_0^{r-1}(1+o(1))\) and the active width is \(\sqrt N\).
Thus \(H_0\) satisfies (2) after increasing \(C_*\).

The uniform Cramér CDF and lower-tail overshoot estimates in
r063/refined63.md yield
\[
 g_0(c):=E_cH_0=1+O(L^2/N)\qquad(c\in I_K).              \tag{5}
\]
In fact,
\[
 E_cH_0=1-F_{N,c}(q_0)+N^{-1}E_c(S-q_0)_++O(N^{-M}),
\]
and the \(d_N/N\) terms cancel; the mean shift is \(O(L)\), the CDF error
is \(O(L^2/N)\), and the lower overshoot is \(O(1)\).

Let \(c_N\) maximize \(g_0\) on \(I_K\). The refined common-curvature
theorem applies to the near-success reference and gives
\[
 g_0''(c)=-p_\mu L+o(L)\qquad(c\in I_K).                 \tag{6}
\]
Hence \(c_N\) is unique for large \(N\), and depends only on the scalar
reference law—not on the received terminal.

## 2. Terminal-to-reference comparison

For \(H\), set \(P=(H-1)_+\), \(D=(1-H)_+\). The two-buffer envelope in
r063/refined63.md, together with the full-window near-level bootstrap of
r062/sharp_response62.md, gives uniformly on \(I_K\)
\[
 \|P-U\|_1+\|D-I_q\|_1=O(L^2/N),                          \tag{7}
\]
where \(U=U(S)\), \(I_q=\mathbf1_{\{S\le q\}}\), and
\(|q-q_0|=O(L)\). The bootstrap starts from the successful action in
(2), and yields \(E_cH\ge1-O(L^2/N)\) throughout the window.

For the reference,
\[
 P_0=(h_0(S)-1)_+,\qquad D_0=\mathbf1_{\{S<q_0\}},
\]
up to the negligible cap tail. The inverse-ramp formula for \(U\), the
central estimate \(U'=N^{-1}(1+O(\sqrt{L/N}))\), and
\(|q-q_0|=O(L)\) imply
\[
 \|U-P_0\|_1+\|I_q-D_0\|_1=O(L^2/N).                      \tag{8}
\]
For the indicators this is a strip of width \(O(L)\) times the critical
density \(O(L/N)\); moving the positive-ramp turn-on costs only \(O(L/N)\).
Thus
\[
 \boxed{\ \|H-H_0\|_{L^1(P_c)}=O(L^2/N),\quad c\in I_K.\ } \tag{9}
\]
This is the efficiency-qualified step: it uses actual critical budget and
a near-unit received mean, rather than a bare angular-envelope norm.

## 3. One action for all critical terminals

Let \(s_c(T)=\partial_c\log f_{T,c}(T)\). The Euler/Fourier block estimate
of r062/sharp_response62.md gives \(\sup_{\mathcal C_B}|s_c|\le CL\).
Off the blockwise central box, the conditional Gamma--Binomial identity
gives \(\|s_c\|_2=O(\sqrt N)\), while the box probability can be made an
arbitrarily high inverse power. Differentiation of the finite
Gamma--Binomial mixture is valid for bounded measurable terminals.
Consequently (9) gives
\[
 \sup_{c\in I_K}|g_H'(c)-g_0'(c)|\le C L^3/N.             \tag{10}
\]
No KKT condition or terminal smoothness occurs here.

The same common-curvature theorem gives, for all large \(N\),
\[
 g_H''(c)\le-\tfrac12p_\mu L\qquad(c\in I_K).             \tag{11}
\]
Comparing (10) to the variational inequality for the constrained maximizer
\(c_N\) of \(g_0\), strong concavity gives the endpoint-safe estimate
\[
 \sup_{c\in I_K}g_H(c)-g_H(c_N)
 \le C\frac{(L^3/N)^2}{L}=O(L^5/N^2).                    \tag{12}
\]
Thus \(c_N\) is universal at the precision relevant here. This robust
statement replaces an unnecessary terminal-by-terminal asymptotic formula
for an optimizer.

### Blockwise action vectors

The conclusion is not restricted to receivers already constrained to use a
common action. Let \(\mathbf c=(c_1,\ldots,c_r)\in I_K^r\) act
independently on the equal blocks, and write
\(g_H(\mathbf c)=E_{\mathbf c}H(T)\). The vector near-level bootstrap and
the pointwise vector-Hessian theorem from r063 give, on the full box,
\[
 D^2g_H(\mathbf c)\preceq-\kappa_rL\,I_r.                \tag{13}
\]
For the reference \(H_0(S)\), the law and terminal are invariant under
block permutations. Its strictly concave vector maximizer is therefore
unique and permutation-invariant:
\[
 \mathbf c_N=(c_N,\ldots,c_N).                            \tag{14}
\]
Its scalar coordinate is the maximizer in (6), because restriction to the
diagonal is the common-parameter reference problem.

For clarity, the near-level input for this vector assertion is uniform:
\[
 \sum_{i=1}^rb\,m(c_i)=N\mu+O(L),\qquad
 \sum_{i=1}^rb\,\nu(c_i)=N\mu+O(N\rho).                   \tag{15}
\]
The fixed-\(r\) triangular-array tilted Fourier/CDF estimate therefore
gives \(E_{\mathbf c}H_0=1+O(L^2/N)\) on \(I_K^r\). Thus the vector
Hessian theorem applies to the reference without assuming an exact
success point.

The same \(L^1\) comparison (9), paired with each individual block score,
gives
\[
 \sup_{\mathbf c\in I_K^r}
 \|\nabla g_H(\mathbf c)-\nabla g_0(\mathbf c)\|_2
 =O_r(L^3/N).                                             \tag{16}
\]
The tail estimate is unchanged because \(r\) is fixed and every central
block score is \(O(L)\). Vector strong concavity and the variational
inequality at \(\mathbf c_N\) consequently yield
\[
 \sup_{\mathbf c\in I_K^r}g_H(\mathbf c)-g_H(\mathbf c_N)
 =O_r(L^5/N^2).                                           \tag{17}
\]
Thus even receiver-dependent block-action vectors in the fixed-block
critical strip can be replaced by the one source-independent common action
\(c_N\).

## 4. Uniform spatial repair

Here is the uniform spatial lemma needed for the repair; it is stronger
than an assertion about the derivative only at \(a=0\).  Fix
\(\eta_0=N^{-1/2}\).  If \(K\) has cap at most \(A_*\), critical budget
with the fixed \(O(L)\) slack, and
\[
 E_{\mathbf c}K(T)\ge1-\eta_0,\qquad\mathbf c\in I_K^r,
\]
then the near-level refined location theorem, evaluated at its origin,
gives
\[
 \frac d{da}E_{\mathbf c}K(T+a\mathbf1)\big|_{a=0}
 =\sum_{i=1}^rG_{d_i}(K,\mathbf c)
 \ge c_r\frac LN.                                        \tag{18}
\]
The statement applies to every translated terminal for which these
hypotheses hold, not just to the original \(K\).

After the vector replacement (17), the loss at \(\mathbf c_N\) is
\(\varepsilon_N=O_r(L^5/N^2)=o(\eta_0)\).  Put
\(K_a(t)=H(t+a\mathbf1)\).  Its orthant budget obeys
\[
 J_+(K_a)\le e^{ra}J_+(H),                                \tag{19}
\]
and, for a preliminary radial offset \(u\),
\[
 J_+\bigl(H(u+a\mathbf1+\cdot)\bigr)
 \le e^{\sum_i u_i+ra}J_{\rm global}(H).                  \tag{20}
\]
Thus every \(K_a\), \(0\le a\le C L^4/N\), remains in the same fixed
critical budget class.

For each \(a\) for which \(K_a\) is near-successful, apply the preceding
origin statement to \(K_a\); this is precisely the derivative of
\(E_{\mathbf c_N}K_a(T)\) with respect to \(a\). Let \(\tau\) be the
first point in this interval where
\(E_{\mathbf c_N}K_\tau=1-\eta_0\), if such a point exists.  Before
\(\tau\), the near-level hypotheses hold, so the repeated estimate (18)
makes the expectation
strictly increasing.  This contradicts a downward first exit from
\(1-\varepsilon_N>1-\eta_0\).  Hence (18) holds throughout the interval.
Choosing a sufficiently large fixed \(C\) and integrating (18) with
\[
 a_N=C\,L^4/N=o(1)                                        \tag{21}
\]
gains at least \(C c_rL^5/N^2\), and therefore repairs every loss in
(12) or (17) at the one action \(c_N\).

On the full global radial source space (the physical realisation used for
the final transfer), the change of variables has no truncated orthant and
the modular identity is exact:
\[
 J_{\rm global}(H(\,\cdot+a_N\mathbf1))
 =e^{ra_N}J_{\rm global}(H)=(1+o(1))J_{\rm global}(H).     \tag{22}
\]
Formula (19), rather than equality, is the corresponding statement for a
terminal represented only on \(\mathbb R_+^r\).

All estimates are uniform for a translated/recentred family with total
radial offset \(|\sum_j u_j|\le C_1L\), after the one enlargement of
\(C_*\). Hence this aligns a whole variable common-action critical strip.

## 5. Fixed-action block-model transfer and critical-strip consequence

The required fixed-action transfer is available for the block model itself.
For one block-action vector \(\mathbf c\), put
\(\mathcal R_{\mathbf c}H(u)=E_{\mathbf c}H(u+T)\). Every bounded global
block profile of finite budget obeys
\[
 \nu\{u:\mathcal R_{\mathbf c}H(u)>1\}\le Q_NJ_{\rm global}(H).
                                                               \tag{23}
\]
For compactly supported \(H\) and a compact receiver set with a uniform
strict margin, use the physical deep profile
\[
 f_D(x)=H(v(x)-D\alpha)\mathbf1_{[-1,1]^N}(x).
\]
At \(v=D\alpha+u\), its fixed-anchor rectangle average is exactly
\(E_{\mathbf c}H(u+T)\). The coarea limits in r026/report.md give
\(\int\Phi_N(f_D)/K_D\to J_{\rm global}(H)\) and
\(|E_D|/K_D\to\nu(E)\), while r028/report.md proves that the corresponding
fixed-global-vector full-input weak constant is \(Q_N\). This proves (23)
in the compact strict case.

For general bounded \(H\), set \(H_M=H\mathbf1_{[-M,M]^r}\). Their
convolutions increase to the continuous limiting convolution; Dini's
theorem gives uniform convergence on each compact strict receiver set.
Exhausting the strict halo and using monotone convergence proves (23).
Thus, within the block model, the preceding common-action replacement gives
the following explicit critical-strip statement.  Let \(E_{\rm crit}\) be
the receivers \(u\) for which \(|\sum_j u_j|\le C_1L\), the local terminal
\(H_u(t)=H(u+t)\) has \(J_+(H_u)\le e^{-q_0+C_*L}\), and some
\(\mathbf c(u)\in I_K^r\) satisfies \(E_{\mathbf c(u)}H_u(T)>1\).  Then,
after the one global repair,
\[
 \nu(E_{\rm crit})\le(1+o(1))Q_NJ_{\rm global}(H).
                                                               \tag{24}
\]
The global source is used there only through the orthant inequality (20).
This is not a realization of a general physical parameter map. The
action-dependent angular physical map remains the general-upper obstacle.
