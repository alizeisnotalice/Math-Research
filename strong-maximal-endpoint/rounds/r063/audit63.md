# Round 63 final audit

## Classifications

| File | Classification | What it establishes |
|---|---|---|
| nonmono63.md | **Rigorous scoped fixed-terminal theorem** | The capped balanced-block Hessian/continuation result for arbitrary bounded measurable, action-independent \(H\). |
| variance63.md | **Rigorous leading-critical zero-variance countertest; not Q-sharp** | Its cap-two ramp has budget \(N\mu-d_N+O(L)\), identical fibre responses, and one common constrained argmax. |
| angular63.md | **Exact kernel identities and a low-efficiency universal-obstruction result** | An unrestricted raw angular terminal has no universal negative curvature. This is not a critical-cost/near-efficient physical counterexample. |

None supplies P1, P2-excess/(7), or a strong-maximal upper bound.

## 4. refined63.md — pass at fixed block count

The corrected refined argument is rigorous for fixed \(r\), exact success,
and the stated capped action-independent terminal model.

The CDF normalization is correct. In the tilted Laplace identity, the
uniform local-CLT error contributes \(C/(\lambda\sqrt N)\) after
normalization, while the Gaussian integral is of order \(1/\lambda\);
hence the relative CDF error is \(O(N^{-1/2})\). Retaining the exact
Gaussian tail gives
\[
 F_{N,c}(N\mu-d_N+w)=\bar\Phi(z_N)
 \{1+O(L^{3/2}/\sqrt N)\},
\quad
 f_{N,c}(N\mu-d_N+w)=
 {\varphi(z_N)\over\sqrt{\mu N}}
 \{1+O(L^{3/2}/\sqrt N)\}.
\]
The density must not be replaced by \(L/N\) with this relative error:
the \(\log L/L\) correction is larger. The consequences used later,
\[
 F(s_N+O(L))=d_N/N+O(L^2/N),\qquad
 f(s_N+O(L))=f(s_N)+o(N^{-1}),
\]
do follow.

The two-buffer construction is essential and correct. It makes \(H\)
small on the lower core at \(q=s_N-(K_0+K_1)L\), while the scalar envelope
uses the separate \(K_0L\) cutoff. Although the CDF mass of this displaced
core differs from that at \(s_N\) by \(O(L^2/N)\), its pairing with the
central location score \(O(\rho)\) is \(o(N^{-1})\). Thus it does not
spoil the \(1/N\) response term.

The overshoot estimate is \(O(1)\), so
\[
 EU=d_N/N+O(L/N).
\]
The remaining denominator correction is \(E\delta^2/N^2=O(L/N)\). This
gives the stated \(L^1\) comparison errors \(O(L^2/N)\). Integration by
parts then yields \(E[U\ell_i]=-EU'(S)=-1/N+o(N^{-1})\), because the
envelope derivative is \(1/N+o(1/N)\) off a lower tail of probability
\(O(\rho)\) and an exponentially rare cap tail. Consequently
\[
 G_{d_i}=f_{N,c}(s_N)+1/N+o(N^{-1}).
\]

For two equal blocks, the plus-residual identity gives
\[
 G_{cc}=-p_\mu\left\{
 {z_N\varphi(z_N)\over\bar\Phi(z_N)}+1\right\}+o(1)
 =-p_\mu\{L-2\log L-\log(2\pi\mu)+2\}+o(1).
\]
The needed near-success tolerance is correctly stronger here:
\(\eta\rho=o(N^{-1})\). This remains a fixed-block auxiliary refinement,
not a physical angular bridge or a global upper-bound result.

## 5. Capped-ramp critical-scale packet

The announced ramp derivation is sound, provided the written section defines
the terminal explicitly as
\[
 H_q(t)=
 \begin{cases}
 0,&t<q,\\
 1+(t-q)/N,&q\le t\le q+N,\\
 2,&t>q+N.
 \end{cases}
\]
It is a fixed cap-two monotone terminal. Its receiving expectation is
continuous and strictly decreasing in \(q\), with derivative
\[
 {d\over dq}EH_q(S)=-f(q)-N^{-1}P(q<S<q+N).
\]
At the canonical lower-tail level, the root equation is
\[
 F(q)=N^{-1}E\min\{(S-q)_+,N\}.
\]
The refined CDF estimate and the \(O(1)\) lower-tail overshoot give a root
\(q=N\mu-d_N+O(L)\); the displayed derivative has size \(\asymp L/N\).

For this terminal,
\[
 J(H_q)=e^{-q}I_N,\quad
 I_N=\int_0^N e^{g_N(u)}\,du+e^{g_N(N)},\quad
 g_N(u)=\log\Phi_N(1+u/N)-u.
\]
The stated \(g_N(0)=g_N'(0)=0\) and
\(-2/N\le g_N''\le-a/N\) imply \(I_N\asymp\sqrt N\). Hence
\[
 -\log J(H_q)=N\mu-d_N+O(L).
\]
Lowering \(q\) by any arbitrarily small positive amount makes reception
strict while preserving that leading-critical budget scale. In the deep
scalar packet the fibre responses are identical, so the capped
common-window theorem gives one common argmax and zero variance.

This is therefore a genuine countertest to an **unqualified P2 over all
leading-critical cap-two inputs**. It is still not a \(Q_N\)-sharp or
\(C_N\) counterexample, and does not defeat an existential P1 that imposes
an additional feedback-excess selection criterion.

## 1. nonmono63.md — pass

The proof removes all use of terminal monotonicity. Its key estimates check:

- \(F_c=\prod_i f_{b,c_i}\le1\), so the actual modular budget makes the
  exceptional-set Lebesgue measure, and hence its \(P_c\)-mass,
  \(O(N^{-K_{\rm big}})\).
- On the lower core, \(\Phi_N\ge H\) and \(e^{-S}\ge
  e^{-s+K_{\rm big}L}\), giving \(E[H I]\le N^{-K_{\rm big}}\), which is
  exactly what is needed for both \(ED\) and \(E[DI]\).
- The scalar \(U(S)\) is monotone even though \(P=(H-1)_+\) is not.
  The identity
  \[
  E|P-U|=EU-EP+2E(P-U)_+
  \]
  has the correct sign, and \(EP\ge ED-\eta\) supplies the necessary
  upper bound on \(EU-EP\).
- The \(\delta\)-moment bounds are valid: with
  \(\delta=(a+S-ES)_+\), \(a=d_N+O(rL)>0\),
  \(E\delta\le a+\operatorname{Var}S/a\) and
  \(E\delta^2\le a^2+\operatorname{Var}S\). The relative \(1/L\) error
  is correct. The central box is used only where its denominator is
  positive; the cap handles its complement.
- Replacing \(P\) by \(U\) in the location response is legitimate:
  \(E[U\ell_i]\le0\) by Stieltjes integration by parts, while the
  \(L^1\) replacement errors pair with the central score supremum and
  Fisher tails at the prior \(o(L/N)\) scale.

Take the already stated \(K_{\rm big}>20\) (or any fixed sufficiently
large value). The propagation and Hessian arguments then retain their
former rates. This is a real fixed-block theorem, not a physical
representation result.

## 2. variance63.md — repaired diagnostic

Sections 3--4 correctly identify that regret-to-variance needs a converse
curvature/boundary condition, and (7) is correctly isolated as the missing
global excess-to-regret transfer.

The repaired depth-limit coarea calculation now supplies the missing fact:
\[
 F_y(c)=E h_N(S_{N,c})=:F(c)\quad\hbox{for almost every }y.
\]
Thus, once the capped scalar common-window theorem applies, strict
concavity gives one common constrained argmax and
\(\operatorname{Var}_\nu(a_y)=0\). This is stronger than a constant
witness-label observation and closes the maximizer-versus-label gap.

The new cap-two ramp in section 5 strengthens this conclusion at the broad
leading-critical scale. Its exact mean equation, derivative estimate, and
Gaussian integral budget calculation are correct:
\[
 q=N\mu-d_N+O(L),\qquad -\log J(H_q)=N\mu-d_N+O(L).
\]
After an arbitrarily small threshold decrease, reception is strict. Thus it
is a genuine zero-argmax-variance countertest to P2 quantified over *all*
cap-two inputs at that sole budget scale.

The note correctly still makes no \(Q_N\)-sharp endpoint claim:
\(\log\mathsf{Sc}=\log Q_N-o(L)\) is not established. It therefore does
not refute an existential P1 or a P2-excess formulation that selects a
feedback-excess packet, but it does rule out unqualified leading-critical
variance forcing.

## 3. angular63.md — valid obstruction, correctly non-efficient

The Jacobian, Gamma--Bernoulli law, conditional angular kernel, and
first/second transport differentiation are consistent. In particular, the
repaired section 5 correctly compares with the frozen-angular radial
response: the full Bernoulli reweighting term is not claimed small.

The conditional example (14) needs its elementary parameter regime stated,
namely \(0<a<c<s/2\); under it the displayed
\(a/[2(s-c)]\) calculation is correct. The unconditioned sign example
also gives positive curvature. These establish only that a universal
negative Hessian for arbitrary bounded raw angular terminals is false.

Crucially, neither construction is shown to have critical modular cost,
near receiving level, endpoint efficiency, or the fixed-terminal
block-budget structure required by nonmono63.md. They must not be promoted
to a counterexample to critical efficient physical curvature. The
conditional transport-seminorm route is a proposed bridge condition, not a
proved rigidity theorem.
