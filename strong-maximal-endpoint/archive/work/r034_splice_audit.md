# r034: audit of the proposed Bellman splice

This note concerns only the scalar `MAX` Bellman operator. It does not turn
the scalar result into a statement about the geometric strong maximal
operator.

Let

\[
 ({\cal T}w)(r)=\sup_{c\geq0}K_cw(r),\qquad
 K_cw(r)={\mathbb E}_c w(r-Z_c).
\]

The supersolution convention in this project is
\({\cal T}W_{t-1}\leq W_t\). The following elementary pasting statement is
fully valid and is the useful safe interface.

**Min-pasting lemma.** If `A_t` and `B_t` are nonnegative, global
supersolutions (for every state and every time), then
\(W_t=A_t\wedge B_t\) is again a global supersolution.

Indeed, positivity gives
\[
 {\cal T}(A_{t-1}\wedge B_{t-1})
 \leq ({\cal T}A_{t-1})\wedge({\cal T}B_{t-1})
 \leq A_t\wedge B_t.
\]
This needs neither a common maximizing action nor regularity at the interface.
It is therefore legitimate to paste the two *global* left/right barriers by a
minimum.

That lemma does **not** justify pasting a barrier which is known only on the
natural tube. Let `O_t` be the global outer minimum and `I_t` be an inner
potential. To use `W_t=O_t\wedge\widetilde I_t`, one needs a global extension
`\widetilde I_t` and, on every state where the inner branch is selected,
\[
  \sup_c K_c( O_{t-1}\wedge\widetilde I_{t-1})\leq\widetilde I_t. \tag{S}
\]
The stronger, convenient condition
\(\sup_cK_c\widetilde I_{t-1}\leq\widetilde I_t\) there is sufficient,
because `O wedge Itilde <= Itilde`. Merely proving the reference identity
for `K_mu I` is not a proof of (S).

For the quadratic Poisson proposal the precise missing comparison is visible
from
\[
 I_t=v_t+U_t,\qquad U_t=K_\mu U_{t-1}+g_t.
\]
At an inner state, an all-action supersolution estimate would require
\[
g_t(r)\ \geq\ \sup_c\{(K_c-K_\mu)v_{t-1}(r)
                    +(K_c-K_\mu)U_{t-1}(r)\}, \tag{P}
\]
up to a separately accounted summable error. The Poisson identity controls
only the `K_mu` part. It is not legitimate to replace the second term in (P)
by an expectation under the reference bridge.

There is a concrete reason that the known one-action argument for `v` does
not automatically extend to `v+U`. In the notation used in r033, the action
root residual for a generic continuation `F` obeys
\[
 D_F'(r,c)=(1+e^{-c})F'(r-c).
\]
For `F=v` this is negative, hence the root is a unique global maximizer.
Nothing in the reference-Poisson construction implies
\((v+U)'<0\) globally, or even on an enlarged tube. Thus an action which
stays wholly in the tube can already invalidate (P); interface-crossing is
not the first problem. A valid splice must either establish this monotonicity
(with enough margin), or prove (P) directly for all `c`.

## Cross-interface actions

There is a useful, but conditional, tail fact. The one-step increment density
satisfies
the uniform estimate
\[
 q_c(z)={e^{-z}(1+{\bf1}_{z\geq c})\over1+e^{-c}}\leq2e^{-z},
 \qquad z>0.
\]
Hence an inner tube and an enlarged tube separated by
\(h_t\asymp\sqrt{t\log t}\) have the uniform one-step exit estimate
\[
 \sup_c\mathbb P_c(Z_c\ge h_t)\le 2e^{-h_t}.
\]
Since the kernel only decreases `r`, the opposite-side displacement due to
the change of tube centre is only `O(1)`. If the extension/outer barrier is
bounded by a fixed exponential of rate strictly below one (or by the actual
linear left piece), its *weighted* exit contribution is at most a constant
times \(e^{-(1-\theta)h_t}\). This is summable in `t`.

It is essential that the last envelope be *proved*, not inferred from the
tail probability. At the critical exponent one has a literal counterexample:
for `Z` exponential of rate one and `F(r)=e^{-r}` on the positive half-line,
for `r>h>0`,
\[
 {\mathbb E}\,[F(r-Z){\bf1}_{\{h\leq Z\leq r\}}]
 =\int_h^r e^{-(r-z)}e^{-z}\,dz=(r-h)e^{-r}.
\]
This elementary integral was also checked symbolically with WolframScript.
The probability of the crossing event is exponentially small in `h`, but
its contribution relative to `F(r)` is not small at all. The terminal payoff
has exactly this critical exponential branch before it turns linear on the
left. Thus a bare ``wide buffer plus q_c tail'' argument is invalid. At
active large times it may be possible to replace this calculation by a
moderate-deviation estimate for the already-smoothed `v_t`, but that is an
additional lemma about the actual continuation and the actual outer branch.
It cannot be supplied by the uniform bound on `q_c` alone.

## Why reference-bridge rarity cannot be used for the outer annulus

There is a direct all-policy counterexample to any claim that the natural
outer annulus is rarely occupied. Put
\(q=N\mu-d_N\), take pure-exponential controls (mean one) for
\[
 k_N=\left\lfloor{d_N+2A\sqrt{N\log N}\over\mu-1}\right\rfloor
\]
elapsed steps, and subsequently use the reference action. At that time,
with `t=N-k_N` remaining and `G_{k_N}` the gamma sum,
\[
 r-t\mu=q-G_{k_N}-t\mu
 =(\mu-1)k_N-d_N-(G_{k_N}-k_N)
 =2A\sqrt{N\log N}+o_{\mathbb P}(\sqrt{N\log N}).
\]
Here `k_N=O(\sqrt{N\log N})` and the gamma fluctuation is only
`O_P(sqrt(k_N))`. Therefore the upper annulus
\([A,3A]\sqrt{t\log t}\) is occupied with probability tending to one under
an admissible policy, whereas it is rare under the reference bridge. A
splice proof must be pointwise/Bellman-valid there; it cannot discard that
region by reference occupation. This observation alone does not disprove a
small *weighted* source, since its density factor may be small.

## Conclusion

The global `min` splice itself is sound. A buffered splice can reduce crossing
terms only after a *subcritical or otherwise weighted* crossing estimate has
been proved; the bare increment tail does not furnish one at the critical
exponential branch. Independently, the first all-action issue for the inner
potential is (P), including actions entirely inside the tube. Any r034 proof
that cites only a reference-Poisson error, a raw crossing probability, or
reference bridge occupation has a blocking gap at exactly these points.

## Final check of r034 report, sections 1--3 and 5--6

Accepted as scoped. The affine/Gamma barrier has the right direction:
`b_N` is nondecreasing and `Z_c >= Exp(1)` under the common coupling, so
`K_c B <= G B`; the affine part supplies the complementary mean inequality.
The displayed Gamma formula is consistent with direct integration (the
last term is `-t F_{t+1}/N`). The minimum rule is only invoked for the two
global barriers. The two natural-boundary estimates use the global barriers,
not an unproved local splice.

The two counterexamples in section 5 have the correct limited force: (13)
refutes a reference-bridge occupancy shortcut but not a density-weighted
source bound, while (14) refutes replacement of a weighted crossing integral
by its probability. Equation (15) correctly remains the unsolved all-action
condition. Finally, (16) follows by writing `v=L^a+w`; its affine difference
is `(m(c)-m(a))/N-delta_a/N=(m(c)-mu)/N`. At `r<=0` the continuation remains
on the affine terminal branch, so the compensated reference value indeed
equals the true MAX value. The matching expansion (17) has the stated sign.
No claim in these sections promotes the candidate splice or matching drift
bookkeeping to a global sharp result.

## Addendum: valid strict-interior weighted crossing module

The preceding critical-exponent warning does not apply to the following
newly specified comparison, which is a real repair in the strict interior.
Suppose the patched potential `W` is nonnegative and globally bounded by the
left screen `O^-_{t-1}=L_{t-1}+B_{t-1}`, while the inner continuation `I` is
nonnegative. Then
\[
 (W-I)_+\leq W\leq O^-_{t-1}.
\]
For `r>=0`, `z>=0`, and `t<=N`, the elementary bound
\(B_{t-1}(r-z)\leq(r-z)_+/N\) gives
\[
 O^-_{t-1}(r-z)
 \leq1+(t-1)\mu/N+z/N\leq1+\mu+z/N. \tag{C}
\]
Together with `q_c(z)<=2e^-z`, uniformly in the action,
\[
 K_c\big((W-I)_+(r-Z_c){\bf1}_{Z_c\ge h}\big)
 \leq2e^{-h}\bigl[1+\mu+(h+1)/N\bigr]. \tag{D}
\]
WolframScript independently simplified the deterministic inequality (C).

On a fixed natural moderate-deviation band one may use the usual lower
bound `f_t(r)>=c_A t^{-K_A}`. Thus, for
`h=delta sqrt(t log t)`, (D) is smaller than `f_t(r)/t^M` for every fixed
`M`, once `t` is sufficiently large. With an inner tube of radius
`(A-delta) sqrt(t log t)` and an outer tube of radius
`A sqrt((t-1)log(t-1))`, a transition from the inner current state to the
left outer complement forces `Z_c>=delta' sqrt(t log t)` after absorbing the
`O(1)` centre displacement; transition to the right complement is impossible
for large `t` because `Z_c>=0`.

This is a valid all-action crossing estimate, unlike the raw-tail shortcut:
the positive discrepancy is controlled by the subcritical affine/Gamma left
screen, not by the critical exponential branch of `I`. It resolves only the
strict-interior interface term. It does not prove the annular comparison,
monotonicity/action localization of `v+U`, or the remaining all-action
residual.
