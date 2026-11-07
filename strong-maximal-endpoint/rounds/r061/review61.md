# r061 independent audit

## Verdict

Sections 1--4 of `block61.md` are correct after the added condition `L>max(W,q_N)`. They obstruct the crude global-to-local budget estimate only, not an endpoint counterexample.

The product-mesh quantile conclusion passes after addition of the completed, self-contained `local_md61.md`; its Fourier/tilting proof has been independently checked. The updated all-feedback L1 extension has the right Young--Fenchel tail mechanism; the cap is needed only for its L2 corollary.

## Checked calculations

In the band example, translated-cone coordinates have Jacobian `1/2` and allowed width `2 min(z,L)`. Hence, for `L>=q_N`,
\[
 J_u/\Phi_N(a)=\int_{q_N}^\infty e^{-z}\min(z,L)dz
 =(q_N+1)e^{-q_N}-e^{-L}.
\]
For `L<q_N` it is `Le^{-q_N}`. This verifies (8), and the Jacobians in (4)--(6) are also right. Chebyshev gives `P(|D|>W)<=mu/(64 log N)`, enough for strict reception.

For the mesh, the upper half of (6) follows from
\(\Phi_N(Ah)\le hA(1+\log(Ac_N))^{N-1}\) and `int h_N=b_N(1-R_N)`. The unchanged-area lower bound is the needed nonexponential denominator control. A local Wolfram calculation confirmed `Var(Y+mu J)=mu` from `mu=1+exp(-mu)`; it is only sanity checking.

## Required repairs

1. **Bootstrap (9)--(10).** This is now closed by the new self-contained `local_md61.md`.  Its uniform tilted local CLT follows from the explicit characteristic-function bound and a three-region Fourier split; exact exponential tilting and the Laplace integral give the density/CDF moderate-deviation formulae.  The fixed-interval argument then correctly chooses `C>kappa sqrt(mu)`, brackets `qtilde_N`, and only then uses MVT.  The required statement is
   \[
   f_m(x)={\log N\over\sqrt\mu N}(1+o(1)),\qquad
   F_m(q_N\pm C\log N)-F_m(q_N)
    =\pm {C\log^2N\over\sqrt\mu N}(1+o(1))
   \]
   with one-sided bracketing constants. Use the CDF mass perturbation
   `A_N^-1(gamma_N^-1-1)=kappa log^2(N)/N(1+o(1))` first to place
   `tilde q_N` in that interval.  Cite `local_md61.md`, rather than the unfinished r016 skeleton.

2. **Actual mixed modular.** Insert
   \[
   \Phi_N(A_N)/\Phi_{N-2}(A_N)=(1+\log A_N)^2=1+o(1).
   \]
   Thus the product denominator is `Phi_N(A_N)D_N(A_N)`, and use of the `N-2` core changes a log score only by `o(1)`.

3. **Receiver scope and finite truncation.** The latest core-section/zero-volume label is correct. Replace “any fixed positive strict margin can be `o(R_N)`” by a positive N-dependent margin. Give a common threshold decrease or sub-cap height increase for strictness, bound its score change by `o(log N)`, then select depth so truncation loss is below that margin.

4. **Fixed blocks.** Conditional on the r046 weighted-tail statement, (14)--(17) prove L1 for unbounded H; only the L2 estimate needs a cap. State the exact proved r046 lemma and location: it must bound `E[H(T)1_{C_a^c}]` from the actual modular budget for arbitrary H, not from diagnostic numerics.

5. **Feedback density/MGF.** The latest text correctly states `c in [0,infty]`, separates \(\sup_zq_c(z)\le2\) from \(q_c(z)\le2e^{-z}\), and supplies the uniform MGF calculation including `c=infty`. The cap-free L1 proof properly uses the r032-style Young bound with
   \(f_{\pi,N}(t)/e^{-t}\le2^Nt^{N-1}/\Gamma(N)\),
   `Lambda=e^s N^(-K-2)`, and
   \(\Phi_N^*(a)\le a\exp(a^{1/(N-1)})\). The selected exact budget exponent and `pC_K/N<=1/2` are now recorded, so the weighted `N^-K` tail closes and cap remains necessary only for L2. Minor typo: after (27), "last term in (25)" should refer to (26).

## Finite certificate

The finished directed run gives `0.112898286568... > 0.1128`, after an explicit tail subtraction. Its interval algebra is coherent: `P=e^(x-u)f`, affine terms use exact Gamma--Bernoulli CDF increments, and the omitted lower coordinate tail is empty because `x>=u>=-15.2`. The code returns the exact value `q=1` on `[1,N]`. Chords between lower endpoint values remain below the true concave q.

The conclusion is one-way: the finite 33-anchor dual lower bound also lower-bounds the primal cost for any continuous constraint system containing those anchors. It neither constructs a feasible continuum input nor proves a halo/endpoint theorem. Keep that distinction when combining it with the separate cone-volume certificate.


Final integration: the new local MD lemma, N-dependent strictness margin and core-sectional scope close the listed repairs. The accepted conclusions remain the universal terminal-profile stability, the fixed-curve cone-witness obstruction and the distinguished-section product calculation; none identifies the general endpoint constant.

Section 7 was separately reviewed and accepted: the pre-reception weighted amplitude tail is uniform, one common truncation/rescaling preserves every previously receiving strategy, and its actual modular log loss is O(N^(1−K)). The step-core calibration is explicitly not Q_N-sharp at the first deficit scale.
