# Independent review of window72

## Verdict

**PASS with one explicit bookkeeping sentence required.**  The wide-window
argument rigorously closes the bounded-input, fixed-reference local
qualification module.  It does not close a general upper bound, and the
document correctly retains the common-source/anchor and moving-action
gaps.

## Verified points

1. For \(m(c)=1+c/(1+e^c)\), the identity \(m'(\mu)=0\) is correct.
   The stated second-derivative estimate gives
   \(Nm(c)-N\mu=O(D^2\log N)\) on
   \(|c-\mu|\le D\sqrt{\log N/N}\).

2. The centered-mgf estimates are valid for \(c\le2\):
   the exponential part is bounded by \(t^2\) on \(|t|\le1/2\), and
   Hoeffding's Bernoulli bound is \(c^2t^2/8\).  Thus the displayed
   \(2e^{-d^2/(6N)}\) tail follows.  The radius
   \(7\sqrt{N\log N}\) leaves a \(6\sqrt{N\log N}\) gap after the
   \(O(\log N)\) mean displacement, proving the uniform \(2N^{-6}\)
   tail.  No Gaussian approximation is being used.

3. The physical \(L^1\) identity is exact:
   the \(e^{-s}\) factor is canceled by the physical radial translation,
   and complete-label Markovness produces no extra label factor.
   Therefore the window has screening exponent \(1/2+o(1)\).

4. The coefficient \(\alpha=(9-\mu)/24\) is used only after specializing
   to the shifted deterministic reference \(c_N\).  It is not assigned to
   arbitrary critical \(c\).

5. The receiverwise profile is legitimate.  It is an auxiliary scalar
   profile, not a receiver-dependent alteration of the common source.
   Jensen gives its windowed Orlicz cost at most the same-source quantity
   \(j_W[\Phi(F)]\).  The tail repair has
   \(\eta_N=2AN^{-6}\), so \(N\eta_N=o(1)\); the cap remains inside a
   fixed enlarged cap.  For the final
   \(\eta_{\rm ref,N}=o(\log N/N)\) extension, the cost multiplier has
   logarithm \(O(N\eta_{\rm ref,N})=o(\log N)\), namely \(N^{o(1)}\).

## Required explicit transfer after rescaling

Equation (11) follows from r069 applied to
\(\widetilde h=h_W/(1-\eta_N)\), whereas Jensen initially controls
\(\Phi(h)\), not \(\Phi(\widetilde h)\).  Insert
\[
 \Phi_N(\widetilde h)
 \le(1-\eta_N)^{-N}\Phi_N(h)
 \le(1+o(1))[K_c^s\Phi_N(F_{u+s})](X,\sigma).
\]
Consequently the captured lower bound for \(\widetilde h\) transfers to
the actual common-source band cost with a \(1-o(1)=N^{o(1)}\) factor.
With this sentence, (11) is justified at the precision claimed.

## Potential note check

The finite-\(c\) extension in potential72 is correct.  The label-summed
radial shell mass is \(ds\), giving \(F_c(A{\bf1})=1/A\).
Exchangeability conditional on \(S\) yields the covariance identity, and
differentiation of \(F_c^{-1}\) cancels the
\(2/(N^2A)\) terms, leaving the strictly positive
conditional-variance integral.  The small-shell uniform-simplex lower
bound is computed correctly.  The deterministic-total-increment argument
is valid: diagonal potential transform \(1/A\) forces total exponent
\(A\), hence deterministic total drift; infinite divisibility then forces
every coordinate variance to vanish.

## Integration resolution

The required rescaling-to-source inequality has been inserted explicitly before (11) in the final window72.md. No additional source alteration is introduced.
