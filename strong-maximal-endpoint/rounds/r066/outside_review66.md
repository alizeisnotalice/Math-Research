# Second independent review: outside_fixed66.md

## Verdict

**PASS for the deliberately narrow fixed-\(r\), equal-block, bounded-cap
radial model. The final text supplies the four stated repairs and explicitly defines the standard \(d_N\) quantile quantities.** I found no
remaining material obstruction in the Cramér, clipped-inverse, or
mean/variance-bootstrap mechanism. The local theorem is a localization result for each predetermined action vector; its broad-grid corollary controls receiver-dependent block actions. it does not establish a physical
action-dependent angular theorem or a general strong-maximal bound.

The needed \(d_N\) facts should be cited or restated where first used:
\[
d_N=\sqrt{\mu N}\,z_N,\qquad
\bar\Phi(z_N)=\sqrt{\mu/N}\,z_N,
\]
so \(d_N\asymp\sqrt{N\log N}\), and for every fixed \(O(L)\) shift of the
threshold,
\[
F_\mu(N\mu-d_N+O(L))=d_N/N+O(L^2/N)\asymp\rho.
\]
Without these facts, (22), (23), (25), and the later assertion that \(q\) is
a fixed multiple of \(d_N\) below the mean do not formally follow from the
scoped statement alone. They are available in the surrounding setup, so
this is a hypothesis/citation issue rather than a counterexample.

## Cramér tilt and density estimate

The fixed-\(\theta\) tilted summand is correctly
\[
\operatorname{Exp}(1-\theta)+c\operatorname{Bernoulli}
((1+e^{(1-\theta)c})^{-1}).
\]
For \(\theta\in[-\theta_0,0]\), all fixed centered moments are uniformly
bounded in \(c\): the only potentially large factors are \(c^k
e^{-(1-\theta)c}\), which are bounded. The exponential component gives a
uniform positive variance lower bound and
\[
|\chi_{c,\theta}(t)|\le
(1-\theta)/\sqrt{(1-\theta)^2+t^2}.
\]
Consequently the Fourier local-CLT route is legitimate for arbitrary
actions, not merely a compact action box.

There are two presentational points that are essential to that proof, both
already identified by root:

1. Choose the saddlepoint \(K'(\theta)=q\) before applying the centered
   density CLT. The lemma is used only after the coarse CDF bound has
   forced \(M-q>0\), and \(M-q=O(\sqrt{NL})\), so this negative tilt exists
   and lies in the stated fixed interval.
2. Write the high-frequency integral explicitly. After rescaling, the
   product version of the displayed characteristic-function bound gives an
   integrable tail of the form \(C\sqrt N e^{-cN}\) (plus the Gaussian
   intermediate tail). This justifies the uniform \(O(N^{-1/2})\)
   normalized density error rather than leaving a non-integrable
   pointwise assertion.

With these supplied, the Laplace identity has error
\(O((\lambda\sqrt N)^{-1})\), while its normal integral is
\(\asymp(1+\lambda)^{-1}\). Thus its *relative* error is
\(O(N^{-1/2})\). The uniform cumulant comparison gives
\(O(x^3/\sqrt N)\); at \(x\asymp\sqrt L\), multiplying by the critical
tail \(\asymp\rho\) gives the required absolute
\(O(L^2/N)\), not the inadequate remainder in the previous draft.

## Clipped scalar budget

This is now quantitatively sound. The joint density is at most one:
every one-block law is a convex mixture of translates of a rate-one Gamma
density, each of supremum at most one. Hence the modular budget makes both
exceptional sets have probability at most \(N^{-K}\). The cap means that
\(P>U(S-q)\) cannot occur on the clipped branch, and off that set it gives
the pointwise inequality \(P\le U(S-q)\). Therefore
\[
F_{\mathbf c}(q)\le E_{\mathbf c}U(S-q)+(A+1)N^{-K},
\]
which is sufficiently sharper than \(o(\rho)\).

The proposed enlargement to \(A\ge2\) removes a harmless degenerate
notation issue. The inverse Taylor bound, uniform exponential tail above
\(\delta N\), and the two moment bounds give \(EU=O(\rho)\) first and,
after mean localization,
\[
EU\le(d_N-\Delta)/N+O(L/N).
\]
No unproved conditional-density or angular-transport assertion is hidden
in this step.

## Bootstrap and full-action sequencing

The coarse bootstrap works as follows. The CDF upper bound forces
\(M-q\gtrsim\sqrt{NL}\); since \(M\le N\mu\), it gives
\(\Delta=O(\sqrt{NL})\). At this point (4) first implies that **every**
\(c_i\) is near \(\mu\), because \(r\) is fixed and each block has size
\(b=N/r\). Only after this observation may one replace the clipped gap by
\(\Delta\gtrsim Na^2\) and assert \(V=N\mu+O(Na)\). The specified repair
of this ordering is necessary and sufficient.

The moderate-deviation comparison then gives
\[
\Delta=O(\sqrt{N/L}),\qquad a=O((NL)^{-1/4}).
\]
At that scale, the Cramér estimate has absolute \(O(L^2/N)\) accuracy and
the normal comparison yields the one-sided lower CDF bound in (27).
Combining it with the sharpened clipped-envelope upper bound gives
\[
a^2\le C\rho a+C\rho^2,
\]
hence \(a=O(\rho)\). Fixed \(r\) converts this RMS conclusion to
\(\max_i|c_i-\mu|=O(\rho)\).

The final broad-grid ordering is also correct: before this localization,
relative block means can be of order \(N\), so the old
\(R=C\sqrt{NL}, W=N\) grid is invalid. The proposed uniform
\(R=CN,W=N^2\) grid has exponentially small discarded probability and
\(1+O(N^{-1})\) random-shift overlap. Localization must occur before
invoking any theorem that assumes the near-\(\mu\) box.

## Scope and residual caveat

This conclusion relies on the cap being fixed, the number of blocks being
fixed, and the radial translation step genuinely producing the scoped
budget/exact-success hypothesis. It does not by itself prove that
translation implication, nor repair a model with within-block
action-dependent angular kernels. Subject to those explicit scope limits,
the repaired proof is closed.
