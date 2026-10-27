# R052 candidate: later theorem activation

The R052 production report records the original conditional status at the time of enumeration. Subsequently the fixed-parameter flexible-budget theorem R051 and the fixed-tuple threshold R053 were independently accepted. For the retained tuple a=73/8, s=293/16, d=1023/2²⁰, R053's accepted theorem gives, for every integer n with L=log n≥2³²,

\[
C_n<F(A_n(t),B_n(t)-\gamma n)
<n\log n-n\log\log n-2^{-16}n,
\quad t=n/[(1-\sigma)L],\quad\sigma=sd^2,
\]

with γ and all finite guards as stated in R053. The R052 c_raw is a scalar surplus; it is not itself the complete operator cost, a continuous-family optimum or a dimension-independent threshold. The accepted R051 theorem applies to fixed admissible tuples with tuple-dependent sufficiently-large-n thresholds, not arbitrary n-dependent tuning.

Separately, R055's accepted tuple a=9, s=19, d=2⁻¹⁰ gives the same displayed coefficient 2⁻¹⁶ for L≥2²⁷. It is not the R052 top candidate and was checked independently of the search. A lower threshold by itself does not establish pointwise ordering of the two full costs.

Source identity for the archive (paths relative to the `cube-dim-7-bounds` project root):

- `upper/round051/flexible_budget/PROOF.md`: SHA-256 `aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6`; independent audit `audit/round051/flexible_budget/proof_review.md`: `66355a1d58d3b6a832e844645b9e61ca9e3c5a2f27bc0958947054aa9144f960`.
- `upper/round053/explicit_candidate_threshold/THRESHOLD_PROOF.md`: `26adb4e31f0026637ea17d8f7943bc2933dce97b8655e121f2604b6ddeea2775`; independent audit `audit/round053/threshold/review.md`: `eef083003c5020cc184e2ee312a0374052c1a4fa962fef7da753e22f0e9fba1b`.
- `upper/round055/simple_tuple_threshold/PROOF.md`: `5d2c61adf0e341b2a11daabec61daedd6e009b65a0139bb468711175f0288e9f`; independent audit `audit/round055/threshold/review.md`: `cd8440a63cd30bd7368db692595a35a309c8050e7921637edc473cc447b8b41b`.

This note records later acceptance; it does not rewrite the status of the original R052 run or substitute the retained-candidate checker for these independent proofs.
