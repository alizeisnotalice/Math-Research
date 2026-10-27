# R032 general-cover audits and R034 threshold

2026-10-28. All three specified R032 theories are **ACCEPTED within their
explicit scopes**. No counterexample to their stated mathematical claims
was found. None of these three documents alone supplies a new explicit
uniform numerical C_n upper bound.

| Source | Verdict | Material scope restriction | Review |
|---|---|---|---|
| geometry/multiscale_theory | ACCEPT | Single-box no-growth obstruction assumes the stated high-dimensional t,tau scale; fractional charging still needs certified regions and one common overlap capacity | `multiscale_review.md` |
| soft/joint_cover_theory | ACCEPT | Exact LP/flow interfaces require a uniform certificate over every finite geometry before yielding an improved C_n | `joint_cover_review.md` |
| parameters/anisotropic_cover | ACCEPT | Budget-family and global-ratio payment obstructions do not bound the actual C_n from below or exclude other covering algorithms | `anisotropic_review.md` |

The anisotropic review was delegated, as explicitly authorized, to one Sol
agent in its independent clone. Its report was copied verbatim and its hash
verified: `5b511db1d2a5efa869cc9b0e4dadf5cbc632a4d0cbcb390bd6ea439a31e1f5c2`.
It notes one harmless endpoint qualification: when delta_i=1 and rho=1,
the strict definition has N_i=1 rather than 2. The asserted N_i=2 statement
is correct almost everywhere and all integral comparisons remain valid.

The source hashes, including each exact reviewed version, are recorded in
`small_receipt.json`. The local exact checker actually passed; there was no
parameter search, n=7 computation or MATLAB. No frozen R033 file was changed.
No commit, push, fetch or rebase was performed by this audit. The explicitly
authorized isolated child clone is documented in its own report.

Separately, `q005_threshold.md` proves that the accepted R034 hand-picked
q=1/20 example gives, for every integer n with L=log n>=2^22,

```
W_n<T+n/(40sqrt L),
C_n<T-n/(40sqrt L),  T=nL-n log L.
```

The proof uses a decreasing explicit function and the exact Fraction base
bound 2474666053/2684354560<1. This is a sufficient threshold, not a claim
that it is the smallest one.
