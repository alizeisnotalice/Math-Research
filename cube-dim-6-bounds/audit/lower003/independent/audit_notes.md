# lower003 local-spacing audit — interim result

This audit registered exactly twelve numerical coordinate directions in
`registration.md`.  They were recomputed from all `4^6` count vectors; no
author candidate list was imported.  The sole improving pilot was

```
a_2 -> (10001/10000) a_2,
```

with numerical ratio `2.06516698104056695845...`, about `3.10e-5` above the
then complete R6 certificate upper endpoint.  It was sent immediately to the
independent exact-certification owner.  The numerical pilot is not itself a
new bound.

The endpoint-choice list and the retained 96 boxes stayed unchanged in all
twelve samples.  That condition is not sufficient to freeze the union
decomposition: the recorded recursive union trace changes in most coordinate
directions, including several where the initial fixed-order jet gives the
wrong finite difference.  Consequently this audit makes no general local
optimality statement and uses the jet only as a diagnostic.  For the `a_2`
positive direction the union trace remains the baseline trace and the jet
increment agrees with the computed increment, which explains why it was a
particularly suitable candidate for exact certification.

The exact verifier subsequently supplied by the main research branch is
reviewed by `replay_main_inward_candidate.py`.  The mathematical interface
being audited is unchanged by the `a_2` rescaling: the rectangular lattice
still has one unit atom per fundamental cell; `64` is the reflected physical
measure of the selected half-cell union, not an assumption that the lattice
covolume is one.  For each selected box the code verifies all coordinate
source endpoints against the same radius and verifies `(2r)^6 <= m`; inward
rounding supplies strict geometric room.  With uniform bumps of support
radius below half the recorded margin, these are finite ordinary Lebesgue
`L1` inputs.  The source collar `[-K-2,K+2]^6` contains every block needed by
receiver cells `[-K,K]^6`, so the displayed finite ratio is a valid lower
bound for the full maximal superlevel set.  Its direct comparisons use the
full frozen R6 upper endpoint and the R5/R4 upper endpoints, rather than
displayed decimals.

This is a review conclusion only; the authoritative current ledger remains
the main project ledger after its owner integrates the independent results.
