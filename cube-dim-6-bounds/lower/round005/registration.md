# Round 005: length-five witness extension

**Frozen comparator.**  The comparison is the complete upper endpoint of
the current main-ledger entry `R6-003-breakpoint`, read and hashed at run
time.  Display decimals are not comparison certificates.

**Mechanism.**  Retain the certified breakpoint spacing
`a2=(u^3+u^4v^3)/4` and the original centered odd/even block witness rule,
but allow block lengths in `{1,...,5}^6` having at least one coordinate 5.
Each block continues to use one radius `r=(product ell_j)^(1/6)/2`; a
surviving center box is checked against every endpoint atom of its chosen
block.  The comparison union is the full union of old and new witness
boxes, never a sum of individual volumes.

**Finite budget discrepancy and deterministic handling.**  The full new
index set has `5^6-4^6=11529`, whereas the received hard limit is 11,429.
This round obeys the smaller number: order new indices by decreasing block
mass `product ell_j`, then lexicographically, and screen the first 11,429.
The remaining 100 are explicitly recorded as untested.  Thus a negative
outcome cannot be stated as exhaustion of all length-five indices.

The old 4-by-...-by-4 witness boxes are rebuilt only as reference boxes for
containment and joint-union comparison; this is not a separate rerun of the
old experimental baseline.  Runtime is limited to 1,200 seconds.  If the
numeric joint union visibly improves, halt promotion and create a separate
exact `Q[u,v]/(u^6-2,v^6-3)` root-isolation/box certificate.  If every
screened new nonempty box is directly contained in an old exact witness box,
produce a finite exact containment certificate instead.  Otherwise report
the numerical and proof status without claiming exclusion.

## Correction sub-batch: omitted 100 and rational union coverage

The intended new-index total is 11,529, not 11,429.  This correction runs
only the 100 entries omitted by the preceding deterministic ordering; it
does not repeat the first 11,429.  Old boxes may be rebuilt solely as
reference inputs for containment acceptance, not as a new baseline search.

For the five numerically non-single-box-contained new boxes, a separate
certificate attempt normalizes to `s_j=2x_j/a_j`.  Rational sixth-root
intervals give *inner* old boxes and *outer* new boxes, preserving exact
shared endpoints 0 and 1.  At `D=10^30` and then `D=10^60`, acceptance is
the exact Fraction residual
`volume(union(old_inner,new_outer))-volume(union(old_inner)) = 0`.
An unresolved or nonzero residual is reported; no 216-dimensional field is
introduced.
