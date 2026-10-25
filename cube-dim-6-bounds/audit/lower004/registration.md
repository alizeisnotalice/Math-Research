# lower004 — two-axis finite endpoint-event diagnostic

Date: 2026-10-25.  Starting point is the authenticated R6-003 breakpoint
spacing.  This round registers a finite numerical screen only: for coordinate
`a1` and for one representative tail coordinate `a3`, enumerate all positive
pair equalities among pre-branch coordinate endpoints arising from count
vectors whose other five coordinate slabs are nonempty.  For each axis retain
at most the eight nearest smaller and eight nearest larger distinct events,
then recompute the full union from all `4^6` count vectors at those event
values.  Thus each axis has at most 16 evaluated event parameters.

Acceptance objects are the event descriptors, raw nonempty and maximal box
counts, full-union numerical value, and comparison with the complete current
R6-003 upper endpoint read at execution time.  Equal box counts are never
treated as a combinatorial-cell certificate.  A numerical improvement is
only sent for exact certification; no new lower bound follows from this
screen.  Stop after these two finite batches or on a diagnostic failure.
