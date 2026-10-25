# Round 004 period-two mass pilot: stopped at geometry budget

The target was a new mechanism: along coordinate 1, use alternating masses
`1+delta,1-delta` over the full cell `[0,2a1] x product [0,aj]`, with the
certified breakpoint spacing retained.  The selected common-radius block
mass was computed from the parity of its first-coordinate start.

An initial implementation clipped a single representative block into the
period cell.  Its `delta=0` calibration produced `0.6114025649`, which is
incompatible with the already certified unmodulated construction.  It was
therefore rejected: it omitted periodic images of center intervals that
cross a cell boundary.

The corrected implementation projected every center interval to the torus
and carried the corresponding source-block period shift in its containment
check.  During the required `delta=0` calibration it generated the 5001st
raw box, above the registered cap of 5,000, before a union was evaluated.
The `delta=1/20` pilot and all remaining axis/delta combinations were not
started.  No numerical ratio, exact certificate, or new lower bound is
claimed.

This only records that the present uncompressed periodic-image expansion
is too large for the round's finite budget.  It does not exclude the
period-two construction.  A future round would need a quotient-cell or
translation-orbit compression with a separately proved equivalence before
resuming this mechanism; the box cap must not be enlarged retroactively.
