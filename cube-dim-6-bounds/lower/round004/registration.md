# Round 004: one-axis period-two alternating-mass pilot

**Frozen comparison.** The comparator is the complete R6-003-breakpoint
upper endpoint in the main ledger, not `2.0669238260108`.

**Mechanism.** Keep the certified breakpoint spacing
`a2=(u^3+u^4 v^3)/4`.  Along coordinate 1 only, assign masses `1+delta`
and `1-delta` to alternate lattice hyperplanes; all other coordinates
retain unit periodicity.  The fundamental cell is
`[0,2a1] x product_{j>1}[0,aj]` and has mass 2.  Every block uses one
common radius derived from its exact selected mass.  Its shifted center
box is intersected with this full period cell and the entire union is
computed; no box volumes are added independently.

**Finite batch.** Run a new-geometry `delta=0` calibration once, then only
axis 1 with `delta=1/20`.  A hard cap is 5,000 raw boxes and 20 minutes.
Only if the latter is geometrically valid and numerically clears the frozen
upper interval by `10^-7` may the follow-up batch consider axes 2--3 and
the two remaining deltas `1/10,1/5`, with at most nine axis/delta pairs in
total.  Any promotion requires a fresh exact union and inward L1 proof.

**Stop condition.** If calibration fails to reproduce the fixed-spacing
geometry at the necessary interface, or the axis-1 pilot does not promote,
report the tested scope and do not run the remaining pairs.

## Sub-batch B: reflection-compressed half-period witness family

The full torus-image expansion stopped before a union could be evaluated.
This is a new representation, not a retroactive enlargement of that budget.
Use the reflection-invariant fundamental cell
`[-a1,a1] x product_{j>1}[-aj/2,aj/2]` and its half cell
`H=[0,a1] x product_{j>1}[0,aj/2]`.  Alternation in coordinate 1 is
preserved by `x1 -> -x1`; the other five unit lattice directions have their
usual reflections.  Thus the reflected witness subfamily has periodic union
measure `64 |U_H|`.  Since one full cell has source mass 2, its lambda-one
ratio is `32 |U_H|`.

For coordinates 2--6 retain precisely the old odd/even centered-block
witness rule.  In coordinate 1 enumerate only consecutive blocks
`{s,...,s+ell1-1}` with `ell1 in {1,...,4}` whose center interval meets
`[0,a1]`; their mass is `ell1` for even `ell1`, and
`ell1 + delta` (respectively `ell1-delta`) for odd `ell1` and even
(respectively odd) `s`.  The common radius is one half the sixth root of
the complete selected block mass.  This is explicitly a witness subfamily:
no assertion that every admissible cube has been enumerated is needed.

**New finite batch and acceptance.**  With the same raw-box cap 5,000 and
20-minute wall-time budget, run delta=0 first.  It must reproduce the
existing fixed-spacing numerical union under this reflected witness rule;
otherwise report only the mismatch and do not evaluate delta=1/20.  If it
does, run the single delta=1/20 case and compare only as numeric discovery
against the frozen complete R6-003 upper endpoint.  Any apparent improvement
requires a separately written exact union certificate and the unchanged
finite-block/smoothing transfer before it is a lower bound.

## Sub-batch C: remaining pre-registered period-two combinations

The half-cell interface has independently passed audit.  Complete the
remaining eight members of the original finite nine-case family without
repeating axis 1, `delta=1/20`: axis 1 with `delta=1/10,1/5`; axis 2 with
`delta=1/20,1/10,1/5`; and representative equal-spacing tail axis 3 with
`delta=1/20,1/10,1/5`.  Coordinate permutation moves the active period-two
axis to the first code coordinate, after which exactly the Sub-batch B
half-cell witness rule and normalization `32|U_H|` apply.  Axes 4--6 are
not separate cases because their spacings and the unmodulated witness rule
are invariant under their coordinate permutations.

There is no new delta-zero calibration and no expansion of the case set.
Each case separately retains the 5,000 raw-box and 1,200-second limits.  A
case is eligible for exact certification only if its numerical value clears
the frozen complete R6-003 upper endpoint with a visible margin; otherwise
it is reported solely as a non-promoted numerical screen for this precisely
defined witness family.
