# Round 003: nearest upper endpoint-order breakpoint

The event-driven batch enumerated raw-box min/max branch switches,
zero-width events, and all pairs of possible second-coordinate endpoints.
It evaluated the nearest sixteen events around the independently certified
`k=16` spacing.  The other five coordinates were held fixed; the affine
argument applies to the full raw-box union, not to a dynamically selected
96-box representation.

The selected upper event has the exact form

```
a2=(u^3+u^4 v^3)/4,
r_(1,1,1,2,2,2) = 2 a2 - r_(2,4,2,3,3,3).
```

The discovery value at that event was about `2.06692382601088`.  A single
primary number-field reconstruction then regenerated 127 nonempty boxes,
kept 96 after containment removal, checked all coordinate endpoint
containments, and computed the full union.  It gave

```
206692382601088/10^14 < R_break < 206692382601089/10^14.
```

The code makes an exact comparison to the full `k=16` interval upper.
Because this primary implementation imports the round-001 arithmetic
layer, it is only a promotion candidate until an independent inward
rational union and finite-L1 certificate pass.  The batch did not launch
the parked period-two mass experiment or any integer continuation grid.
