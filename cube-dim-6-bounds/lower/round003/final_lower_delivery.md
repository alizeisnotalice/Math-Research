# Round 003 lower delivery: breakpoint refinement

## Result

The current certified lower result of this round is

```
C6 > 5167309565027/2500000000000 = 2.0669238260108.
```

It is an internal, strict Lebesgue lower certificate for centered,
axis-parallel cubes with one common radius and ordinary finite smooth
`L1` inputs.  It is not an upper bound, a claim of optimality for the
rectangular family, or a new geometric mechanism.

The primary number-field union calculation gives the narrower interval

```
206692382601088/10^14 < R < 206692382601089/10^14.
```

The independent inward rational certificate has a finite L1 witness with
receiver radius `2^50`, two-cell source collar, a positive uniform bump
halfwidth, and a strict threshold.  It compares its full rational lower
directly with the complete R6-003-k16 interval upper, rather than with a
displayed decimal.

## Construction and proof interface

Start with the round-001 rectangular lattice and retain its first,
third, fourth, fifth and sixth spacings.  Put

```
u^6=2, v^6=3, u,v>0,
a2=(u^3+u^4 v^3)/4.
```

This is the nearest upper endpoint-order event after the certified local
`k=16` line point.  It is characterized by

```
r_(1,1,1,2,2,2)=2a2-r_(2,4,2,3,3,3).
```

The author certificate rebuilds all 127 raw nonempty boxes, verifies each
block with a common radius `r=(prod ell)^(1/6)/2`, retains 96 noncontained
boxes only as a computation reduction, and computes the full union in
`Q(u,v)`.  The independent certificate instead inward-rationalizes raw
boxes and uses an integer union; it does not import author arithmetic.
The finite-block, inner-box, same-mass smoothing proof in
`second_spacing_transfer.md` gives the passage to general L1 inputs and
all scales.

## Event batch and stop condition

The breakpoint screen treated the full raw-box union.  With only `a2`
variable, its other five coordinate endpoints are fixed.  When all raw
`a2` endpoint branches, zero-width statuses and endpoint orders are fixed,
finite five-dimensional cell decomposition makes the union affine in
`a2`.  It enumerated 1,461 numerical event values and evaluated the 16
nearest to the certified `k=16` point; no integer continuation grid or
period-two mass batch was started.  The selected event was certified once.
This round now stops lower candidate generation.  A later round may use
other axes or the parked period-two construction after a fresh bounded
registration.

## Minimal reproduction and source locks

The following are the required integration paths relative to the main
project root (the author-produced copies currently reside in this
agent's independent clone until the main agent stages them):

1. `lower/round001_rectangular_certificate.py` — initial exact arithmetic
   layer, SHA-256 `9072275d47bd66c3ba9ec7de0c435efdba0785c39a02b122ea1ebb960403ec82`.
2. `lower/round003/certify_second_spacing_breakpoint.py` — primary event
   reconstruction, SHA-256 `0a0d61d1e6fb05ef58cb016709b93eed195944bc94cb9a4a09d4b11710ab2db2`.
3. `lower/round003/second_spacing_breakpoint_certificate_results.json` —
   primary output, SHA-256 `bd62b466f22c3c106a42cb9cfbceb8e212b6338ac5bf80ff976ef4f29b4613d3`.
4. `audit/lower003/rational_inward_breakpoint.py` — independent inward
   verifier, SHA-256 `5048611df7f55070baac57138d34373f3bd2b1979dedfd15143880ebc13fff5a`.
5. `audit/lower003/rational_inward_breakpoint_results.json` — independent
   result, SHA-256 `fe9ef75a71ade0397671cdbda502169a03b33f7b91f741b75e2a220d2dd289fe`.

Actual author-side runs: the 15-point local line finished exit 0 in
0.493 seconds; the single promoted breakpoint exact computation finished
exit 0 in 10.623 seconds.  No MATLAB session or background process was
started for this lower batch.
