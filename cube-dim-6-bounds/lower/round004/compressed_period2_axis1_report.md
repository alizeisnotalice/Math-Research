# Round 004B: reflection-compressed axis-one period-two pilot

## Registered object and normalization

Keep the R6-003 breakpoint spacing
`a=(1+u^4/2,(u^3+u^4v^3)/4,(1+u)/2,(1+u)/2,(1+u)/2,(1+u)/2)`, where
`u^6=2` and `v^6=3`.  On lattice site index `k` in coordinate 1 use
weight `1+delta` for even `k` and `1-delta` for odd `k`; the remaining
coordinates have unit weights.  The full periodic cell is
`[-a1,a1] x product_{j>1}[-aj/2,aj/2]` and its source mass is exactly 2.

For a consecutive first-coordinate block `{s,...,s+ell1-1}` and unchanged
coordinate block lengths `ell2,...,ell6`, the selected mass is

`m=(ell2...ell6) ell1` if `ell1` is even, and
`m=(ell2...ell6)(ell1+delta)` or
`(ell2...ell6)(ell1-delta)` if `ell1` is odd and `s` is even or odd.

The cube radius is `r=m^(1/6)/2`, so that `(2r)^6=m`.  Every retained
center box is checked against both endpoint atoms in every coordinate;
therefore its cube contains the selected block (it may contain extra
atoms, which only strengthens the lambda-one witness inequality).

The selected center boxes lie in
`H=[0,a1] x product_{j>1}[0,aj/2]`.  Reflection in coordinates 2--6 is
the original unit-lattice symmetry.  In coordinate 1, reflection sends a
block start `s` to `-(s+ell1-1)`; for odd `ell1` it has the same parity as
`s`, and for even `ell1` the selected mass is parity-independent.  Hence
the mass rule and radius are preserved.  Except for null cell boundaries,
the 64 reflected half cells have equal witness-union measure.  Thus this
subfamily gives the periodic atomic ratio

`|U_full| / 2 = 64 |U_H| / 2 = 32 |U_H|`.

This is an admissible *witness subfamily* calculation.  It does not claim
to enumerate every cube that can witness the alternating-mass measure.
No finite-block/smoothing transfer is invoked here because no candidate is
promoted to a lower-bound certificate.

## Actual finite batch

The old full-torus representation had stopped at its 5,000-box cap.  This
new run did not increase that cap: it enumerated only starts `s` satisfying
the exact intersection condition

`(s+ell1-1)a1-r < a1` and `s*a1+r > 0`.

The finite outer range is derived with the maximum allowed radius and each
candidate is then retested with its parity-dependent actual radius.  Block
lengths were `ellj in {1,2,3,4}`.  The run was limited to 5,000 raw boxes
and 1,200 seconds, first at `delta=0` and then, only after calibration, at
`delta=1/20`.

| case | raw / maximal boxes | `32|U_H|` | status |
|---|---:|---:|---|
| `delta=0` | 181 / 144 | 2.0669238260108806 | matches existing numerical R6 witness union to `4.44e-16` |
| `delta=1/20` | 182 / 145 | 2.051499464558528 | noncompetitive numerical witness result |

Elapsed wall time was 0.141 seconds; the caps were not approached.  The
frozen complete R6-003 upper endpoint has decimal value
`2.066923826010880910097778283840054...`.  The second row is lower by
`0.015424361452352910097778283840054...`; it is not sent to exact union
certification.

## Evidence status and scope

This is numerical discovery evidence only.  The delta-zero agreement is a
geometric calibration of the half-cell/reflection implementation, not an
independent exact re-certification of R6-003.  The sole explicitly checked
new pattern (axis 1, `delta=1/20`, with the stated block witness rule) did
not improve the current R6 lower construction.  It neither rules out other
period-two witnesses nor disproves the alternating-mass mechanism.

Reproduce with:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/lower/round004/compressed_period2_axis1_screen.py
```

The code reads the integration clone's ledger only to freeze and hash its
comparison endpoint; it never writes outside this lower-bound work clone.
The input ledger SHA-256 is
`de0b0a2eb842bffee5dbf7dab9a58bd658cc155c799ffacff0976e1df4a2ec1b`.
At the recorded run the source SHA-256 was
`c10e5277c2c3bc0f2ef991c9009f18f6e7dd43515f2ce16c6d5643a01f91f465`,
and the results SHA-256 was
`36652b1ccd3c8e149b2b39e3fa63cd2342555e289b1754ad2b3e8423e3f0a355`.
