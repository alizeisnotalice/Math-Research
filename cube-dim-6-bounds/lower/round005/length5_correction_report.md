# Round 005 correction: all 11,529 cap-five indices and union-cover attempt

## Corrected finite scope

The prior cap was a transcription error: the exact new-shell size is
`5^6-4^6=11,529`.  This correction ran **only** the 100 vectors omitted by
the previous decreasing-product order; it did not rerun the first 11,429.
All 100 produced an empty numerical canonical centered half-cell witness
box.  Hence all 17 numerically nonempty cap-five new boxes had already
appeared in the first screen.

This completes the numerical enumeration of the stated standard centered
half-cell length-at-most-five family.  It does not concern shifted blocks,
period-two weights, length six, or any other mechanism.

## Rational union-containment attempt for the five boundary boxes

The five cap-five boxes not numerically contained in one old maximal box
were

```
(3,5,5,5,5,5)
(3,5,4,5,5,5)
(3,5,5,4,5,5)
(3,5,5,5,4,5)
(3,5,5,5,5,4).
```

For a rigorous sufficient cover test, coordinates were normalized as
`s_j=2x_j/a_j`.  For old `ell_j<=4` boxes, rational sixth-root intervals
gave inner intervals

```
odd ell:  [0, min(1, q_lower-ell+1)]
even ell: [max(0, ell-q_lower), 1],  q=2r/a_j,
```

with the direction adjusted in the even case so the old interval is an
inner one; the new five boxes used the corresponding outer endpoint
bounds.  Shared endpoints 0 and 1 were represented exactly.  All interval
and union computations used `Fraction`, with roots of 2, 3, and 5 enclosed
by rational sixth-root intervals.

The required sufficient equality

`volume(U_old_inner union U_new_outer) = volume(U_old_inner)`

did **not** hold.  The exactly positive residual was approximately
`5.99970719705012078045474456773e-35` at `D=10^30` and
`5.99970719705012078045474453560e-65` at `D=10^60`.  Its decay is
consistent with a boundary-touching situation, but that is only an
inference: this inward/outward test does not certify union containment.

Accordingly, no exact containment conclusion is drawn.  In particular, the
earlier floating equality of the joint and old unions remains numerical
evidence only.  No 216-dimensional number-field calculation was forced,
as required by the round instructions.

## Evidence files

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/lower/round005/length5_remaining100_screen.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/lower/round005/rational_union_cover.py
```

Actual results: the remaining-100 result hash is
`e0d5308d67fb5e371143c568c289ffbf388c0624327c645cdfdc8b2524a3ac65`;
the rational-cover result hash is
`be4e6f5e557a13923b95d7141cc9b844332d7831fad93c2f119bc68873fdaf0e`.
The corresponding source hashes are `98445c5c1c6380025c3aa1fd9263847aee4e7b4c8692914c5eb0a63293b0ab40`
and `7e5f21ef4fb062111016c6f200a188309043b12da1549deb50f8c4878ca673df`.
