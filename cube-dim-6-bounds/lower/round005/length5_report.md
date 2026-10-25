# Round 005: capped length-five witness extension

## Finite scope and construction

The certified R6-003 breakpoint spacing and its centered odd/even witness
rule were retained.  For a block length vector `ell`, the selected source
mass is `m=product ell_j` and the single common radius is
`r=m^(1/6)/2`; hence a checked block contributes a lambda-one atomic
witness.  For each nonempty center interval, both endpoints were checked
against both endpoint atoms in every coordinate.  The value below is always
computed as the union of all old and newly retained boxes, not as a sum of
box volumes.

There are exactly `5^6-4^6=11529` new vectors with at least one entry 5.
The received finite cap was 11,429, so the run obeyed that smaller cap.  It
processed vectors in decreasing `product ell_j`, breaking ties
lexicographically; 100 vectors remain untested (their ordered-list
SHA-256 is `13b37c46ac6f31e14009742b69daf28d6a65983e8b4c2cf4624fc647add7e6fe`).
Consequently this report does **not** claim to have exhausted the complete
length-five family.

The prior 127 nonempty / 96 maximal length-at-most-four boxes were rebuilt
only as a reference collection for containment and a joint-union
calculation.  This was not a separate baseline search.

## Actual run

The 11,429 screened new vectors produced only 17 numerically nonempty
witness boxes.  Twelve were numerically contained in one of the 96 old
maximal boxes.  Five were not directly contained in a single old maximal
box:

```
(3,5,5,5,5,5)
(3,5,4,5,5,5)
(3,5,5,4,5,5)
(3,5,5,5,4,5)
(3,5,5,5,5,4)
```

After joint containment reduction there were 101 maximal boxes.  The
floating joint union and old-reference union were, respectively,

```
64 |U_old ∪ U_new| = 2.0669238260108800
64 |U_old|         = 2.0669238260108806
```

The displayed difference is `-4.44e-16`, ordinary floating roundoff rather
than evidence of a decrease.  There is no visible numerical increase over
the frozen R6-003 construction, so no candidate was promoted to exact
certification.  Runtime was 0.081 seconds, well below the 1,200-second
budget.

## Certificate status and limitation

The direct-containment condition needed for a finite exact containment
certificate did not hold: five boxes were not inside one old maximal box.
Their apparent coverage by the *union* of old boxes is only a numerical
observation in this run.  Exact treatment of these new radius values would
involve `5^(1/6)` in addition to `u=2^(1/6),v=3^(1/6)`, i.e. the 216-degree
field explicitly excluded from forced expansion in this batch.  No exact
containment, exact union equality, new lower bound, or general-L1 transfer
claim is made.

Thus the only excluded object is a numerical promotion for the stated
11,429-index priority subset under this specific witness rule.  The 100
omitted indices and all other length-five witness choices remain open.

Reproduce with:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/lower/round005/length5_screen.py
```

The main-ledger input had SHA-256
`da73e41e9987caa273c2438ebf9d73e582c699735ee16f688a83daff04190fd9`.
At execution the source SHA-256 was
`2b6ce7f5c101b2227dda570b7cb4115699facc1cc8716da92d32be82bde4ac59`;
the machine-readable result file is included separately.
