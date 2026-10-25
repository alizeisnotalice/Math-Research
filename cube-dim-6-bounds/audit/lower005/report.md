# lower005: rigorous length pruning at the R6-003 spacing

This audit concerns only the standard centered consecutive-block, positive
half-cell family at the fixed R6-003 breakpoint spacing.  It does not address
arbitrary shifted blocks, the separate period-two construction, or a new
maximal-operator bound.

For a length vector `ell`, write `P=prod ell`.  A common-radius center box
can be nonempty only if, for each coordinate,

```
(ell_j-1) a_j < P^(1/6).                                      (1)
```

Multiplying (1) yields the additional necessary geometric-mean condition

```
prod a_j * prod ((ell_j-1)/ell_j) < 1.                         (2)
```

`span_prune.py` evaluates both conditions using rational isolating intervals
of width `10^-35` for every required sixth root and for all spacing terms.
There is no floating comparison in a classification.

| length cap | vectors | excluded by (1) | extra exclusions by (2) | not excluded by (1)--(2) | new-shell survivors |
| --- | ---: | ---: | ---: | ---: | ---: |
| 5 | 15,625 | 15,481 | 0 | 144 | 17 |
| 6 | 46,656 | 46,505 | 0 | 151 | 7 |

Thus the product condition is mathematically valid but contributed no
additional exclusions in these two finite boxes; it must not be credited
with a pruning it did not perform.  The reason the apparent all-long tail
patterns survive (2) is that a short first coordinate contributes a small
factor `(ell_1-1)/ell_1`.

The initial containment check in this round had a factor-two radius error;
its containment conclusion is withdrawn.  See
`correction_radius_factor.md`.  The corrected audit, which explicitly uses
`r=(prod ell)^(1/6)/2`, gives the following valid direct-single-box results:

| shell | span survivors | certified contained in one old `<=4` box | certified not contained in any one old `<=4` box |
| --- | ---: | ---: | ---: |
| `max ell=5` | 17 | 12 | 5 |
| `max ell=6` | 7 | 7 | 0 |

The five cap-5 vectors not directly contained are the four vectors obtained
by placing the sole `4` in one tail coordinate of `(3,5,5,5,5,5)`, together
with `(3,5,5,5,5,5)`.  They still may be covered by the
*union* of multiple old boxes, or they may enlarge the union; this audit does
not decide that question.  Thus cap 6's seven new vectors can be discarded
by direct containment in this standard family, while cap 5 requires a full
union-coverage or new-union analysis for its five remaining vectors.

Actual commands used the existing Python 3.15 environment and completed in
under four seconds total:

```
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/audit/lower005/span_prune.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/audit/lower005/containment_check.py
```
