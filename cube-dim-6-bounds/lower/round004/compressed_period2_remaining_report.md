# Round 004C: remaining eight period-two witness screens

## Exact finite family screened

This completes, without repetition, the eight members remaining after
Round 004B had evaluated axis 1 with `delta=1/20`.  The active alternating
axis is moved to the first computational coordinate by a coordinate
permutation; the period cell then has source mass 2 and the independently
audited reflection calculation gives the numeric witness ratio `32|U_H|`.
Block lengths in every coordinate were in `{1,2,3,4}`.  The active-axis
starts are precisely those whose block center interval meets its half-cell,
and the other five coordinates use the original odd/even centered witness
rule.  Every box is endpoint-checked against the selected atoms and uses
one common radius from the entire selected block mass.

The three tails have equal spacing and otherwise symmetric witness rules,
so axis 3 is a representative for axes 3--6; no unregistered duplicate
tail cases were run.  Each screen independently had the registered 5,000
raw-box cap and 1,200-second limit.  The comparison was to the then-current
complete R6-003 exact upper endpoint, whose input ledger SHA-256 was
`c47ca8c7e14684f106f501d94f1dbdb6c34c843ff54599932f8d7d2f41de39d4`.

| active axis | delta | raw / maximal boxes | numeric `32|U_H|` | margin vs frozen upper |
|---:|---:|---:|---:|---:|
| 1 | 1/20 | 182 / 145 | 2.051499464558528 | -0.015424361452353 |
| 1 | 1/10 | 182 / 150 | 2.032080368286628 | -0.034843457724254 |
| 1 | 1/5 | 177 / 141 | 1.970546840709035 | -0.096376985301846 |
| 2 | 1/20 | 256 / 145 | 2.053230503211683 | -0.013693322799198 |
| 2 | 1/10 | 256 / 146 | 2.038312778346952 | -0.028611047663929 |
| 2 | 1/5 | 250 / 146 | 2.004010038348869 | -0.062913787662013 |
| 3 (tail representative) | 1/20 | 256 / 145 | 2.055080022377521 | -0.011843803633360 |
| 3 (tail representative) | 1/10 | 256 / 146 | 2.042376275643958 | -0.024547550366923 |
| 3 (tail representative) | 1/5 | 249 / 145 | 2.005867992821707 | -0.061055833189175 |

The eight new screen evaluations together took 2.34 seconds.  The best new
one is active axis 3 at `delta=1/20`; it remains visibly below the frozen
endpoint.  Consequently no exact union computation, exact comparison, or
finite-block/smoothing transfer was started for any member.

## Evidence boundary

All entries are numerical witness-subfamily screens, not lower-bound
certificates.  They exclude only these nine fully specified combinations
under this block-witness rule as an immediately competitive route.  They do
not exclude alternative cubes for the same period-two measure, different
weights, another period, or another construction family.

Run the eight new cases with:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/lower/round004/compressed_period2_remaining_screen.py
```

The source code SHA-256 at execution was
`891e7158bd24fd720170d544a7b863b7b55499698d0aade21fcf1da84234dbd7`.
The machine-readable results are
`compressed_period2_remaining_results.json`; its SHA-256 is recorded with
the round artefacts.
