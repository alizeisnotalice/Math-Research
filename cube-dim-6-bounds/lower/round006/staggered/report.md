# Round 006 D: staggered unit-mass columns

The new finite batch is complete with no promotion. The best screened ratio
is `2.037927812469202`, at ordered axes `(i,j)=(1,3)`, `delta=1/32`.
It is approximately `0.028996013541679` below the current complete R6-003
upper endpoint. These floating calculations produce no new C6 bound and do
not exclude the full staggered-lattice family.

## Geometry and normalization proved for this family

Keep the exact R6-003 spacings
`a=(1+u^4/2,(u^3+u^4*v^3)/4,(1+u)/2,...,(1+u)/2)`, with positive
`u^6=2,v^6=3`. Put unit atoms at

```
x_i=k_i*a_i,
x_j=a_j*(k_j+delta*(-1)^k_i),
x_k=k_k*a_k     (k other than i,j).
```

The period in axis i is `2ai`, and all other periods are `ak`. A half-open
period cell contains two atoms, counting the two i parities, so its mass
is exactly 2. The phase shift in j does not change this mass.

A contiguous **index** block with starts q and lengths ell has exactly
`M=prod(ell)` atoms. It is not generally a Cartesian product of spatial
coordinates, but it has explicit coordinate extrema. The extrema in i and
the four unshifted axes are `qk*ak` and `(qk+ellk-1)*ak`. Define the offset
set `D={delta*aj*(-1)^ki: qi<=ki<qi+elli}`. The j extrema are

```
L_j=qj*aj+min(D),
U_j=(qj+ellj-1)*aj+max(D).
```

For `elli=1`, D has one offset; for `elli>=2`, it contains both offsets.
Taking the one common cube radius `r=M^(1/6)/2`, the center box is the
product of intervals `[U_k-r,L_k+r]`, clipped to

```
H=[0,ai] x [-aj/2,aj/2] x prod_other[0,ak/2].
```

Every center in that box is within r of every selected atom in every
coordinate, so its cube contains all M selected unit masses. Its volume is
`(2r)^6=M`. This verifies the common-radius witness rule directly.

Reflecting x_i sends k_i to -k_i, preserving parity and therefore the
j offset. Reflecting any of the four unshifted coordinates also preserves
the source. **No j reflection is used for positive delta.** The 32 sign
images in precisely these five axes give half cells with disjoint interiors
that fill one full period rectangle. Their reflected witness unions have
equal measures and disjoint interiors; hence the reflected union measure
is exactly `32|U_H|`, and its periodic atomic witness ratio is

```
32|U_H| / 2 = 16|U_H|.
```

This is a witness subfamily ratio, not a claim of exhaustive maximal-level
coverage. Average at least one gives weak-type lower witnesses by taking
levels increasing to one. No new lower endpoint was promoted, so this batch
does not claim finite-source truncation, smoothing or general-L1 certification.

The negative-delta source is the j reflection of the positive-delta source;
the symmetric full j interval of H is preserved. It therefore has the same
full level-set volume and the same witness-family union volume.

## Finite range and actual results

Lengths in every axis are 1--4. All starts in i and j that can meet H are
enumerated inside a radius-derived finite outer range and then filtered by
actual interval intersection. The four other starts are the old central
witnesses. The j outer bound is deliberately enlarged by `abs(delta)*aj`
and half a j period; it covers every actual shifted intersection. Every
retained box passes source-extrema/common-radius containment checks.

There was one new-geometry calibration: `(1,2),delta=0` returned
`2.0669238260108806`, with 348 raw / 218 maximal boxes. Its floating error
against the complete current endpoint is `-4.440892098500626e-16`.
The first actual candidate `(1,2),delta=1/16` was run once; the remaining
11 combinations were then run without repeating it.

| ordered axes (i,j) | delta | raw / maximal boxes | `16|U_H|`, float only |
|---|---|---:|---:|
| (1,2) | 1/16 (pilot) | 307 / 205 | 2.004458240470025 |
| (1,2) | 1/32 | 349 / 217 | 2.031549368285088 |
| (1,2) | 1/8 | 279 / 201 | 1.953005585237 |
| (2,1) | 1/32 | 486 / 213 | 1.9590027051470782 |
| (2,1) | 1/16 | 486 / 209 | 1.8276055469822605 |
| (2,1) | 1/8 | 329 / 189 | 1.7082232120323528 |
| (1,3) | 1/32 | 307 / 205 | 2.037927812469202 |
| (1,3) | 1/16 | 307 / 205 | 2.0095155230914616 |
| (1,3) | 1/8 | 269 / 193 | 1.958635746112176 |
| (3,1) | 1/32 | 486 / 213 | 1.9608412496403125 |
| (3,1) | 1/16 | 486 / 209 | 1.8288839799689787 |
| (3,1) | 1/8 | 329 / 189 | 1.7059844694278803 |

Union volume is computed using recursive coordinate slabs and merged
one-dimensional intervals. Overlapping boxes are never summed directly.
Containment pruning uses actual floating inclusion, with no tolerance
extension. The comparison endpoint is the exact complete current R6-003
upper endpoint, preserved in the result JSON, rather than the displayed
simple lower bound. All 12 candidates were below it. No candidate was sent
to E for inward certification and no additional scan was started.

The additional exact Python Fraction check in `check_geometry.py` compares
direct atom lists with the derived extrema and reflection formulas for
2,352 block cases (starts -3 through 3, lengths 1 through 4, three deltas).
It checks the source geometry identities, not numerical union volumes.

## Execution and reproducibility

Registered caps: 12 candidates, 5,000 raw boxes per candidate, 1,200 seconds
globally, with stops on calibration mismatch, resource cap, or numerical
promotion. Maximum actual raw count was 486. The screen ran as PID **13029**,
execution session **6364**, and finished with exit status **0** after
**1.710256583 seconds**. No calculation process or MATLAB call survives.
Only the assigned staggered directory was written; no staging, commit or push.

Input ledger SHA-256:
`096aaff33ea186359046528cfccdde86367610f7fa081d0bcf7d9b273f69ad6b`.
Frozen thresholds SHA-256:
`e0528774b908a5135e24337dacfd6c15c5cc3e156da4578cd694fdece5eecd6d`.
Screen source SHA-256:
`3baa94ca9e194332d8b5223947f71d14c3b8e6b461830dfc263963053e1fd9d5`.

From this directory:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python staggered_screen.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_geometry.py
```

The script records input/source hashes, exact comparator strings, PID,
counts, times and all candidate results in `staggered_results.json`, created
before computation and checkpointed after every candidate. This closes only
the registered finite batch; no bound or complete-family exclusion results.
