# Round 006 D: two-axis correlated parity weights

The mechanism is mathematically admissible, but the registered 12-candidate
floating screen produced no promotion. The best screened ratio is
`2.056951307408949`, for axes `(3,4)` and `delta=1/20`. It is approximately
`0.009972518601932` below the current complete R6-003 upper endpoint. This is
numerical discovery evidence, neither a new C6 bound nor a strict exclusion
of the full correlated-parity family.

## Exact geometry and normalization

Let `u=2^(1/6)`, `v=3^(1/6)`, and keep the R6-003 spacings

```
a=(1+u^4/2, (u^3+u^4*v^3)/4, (1+u)/2, (1+u)/2, (1+u)/2, (1+u)/2).
```

At lattice site `(a1*k1,...,a6*k6)`, put the positive weight
`w(k)=1+delta*(-1)^(ki+kj)`, with `0<delta<1`. The period lengths are
`2ai,2aj` in the special axes and `ak` otherwise. The four masses in a
period cell sum to `4`; summing over either special parity alone gives
`(1+delta)+(1-delta)=2`. Thus the average lattice-site mass is exactly 1.

For a Cartesian block of consecutive indices with lengths `ell` and starts
`q`, write

```
S(q,l)=sum_{k=q}^{q+l-1}(-1)^k
      =0                  if l is even,
      =(-1)^q             if l is odd.
M=prod(ell)+delta*S(qi,elli)*S(qj,ellj)*prod(other ell).
```

This follows by factoring the parity term in the finite sum, so it is valid
for blocks crossing period boundaries as well. Use the **one common radius**
`r=M^(1/6)/2` in all six directions. Every center in

```
B(q,ell)=prod_k([(qk+ellk-1)*ak-r, qk*ak+r]) intersect H,
H=[0,ai] x [0,aj] x prod_other[0,ak/2]
```

has a cube containing the entire selected block. Its mass is at least M
and its volume is `(2r)^6=M`; extra atoms can only help. Empty intervals
are discarded. In other axes we use the old central starts
`q=-floor(ell/2)` for odd ell and `q=1-ell/2` for even ell. In the two
special axes every start meeting H is enumerated, with the finite outer
range derived from `r<=rmax` and subsequently tested at its actual radius.

Coordinate sign reflections preserve every source weight because
`(-1)^(-k)=(-1)^k`. They also preserve each block's selected mass and its
common radius. The 64 sign images of H have disjoint interiors and fill
the full period rectangle
`[-ai,ai] x [-aj,aj] x prod_other[-ak/2,ak/2]`. Reflect the witness union
`U_H` with them. These copies have equal measure and their interiors lie
in distinct half cells; hence their union has measure **exactly**
`64|U_H|`. Division by the period's mass 4 gives the witness ratio

```
64|U_H|/4 = 16|U_H|.
```

This equality is for the explicitly reflected witness union. It does not
assert that these block witnesses exhaust the periodic maximal level set.
At level one their averages are at least one; a weak-type strict-level
statement is obtained by taking levels below one and passing to the limit.
No promoted lower endpoint, finite-source truncation, smoothing certificate,
or general-L1 transfer was needed or claimed by this failed numerical batch.

For negative delta, translation by one lattice spacing in either special
axis exchanges `delta` and `-delta`; all maximal-level volumes are equal.
For this particular half-cell witness family the same conclusion holds:
the map `x_i -> ai-x_i` preserves H and sends integer indices `k_i` to
`1-k_i`, changing parity sign. It sends a start `q_i` and length `ell_i`
to `2-q_i-ell_i`, with the same length. Since all intersecting starts in
the special axes are included, this bijects the positive- and negative-
delta witness boxes. No extra negative-delta screen is required.

## Registered finite batch and actual results

Each length lies in `{1,2,3,4}`. The four axis-pair representatives cover
the two exceptional spacings and the four equal tail spacings. Caps were
12 positive-delta candidates, 5,000 raw boxes per candidate, and 1,200
seconds globally. A separate required delta-zero calibration was run once.
The code unions boxes by recursive coordinate slabs, with an interval-union
base case; it never sums overlapping box volumes. Containment pruning uses
actual floating inequalities without tolerance enlargement.

Calibration `(1,2),delta=0` returned `2.0669238260108806`, 347 raw / 219
maximal boxes, within `4.440892098500626e-16` of the current R6 upper endpoint
converted to floating point. All source endpoint/common-radius checks passed.

| axes | delta | raw / maximal boxes | `16|U_H|` (float only) |
|---|---|---:|---:|
| (1,2) | 1/20 | 347 / 219 | 2.0530466606525697 |
| (1,2) | 1/10 | 347 / 219 | 2.037516473171338 |
| (1,2) | 1/5 | 351 / 221 | 2.000587280853158 |
| (1,3) | 1/20 | 347 / 219 | 2.054870841276996 |
| (1,3) | 1/10 | 347 / 219 | 2.041397444079832 |
| (1,3) | 1/5 | 351 / 221 | 2.0014392937895633 |
| (2,3) | 1/20 | 544 / 216 | 2.055155970918871 |
| (2,3) | 1/10 | 544 / 216 | 2.0423420429596844 |
| (2,3) | 1/5 | 544 / 216 | 2.0128509002970767 |
| (3,4) | 1/20 | 544 / 216 | 2.056951307408949 |
| (3,4) | 1/10 | 544 / 216 | 2.0460958824228865 |
| (3,4) | 1/5 | 544 / 216 | 2.013565294643721 |

The numerical comparison used the complete exact upper endpoint stored in
the ledger and frozen phase thresholds, rather than the displayed simple
lower bound `2.0669238260108`. The exact endpoint string is preserved in
`correlated_results.json`. It starts
`2.066923826010880910097778283840054...`. No candidate exceeded it, so
no candidate was sent for rational inward certification and no blind scan
was added. The finite batch does not rule out different spacings, blocks,
more general joint weights, or other witnesses.

`check_mass_symmetries.py` additionally used Python Fraction arithmetic to
check the direct finite mass sum, factorized formula, sign reflection, and
negative-delta half-cell equivalence for all starts -3 through 3, lengths
1 through 4, and three deltas: 2,352 cases passed. This supplements the
algebraic proof; it is not a volume certificate.

## Actual run, hashes and reproduction

The screen ran as PID **12079**, execution session **16336**, and finished
with exit status **0** in **3.558591458 seconds**. Its actual state is
finished; no child process, MATLAB call, or surviving calculation handle
exists. All writes are confined to this owned directory; no files were
staged, committed or pushed. The source was registered before execution.

Input ledger SHA-256:
`096aaff33ea186359046528cfccdde86367610f7fa081d0bcf7d9b273f69ad6b`.
Frozen threshold file SHA-256:
`e0528774b908a5135e24337dacfd6c15c5cc3e156da4578cd694fdece5eecd6d`.
Screen source SHA-256:
`8d827bfb2ce89a01411fed0ae7df8ee63a1041e49b0076190541d5cd532d2033`.
Exact identity-check source SHA-256:
`5a5a650c08237eeb98b18f781ed392e3b30923bc44b2a929e522202fbd4580a0`.

From this directory, reproduce using:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python correlated_screen.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_mass_symmetries.py
```

The screen reads the main ledger at execution time, hashes it and stops
on calibration mismatch, raw-box/time cap, or any promotion trigger. Its
checkpoint/result JSON records the exact input endpoint, source hash,
per-candidate times, counts, domain normalization and explicit no-C6-claim
scope. This report closes the finite batch without changing any bound.
