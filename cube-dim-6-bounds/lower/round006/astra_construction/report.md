# Round 006 D Astra: joint gap/mass construction

The registered batch completed with **no new lower-bound candidate**.
Its best floating witness ratio was `2.0073678363809004`, at active axis 3,
`epsilon=1/8`, `delta=0`. It lies approximately `0.059555989629980655`
below the complete R6-003 upper endpoint, whose decimal expansion begins
`2.066923826010880910097778283840054...`. The exact rational comparator,
not the shortened display lower bound, is stored in `dimer_results.json`.
No result was sent to E for certification and no C6 bound was changed.

## New mechanism and proved geometry

Let `a` be the R6-003 spacing vector, with positive `u^6=2,v^6=3`:

```
a=(1+u^4/2, (u^3+u^4*v^3)/4, (1+u)/2, (1+u)/2, (1+u)/2, (1+u)/2).
```

Select an active axis i. For `0<=e<1`, set `s=(1+e)^(1/6)`,
`b_i=a_i*s/(1+e)` and `b_j=a_j*s` for the other axes. Along axis i use
`p(k)=b_i*(k+e*(k mod 2))` with weight `1+d*(-1)^k`, `0<=d<1`;
the five remaining coordinates form the product lattice with spacings b.
The source period has two atoms in its active motif, total mass 2, and
volume `2*prod(b)=2*prod(a)`. Thus the density is constant throughout the
path, even when the two weights differ.

For an active consecutive block with start q and length l, its span is

```
b_i*(l-1)                              if l is odd,
b_i*(l-1+e*(-1)^q)                     if l is even,
```

and its mass factor is `l+d*(-1)^q` for odd l and l for even l.
Hence geometry acts on even blocks while mass acts on odd blocks. This
joint mechanism is absent from the earlier fixed-spacing parity screens;
it also differs from the earlier cross-axis shifted-column source.

For each block, multiply its active mass factor by the five inactive
lengths to obtain M. The **same radius in all axes** is `r=M^(1/6)/2`.
The center witness box has coordinate intervals
`[max_source_coordinate-r,min_source_coordinate+r]`, intersected with

```
H=[0,2b_i] x product_{j!=i}[0,b_j/2].
```

Each center cube then contains the entire selected block and has volume M.
All active starts whose boxes meet H are enumerated; in inactive axes the
old central odd/even starts are used. Reflection in the five inactive
coordinates preserves the source and yields 32 copies with disjoint
interiors in a full period cell. The resulting witness ratio is
`32|U_H|/2=16|U_H|`. No active-axis reflection is assumed for the computation.

For completeness, reflection about the midpoint of the first active gap,
`x -> b_i*(1+e)-x`, maps index k to `1-k` and exchanges the two weights.
Thus changing the sign of d gives the same full-period witness volume.
Reflection about zero exchanges the two signed gap displacements at fixed
b. These are source symmetries, not a claim of optimality.

At e=0,d=0 this is the original lattice. At e=1, adjacent opposite-parity
atoms coalesce: an even active block of length 2m starting at odd q
contains m double atoms, each of total mass 2 independently of d.
Its effective active spacing is `2b_i=u*a_i`, and every inactive spacing
is `u*a_j`. Global spatial dilation by u and atom mass scaling by 2
recover exactly the original witness ratio. The initial and coalesced
endpoints therefore give the same baseline; the batch tested nontrivial
interiors of this path.

## Actual finite results

The 24 cases were the product of active axes 1,2,3, `e` in
`{1/8,1/3,2/3,7/8}`, and `d` in `{0,1/5}`. Axis 3 represents the four
permutation-equivalent tail spacings. Active lengths were 1..8 and
inactive lengths 1..4. Active length 8 includes four pairs near the
coalesced endpoint; these are not repetitions of the old fixed-spacing
length-5/6 extension calculation.

| Active axis | e | d=0 ratio | d=1/5 ratio |
|---|---:|---:|---:|
| 1 | 1/8 | 1.9681855036344884 | 1.8824848971936943 |
| 1 | 1/3 | 1.6949651304096690 | 1.6923104507546871 |
| 1 | 2/3 | 1.6886250515019148 | 1.6987230804632490 |
| 1 | 7/8 | 1.8423496296596888 | 1.8511562714210745 |
| 2 | 1/8 | 2.0045304784553480 | 1.9760334130603840 |
| 2 | 1/3 | 1.9609484553342167 | 1.9664237751426680 |
| 2 | 2/3 | 1.9746563359850722 | 1.9771482588291970 |
| 2 | 7/8 | 1.9997647549749630 | 2.0014177731921620 |
| 3 | 1/8 | 2.0073678363809004 | 1.9716647570926247 |
| 3 | 1/3 | 1.9468841659259386 | 1.9518128955497910 |
| 3 | 2/3 | 1.9682486135385802 | 1.9700867130921484 |
| 3 | 7/8 | 1.9954678492172837 | 1.9969434078589519 |

The independent prescribed endpoint check `(i,e,d)=(3,1,1/5)` returned
`2.0669238260108806`; the initial `(1,0,0)` calibration gave the same value.
Each differs from the complete old endpoint converted to float by
`-4.440892098500626e-16`. This checks the normalization and paired endpoint
implementation; it is not a new candidate or an exact volume certificate.

The geometry verifier directly enumerated 1,056 rational motif blocks,
checking mass, extrema, span and both reflection identities with Python
`Fraction`, and 24 exact double-atom cases. All passed. Every retained
floating box also passed its common-radius/source-endpoint checks.
The union routine recursively partitions into slabs and merges intervals;
it does not sum overlapping raw box volumes. Containment pruning uses no
tolerance enlargement.

## What this failure establishes, and the concrete remaining obstacle

Only the 24 listed numerical values for this precise witness family have
been tested. The batch does not exclude other e,d, arbitrary six spacings,
larger block families, non-product motifs, or the full dimer family.
In particular, no exact upper bound on this family was computed.

Two structural facts explain why further grid refinement along this path
is not presently justified. First, for any fixed odd l and inactive mass P,
the two radii satisfy, strictly for `d!=0`,

```
(P*(l+d))^(1/6)/2 + (P*(l-d))^(1/6)/2 < (P*l)^(1/6).
```

This is strict concavity of the sixth root. It is a genuine penalty in the
sum of corresponding radii, though overlap redistribution prevents it
from proving a union-volume bound. Second, the chosen constant-density
path forces all five inactive spacings to expand by s, so the gain from
short active gaps competes against lost overlap in five directions.
The data show a valley between two equivalent endpoints rather than a
numerically promising shoulder. This interpretation is not a proof of
local or global maximality.

The next mathematically specific step, if this mechanism is continued,
is to free the spacing tangent instead of adding another e,d grid. Write
`b_j(t)=a_j*(1+t*h_j)`, `e(t)=t*E`, `d(t)=t*D`, initially with
`sum h_j=0` if constant density is retained. The tested path uses only
`h_i=-5E/6` and `h_j=E/6` for j other than i. For a raw witness at t=0:

* `r'(0)=r(0)*D*(-1)^q/(6*l)` for odd active l, and zero for even l;
* `p_q'(0)=a_i*(q*h_i+E*(q mod 2))` in the active coordinate;
* inactive source endpoints have derivative `a_j*q_j*h_j`;
* clipped face derivatives use min/max of these derivatives at tied faces.

Thus the one-sided first variation of the **whole union**, including
zero-width newborn slabs and endpoint ties, is piecewise linear on the
finite endpoint-order tangent fan. An increasing feasible tangent would
select a defensible new candidate and scale; nonpositive maxima on that
fan would give a genuine local obstruction for the specified raw family.
Neither conclusion has yet been computed. Handling these tied endpoints,
rather than merely differentiating the currently retained boxes, is the
precise unresolved step. It changes the next action from a blind interior
grid to an event/shape-derivative calculation.

## Execution, limits and reproducibility

The registered 24-case screen used PID **14829**, exec session **62042**,
and exited **0** after **9.760100999992574 seconds**. It visited **666451**
union states, including calibration, and at most **511** raw boxes in a
candidate. Limits were 24 candidates, 5000 raw boxes per candidate,
1,000,000 union states and 1200 seconds globally. No limit was reached.
The separate identity/endpoint verifier used PID **15076**, exited **0**,
and took **0.1263260419946164 seconds**. All processes have finished;
there was no MATLAB call and no child agent or calculation process.

Relevant SHA-256 values:

* Main input ledger: `096aaff33ea186359046528cfccdde86367610f7fa081d0bcf7d9b273f69ad6b`.
* Frozen thresholds: `e0528774b908a5135e24337dacfd6c15c5cc3e156da4578cd694fdece5eecd6d`.
* Registered screen source: `9b8cdc0cc12c99a6ee4645a44274b5495dc254cab2542a2b5049886559d29b1c`.
* Screen result: `180dcb5d4b4d61657b0147d4c2efda3e74dce28fc67fe6a5db9dd47e25f04110`.
* Geometry source: `859935904c629b9d20576d0edcb8eee9ba036245c17a5a14adf1ca6c83ba8e0c`.

The result JSON also records the protocol, roles, fifth-round extension
certificate, and both old-D report hashes. From this directory reproduce:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B dimer_screen.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B verify_geometry.py
```

Only this owned directory was written. There were no changes to the old D
outputs, main ledger, other projects, Git index, commits, or remote state.
