# A038: three-point diagnostic of the Green floor 1/65000

The three authorized new full G+six-gradient calls completed and stopped
at the call limit. All three values strictly decrease from the saved B032
seed, but all three are strictly above1/65000. This batch therefore does
not refute that conditional global floor, and does not prove it either.
The actual I/g model remains unresolved; C6 is unchanged at64. No H, K,
p, spatial grid, optimizer or MATLAB computation was run.

## Exact inputs and results

The B032 seed was read from its complete saved native endpoints and was
never recomputed. Its report SHA256 was checked against the supplied
`d7b7c90e8f4ab09bae84caece1b92529ce793067d09556bd9f311a843ff8b723`.
The seed is

    (100000000000000,99999998999999,99999997999988,
     96002619105455,55264180466436,95885707114)/10^14.

All three new points have denominator10^14 and the same first three
numerators100000000000000,99999998999999,99999997999988.

| New point | Fourth numerator | Fifth numerator | Sixth numerator | G, display only |
|---|---:|---:|---:|---:|
| 1 | 99232423642001 | 55975441278582 | 74187318416 | 1.6155261376420057e-5 |
| 2 | 99999996999867 | 56645900209181 | 57057470436 | 1.6019870773260950e-5 |
| 3 | 99999996999867 | 57255226550170 | 43659674865 | 1.6005426481742880e-5 |

The proof uses the original Arb512 dyadic endpoints, not those decimal
displays. Independent exact postprocessing gives the following rational
outer intervals (strict inequalities at both ends):

| Point | Lower rational bound | Upper rational bound |
|---|---|---|
| 1 | 323105227528401/20000000000000000000 | 807763068821003/50000000000000000000 |
| 2 | 320397415465219/20000000000000000000 | 100124192332881/6250000000000000000 |
| 3 | 1600542648174287/100000000000000000000 | 100033915510893/6250000000000000000 |

Each displayed lower rational is strictly greater than1/65000. The full
saved endpoint comparisons establish

    1/65000 < G(point3) < G(point2) < G(point1) < G(B032 seed).

No fourth point was generated or evaluated. The first-counterexample stop
condition never fired; the three-call cap did.

## Fixed rule and boundary handling

Every proposal uses the exact midpoints of all six preceding saved
gradient intervals; the five free components determine the negative
direction normalized by their maximum absolute value. The step is always
1/20. The frozen rule then clips sequentially to the ordered face with
coordinate gaps at least1e-8, rounds exactly to nearest-even on the10^14
lattice, and applies only the previously registered gap repair if needed.
Every candidate file records the preceding complete native file path/hash,
gradient endpoints' exact midpoints, direction, raw proposal, clipping,
quantization error and repair.

At point1, clipping affected x2 and x3. At points2 and3 it affected x2,x3,
and x4. The unmodified quantized candidates had respectively1,3,3 zero
Fourier forms. No kernel was called on those unsupported inputs. In all
three cases the first prescribed repair k=1 passed: add11^j lattice units
to every gap except the largest, and subtract their sum from that largest
gap. All final inputs remain strictly ordered and inside the cube. The
boundary approach therefore caused clipping but no interface failure;
there was no new zero-node rule, changed step or alternate seed.

The full primitive signed/permuted Fourier-form set has9308 forms. Each
final input has zero zero-forms and minimum absolute integer form value10.
The conservative raw node numerator bound is1.2e15 and denominator bound
1.6e16, both strictly within int64. Initial and repaired preflights are
retained. Exact duplicate checks used159 prior Green inputs and then each
new point successively; no previously certified point was recomputed.

## Kernel identity, saved endpoints and checks

This batch reused the frozen C013 generic nonzero-node body32 Green
engine. It did not implement, compile or alter a Green kernel. The locked
engine source, binary and fan hashes are:

| Object | SHA256 |
|---|---|
| C013 generic_green.cpp | 516d69876325810868728ff450db7a8844dcd6f7258223e6495db819b9bda7e2 |
| C013 native binary | 9e33f438630d178af864a4720572dee403548a79656b6a7978db3fe0e0ed4a88 |
| Original body32 fan | 7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812 |
| Batch registration | e31eb3743040c5a53c1ae0698d7f2790bb189d3e13371e9a9354b4b68f8c3246 |

Fifty inherited/current source paths, including the B032 frozen records,
engine and runtime libraries, were hash-checked before execution and
again in postchecking. The frozen B032 whitelist was verified in full.
The reused rule's G0 comparison constant is set to1/65000 only in the new
driver's memory; no old source file is written.

Each native output reports967680 full logical terms,483840 computed sign
pairs,9780 unique node values, zero_node_terms=0, confluent_terms=0 and
precision_bits=512. Full lower/upper mantissas and common binary exponents
for G and all six derivatives are retained in `native_outputs/01..03.stdout.json`.
All seven interval widths pass the2^-256 gate. The native Euler interval
and an independent reconstruction of x dot grad G +5G both contain0.
Euler is an identity check, not a replacement for the previously reviewed
full-group formula and its gradient chain rule.

`check_saved.py` independently reconstructs exact dyadic values, the
normalized direction, clipping, nearest-even rounding and k1 gap repair;
it checks each saved input, all9308 forms, hashes, int64 bounds, full term
counts, Euler, strict descent and strict survival above1/65000. It never
invokes a native kernel. Its first run passed in0.048629 seconds.

## Essential-floor interpretation and remaining gap

The three exact inputs lie off the genuine body32 S2: the saved rigorous
lower bounds on their l-infinity distance to S2 are respectively

    22012279360709/100000000000000,
    43354099790819/200000000000000,
    4274477344983/20000000000000.

These follow because a body32 pair span has at most two absolute
coordinate values; among any three distinct absolute coordinates, one
coordinate approximation error must reach half the smallest pair gap.
Taking the maximum over triples gives the saved positive bounds.

Green is continuous off S2; the nonzero-node local finite formula also
supplies continuity at these points. Thus these data further lower the
known pointwise upper obstruction on the essential cube infimum:

    g <= G(point3) < 100033915510893/6250000000000000000.

That upper obstruction remains strictly above1/65000. It supplies no
lower certificate for g, no face integral I bound, and no verification of
the proposed Hnorm<=1/14000 condition. The finite descent has not located
or certified a local or global minimum. In particular, clipping several
coordinates does not justify concluding that no better point exists.

## Actual resources, files and freeze

Started UTC2026-10-25 20:50:39 with a600-second total cap. The driver
PID44671 completed all three calls in3.760972 seconds. Native PIDs44679,
44680,44681 returned exit0, with wall times1.229932,1.233048,1.243438 seconds.
The exec session63905 closed, with launch chunkc5e0de and final chunk22670a.
The checker PID44842 returned exit0 (chunkda0ff2). All child handles are
closed. There were no numerical failures, retries, extra points or extra
kernel calls. One initial read-only directory lookup found no C013 copy in
the upper-astra clone; the already recorded audit-clone path was then used.

The necessary package comprises registration.md/json, run_batch.py,
candidate_01..03.json, three raw native outputs/stderr files, run_receipt.json,
results.json, check_saved.py, checked_results.json, this report and the
final whitelist manifest. Original native outputs are retained because
their complete endpoints are irreducible evidence for the proposal chain.

Recheck saved evidence without any kernel call:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round038/green_floor_diagnostic/check_saved.py
```

The driver refuses a second batch when its run receipt exists. No further
candidate or changed target is authorized or launched here. The parent
retains final proof review and integration responsibility.
