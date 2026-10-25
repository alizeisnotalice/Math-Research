# B045: restricted normalization descends once, then the fixed step returns

The authorized batch stopped after two new complete G+six-gradient calls.
The first point is strictly below the saved A038 seed. At the second
point G is strictly above the first, so the registered no-descent stop
applied immediately. No third proposal, smaller step, backtracking,
replacement seed or further kernel was generated.

Neither point refutes g>=1/65000. This is a fixed-step diagnostic,
not a restricted-face minimum, a global minimum, a verified global
floor, or an update of C6. Public and internal upper bounds remain 64.

## Registered inputs and the new rule

The read-only primary round045 registration/results identified that
normalizing by coordinates2-6 suppresses the feasible last-two movement.
The saved exact ratio of last-two to old normalization is approximately
0.10659867078859825. This motivated the newly authorized rule, rather
than repeating the old normalization. This batch does not treat that
display number as a proof of improvement or convergence.

The only initial seed is A038 point3, with denominator10^14 and numerators

    (100000000000000, 99999998999999, 99999997999988,
     99999996999867, 57255226550170, 43659674865).

Its full native G/gradient output SHA256 is
`7bd6d46cac3cf1d7daf4968ebd012e1108a5f558f9c51eb51f262cdf6fcd4c1f`.
It was read and validated, never recomputed. The first four numerators
remain exactly fixed throughout; they are not replaced by four equal
coordinates. The new preparation locks73 external sources, including
the primary diagnostic, all A038 frozen files, the old engine/runtime
chain, and the historical point records.

The no-repeat set contains the inherited159 pre-A038 exact keys plus the
three A038 points. Conservative supplemental keys from all five A008
diagnostic/exact points and the symbolic x0 give167 unique keys after
deduplication. Old diagnostic values are used only to forbid repeats,
never as numerical certificates. Each new input is checked against
this set and earlier inputs in the current batch.

For each saved full gradient, take its exact dyadic midpoint g_j and
M=max(|g_5|,|g_6|). The proposal has direction zero in coordinates1-4
and -g_j/M in coordinates5,6, with step1/20. If M=0 the batch stops.
For gamma=10^-8, the exact ordered clipping rule is

    x5'=max(2gamma,min(x4-gamma,raw5)),
    x6'=max(gamma,min(x5'-gamma,raw6)).

Nearest-even rational rounding uses denominator10^14. Full preflight
checks all9308 exact primitive Fourier zero forms, ordered feasibility
and conservative int64 numerator/denominator bounds. On a zero form
only, the registered repair tries adding k to the sixth numerator,
k=1,...,16, always based on the same quantized input; the first fully
feasible generic point is used. No other coordinate can change in a
repair. Both actual proposals were generic without any repair.

Registration SHA256:
`c679a535638908452ac2f7701ab35151fa035451ff9beba39d863fcba70f84fe`.
The batch has a600-second cap, at most four complete native calls,
and120 seconds per call. Before each call, the exact candidate and its
hash are written to the run receipt. The first strict target rejection,
failure to certify strict descent, or unsupported/failed computation
stops the batch; there are no retry or backtracking rules.

## The two exact points and certified comparisons

Both new points have the same first four numerators as the seed and
the same denominator10^14. Their last two numerator pairs are

    point1: (62255226550170,1000000),
    point2: (57255226550170,1000000).

Thus x6=10^-8 at both points. Point2 shares the seed's x5 but not its x6;
it is a new exact input, not a repeated seed evaluation. Both points
passed exact no-repeat, range and generic checks. Their minimum
absolute primitive zero-form numerator is1, so their near-wall nature
was retained exactly rather than handled with a floating tolerance.

The complete raw Arb512 mantissa/exponent endpoints are retained in
`native_outputs/01.stdout.json` and `02.stdout.json`. The following
strictly exterior rational bounds were reconstructed from them:

    1596157362896397/10^20 < G(point1)
                              < 798078681448199/(5*10^19),
    1600541724642269/10^20 < G(point2)
                              < 160054172464227/10^19.

All lower endpoints are strictly above1/65000. Displays only:
G(point1)=1.5961573628963972e-5 and
G(point2)=1.60054172464227e-5. Exact dyadic differences establish

    G(point1)_upper < G(A038 point3)_lower,
    G(point2)_lower > G(point1)_upper.

The second inequality is stronger than merely failing to certify a
decrease. It triggers `NO_CERTIFIED_STRICT_DESCENT_IMMEDIATE_STOP`.
The exact difference intervals are in `results.json` and independently
reconstructed in `checked_results.json`.

The first new gradient has positive coordinate5 derivative, approximately
5.211354122597271e-7. The registered negative-gradient rule consequently
moves x5 back by exactly1/20 on step2; its new derivative is negative,
approximately -2.2125198716017544e-6. The sixth derivatives are positive
at both points (about9.45e-13 and9.69e-13), and the prescribed lower
coordinate bound clips the negative sixth movement. These statements
come from saved full intervals; no new scalar kernel was computed.
They explain the fixed-step return but do not certify the location,
value, uniqueness or existence of a two-coordinate minimum. This batch
does not select a new step or continue a local search.

## Full-group and essential-infimum scope

Both calls use the original body32 C013 engine, without source changes:

| Input | SHA256 |
|---|---|
| C013 source | 516d69876325810868728ff450db7a8844dcd6f7258223e6495db819b9bda7e2 |
| C013 binary | 9e33f438630d178af864a4720572dee403548a79656b6a7978db3fe0e0ed4a88 |
| Body32 fan | 7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812 |
| New driver | 560ad259877eb23476263201d86f2a069073d4a0d8f2cffd5dadde3e3b66d8f8 |

Runtime FLINT/GMP/MPFR library paths and hashes are inherited explicitly
in registration and checked before launch and in saved-data verification.
Both outputs report967680 full logical terms,483840 sign-paired terms,
512-bit precision, zero_node_terms=0 and confluent_terms=0. No zero or
repeated node is silently skipped. Each G and derivative enclosure passes
the2^-256 width gate. Both the native Euler interval and an independent
exact reconstruction of x.gradG+5G contain zero. These are integrity
checks in addition to, not replacements for, the previously reviewed
C013 full-group formula and derivative chain rule.

The exact lower bounds on distance to the true body32 S2 in l-infinity
norm are respectively

    3774477344983/20000000000000,
    4274477344983/20000000000000.

The checker derives these by taking half the smallest pairwise gap in
each three-coordinate selection, then the maximum over such selections.
A point of S2 has at most two distinct absolute coordinate values;
approximating three separated nonnegative values with two levels forces
at least one error of half their smallest separation. Thus both bounds
are rigorous and strictly positive. Green is continuous off S2, so the
new first point gives a valid pointwise upper obstruction on the
essential cube infimum:

    essinf_cube G <= G(point1) < 798078681448199/(5*10^19).

That upper obstruction is still above the proposed floor1/65000. It is
not a lower bound on G, not an integral I certificate, and not a proof
of either conditional premise needed for a new weak-type constant.

## Independent saved-data check and actual resources

`check_saved.py` does not import the new driver or run any kernel. It
independently rebuilds dyadic midpoints, last-two normalization, exact
clipping, nearest-even rounding, every allowed repair if present,
preflight, all six derivative endpoints, Euler, true S2 separation,
both strict comparisons and the first-stop logic. The frozen preflight
helper supplies the full9308-form orbit; its source and fan are locked.
The checker verifies that no third proposal/output exists.

The sole driver PID48625 ran for2.501291 seconds. Native PID48633 exited0
after1.234056 seconds; native PID48636 exited0 after1.230259 seconds.
Tool session50048 closed with exit0, launch chunkcd8280 and final
chunk7de585. There were no native failures, retries, seed calls, repairs,
extra points or orphan handles. The independent pure checker PID48715
passed in0.038634 seconds, exit0, tool chunkddad9d.

Only the new owned directory was written. There were zero H/K/p calls,
zero grid/optimizer/MC/MATLAB calls, no old-source edits, no agent
expansion, no commit and no push. All exact proposals and raw native
outputs are preserved. The finite batch is complete and stopped; the
parent retains final proof integration.

Read-only reproduction without a kernel:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round045/restricted_face_diagnostic/check_saved.py
```

The driver refuses a second run once its receipt exists. No continuation
is implied by the unused two-call budget.
