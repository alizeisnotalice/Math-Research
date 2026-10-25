# Round 007 D: a positive whole-union shape direction

The exact direction batch found a positive right derivative on its third
direction and stopped. Increasing only the second lattice spacing gives

```
312170632105328/10^15 < R'(0+) < 312170632105329/10^15.
```

The single selected step `t=1/1000000` passed an exact endpoint-order
check. The resulting **primary candidate** satisfies

```
206692413818151/10^14 < R_new < 206692413818152/10^14,
R_new approximately 2.066924138181513015... .
```

Its complete exact lower endpoint is strictly above the complete current
R6-003 upper endpoint, by more than `3.12170632105e-7`. The complete
intervals, exact coefficient vector and rational comparison are in
`positive_candidate.json`. This D report does not mark a final new C6
bound: E has received the candidate for independent inward-volume and
finite-L1 certification; main retains final acceptance responsibility.

The candidate source remains the existing rectangular unit lattice.
The new element here is the whole-union first-variation method and its
exact finite-step verification, not a claim of a novel source family or
global optimality. The earlier a2 breakpoint was not claimed to be a
local maximum; the positive right derivative now proves it is not a
local maximum within this explicitly retained witness family.

## Full raw family and exact arithmetic

The active representation axis is 3. Active lengths 1..8 and inactive
lengths 1..4 were registered. Let `a_min=min_j a_j`; exact Q(u,v)
arithmetic proves `a_min>1` and

```
((l-1)*a_min)^6-l*4^5 > 0, l=5,6,7,8.
```

Those active blocks have strictly negative available width at baseline,
so continuity excludes them in a neighborhood of t=0 for every fixed
registered direction. No fifth/seventh-root extension is needed. All
remaining active starts q=-4..6 are included. For ell<=4, r<=2 and
a_min>1; q outside this range has a strict gap from the receiver interval
[0,2a_i], and cannot enter for sufficiently small t.

The baseline builder discards only strictly negative-width boxes. It
retains exactly **382 raw boxes** in the full-active/half-inactive cell,
without containment pruning. In this particular instance all 382 have
positive baseline widths; endpoint contacts between different boxes still
matter. The algorithm itself preserves any zero-baseline/positive-linear
width, as checked by the analytic toy cases. No dynamically chosen
96-box representation is differentiated.

Endpoint jets `(F0,F1)` lie in the original 36-coefficient field
representation with u^6=2,v^6=3. The union routine sorts coordinate
endpoints lexicographically, represents interval positions by integer
ranks, and recursively partitions complete slabs. Slab widths and child
volumes are multiplied modulo t^2. Overlapping boxes are never summed
as independent volumes.

`theory.md` proves correctness when distinct true endpoint jets separate
to first order. It then proves a finite-union Lipschitz extension: replacing
all C1-right endpoints by their affine jets changes union volume by o(t).
Thus equal jets with different higher-order germs may be merged for the
first-order **volume** calculation, while their true endpoint order and
finite-step behavior remain unresolved. The code records such unresolved
germ groups instead of calling the true endpoints identical. E separately
reviewed this interface and the source. In the three directions actually
evaluated, the endpoint functions are affine and the unresolved-group
count is zero, so the positive result does not rely on a nonlinear tie.

The four registered self-checks all passed:

| Check | Constant volume | Right derivative |
|---|---:|---:|
| Separating interval contact | 2 | -1 |
| Overlapping interval contact | 2 | 0 |
| One newborn thin direction | 0 | 6 |
| Two newborn thin directions | 0 | 0 |

## Three actual directions, then the required stop

Spacing directions mean `b_j(t)=a_j*(1+t*h_j)`; all masses stay one.
Each exact constant volume coefficient was identical, coefficient by
coefficient, to the old R6-003 certificate.

| Direction | Right derivative, float display | Certified sign | Union states |
|---|---:|---:|---:|
| h1=+1 | -0.1811256459374387 | negative | 3643 |
| h1=-1 | -0.11434580137695732 | negative | 3726 |
| h2=+1 | 0.3121706321053281 | positive | 3571 |

The remaining 13 registered directions were not run. This is neither a
complete tangent-fan optimization nor an audit of all local directions.
It is a rigorous positive direction sufficient to select one new input.

## The single finite input delivered to E

For positive u^6=2,v^6=3, the new spacings are

```
a1=1+u^4/2,
a2=(1000001/1000000)*(u^3+u^4*v^3)/4,
a3=a4=a5=a6=(1+u)/2.
```

All source masses remain one, and each witness uses the same radius
`(prod ell)^(1/6)/2` in every coordinate. This is within the fixed
Lebesgue/centered/axis-parallel/common-radius model.

Only coordinate 2 moves, so all its unclipped endpoint functions and
clipping faces are affine in t, while all other coordinates stay fixed.
`confirm_step.py` collects **57 distinct affine endpoints** for that
coordinate from all 382 retained raw boxes and checks their **56 adjacent
lexicographic orders** at exactly t=1/1000000. Every separation has a
strictly positive rational lower bound. Linear interpolation then
preserves the entire order for all 0<t<=1/1000000. This includes branch
changes in clipping and inter-box endpoint contacts.

Consequently this retained witness union has exactly affine volume on
that interval, and `R_new=R_old+R'(0+)/1000000`. Its field coefficients
and enclosing rational interval are obtained by that identity, with no
additional union call. Strictly empty baseline boxes that might enter
later are not needed for this lower witness. No step-size grid was run,
and no second step was selected.

E may independently use the standard all-ell<=4 positive-half-cell
family at these spacings; it represents the same unit lattice and can
be checked without using the jet implementation. Finite truncation,
smoothing, strict threshold choice and the final L1 lower endpoint are
E/main acceptance tasks, not inferred from the derivative alone.

## Execution and provenance

The main direction batch ran as PID **17473**, exec session **73571**,
and exited **0**. Total time was **5.27727879100712 seconds**, with
**10946 union states** including four self-checks. Limits were 16
directions, 1,000,000 states and 1200 seconds; the positive-direction
stop was reached, not a resource limit.

The one-step check ran as PID **17783**, exited **0**, and took
**0.7826590829936322 seconds**, with **zero additional union calls**.
Both computations have ended. No MATLAB, subagents, commits or pushes
were used, and round006 files were not modified.

The input ledger hash is
`096aaff33ea186359046528cfccdde86367610f7fa081d0bcf7d9b273f69ad6b`,
and the full frozen-threshold hash is
`e0528774b908a5135e24337dacfd6c15c5cc3e156da4578cd694fdece5eecd6d`.
The executed jet source hash is
`4c559f0d5e4610d9417139fda70253d76d4b0827746e10f47587587e764a40f4`.
Both output JSON files record their remaining input/source hashes,
including the exact baseline certificate and old arithmetic layer.

Reproduce from this directory with the existing local Python:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B shape_jets.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B confirm_step.py
```

The archival payload comprises this report, two registrations, `theory.md`,
the two source scripts and their two result JSON files. Their full exact
integers and coefficients are essential for review; no intermediate grid,
Monte Carlo samples or unbounded process is part of the payload.
