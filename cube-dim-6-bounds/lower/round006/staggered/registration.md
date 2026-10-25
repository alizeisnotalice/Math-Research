# Round 006 D — staggered columns discovery registration

Registered before first execution, 2026-10-25. This new mechanism was assigned
after the separate correlated-parity batch completed without promotion.
Only this directory is written; no staging/commit/push/MATLAB/agents by D.

Keep R6-003 breakpoint spacings. Put unit source mass at
`x_i=k_i*a_i`, `x_j=k_j*a_j+delta*a_j*(-1)^ki`, other `x_k=k_k*a_k`.
Period i is `2ai`, all others `ak`, so period mass is 2. Reflection in i
preserves parity, and the four other unshifted axes may be reflected. Axis j
is not reflected. On
`H=[0,ai] x [-aj/2,aj/2] x prod_other[0,ak/2]`, 32 reflection copies tile
the full period cell, giving the candidate ratio `16*union_volume(H)`.

Use contiguous index block lengths 1--4 in every axis. The source mass is
`prod(ell)` and the one common cube radius is `prod(ell)^(1/6)/2`. In axis j,
source minimum/maximum include the minimum/maximum offsets among the selected
i columns. Both special-axis starts are exhaustively enumerated using a
finite radius-based outer range, then actual intersection with H is imposed.
Other starts are the old central witnesses. Coordinate endpoint containment
is checked, and union uses coordinate-slab recursion without overlap addition.

One new-geometry delta-zero calibration is required. First candidate is
ordered pair `(1,2),delta=1/16`. If calibration and pilot complete, use the
remaining combinations of ordered pairs `(1,2),(2,1),(1,3),(3,1)` with
`delta=1/32,1/16,1/8`, excluding the already-run pilot: total 12 candidates.
Raw-box cap 5,000 per candidate, global wall cap 1,200 seconds. Stop on any
calibration mismatch, cap, or candidate above the current complete R6 upper
endpoint (floating trigger guard 1e-10). Promotion exports source-block labels
for E and stops all further scans. No floating result is a C6 bound.

`staggered_results.json` is created before computation and checkpointed after
each candidate, recording PID, input/source hashes, exact threshold strings,
finite scope, counts, times, and actual final process status.
