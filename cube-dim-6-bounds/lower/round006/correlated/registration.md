# Round 006 D — two-axis correlated parity finite discovery batch

Registered before first execution, 2026-10-25. Owner D writes only this directory;
the main ledger and thresholds are read-only. No commit, push, staging, MATLAB,
or agents are used by D.

The fixed R6-003 breakpoint spacings are used. Weights are
`w(k)=1+delta*(-1)^(ki+kj)` on the lattice, where `delta in {1/20,1/10,1/5}`.
The 12 finite candidates use representatives `(i,j)=(1,2),(1,3),(2,3),(3,4)`;
axes 3--6 have equal spacings. A single required `delta=0,(1,2)` calibration
precedes these candidates. The first candidate is `(1,2),delta=1/20`.

Lengths of every contiguous block are in `{1,2,3,4}`. Starts in the two special
axes are exhaustively enumerated within the radius-bound intersection range;
other axes retain the old central witnesses. The selected mass is
`prod(ell)+delta*S(si,elli)*S(sj,ellj)*prod(other ell)`, where `S` is 0 for
even length and `(-1)^start` for odd length. Each box uses the single common
radius `mass^(1/6)/2`; source endpoint containment is checked in all six axes.

The half cell is `H=[0,ai] x [0,aj] x prod(other [0,ak/2])`. Coordinate
reflections preserve the weights. The full period cell has mass 4, so the
reflected union ratio is `64|U_H|/4=16|U_H|`; the report supplies the proof.
Union volume uses a recursive slab union, never sum of overlapping boxes.

Global wall budget 1,200 seconds, raw box cap 5,000 per candidate, 12 positive
delta candidates. Calibration mismatch stops immediately. A numerical result
above the current complete exact R6 upper endpoint (with a `1e-10` floating
trigger guard) stops further screening and exports all source labels/boxes for
independent rational inward certification. The script snapshots input hashes,
exact comparison endpoints, PID, elapsed time and checkpoints to
`correlated_results.json` before and after each candidate. The completed process
status is recorded; numerical results are never called new C6 bounds.

Run with the installed Python:

```
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python correlated_screen.py
```
