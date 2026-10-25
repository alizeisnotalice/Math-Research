# Round009 completion: all current raw a1 ordering events

The task resumed from the locked8/10 checkpoint and executed only the last
two streams, registering and saving each before the next. No old stream was
rerun. The work concerns only the current1024 candidates and433 raw affine
endpoints, at fixed tail spacing and angle, with a1 in[3/2,2].

## Slopes 0 and 2: 11130 pairs

Exact outward intervals classify9859 pairs left of the window and1271
inside it. Deduplication leaves1072 window expressions, all already present
in completed streams. No new unique event or unresolved cluster appears.

Actual process PID38293, nice10,1 thread,
peak RSS34029568 bytes,0.143286042
seconds,11130 nodes, exit0. It exited before the last registration.

## Slopes -1 and 0: 16218 pairs

Exact classification gives9110 left,1018 window and6090 right pairs. The
687 distinct window expressions are all already known. There is no new
unique event or unresolved cluster.

Actual process PID38632, nice10,1 thread,
peak RSS31473664 bytes,0.188155542
seconds,16218 nodes, exit0. Both stream jobs used all six numerical-library
thread variables fixed to1. No pool, window or concurrent lower job was
created. Outside pairs are recorded by exact aggregate ranges and deterministic
pair-trace hashes, not large redundant arrays.

## Complete finite event coverage

All ten different-slope group pairs are accounted for, covering68351 raw
function pairs. The five slope groups are-1,0,1/2,1,2. The compiler also
verifies the constant order within each parallel group exactly: the fifth
power of an offset c*d^(1/5) is the rational c^5*d, and the real fifth-power
map is strictly increasing. No parallel pair can introduce a missing a1 event.

The saved result contains1947 unique formal potential event positions,
all inside the closed window, with disjoint exact isolating intervals and
no unresolved cluster. Boundary events are exactly3/2 and2. These critical
positions induce1946 open parameter components on which the raw endpoint
order cannot change. Explicit433-function permutations in every component
were not evaluated; only the ordering-event arrangement is now complete.

The certificate compiler only reads saved stream outputs. It used PID
39017, nice10,1 thread, peak RSS
29097984 bytes,0.040041250 seconds
and4415 validation nodes. It did not recompute an event stream or union.

Resume from`complete_events_10of10.json` for the complete sorted event catalog;
`checkpoint_10of10.json` preserves the last execution state. The9/10 checkpoint
is retained because it is the exact registered input to the final stream.
Full verbose execution outputs remain locally available and hash-identified;
rebuildable linked-case lists and duplicate singleton-cluster arrays are
excluded from the append whitelist. Every event value and raw-function pair
link in the saved compact checkpoints is preserved.

## What remains outside this result

There are no unprocessed streams for the present433 functions. The17 unresolved
planar capture cases remain entirely untouched:

    (8,2), (27,10..16), (36,10..14), (48,8..11).

The original18 materialized-packet queue and the prior analytic whole-plane
qualification(8,1) remain recorded; no fixed-ID partition was added.

This establishes completeness only of the a1 ordering events for the current
finite family. It does not evaluate a union at these events, solve the a1
optimization, establish uniform all-scale parameter completeness, or yield
a C5 claim. No parameter grid, angle search, plane calculation, commit or
formal-ledger update was made. Old8-flow locks and earlier packages were
hash-verified unchanged. All jobs have exited; no live exec handle remains.

## Portable replay

Restore the five dependency manifests listed in the append whitelist, together
with this append's ten files, preserving their relative project paths. The
saved-certificate compiler reads only those registered saved inputs; it does
not rerun an event stream. Run from the project clone:

```sh
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round017/append_last_two/compile_complete_events.py
```

Use a separate restored replay copy because the compiler writes its output.
Runtime PID, timing and RSS fields change on replay; compare the exact event
expressions, rational intervals, coverage and registered source/input hashes.
The two individual streams can also be replayed there with
`one_window_stream.py zero_two` and then `one_window_stream.py minus1_zero`,
under the same environment. Their verbose outputs are reproducible and are
not required by the compiler, which uses the compact checkpoints.
