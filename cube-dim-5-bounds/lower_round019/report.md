# Lower round019: bounded planar capture diagnostics

Only the remaining16 finite (K,n) cases were registered. Their common
physical radius is r=(K*n)^(1/5)/2 and s=(1+u)/2, u^5=2. The normalized
radius rho=r/s=(K*n)^(1/5)/(1+u) is independent of a1. Fourteen distinct
masses occur: N=432 is shared by three cases, so only28 arrangements were
executed. Every inner/outer rational pair has width1/100000, below the
registered limit1/10000. The registration preserves exact root enclosures.

## Result and exact scope

Eight cases are rigorously empty: (27,11..16), (36,13..14). The outer
closed-period maximum is smaller than n in each case, and the actual
capture set is contained in that outer set. One case, (48,8), captures
at least eight points throughout its inner period, so its actual capture
set is the entire plane and its period area is exactly25s^2.

The other seven cases have only rigorous local area bounds, listed below.
They remain unresolved as exact capture cases. Areas in this table are
normalized by s^2; physical area bounds multiply the lower endpoint by the
saved lower s bound squared and the upper endpoint by the upper bound squared.

| K | n | N | exact actual min/max | normalized capture area lower | upper |
|---:|---:|---:|:---:|---:|---:|
| 27 | 10 | 270 | 7/10 | 104976/390625 | 6723649/25000000 |
| 36 | 10 | 360 | 8/12 | 604173893/100000000 | 151091917/25000000 |
| 36 | 11 | 396 | 8/12 | 736049/312500 | 1884921/800000 |
| 36 | 12 | 432 | 8/12 | 11042329/25000000 | 44182609/100000000 |
| 48 | 9 | 432 | 8/12 | 301302671/12500000 | 1205227391/50000000 |
| 48 | 10 | 480 | 8/12 | 524369757/25000000 | 2097619133/100000000 |
| 48 | 11 | 528 | 9/13 | 716120279/50000000 | 179051039/12500000 |

Inner and outer minima/maxima agree for all fourteen masses. Monotonicity
therefore certifies those extrema at each actual algebraic radius as well.
This does not make the seven partial area intervals exact areas.

## Rational endpoint and lattice guards

The normalized planar lattice is R Z^2, with R=(4,-3;3,4)/5. It has25
distinct representatives (U/5,V/5) in[0,5)^2. The exact integer preimage
is ((4U+3V)/25,(-3U+4V)/25); the divisibility tests and the invertibility of
R prove distinctness. Translation by(5h,5k) corresponds to the integer shift
(4h+3k,-3h+4k). The225 points for h,k in{-1,0,1} include every point whose
radius<2 window can meet the period. All omitted translates are farther
than the window radius in at least one coordinate.

Each window has four rational endpoints, clipped to[0,5]^2. Both endpoint
sets form51 by51 open cells. Integer bitmask intersections give their exact
counts and area histogram. A separate integer y-event sweep independently
recomputes the complete histogram, and agrees for all28 arrangements. Each
histogram sums to25, and every cell-mask trace has a saved SHA256.

All chosen radii satisfy2rho not in(1/5)Z. Thus starts and ends never
coincide along a coordinate. At any closed event stratum all touching
windows are captured in one neighboring open quadrant; the periodic seam
has no genuine event. Consequently the closed-period maximum equals the
open-cell maximum. The closed-period minimum cannot drop below an adjacent
open-cell minimum, and is attained in an open cell. The reported extrema
therefore also cover event lines, corners and the periodic boundary.

Integer fifth-power bisection at denominator2^80 encloses both u and N^(1/5).
Dividing the lower numerator root by1+the upper u bound gives a lower rho
bound; the reverse endpoints give an upper bound. The registered rational
endpoints strictly surround that interval. Hence their capture sets bracket
the actual set in the required direction, with no floating-point intervals.

Any selected planar n-tuple consists of distinct integer indices. Its
Cartesian product with K distinct unrotated points consists of Kn distinct
five-dimensional atoms, all at the same r. This diagnostic supplies the
four planar endpoint conditions only; a downstream five-dimensional witness
must also satisfy the six unrotated endpoint conditions at this same r.
Those unrotated conditions and a1 objectives were not newly evaluated here.
The planar cases introduce no a1 dependence or new a1 event.

## Historical evidence and resource limits

The old completed profile starts at physical side4. Every newly registered
side N^(1/5) is strictly below4. No old profile interval was repeated.
The same-radius old cases are cited in the registration, including the
N384,n12 historical capture; the current rational endpoint pairs are new,
and the requested threshold for N384 is n8. Duplicate mass432 is computed
only once per endpoint and reused for thresholds16,12,9.

The single actual process was PID42396, nice10, one thread; all six numeric
library thread variables were1. Each endpoint was saved before the next.
There were72828 elementary cells and119740 total cell/sweep-event nodes.
Every arrangement used2601 cells and4183..4387 nodes, below10000; every
arrangement took less than0.003 seconds, below30 seconds. The sum of actual
arrangement times was0.060146458040 seconds. Total process wall time was
7.574789000006 seconds, including the short pre-computation handoff.
Peak RSS was20905984 bytes. Session99179 exited with code0; ps
confirmed PID42396 absent. No process, pool or live handle remains.

Only the registered batch was completed. No a1 union,1946-interval
optimization, new angle, new lower bound, full-scale optimality claim,
formal-ledger change or commit was made. All old017/018 input hashes were
verified unchanged. Historical queues stay intact; this019 result is a new
append-only qualification/diagnostic layer.

## Compact evidence and replay

The archive contains one lossless arrangement_certificates.json container
rather than28 loose checkpoint files. Every original filename, full result,
trace, runtime and SHA256 is preserved; the final checkpoint is included.
The actual loose files remain locally available. This does not discard any
saved rational count/area certificate.

Restore the registered input dependencies and the seven archived files into
a separate replay clone with the same relative project paths. Run:

```sh
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round019/planar_bracket_diagnostic.py
```

Do not use --prepare or set LOWER019_HANDOFF for replay. With no loose
endpoint files this reconstructs only the registered28 arrangements and
writes new checkpoints. Runtime fields vary; compare exact radius parameters,
counts, area histograms, trace hashes and case areas. Alternatively, unpack
each container record with the serialization rule saved in the container
and verify its original SHA before resuming. No external numerical package
is required. Process-control start.flag is excluded from the archive.
