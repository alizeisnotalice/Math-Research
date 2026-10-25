# Lower round018: one exact planar two-cover qualification

This batch resolves only the previously pending planar case (K,n)=(8,2).
The tail spacing is s=(1+u)/2 with u^5=2. The planar lattice is
s R Z^2, R=(4,-3;3,4)/5. All 16 selected atoms have the common radius
r=16^(1/5)/2, so (2r)^5=16. The planar result is independent of a1.

## Exact rational coverage certificate

Normalize physical plane coordinates by s. The lattice R Z^2 has a square
period 5Z^2 and exactly 25 distinct representatives in [0,5)^2. Write each
representative as (U/5,V/5), 0<=U,V<25; its integer preimage is
((4U+3V)/25,(-3U+4V)/25). The two displayed divisibility conditions select
exactly 25 representatives. Adding (5h,5k) changes the integer preimage by
(4h+3k,-3h+4k), so every representative and translate is a genuine lattice
atom, without duplicates.

The 225 atoms from h,k in {-1,0,1} include every atom whose closed radius
4/5 window can meet [0,5]^2. Omitted translates have at least one coordinate
distance greater than 4/5 from the domain. Their clipped window endpoints
form 25 by 25 elementary cells. Every cell has an explicit pair of distinct
integer atom indices, with all four planar endpoint inequalities checked
using exact fractions. Closed cells tile the entire closed period, including
boundaries. The exact count-area histogram is {2:12,3:12,4:1}; its area is25.
Thus every point is captured by at least two distinct atoms at radius4/5.

Integer bisection of fifth powers gives the strict physical radius margin

    r-(4/5)s >= 133842852809855947736597/12089258196146291747061760 >0.

Hence the actual physical two-capture set is the entire [0,5s]^2 period and
its translates, and its area is exactly25s^2. This is a solved qualification,
not a candidate. The saved 625-cell certificate also provides a finite local
fixed-pair selector; it is not a new global five-dimensional union partition.

## Distinct mass and the ten geometric endpoint guards

The surviving raw unrotated shape for this K,n in a1 in[3/2,2] is (2,2,2);
the other K=8 shapes are already uniformly excluded by their first/tail span.
Its unrotated index set is {0,1}^3. Combining its eight distinct indices with
each certified pair of distinct planar indices gives exactly16 distinct
five-dimensional atoms. Three explicit product examples are saved; in general,
injectivity follows from the Cartesian product and the invertibility of R.

For the first axis the permissible center interval is

    [max(0,a1-r), min(a1/2,r)].

For each of the other two unrotated axes it is

    [max(0,s-r), min(s/2,r)].

These formulas are the intersection of the half-period with
[atom_max-r, atom_min+r], and therefore prove the six unrotated endpoint
inequalities. Each planar cell supplies the other four endpoint inequalities
at radius(4/5)s, strictly smaller than the same common r. No atom receives an
individually chosen radius. For positive-volume centers, an interior collar
allows a common radius decrease; no full transfer theorem is newly asserted.

For shape(2,2,2), positive five-dimensional volume occurs precisely when
3/2<=a1<2r=16^(1/5). At equality the first-axis interval is a singleton,
so its five-dimensional volume is zero; the interval is empty only above
2r. The originally printed result used the shorthand "empty at and above";
the archived JSON text clarifies this endpoint without changing arithmetic.
This activation endpoint is already present in the completed017 event
catalog, as the formal expression [[16,"1"]]. It introduces no missing event.

## Actual execution and remaining scope

Preregistered limits were30 seconds and10000 candidate cells. The sole
calculation used PID41466, nice10, one thread, with all six numerical-library
thread variables set to1. It checked625 cells in0.008406667 seconds and had
peak RSS21659648 bytes. Session56654 exited with code0; no live handle remains.

The unresolved planar list decreases from17 to16:
(27,10..16), (36,10..14), (48,8..11). No additional case was started.
The old017 checkpoint remains unchanged and continues to record its historical
17-case queue; the new018 result is the append-only override for this one case.
No1946-component union, angle grid, parameter optimization, new lower bound,
formal-ledger change or commit was performed.

## Replay

Restore the registered017 inputs at the recorded relative paths and run in
a separate replay copy:

```sh
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round018/qualify_K8_n2.py
```

Do not set LOWER018_HANDOFF for replay. The short handoff files are local
process-control artifacts and are excluded from the archive. Compare exact
certificate fields, allowing runtime fields and the documented endpoint
wording clarification to differ. All arithmetic uses the Python standard
library; there is no floating-point interval-library dependency.
