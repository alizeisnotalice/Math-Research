# Lower round020: three symbolic planar strata and seven exact areas

The new calculation is a symbolic capture-area construction on three open
normalized-radius strata, not a denser version of the019 endpoint grid. Only
three arrangements were executed. The earlier cases with known maximum counts
are not claimed as first exclusions.

## Geometry and a finite topology proof

Let u>0 satisfy u^5=2, s=(1+u)/2, and R=(4,-3;3,4)/5. The normalized planar
lattice is R Z^2. After scaling normalized coordinates by5, its points are
z=(4m-3n,3m+4n), with square period25Z^2. The225 integer-index atoms are the
same distinct25 period representatives and h,k in{-1,0,1} translates as019.
All relevant rho<1.7, so omitted translates cannot meet the full period.

Each uncut window endpoint is c-5rho or c+5rho with integer c; clipping
adds the constants0 and25. Each axis has452 raw endpoint descriptions and
152 unique affine functions. Exact identities are removed before ordering.
Across the union of the two identical function sets,11476 nonidentical pairs
are checked per stratum. Slopes are-5,0,5 and constants are integers, hence
every finite zero of a difference belongs to(1/10)Z. No such zero lies in any
registered open stratum. Distinct parallel pairs have no zero.

The three single order samples are29/20 in(14/10,15/10),31/20 in(15/10,16/10),
and33/20 in(16/10,17/10). They determine only the endpoint order and atom
bitmasks. All cell widths remain exact affine functions of rho. In each
axis,52 clipped endpoint functions give51 intervals. Every width is
nonnegative at both closed stratum endpoints and positive at the midpoint;
by affinity it is strictly positive throughout the open stratum. All window
clipping/emptiness choices and cell masks therefore hold throughout that
stratum. Each result preserves the full51-by51 cell mask array and widths.

Multiplying the two affine widths gives a quadratic scaled-area polynomial.
Division by25 converts to normalized plane area. Summing by exact captured
point count gives the following coefficients, where each row means
a*rho^2+b*rho+c. Boundaries have zero area; the formulas are exact on each
open stratum. No topology claim is made at its excluded radius endpoints.

## Complete count-area histogram polynomials

| m, rho in(m/10,(m+1)/10) | count | a | b | c |
|---:|---:|---:|---:|---:|
| 14 | 7 | 400 | -1200 | 900 |
| 14 | 8 | -500 | 1280 | -791 |
| 14 | 9 | -300 | 1040 | -868 |
| 14 | 10 | 400 | -1120 | 784 |
| 15 | 8 | 200 | -660 | 544 |
| 15 | 9 | 100 | -440 | 452 |
| 15 | 10 | -800 | 2560 | -2036 |
| 15 | 11 | 400 | -1160 | 840 |
| 15 | 12 | 100 | -300 | 225 |
| 16 | 9 | 400 | -1360 | 1156 |
| 16 | 10 | -200 | 540 | -340 |
| 16 | 11 | -800 | 2680 | -2232 |
| 16 | 12 | 500 | -1540 | 1185 |
| 16 | 13 | 100 | -320 | 256 |

For each of the three strata, summing every count row gives the exact
coefficient identity(0,0,25), so the complete normalized period has area25.
For every threshold n=1..16, the result also saves the exact polynomial
obtained by summing rows with count>=n.

## The seven previously partial cases

Here A is normalized capture area, so physical capture area is s^2*A.
Each actual radius is strictly inside its assigned stratum by the already
saved019 rational brackets. The exact threshold-area coefficients are:

| K | n | N=Kn | m | a | b | c |
|---:|---:|---:|---:|---:|---:|---:|
| 27 | 10 | 270 | 14 | 400 | -1120 | 784 |
| 36 | 10 | 360 | 15 | -300 | 1100 | -971 |
| 36 | 11 | 396 | 15 | 500 | -1460 | 1065 |
| 36 | 12 | 432 | 15 | 100 | -300 | 225 |
| 48 | 9 | 432 | 15 | -200 | 660 | -519 |
| 48 | 10 | 480 | 15 | -300 | 1100 | -971 |
| 48 | 11 | 528 | 16 | -200 | 820 | -791 |

For each row take the positive real root w with w^5=N. Then

    rho = w/(1+u),
    physical area = [a*w^2+b*w*(1+u)+c*(1+u)^2]/4.

This is a direct exact algebraic expression. It is not inferred from decimal
display or an interval midpoint. All seven previous partial exact-area
questions are resolved. The archived result also gives the exact formulas
for the other nine019 cases, agreeing with their existing zero/full results.

Every saved019 inner and outer endpoint was substituted into the applicable
polynomial. All count-histogram entries for28 endpoint records match exactly
as rational numbers, and all16 threshold-area pairs agree exactly with019.
This stage evaluated polynomials only and did not rerun an019 arrangement.

## Distinct atoms, common radius and scope

At a five-dimensional witness, n selected planar lattice indices are distinct.
Their Cartesian product with K distinct unrotated indices has Kn distinct
unit atoms. The same common radius is r=(K*n)^(1/5)/2 on every coordinate,
and rho=r/s. No atom has an individually chosen radius. The full masks allow
a deterministic choice of n distinct planar atoms wherever count>=n.

The current result supplies exact planar capture measures. It does not
check all six unrotated endpoint conditions or construct a five-dimensional
union. Planar masks depend on rho but not on a1; no a1 interval objective or
additional a1 event was evaluated. Areas of different cases must not be added
as a union area. There is no new C5 bound, parameter optimum or all-scale
optimality conclusion.

## Actual execution and reproducibility

Only the three registered arrangements ran, each with2601 cells, for7803
total cells. Per-stratum time was
0.036978250020, 0.034298166982, and 0.036814999999 seconds;
their sum was0.108091417002 seconds. Each was below the30-second and10000-cell
checkpoint limits and was saved before the next. The sole PID was46755,
nice10, one thread, with all six numerical-library thread variables set to1.
Total process wall time was9.824672707997 seconds including the short
handoff; peak RSS was22642688 bytes. Session68500 exited with
code0 and ps confirmed PID46755 absent. No process or live handle remains.

The nine-file whitelist retains the actual source and preregistration, scaled
integer atom manifest, all three full cell/mask certificates, exact results,
this proof report and the completion receipt. Reconstructable checkpoint and
process-control start.flag are local only. All019 input hashes were verified
unchanged; no017/018/019 lock, formal ledger or main file was modified.

Restore the four registered019 dependency inputs at their relative project
paths, then run in a separate replay clone:

```sh
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round020/symbolic_strata.py
```

Do not set LOWER020_HANDOFF for replay. To rerun the three arrangements rather
than load saved certificates, omit stratum_m14.json through stratum_m16.json
in that separate replay copy. Runtime fields vary; compare the affine
functions, masks, exact polynomial coefficients and019 rational evaluations.
Python standard-library exact integers and fractions suffice. No new task
was started after the registered batch, and nothing was committed.
