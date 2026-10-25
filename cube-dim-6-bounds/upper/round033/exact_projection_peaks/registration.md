# A033 exact projected Cauchy peaks

Start UTC2026-10-25 20:30:57; deadline1200 seconds. Only this directory is
writable. No subagents, MATLAB, native G/H calls, optimization, spatial grids,
stage, commit or push. Read old fan geometry and B008 only; do not rerun them.

Before the full pass: prove the one-sided pushforward trace formula, audit
integer fan rates/determinants against the exact geometric certificate, and
classify all4960 unordered body triples into the five proposed orbits.
Stop on any failure. The sole full enumeration has nine projection cases,
at most1e6 unique pattern/cone items total and9*967680 logical terms.
The actual scheme has701568 pattern/cone items. No global-sign pairing.

Normals: r1 all-ones; r2 Hamming k=1,2,3; r3 four projective-class count
partitions(4,1,1,0),(3,2,1,0),(2,2,2,0),(3,1,1,1),(2,2,1,1).
Trace directions v=(1),(1009,1),(1009^2,1009,1). All actual normal facets
must have nonzero v pairing. Preserve every selected m=r contribution with
its permutation/sign/cone index, normal determinant and exact weight.
The minimal checker rechecks selected terms and their exact sums, without
re-enumerating all patterns. Nine Arb comparisons with the already stored
B008 Q_d only; no old singular-tube or whole-face budget computation.

Output: exact rational R, q_d=sqrt(Gram)R/(2pi)^d, kappa=R/(2pi)^d;
q_d is a projection-density peak and scales as t^-d. This is no C6 bound.
