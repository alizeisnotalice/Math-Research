# B050: complete exact sparse4 fan, with all remaining parents certified

The body32+(1/80)sparse4 frame now has a complete exact conic fan on the
ordered positive Fourier chamber. The final object has478 primitive
rays and13,229 full-dimensional simplicial cones, refining4,266 full
arrangement cells over all21 old parents. All20 previously unrefined
parents were completed and checkpointed; the uncovered list is empty.
The already audited first parent was inherited without solving or
re-enumerating any of its462 active systems.

The separate saved-certificate checker passed the entire new active
classification ledger, complete vertex pools, facet-adjacency closure,
pulling recursion, true152-column signs/rates, determinants and shared
parent faces. This establishes a geometric new-frame kernel interface.
It computes no kernel value, Green value, time average, tube constant,
g/I budget or C6 bound. An external full-fan review and any numerical
engine domain review remain separate; F046's independent PASS applied
only to the inherited first parent.

## Locked inputs and actual limits

The frame is the exact152-column object from B044, SHA256

    67a12cd67961baf5ee554c138c82d6b4b627701c66dc691cb6f8b2a823c9be1c.

The final A046 manifest pin is
`e0b943a42323ddb0f9cf2c937141bb252ee17c93e55927fd3de343cba6f74ca7`;
its inherited local fan is
`bae26cd08d31dac644e5ae0e77483ae912de0f687a03d6bc3d7b98a4a35a78b6`.
The F046 final report/manifest/results pins are respectively

    07dc3074ba73c184a8eea7f73f4e3a55bbfec55f583554a9f38b8570b448cd7d
    6b347c5b9aa900f11893afbd0a7e706fda066c7891496d119eed7226bdd59edd
    6bec898cb72118c8526b4363acd94cc7513aa94a95fad24e0ec0b61ed32eb1ad.

These were checked before execution. The explicit parent release after
F046 PASS is recorded in registration. The final registration SHA is
`d119907131a2910cbae47796b1d2d54c79f2d5222c6bb3430a3819197daaa649`.
It locks75 source/runtime files, including all installed Python-FLINT
extension binaries and its FLINT/GMP/MPFR libraries. No environment was
installed or changed. Python3.15.0, python-flint0.9.0 and FLINT3.6.0 were
used for exact arithmetic.

The limits were20 new parents,5,000,000 active-system/matrix operations,
200,000 full cells,1,000,000 full simplices, and1,200 seconds of actual
construction plus the initial saved verification. There was no artificial
CPU or memory shutdown setting. The exact planned new active count was
783,994. The saved checker adds783,994 determinant/witness checks, without
a solve or new vertex search, giving1,567,988 combined matrix operations.
All counts and times remain below the registered caps.

The original preparation registration was amended before any active
calculation to finish the saved checker, pin the compiled dependencies,
and serialize the actual-start clock. The old and new preparation hashes
are preserved in `preparation_history.json`. One text patch context
mismatch was corrected during preparation; it performed no computation
or partial file edit. There were no mathematical failures, active-solve
retries, batch restarts, checkpoint resets or budget extensions.

## Complete vertices without a sign-pattern explosion

For a parent with six ray columns R, use the slice
P={u>=0,sum u=1}. For each sparse4 form w, its pullback is wR. Only forms
with strict positive and negative ray values need cuts. With m such
cuts, the boundary list consists of six coordinate hyperplanes and m
cut hyperplanes. A full-dimensional vertex requires five independent
active boundaries in the affine5-space. The producer therefore solves
each of the binom(6+m,5) equality systems once, with no sign assignments
in the solve. Nonnegative solutions are deduplicated into a global
vertex pool; each keeps its first independent active witness.

All783,994 ordered systems were classified exactly:

| Classification | Count |
|---|---:|
| Singular | 247,176 |
| Nonsingular, outside the parent simplex | 319,333 |
| Nonsingular, inside the parent simplex, before vertex deduplication | 217,485 |

Each parent's compressed active ledger stores every combination in order,
its classification, the full rational outside solution when applicable,
or an inside-vertex enumeration ID. Its uncompressed digest and compressed
file hash are retained. The saved checker recomputed each determinant:
singular records must have determinant zero; outside records must have
nonzero determinant, satisfy the six equations and contain a negative
coordinate; inside records must have nonzero determinant and their saved
vertex must satisfy all equations. It also checked first-appearance IDs,
active witnesses, the full lexicographic combination stream, and exact
exhaustion. Thus vertex completeness is checked directly, not inferred
from BFS or a volume sum. The checker did not solve these systems again.

All subsequent sign-cell tests filter this complete vertex pool against
the chosen signed halfspaces. Every feasible pool point is a vertex of
that cell, since its five independent active equalities are still
boundaries of the cell. Conversely every cell vertex occurs in the pool.

## Full-cell BFS and boundary coverage

A finite exact generic interior seed is obtained from positive u
proportional to(1,z,...,z^5). Each cut is a nonzero polynomial of degree
at most5, so at least one of5m+1 distinct positive integer z values
avoids all zeros. No spatial kernel or optimization uses this seed.

For each full cell, a cut supports an interior facet precisely when
the cut-zero subset of its complete vertices has affine dimension4.
The adjacent full cell is obtained by flipping just that cut's sign.
Distinct cuts cannot coincide, and a crossing cut facet cannot lie
entirely in a parent boundary facet. The producer checks the adjacent
cell's dimension5 and exact common facet before enqueueing it.

The full-cell adjacency graph in a simplex is connected: a generic
polygonal path between strict interior points crosses cut hyperplanes
one at a time away from their codimension-two intersections. Thus BFS
queue exhaustion gives every full cell. The saved checker independently
filters all saved cells, reconstructs every cut facet and forced neighbor,
and verifies the recorded graph is connected and closed under all such
facet flips. It performs no new sign-cell search. No2^m pattern loop was
used. Lower-dimensional zero and tie regions are retained by closure:
strict full-cell points are dense in the parent and the finite union of
closed full cells is closed.

These completeness arguments, together with full vertex classification,
are the reason BFS is valid here. Neither BFS alone nor exact volume
equality would establish the necessary coverage.

## Pulling, true frame rates and cross-parent compatibility

All parents use one physical-ray order: the twelve inherited A046 rays
first, in their frozen ID order, then all other primitive physical rays
in lexicographic order. This order restricts to the original A046 order;
the first parent's14 simplices remain unchanged. No first-parent vertex,
active-system or local-fan construction was rerun.

Every full cell and recursively used face has a complete saved vertex
set. The producer obtains facets as proper boundary-zero intersections
of dimension one less. The checker independently takes the inclusion-
maximal proper zero faces and verifies those dimensions and facet sets.
The least vertex is pulled, and facets not containing it are recursively
triangulated. The saved face DAG includes each apex, facet and simplex.
Coning from the apex over the exit facets proves coverage and disjoint
interiors by induction. There are20,888 new face-DAG nodes; the inherited
first parent's proof DAG remains its separately audited object.

For each rational u vertex, Ru=s r is recorded with primitive integer
r and exact positive rational s. The true rate is recomputed from all
152 column numerators with common denominator80:

    q(r)=(1/80) sum_j |column_numerator_j.r|.

Every full cell and every output simplex separately passes all152 weak
sign tests. Its exact slope satisfies slope.r=q(r) at every ray. Old
rates q are explicitly discarded when parsing the parent fan.

For every child, both its primitive physical determinant and rational
slice determinant are checked, with

    |det R_child| product_i s_i = |det R_parent| |det U_child|.

Every determinant is nonzero and every parent's relative volume sums
exactly to1. These are redundant checks of the previously established
partition and pulling coverage. The global xi1=1 slice has exact volume
1/120, computed by summing |det R_child|/(120 product_i r_i1). This is
the ordinary ordered chamber slice, not a psi=1 section; no heat-kernel
peak integral was evaluated.

All126 old parent facets are accounted for. The28 paired interior facets
have identical induced triangulations on both sides. The70 unpaired
facets each have an exact ordered-chamber wall witness. These checks
include the inherited first parent's saved boundary subdivisions, solely
as read-only combinatorial data. Common global pulling order and complete
arrangement vertices give compatibility on lower-dimensional shared
faces as well. Positive rescaling between parent slices preserves the
same ray/face incidence relations.

Across signed-permutation copies of the ordered chamber, a shared wall
is fixed pointwise by its adjacent simple reflection (a tied-coordinate
swap or final-coordinate sign reversal). Therefore the induced boundary
subdivision also matches across chambers. Full-dimensional interiors
have unique signed-permutation chamber representatives. This justifies
the full2^6*6!=46,080 multiplicity without a boundary exception or an
extra factor2. Any later global-negative pairing is a separate algebraic
rewriting, not an extra weight in this fan.

## Per-parent settlement

Parent1 is inherited:12 vertices,12 full cells,14 simplices, zero new
active systems. All other parent certificates were produced once:

| Parent | Active systems | Vertices | Full cells | Simplices |
|---:|---:|---:|---:|---:|
| 2 | 6,188 | 32 | 63 | 106 |
| 3 | 3,003 | 19 | 36 | 38 |
| 4 | 11,628 | 54 | 138 | 275 |
| 5 | 6,188 | 32 | 63 | 105 |
| 6 | 462 | 12 | 12 | 14 |
| 7 | 53,130 | 84 | 225 | 698 |
| 8 | 33,649 | 64 | 142 | 320 |
| 9 | 80,730 | 204 | 575 | 2,819 |
| 10 | 53,130 | 117 | 272 | 1,037 |
| 11 | 11,628 | 35 | 49 | 102 |
| 12 | 42,504 | 44 | 175 | 215 |
| 13 | 65,780 | 99 | 286 | 891 |
| 14 | 142,506 | 215 | 813 | 2,860 |
| 15 | 98,280 | 164 | 515 | 1,638 |
| 16 | 26,334 | 56 | 112 | 227 |
| 17 | 98,280 | 109 | 397 | 907 |
| 18 | 26,334 | 65 | 165 | 389 |
| 19 | 15,504 | 51 | 132 | 384 |
| 20 | 4,368 | 28 | 42 | 133 |
| 21 | 4,368 | 21 | 42 | 57 |

Local vertex counts include shared rays; their global deduplication gives
478 primitive rays. Every parent receipt has relative volume exactly1.

## The new interface and its computational limit

For any t>0, the proved frame symbol is linear with rates q_j on each
new simplicial cone with columns r_j. Its exact Fourier-cone integral is

    |det R| product_(j=1)^6 (t q_j-i r_j.x)^-1.

Summing all13,229 ordered-chamber cones and all46,080 group copies, then
multiplying by(2pi)^-6, gives the true p_t Fourier integral. The real
parts of tq_j are positive, so this positive-time cone integration is
absolutely valid. No value of t or x was substituted in this task.

A future direct full-group kernel would have609,592,320 logical terms,
an exact cone-count factor13229/21 over the old967,680-term body32
formula. This is a concrete structural cost increase, not a measured
runtime or an automatic improvement of a weak-type constant. Rational
rates, integer scaling/range, zero Fourier nodes, repeated nodes,
precision/cancellation and any Green time integration still need their
own engine review. In particular integrating individual cone expressions
down to time0 is not asserted safe merely because the positive-time
fan identity is valid. No existing C013/C014 numeric engine is declared
compatible automatically.

The old body32 tubes, Green data and all-time obstruction do not become
new-frame statements. No weak C6 conclusion follows from a new fan alone.
Public and internal upper bounds remain64.

## Artifacts, clocks and reproduction

The final rational fan is `global_fan.json.gz`, SHA256

    2405ccad9b2b0c344167539e769e3c64e3275e6e4ead208fbedc3167be47c948.

It stores478 primitive rays with exact rational rates and13,229 cones
with exact determinant weights. Parent sources are identified on each
cone. `cross_parent_faces.json.gz` stores the complete shared/boundary
face checks, SHA256

    4a02fe422405787deac463836b9f29814e908631407d613929b441f9f3ec0c02.

The saved verification result SHA256 is
`9dfb39f4eb8e20e2a1b96d1acfbff534c9470984c026050fe3f318ddcff3c752`.
For parents02-21, `parents/NN.fan.json.gz` contains complete vertices,
full cells/adjacency, face DAGs and child simplices; `NN.active.jsonl.gz`
is the complete classification ledger. These compressed ledgers are
needed for the independent vertex-completeness check, not incidental
logs. The main receipt and checkpoint preserve every settled parent
and the now-empty uncovered/global-obligation lists.

Preparation and gate handling were separate from the actual clock.
The first preparation registration was at2026-10-25 21:54:33.878 UTC;
the final pre-run preparation at21:58:24.433 UTC. Earlier source reading
and theory work were not counted as solver time. The actual clock began
at21:58:50.482 UTC. Constructor PID50756 completed in18.372790 seconds,
exit0; tool session96654 closed (launch8e337f, completion89e04a).
The saved checker PID50800 completed in11.583378 seconds, exit0;
session29446 closed (launch8b4659, completion7f0240). From actual batch
start through recorded checker completion,45.593994 seconds elapsed,
including the interval between calls, within the common1,200-second cap.

The constructor is single-use and refuses to reset an existing receipt.
The recorded checker added no solve, new vertex, sign-cell search or
first-parent enumeration. All handles are closed. No kernel/G/H/K/p,
heat-time parameter, native kernel subprocess, spatial grid, optimizer,
MATLAB, installation, old-source write, agent expansion, commit or push
was performed. FLINT exact linear algebra was used as registered.

The full read-only saved check can be reproduced with:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round050/sparse4_global_fan/check_saved.py
```

The default replay has its own finite20-minute limit and launches no
kernel or solver search. The original `--record` check shared the
construction clock. The task is now stopped and frozen; no numerical
kernel or new time experiment is authorized by the unused limits.
The parent retains final mathematical review and integration.
