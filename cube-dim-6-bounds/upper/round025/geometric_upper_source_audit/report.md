# C025: geometric source audit and a support-conditioned strict gap

**Outcome:** no new global upper bound for C6. The current bound remains
64. This batch proves a six-dimensional support-conditioned estimate,
refutes a uniform-overlap shortcut, and isolates a concrete selection
lemma that would be sufficient for a global improvement.

## Proven mathematical output

For f>=0 supported in a cube Q, assume lambda|Q|<=integral f. The full
centered, common-radius cube maximal operator satisfies

    lambda |{M_c f>lambda}| <= (64-43/15625) integral f.

This holds for all radii but only at heights compatible with the mass and
one cube containing the support. It is **not** the general-L1/all-height
estimate defining C6. The proof first treats radii no greater than the
support radius: of 64 possible corner regions, at most 21 can carry the
mass necessary to produce a bad center, leaving 43 empty regions. For
the all-radius extension, enlarge the containing cube until its volume
is exactly integral f/lambda; larger radii cannot be bad.

The analytic proof is in theorem.md, Section 2. Eight exact Fraction
checks pass in geometry_results.json, including 22*(3/5)^6>1 and the
exact constant. These checks are arithmetic sanity checks, not a finite
grid standing in for the proof. No kernel, MATLAB, optimization, spatial
grid, or new time-parameter evaluations were made.

Section 1 constructs 64 radius-one cubes with an interior common point
and pairwise center exclusion, while every pairwise overlap ratio and the
total volume redundancy tend to zero. Therefore a uniform geometric
intersection bonus is false. This does not construct a near-extremizing
L1 function or establish fixed-Lebesgue sharpness.

## Source finding and its scope

Four web queries and three primary author documents were actually read
and pinned; all URLs, reading locations and hashes are in search-log.md
and source_receipt.json. The current Lin README is byte-identical to the
old cache. Its general cube result remains 2^d, its 3.879 estimate is
planar, and the displayed higher-dimensional improvements use Euclidean
balls. The newly inspected UpperBound.lean locates the covering proof's
volume and input-mass losses. The author's balayage article still asks
for a positive admissible kernel and does not provide a six-dimensional
covering refinement in the sections inspected.

These are source/scope findings, not a fresh verification of the entire
Lean development or balayage theory. Search snippets did not support
any theorem. The search was bounded and does not prove absence of an
explicit refinement elsewhere. Previously verified Aldaz general-measure
sharpness is retained only as background and is not a Lebesgue obstruction.
The three current documents have no verified publication/release date;
the bibliography therefore carries only the relevant historical verified
Aldaz item, with that limitation made explicit.

## Concrete remaining acceptance object

Section 3 proves the exact loss identity for disjoint selected cubes:

    a M-lambda|E| = a(B+O)+lambda D,   a=(1+tau)^6,

where B is witness mass excess, O is unselected mass, and D is enlarged
volume minus bad-set volume. None of these terms alone has acquired a
uniform positive lower bound from the literature.

Section 4 states (L6): select disjoint cubes whose averages reach lambda,
so that the part of {M_c f>lambda} not detected by any isolated component
at height (1-43/4000000)lambda has measure at most
(43/62500) integral f/lambda. This is **unproved** and could be false.
If established for all L1 inputs and heights, the support-conditioned
estimate would yield C6<=64-43/31250. That number is conditional only.
The next acceptance task is to prove or refute this explicit foreign and
unselected mass estimate. The support-conditioned result already handles
all radii within an isolated component; general cross-component coupling
is the missing step. No follow-up calculation was started.

## Provenance and resource receipt

Registration and old-source hashes preceded every new query and exact
geometry check. The geometry sub-plan was written before its single
execution; an incorrectly typed manual timestamp was corrected from the
filesystem birth time afterward, with the original input hash retained in
geometry_results.json. No mathematical input or check was changed and no
check was rerun.

Owned directory only: upper/round025/geometric_upper_source_audit.
All previous rounds remain read-only. Counts: 4 queries, 3 primary
current documents, 8 named exact checks, 0 kernel passes, 0 MATLAB calls,
0 subagents. No stage, commit, or push. The final receipt records elapsed
time, the file whitelist and hashes. Frozen at the end of this batch.
