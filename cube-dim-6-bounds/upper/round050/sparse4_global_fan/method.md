# B050 registered exact arrangement method

This task concerns the B044 frame body32+(1/80)sparse4 only. Parent1's
A046/F046 certificate is inherited read-only. There is no kernel,
heat-time parameter, new spatial sample, old-tube transfer, or C6 claim.

For each other parent let R be its six independent integer rays and
P={u>=0,sum u=1}. Cuts are precisely the new physical forms w for which
wR has both signs on P's vertices. All other new forms have fixed weak
sign; the old body forms are already sign-constant on the old parent.
The six coordinate boundaries and m crossing hyperplanes form H.

Every vertex of every sign cell has five independent active boundary
normals within sum u=1. Exhausting the binom(6+m,5) equality systems once
therefore gives all arrangement vertices in P, independent of signs.
Conversely an independent active-system solution in P is a vertex of
every sign cell in which it is feasible. All solves and feasibility
checks are exact rational operations. The ordered combination stream is
classified completely and saved; its counts and digest accompany the
complete vertex pool and each vertex's active witness. Time/cap failure
before exhaustion means that parent is uncovered, regardless of any
already discovered vertices. No inference from a partial pool is allowed.

For a sign vector sigma, filter this complete pool against sigma_j h_j>=0.
The retained points are exactly the complete vertices of the bounded
cell P_sigma. Its affine dimension is rank of the homogeneous vertex
rows minus one. A strict interior seed can be chosen with positive
coordinates proportional to (1,z,...,z^5): each nonzero cut is a
polynomial of degree<=5, so among 5m+1 distinct positive integer z
values at least one avoids all cut zeros. This is a finite exact
combinatorial seed construction, not a numerical optimization.

For a full five-dimensional cell, a cut j supports a facet exactly when
its zero subset of the complete cell vertices has affine dimension4.
Such a cut facet cannot be contained in a parent boundary facet: their
distinct hyperplanes would otherwise coincide on an affine4-space,
which would force proportional homogeneous forms, incompatible with
the cut crossing the parent. At a relative-interior point of the facet,
all other cuts and parent coordinates are strict. Crossing it flips
only sign j and gives a unique neighboring full cell. Distinct crossing
forms cannot define the same hyperplane because R is invertible and
the sparse4 orbit has no repeated unoriented lines.

The adjacency graph of full arrangement cells inside the simplex is
connected. Join two strict points by a polygonal path in the simplex
interior, with a generic intermediate point chosen to avoid the finitely
many codimension-two intersections and parallel degeneracies. It
crosses finitely many cut hyperplanes, one at a time, through relative
facet interiors. Hence BFS over exactly the facet flips visits every
full cell. The implementation verifies every proposed neighboring cell
has dimension5 and that the common vertex set is the recorded zero
facet. Queue exhaustion, not a numerical volume, is the BFS stopping
certificate. All remaining lower-dimensional sign cells belong to
closures of full cells: full-cell points are dense and there are only
finitely many closed full cells. No enumeration of 2^m sign patterns is
required and no zero face is discarded.

The saved checker independently verifies the entire ordered active ledger:
singular records have determinant zero; every outside record has a
nonzero determinant and its saved rational solution satisfies all six
equalities with a negative coordinate; every inside record has a
nonzero determinant and matches a saved nonnegative vertex satisfying
the equalities. First-appearance IDs and active witnesses must agree,
and the complete lexicographic combination stream must be exhausted.
This proves completeness of the saved vertex pool without solving for
new points. It adds at most one determinant/witness check per classified
system, so the planned construction and verification combined remain
within5000000 matrix operations. The checker also verifies facet-neighbor
closure/connectedness, maximal proper zero-face facets, pulling DAG,
true rates/signs/determinants and cross-parent faces. The initial recorded
verification shares the construction's1200-second actual-time window;
later explicit read-only replays have their own same finite limits.

Within every full cell and recursively every face, intersect its complete
vertices with each boundary hyperplane. The distinct proper intersections
of dimension d-1 are exactly its facets. Pull the least vertex and cone
over the recursively triangulated facets not containing it. Ray exit
from that vertex proves coverage and disjoint relative interiors by
dimension induction. Save the face DAG, facets, apex and simplices.

Use a single global order of primitive physical rays: the twelve A046
rays come first, in their frozen local ID order, then all other rays
in lexicographic integer order. This restricts to the existing A046
order, so parent1 need not be retriangulated. On shared faces the complete
arrangement rays are identical and pulling with a common order restricts
consistently. Positive rescalings between different parent affine slices
preserve these face/ray incidence statements. The implementation also
compares induced triangulations on each old parent facet: interior old
facets occur in pairs and must match exactly; unpaired facets must lie
on a wall of the ordered positive chamber. This is a redundant finite
compatibility check, not a replacement for the restriction proof.

For every vertex write Ru=s r, with r primitive integer and s>0 rational.
Compute q(r)=sum over all152 columns |a.r| using the true common denominator80.
Every full cell and every output simplex is checked for weak sign
constancy of all152 forms and the exact corresponding linear slope.
For a six-ray child verify

    |det R_child| product s_i = |det R_parent| |det U_child|.

Every child is nondegenerate. Relative simplex volumes |det U_child|
must sum exactly to1 separately for each parent. This is redundant to
the exhaustive vertex/halfspace/BFS/pulling proof. It cannot repair an
incomplete active-system enumeration or an unvisited cell.

The old21 parents are the already audited full cells of the original
body-form arrangement. Their signed-permutation copies partition R6
up to lower-dimensional walls. Combining their complete refinements
therefore gives the new exact full conic decomposition. The full group
has2^6*6! elements; determinants are absolute, no extra negative-ray
factor is inserted. A later numerical engine must separately handle
zero/repeated Fourier nodes and verify its rational-rate/int64 domain.
No numeric kernel or automatic old-engine compatibility follows here.

Limits: at most20 new parents,5000000 active systems,200000 full cells,
1000000 full simplices,1200 seconds from actual batch start. Each parent
is settled atomically and checkpointed. A cap/timeout/failure records
completed parents and an explicit uncovered list, retains the current
partial ledger, and stops without a rerun, raised limit, or reset.
Preparation/gate wait is recorded separately. There is no artificial
CPU or memory shutdown setting. Only a complete21-parent result may be
called a full frame fan; otherwise the result is explicitly partial.
