# Round 039: no positive nonsmooth first variation in the inherited127 family

For the exact round019 subfamily, this batch proves

\[
R(\varepsilon)=R(0)+c|\varepsilon|+O(\varepsilon^2),\qquad c\leq0.
\]

Thus a positive nonsmooth first-order gain is excluded. We do **not**
claim c=0, a sign for the quadratic term, or exclusion of every finite
epsilon. There is no new C6 bound; the accepted008 bound is unchanged.
The proof uses exposed-facet structure and exact domination, rather than
evenness or radius concavity alone. No new union was evaluated.

## Locked interface and exact velocities

Use the accepted008 spacings, positive u^6=2,v^6=3:

```
a1=1+u^4/2,
a2=u^2*v^2/3+u^2*v^4/6,
a3=...=a6=(1+u)/2.
```

Keep **all127** exported raw length labels, their half-cell source starts,
coordinate reflections, and receiver translations. No old96-box pruning
is assumed. No previously empty labels or larger lengths are introduced.
Put w(k)=1+epsilon*(-1)^k1, |epsilon|<1, on the unchanged lattice. The
true double-cell source mass is two. Write V(epsilon) for the union of
the127 boxes in H=product[0,aj/2] in the even receiver cell. The neighboring
cell uses V(-epsilon); each has64 disjoint reflected half-cell interiors.
As established in019,

\[
R(\varepsilon)=32[V(\varepsilon)+V(-\varepsilon)].
\]

For length ell let q1=-floor(ell1/2) for odd ell1 and q1=1-ell1/2
for even ell1. With m=product(ell), S=0 for even ell1 and S=(-1)^q1
otherwise, the selected mass and common radius at level one are

\[
M(\varepsilon)=m(1+\varepsilon S/\ell_1),\quad
r(\varepsilon)=\tfrac12M(\varepsilon)^{1/6},\quad
r'(0)=\frac{r(0)S}{6\ell_1}.
\]

The same radius is used in all six coordinates. Each coordinate interval
is the source-containment slab clipped to [0,aj/2]. At a **noncell** face,
the lower endpoint has velocity -r'(0) and the upper endpoint has velocity
+r'(0). In this family ell1=1 has S=+1, ell1=3 has S=-1, and even ell1
has velocity zero. At cell faces, clipping can give different one-sided
velocities, but it can never add volume outside the fixed H.

All radius functions are analytic near zero. Their finitely many endpoint
comparisons have a stable order on each sufficiently small punctured side
unless they agree identically. Finite-box union volume is consequently
analytic separately on these two sides, with bounded first and second
derivatives up to zero. This justifies the displayed first-order expansion.
In particular

\[
c=32\big(V'_+(0)-V'_-(0)\big).
\]

## Local facet criterion: where a positive cusp could occur

Partition the baseline box-face hyperplanes by all transverse endpoints.
Ignore boundaries of these patches, which have codimension at least two.
A sufficiently small fixed neighborhood of an interior patch has one of
the following forms:

1. A box already covers a neighborhood on both sides of the hyperplane.
   Small perturbations do not change occupancy there at first order.
2. Only boxes with the **same orientation** end on the hyperplane, and
   there is no covering box. This is a potentially exposed outer face.
   The moving union face follows the outermost endpoint velocity. The
   two-sided slope jump can be positive only if coincident participating
   faces have different velocities on a set of positive five-dimensional
   area. If their velocities agree, their first-order displacements agree,
   so their even first variation is zero.
3. Boxes end there with both orientations. The baseline fills a whole
   neighborhood of the hyperplane up to a null set. Perturbation can create
   a thin gap, but cannot increase the already-full neighborhood volume.
   Its even first-order contribution is nonpositive.
4. The hyperplane is a face of the fixed half-cell H. A baseline covered
   inner neighborhood is already full within H; outward motion is clipped
   and inward motion can only cause loss. Its even contribution is again
   nonpositive.

The discarded edge/corner neighborhoods have volume O(epsilon^2), because
all endpoints move O(epsilon) and there are finitely many boxes. Thus to
rule out c>0 it suffices to exclude **coexposed, same-oriented, unequal-
velocity facets of positive five-dimensional area**. This criterion does
not require computing the union volume or the magnitudes of the remaining
nonpositive contributions.

## Complete bounded inventory using the accepted base ranks

The existing E008 result already certifies all120 adjacent endpoint gaps
and their exact equality status at t=1. There is one boundary tie in axis
two, between ranks11/12, and no other distinct-group coincidence. All
labels referring to an identical exported endpoint ID also have the same
baseline coordinate. These inherited certificates give a complete baseline
coincidence partition; **no new field sign comparison or baseline union**
is needed.

`facet_inventory.py` reads the same127 distinct labels and processes their
six coordinates. It groups face occurrences by axis, accepted baseline
coincidence rank, and lower/upper orientation. It computes exact symbolic
velocities from the formula above. Cell-face occurrences are excluded from
the positive-cusp inventory by case4 of the local criterion.

There are exactly **10** noncell same-orientation groups with unequal
velocities: two upper-face levels in each of axes2 through6. Their data are

| mass | velocity of ell1=1 upper face | velocity of ell1=2 upper face | occurrences per axis |
|---:|---|---|---:|
|108|u^2*v^3/12|0|18|
|162|u*v^4/12|0|8|

All other noncell same-oriented coincidence groups have a single velocity.
No ell1=3 face occurs in an unequal-velocity same-orientation group.
The known324/2304 opposite-oriented endpoint tie is not mistaken for a
positive same-oriented max-envelope mechanism.

The full list of these10 groups, including source length labels, endpoint
IDs and exact velocities, is retained in `facet_inventory_results.json`.
This is a facet-label inventory, not a spatial partition/union computation.

## The25 static faces are hidden by15 existing static boxes

The zero-velocity labels in the10 groups consist of exactly these25 boxes:

- 20 labels B108 with ell1=2 and other lengths a permutation of
 [1,2,3,3,3];
- 5 labels B162 with ell1=2 and other lengths a permutation of
 [1,3,3,3,3].

For each B replace its unique other-axis length1 by2, leaving the first
length and every other length unchanged. Call this label C. The resulting
15 distinct C labels **already belong to the original127**:

```
B108 -> C216=(2, permutation[2,2,3,3,3]),
B162 -> C324=(2, permutation[2,3,3,3,3]).
```

Both B and C have even ell1, so their selected masses, radii and geometry
are independent of epsilon. Their radii satisfy rC>rB. In every unchanged
coordinate the containment slab grows monotonically with radius. In the
changed coordinate k>1, the new length2 interval is **all** [0,ak/2],
because rC>ak. It therefore contains B's length1 interval. Consequently

\[
B(\varepsilon)\subseteq C(\varepsilon)=C(0)
\quad\text{for every }|\varepsilon|<1.
\]

At any of the listed noncell upper faces, that coordinate has length3.
Since rC>rB, C's upper endpoint lies strictly beyond B's upper endpoint:
both the clipping cap and the unclipped endpoint of C exceed the noncell
endpoint of B. Except on transverse box-face boundaries of zero
five-dimensional measure, C covers a normal neighborhood of B's entire
candidate face. Thus the static face cannot be coexposed with the moving
face on positive five-dimensional area. It is harmless to retain all127
labels; the domination proof shows exactly why these25 cannot contribute
the suspected velocity split.

Here is a small exact bound proving the needed rC>ak simultaneously for
every changed coordinate. Set x=u^2,y=v^2. Integer-power comparisons give

```
u<6/5,   x<4/3,   y<3/2,
a_tail<11/10<6/5,
a2=x*y/3+x*y^2/6<7/6<6/5,
r216=sqrt(6)/2>6/5,
r324>r216.
```

`dominance_relations.py` verifies exactly12 bounded relation certificates:
six rational power/order bounds, three integer mass-order relations, two
label-family classifications, and one inclusion of the15 C labels in the
accepted127-label set. The25 explicit B->C maps are saved, making the
domination audit reproducible without a union routine.

## First-order conclusion and remaining scope

After the preceding exact domination, every potentially exposed noncell
same-oriented coincidence group has only one participating velocity.
Its even first variation is zero. Any opposite-oriented or fixed-cell
contribution is nonpositive, while edge/corner terms are second order.
Summing the local contributions proves

\[
V'_+(0)-V'_-(0)\leq0,\qquad R'_+(0)=c\leq0.
\]

This closes019's specific unresolved **positive first-order cusp** question
for the stated127-label family. It does not calculate c, prove strict
negativity, or determine a possible positive quadratic term when c=0.
A claim of full local optimality or exclusion at a nonzero epsilon would
need additional information and is not made. No epsilon=1/1024 candidate,
new block, or accurate union is evaluated or automatically authorized.

## Execution, dependencies and stop

Registration SHA-256:
`8363530e146d6200c97b3571d2cccc3943de02a893bb6caf61ffc929446c8709`.
Accepted interface SHA-256:
`99e762ad8fb63303f90b56c5d7f90fb855a3335903cae21650a23bfb07fea933`.
Accepted E008 result SHA-256:
`6849d3774ab3db1524969715970bbb3aecd07de4467faee5d77b69775cf5ed3d`.

| actual run | PID | exit | CPU seconds |
|---|---:|---:|---:|
|facet inventory, exec752908|45172|0|0.001490|
|12 relation certificates, execae26fa|45347|0|0.000900|

Total measured mathematical CPU is0.002390 seconds, below the2-second
target. The distinct inherited label set remains127; the second program
only checks membership against that same frozen set, not a new family.
The scripts use Python standard-library integer/Fraction arithmetic and
JSON/hash utilities. There are zero numerical root evaluations, new field
sign tests, union calls, baseline replays, grids, optimizers, nonzero
epsilon candidates, MATLAB calls or delegated agents.

`manifest.json` records source/result hashes, execution exits and bounds.
All writes are confined to this owned directory. No stage, commit or push
was performed. Both processes have exited. **Stopped; upper-bound work
retains absolute priority.**
