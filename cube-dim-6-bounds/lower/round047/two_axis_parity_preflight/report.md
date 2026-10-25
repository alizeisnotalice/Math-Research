# Round047: two-axis correlated-parity mass preflight

The new model has a valid four-phase/common-radius interface. It produces
new unequal-velocity coincident facets outside the round039 candidate list.
However, a new central-block merging lemma hides every such same-oriented
cofacet with positive five-dimensional overlap. Consequently this exact
127-label subfamily satisfies

\[
R(\varepsilon)=R(0)+c|\varepsilon|+O(\varepsilon^2),\qquad c\leq0.
\]

This excludes positive nonsmooth **first-order** gain in the new model.
It does not calculate c or exclude second-order/nonzero-epsilon gain.
There is no new C6 bound, no new union evaluation, and no change to either
accepted threshold. The argument does not import039's exclusion conclusion.

## 1. Exact new source and four phases

Use the locked008 spacings, positive u^6=2,v^6=3:

```
a1=1+u^4/2,
a2=u^2*v^2/3+u^2*v^4/6,
a3=...=a6=(1+u)/2.
```

At site a*k use positive mass

\[
w_\varepsilon(k)=1+\varepsilon(-1)^{k_1+k_2},\qquad |\varepsilon|<1.
\]

The periods are2a1,2a2,a3,...,a6. A half-open2-by2 source cell contains
four sites with masses1+epsilon,1-epsilon,1-epsilon,1+epsilon: **total4**.
The receiver-cell phases s=(s1,s2) have character chi_s=(-1)^(s1+s2),
with two positive and two negative phases.

Keep all127 inherited central raw length labels, their source starts,
reflections and translations. At the fixed008 candidate these are exactly
the nonempty standard length-1-through4 boxes, as certified in the accepted
008 family. No96-box pruning, length extension, previously empty label or
new spatial grid is introduced.

For each coordinate define the central source index block

```
J_l={q(l),...,q(l)+l-1},
q(l)=-floor(l/2) if l odd, and 1-l/2 if l even.
S(l)=sum_{k in J_l}(-1)^k
    =0 if l even, +1 if l=1, -1 if l=3.
```

For ell=(ell1,...,ell6), P=product_{j>=3}ellj and m=product_j ellj,
factoring the finite correlated sum gives the **exact selected mass**

\[
M_s(\ell,\varepsilon)
=P[\ell_1\ell_2+\varepsilon\chi_s S(\ell_1)S(\ell_2)]
=m[1+\varepsilon\chi_s\rho_\ell],
\quad \rho_\ell=\frac{S(\ell_1)S(\ell_2)}{\ell_1\ell_2}.
\]

Thus either even length makes the mass constant. For the odd cases,

| (ell1,ell2) | rho |
|---|---:|
|(1,1)|1|
|(1,3) or (3,1)|-1/3|
|(3,3)|1/9|

This is the full block count, including blocks crossing a period boundary;
fixed total period mass alone would not determine it.

## 2. Shared radius and ratio symmetry

For a chosen level lambda>0, use **one radius in all six coordinates**,

\[
r_s(\ell,\varepsilon,\lambda)=\tfrac12(M_s/\lambda)^{1/6}.
\]

At lambda=1 write r0=m^(1/6)/2. Its exact first derivative is

\[
r'_s(0)=\chi_s\rho_\ell r_0/6.
\]

The positive-half-cell witness box is

```
B_s(ell,epsilon)=product_j([(qj+ellj-1)*aj-r_s, qj*aj+r_s]) intersect H,
H=product_j[0,aj/2].
```

It contains exactly the centers whose common-radius cubes contain the
selected block, after clipping to H. Extra positive source masses help.
Noncell lower/upper endpoint velocities are respectively -r'_s(0) and
+r'_s(0). At a cell face clipping is handled by one-sided velocities;
outward motion cannot add volume outside fixed H.

Local coordinate reflection about each receiver-cell center preserves
(-1)^(s1+s2+h1+h2), since changing any h_j to -h_j preserves its parity.
It therefore preserves block mass and common radius. Each phase supplies
64 equal-volume reflected half-cell copies with disjoint interiors.
Let V(epsilon) be the union volume in H for chi=+1. The two chi=-1 phases
use V(-epsilon). Division by period mass4 yields

\[
R(\varepsilon)=\frac{64[2V(\varepsilon)+2V(-\varepsilon)]}{4}
=32[V(\varepsilon)+V(-\varepsilon)].
\]

Translation by one spacing in either first or second coordinate exchanges
the signs. R is even, but that fact alone does **not** imply differentiability
or a zero right derivative. As for any fixed finite box family with analytic
radii and min/max clipping, V is analytic on each sufficiently small side
of zero; hence c=32(V'_+(0)-V'_-(0)). The sign needs genuine facet analysis.

## 3. New facet splitting, using the accepted exact base partition

The inherited E008 certificate gives the complete endpoint rank/equality
partition at the fixed candidate:120 adjacent gaps and the single allowed
right-boundary tie in axis2. Identical exported endpoint IDs also give
identical baseline positions. Reusing this **base** certificate does not
reuse the old mass-direction velocities. The new script assigns the new
rho-dependent velocity to each occurrence of each endpoint and processes
all six coordinates of the same127 labels.

The indexing pass finds **31** unequal-velocity, same-oriented, noncell
groups. This is not31 separate spatial exposure calculations: it is one
finite label/rank inventory, followed by two class-wide bookkeeping checks.
Every group has one common baseline mass and one common normal-axis source
length. The groups are:

* axis1 upper faces at masses2,4,8,16,2304:5 groups;
* axes3--6 lower faces at masses4,8,16:12 groups;
* axes2--6 upper faces at masses108,162:10 groups;
* axes3--6 lower faces at mass2304:4 groups.

The actual endpoint IDs, source labels and exact velocities are retained
in `preflight_results.json`. There are21 new groups beyond039's ten-group
list. The old25-static-box dominance argument alone does not cover them.

A minimal explicit new pair is

```
A=(1,2,1,1,1,1),  mA=2, rhoA=0,
B=(1,1,2,1,1,1),  mB=2, rhoB=1.
```

Their axis1 upper faces coincide at r2=2^(1/6)/2, but have velocities0
and r2/6. Both first lengths are1, so this pair escapes039's specific
ell1=2 static-label deletion. Their transverse facet intersection has
positive five-dimensional area: in coordinates2 and3 a length2 interval
of one box lies inside the other's length1 half-cell interval, while the
other coordinates have identical positive intervals. Positivity follows
from the accepted nonempty A/B geometry and the central length1/2 formulas.
This is a real possible max-envelope split, not just a formal total-mass
argument. Its exposure must still be checked.

## 4. New merge lemma: all equal-mass cofacet overlaps are hidden

The central source blocks satisfy the exact chain

```
J1={0} subset J2={0,1} subset J3={-1,0,1} subset J4={-1,0,1,2}.
```

**Lemma.** Let two distinct inherited labels A and B have the same baseline
mass m and hence the same radius r=m^(1/6)/2. Suppose their coincident
noncell facets have a transverse intersection of positive five-dimensional
area. Define C by the coordinatewise maximum of their lengths. Then the
overlapping facet interiors are covered by a strict neighborhood of the
ordinary baseline box for C. Moreover C is already a member of the fixed
127-label family.

**Proof.** At a point x of that facet intersection, the common radius-r
cube contains the A source block and the B source block. Because all J_l
are nested, their coordinatewise source union has lengths C=max(A,B).
An axis-parallel cube containing both original Cartesian blocks contains
their coordinatewise Cartesian hull as well: in each coordinate it contains
both endpoint sets and hence the endpoints of J_Cj. Thus it contains the
entire C block at radius r.

Distinct labels with equal products cannot be componentwise ordered.
Consequently product(C)>m and rC=product(C)^(1/6)/2>r. Increasing the radius
from r to rC gives **strict slack in every source-containment slab**, including
the facet's normal direction. Away from transverse half-cell boundaries,
x therefore lies in an open neighborhood inside the clipped C box.
This also proves that C is nonempty. Since Cj<=4 and the accepted family
contains every nonempty standard central length-1-through4 box at these
spacings, C belongs to the original127 labels. No new label is added.
The positive radius slack persists under sufficiently small perturbation,
even if C's own correlated mass and radius vary. This proves the lemma.

For the explicit new pair above,

```
C=(1,2,2,1,1,1),  mC=4, rC>r2.
```

Its source label is already present (indices A=16,B=8,C=24 in the locked
interface). Its even second length makes it static in this model, though
staticity is not required by the lemma. It covers the entire A/B facet
overlap with positive normal slack. Thus the new positive-area overlap is
**not exposed**.

The class-wide inventory shows that every unequal-velocity same-oriented
group has a common baseline mass. Different velocities imply different
labels. Applying the lemma to any distinct pair in such a group excludes
positive-area coexposure. If the projected intersection has zero area,
it cannot contribute a first-order volume term anyway. This handles all31
groups without enumerating a new union or checking31 spatial coverage sets.

## 5. First-order acceptance and limits

At a regular exposed outer-facet patch, only unequal velocities among
coincident same-oriented participants could give an upward slope jump.
The merge lemma removes all such coexposed overlaps. Remaining same-oriented
participants have equal first-order velocity and give zero even first-order
contribution. If both facet orientations occur at a patch, the baseline
fills a neighborhood and perturbation can only create a gap, giving a
nonpositive even contribution. At a fixed H face outward gains are clipped
and inward losses are nonpositive. Codimension-two edge/corner effects
are O(epsilon^2). There are finitely many patches and positive radius
margins, so the preceding local classification proves c<=0.

This is a new-model first-order exclusion. It is not an inherited039
conclusion and not a general no-improvement theorem for alternating masses.
The exact coefficient c, its possible strict negativity, the second-order
term, and finite-epsilon improvement remain uncomputed. No nonzero epsilon
candidate, accurate union, optimization or automatic lower-bound certificate
is authorized by this preflight. No Cn bound changes.

## 6. Exact finite checks, execution and stop

The script executes12 bounded certificate predicates:

1. exactly127 distinct inherited labels/boxes;
2. accepted complete120-gap base certificate and its one right tie;
3. exact four-phase mass4 (epsilon=1/7 solely as a rational count check);
4. phase-character multiplicities two plus/two minus;
5--7. the three strict central-source nesting relations;
8. uniform baseline mass in every unequal-velocity noncell same-oriented group;
9. common normal-axis source length in those groups;
10. the explicit new A/B labels and their componentwise maximum;
11. their masses2,2,4 and correlated responses0,1;
12. the existing C label's membership in the locked family.

No numerical root or new field-sign evaluation is used. The count check
at1/7 is not a geometric epsilon candidate; no new radius or volume is
evaluated at that value. The pairwise merge proof is symbolic, not an
unrecorded numerical geometry scan.

Registration SHA-256:
`597e5d9b5c172b85cffba3ed8be5d4e1b89cb24477a9abed40fd2bdd885853d8`.
Accepted interface SHA-256:
`99e762ad8fb63303f90b56c5d7f90fb855a3335903cae21650a23bfb07fea933`.
Accepted E008 result SHA-256:
`6849d3774ab3db1524969715970bbb3aecd07de4467faee5d77b69775cf5ed3d`.

Actual run: PID49662, exec chunkbf6379, exit0, measured mathematical CPU
0.001540 seconds. It reads the127 labels once and processes their six
coordinates. Python standard-library integer/Fraction arithmetic suffices.
All code/result hashes, source dependencies and exact exits are recorded
in `manifest.json`. The only writes are the five files in this owned
directory. No union, baseline replay, length extension, grid, optimizer,
MATLAB call, delegation, staging, commit or push was performed.

**Finished and stopped. Upper-bound work retains absolute priority.**
