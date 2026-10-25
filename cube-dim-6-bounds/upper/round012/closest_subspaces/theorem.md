# Exact closest body and augmented subspaces in six dimensions

Let V={(1,s2,...,s6):sj=+/-1}, interpreted as32 unoriented lines. Let
S1(V) be their union and S2(V) the union of all spans of two distinct lines.
The latter already contains S1. For x in R6 put z=sort_desc(|x|), and define

    W(a1,...,am)=sum a_i^2-(sum a_i)^2/m.

Every list below has positive length. Define for k=1,...,5

    E_k(z)=W(z1,...,zk)+W(z_(k+1),...,z6).

Then on all of R6, with no strict-order or nonzero-coordinate hypothesis,

    dist(x,S1(V))^2=W(z1,...,z6),
    dist(x,S2(V))^2=min_(k=1,...,5) E_k(z).                    (1)

## 1. Full sign and permutation reduction

An orthogonal signed-coordinate permutation carries a cube-body line to
another cube-body line. If the first coordinate of a resulting representative
is negative, reverse the entire representative; its line and every pair
span are unchanged. Thus S1 and S2 are invariant under the full signed
permutation group, and Euclidean distances are preserved. It is enough
to prove (1) for nonnegative sorted z.

For a body line v, the squared projection residual is

    ||z||^2-(sum v_i z_i)^2/6.

The absolute signed sum is at most sum z_i, achieved by the all-positive
direction. This proves the S1 formula, including zeros and ties.

## 2. Every body-pair plane is a signed two-block plane

For distinct v,w let A={i:vi=wi} and B={i:vi=-wi}. Both blocks are
nonempty. Their pair plane is spanned by the two orthogonal vectors
v*1_A and v*1_B. Its squared distance is

    ||z||^2-(sum_A v_i z_i)^2/|A|-(sum_B v_i z_i)^2/|B|.

For a fixed partition, each absolute signed block sum is at most its
unsigned sum. The residual is therefore at least W(z_A)+W(z_B).
Conversely every nonempty partition is realized with both unsigned sums:
take v=(1,...,1) and w=1_A-1_B, reversing w if necessary so its first
coordinate is positive. The two representatives are distinct. Hence

    dist(z,S2(V))^2=min_(nonempty A,B partition)[W(z_A)+W(z_B)].

This proves the reduction from signed planes to unsigned partitions, not
merely an enumeration of sample vectors.

## 3. Contiguity without an empty-cluster or tie exception

Fix |A|=k and let S=sum z_i and s=sum_A z_i. The objective is

    sum z_i^2 - f_k(s),
    f_k(s)=s^2/k+(S-s)^2/(6-k).

The function f_k is convex because f_k''=2/k+2/(6-k)>0. The subset sum s
lies between the sum of the k smallest values and the sum of the k largest
values. Convexity gives

    f_k(s)<=max(f_k(smallest-k sum),f_k(largest-k sum)).

Thus a best partition of that cardinality may be taken as the top k or
the bottom k. A bottom-k partition is the complement of the top-(6-k)
partition, and its objective is the same. Taking k=1,...,5 proves (1).

This argument keeps both clusters nonempty and does not rely on a nearest-
center reassignment that might empty one cluster. If values tie, the extreme
subset sums are still well defined and unchanged by tie ordering. If all
values coincide or some are zero, every displayed identity and convexity
inequality remains valid. Nearest planes can be nonunique; the minimum
value is nevertheless exactly the five-term expression.

## 4. Augmented frame: the proposed eight terms are valid, two are redundant

Let V_beta=V union {beta e_i:i=1,...,6}, beta>0. Nonzero scaling does not
change any line or pair span. Every pair is body-body, body-axis, or axis-axis.

For span(v,e_i), coordinate i is free; the remaining five coordinates lie
on a signed constant line. As above, unsigned signs maximize the block sum,
so its best squared residual after varying v is W(z with coordinate i removed).
Write H=sum z_i^2 and S=sum z_i. As a function of the removed value t,

    L(t)=H-t^2-(S-t)^2/5,   L''(t)=-12/5<0.

The minimum over any observed value between min z and max z occurs at an
endpoint of that interval. Thus the proposed two leave-one-out terms are

    W(z2,...,z6) and W(z1,...,z5).

They are exactly E_1 and E_5 already in (1). There is a direct geometric
reason: if w is obtained from v by flipping coordinate i, its body line is
also in V, and

    w=v-2v_i e_i,
    span(v,e_i)=span(v,w).

So body-axis planes add no new subspace to base S2.

For span(e_i,e_j), the error is the sum of squares of all remaining
coordinates, minimized by retaining the two largest absolute coordinates.
Consequently the proposed eight-term minimum is correct and simplifies to
the six-term formula

    dist(x,S2(V_beta))^2
      =min(E_1(z),...,E_5(z),sum_(i=3)^6 z_i^2).                (2)

The eight-term version is (2) with E_1 and E_5 repeated. This holds for
every beta>0; beta does not enter the span geometry. Similarly,

    dist(x,S1(V_beta))^2
      =min(W(z1,...,z6),sum_(i=2)^6 z_i^2).

The augmented703 distinct direction pairs consist of496 body-body,192
body-axis and15 axis-axis pairs. Each body-axis plane duplicates a base
plane. The15 coordinate planes are distinct and none contains a body
direction, so none equals a base pair plane. The geometrically distinct
augmented pair-plane union is therefore the496 old planes plus15 new
coordinate planes. This counts511 planes, despite703 direction pairs.
For completeness, each base pair plane contains exactly two body lines:
in its nonempty blocks a body vector must have each block multiplier equal
to +1 or -1, leaving only the four vectors +/-v,+/-w. Thus different base
unordered pairs give different planes. There are96 distinct body-axis
planes, each counted twice by its192 body-axis pairs, all already in base S2.

## 5. Limits of the simplification and rank-dependent selections

Equations (1)-(2) concern nearest distances. Within the sorted absolute
chamber the base formula uses five candidate block planes; globally the
physical representatives depend on the signed permutation used to sort x.
The physical union is not a union of five fixed global planes, and its
tubes are not disjoint. Ties can create several nearest partitions, and
original pair-plane tubes can overlap even away from exact ties. A fixed
tie-breaking rule may choose an argmin for an algorithm, but does not turn
the original tubes into a disjoint geometric decomposition.

One membership consequence is S2(V)={x:|x_i| have at most two distinct values}.
The augmented set additionally contains every vector supported on at most
two coordinates. These are span-membership statements only; no density,
Green singularity, covering constant or operator bound is proved here.

Three distinct augmented directions need not be independent: take
v=(1,...,1), w=(1,-1,1,1,1,1), and beta e2. They satisfy
beta v-beta w-2(beta e2)=0 and have rank2. Thus an argument based only on
the number of large coordinates cannot assume any three directions are
independent. For a future large-set decomposition, choose a maximal
linearly independent subset by rank-greedy selection. It spans the whole
large-direction set and its size is the actual rank; remaining large
directions must stay in that same span. Analytic use of this selection
requires its own remainder and density estimates, absent from this task.
