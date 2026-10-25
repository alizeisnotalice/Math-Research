# True pair-frame S2 increment and rank safeguards

Let S2 be the union of spans of two distinct unoriented columns of
V union beta P. Scaling beta>0 does not affect this geometry. There are
1891 direction pairs. Their geometrically distinct planes number1251:
the496 old body planes plus **755 new planes**. Plane intersections are
preserved; these are counts of plane objects, never disjoint tube volumes.

## Body-body and body-pair classes

Every old body-body plane splits coordinates into two signed constant
blocks. It contains exactly two body lines: a further body vector must
have both block multipliers in{+1,-1}. Thus the old496 pair planes are
distinct.

For p=e_i+s e_j, there are16 body lines v with v_j=s v_i. Then p is
the signed restriction of v to{i,j}, so span(v,p) is the old body plane
with a block of size2. The other16 body lines satisfy v_j=-s v_i;
then v.p=0 and span(v,p) is new.

There are30*16=480 matching body-pair pairs. They duplicate240 old
size2/4 block planes, each with its two body lines and one pair line.
There are another480 orthogonal body-pair pairs, all distinct new planes.
To prove uniqueness, a body vector in span(v,p) must agree with +/-v
outside the two support coordinates; the conditions on the two remaining
entries force its p coefficient to zero. A pair vector in that span must
have v coefficient zero because the other four coordinates would otherwise
all be nonzero. Thus such a plane has exactly one body and one pair line.

## Pair-pair classes

Two different support2 pair lines have three possibilities:

* Same support: the plus/minus pair spans a coordinate2-plane. There
  are15 pairs/planes.
* One shared coordinate: their union has three coordinates. In normalized
  signs they are e_i+s e_j and e_i+t e_k, and their difference is the
  third pair line s e_j-t e_k. The plane has a normal with three+/-1
  entries and has exactly three pair lines. For each of the20 coordinate
  triples there are4 unoriented normals, hence80 planes. Its240 direction
  pairs count each plane three times.
* Disjoint supports: choose four coordinates, one of their three pairings
  and one of four sign choices. This gives15*3*4=180 planes. Each plane
  contains exactly its two pair lines, because a combination with both
  coefficients nonzero has four nonzero coordinates.

These classes have support unions of size2,3,4, so cannot coincide across
classes. None contains a body direction, whose six coordinates are nonzero,
so none equals any plane of the two body-containing classes. Therefore

    pair-pair direction pairs=15+240+180=435,
    new pair-pair planes=15+80+180=275,
    total planes=496+480+275=1251.

The complete pair accounting is496 body-body+960 body-pair+435 pair-pair
=1891. As another consistency identity,240 old planes and80 triangle
planes have three frame lines; all other planes have two. Thus there are
exactly320 dependent triples of distinct frame directions, derived from
the classification without enumerating all triples.

## Why the geometric S2 is the true Green singular set

This statement can be justified afresh for the new frame; it does not
inherit old quantitative barriers. Removing any two new-frame columns
leaves at least30 body columns. Their symbol is bounded below by
(sqrt(32)-2sqrt(6))||xi||>0, with all remaining pair contributions
nonnegative. The remaining sum Y therefore has a continuous full density.
Its orthogonal4-dimensional quotient at a selected pair plane has the
full new symbol restricted to that quotient, hence density q_t with
q_t(0)=t^-4 q_1(0), q_1(0)>0.

Use positive coarea, not conditioning on a null fiber. At a fixed point
x of the plane, retain ||P_UY||<=L. Scaling gives the retained quotient
density t^-4 integral_(||u||<=L/t)r_1(u,0)du. By monotone convergence this
is at least q_1(0)t^-4/2 for all sufficiently small t. The coarea fiber
identity is legitimate by solving an independent four-column quotient
basis and Tonelli; it gives the continuous quotient density version.
On the retained region the two solved pair parameters are bounded, so
their Cauchy factors have product at least c t^2. Retaining the exact
pair Jacobian proves p_t(x)>=c' t^-2 for small t. Thus every K_T and G
is infinite on S2, including S1 and all intersections. No numeric value
of this new barrier constant is asserted.

For a compact set at distance at least delta>0 from S2, let
M=sum_columns||a_i|| and eta=delta/(2M). If all variables larger than
eta had rank at most2, their contribution would lie in an S2 plane,
while the small-variable remainder has norm at most M eta=delta/2.
This is impossible at such an x. The large directions therefore have
rank at least3. Select three independent ones greedily; do not assume
every triple is independent.

For each selected independent triple I, its three Cauchy factors on
|Z_i|>eta are at most t/(pi eta^2). Positive coarea deletion contributes
1/sqrt(det Gram(I)), and the remaining orthogonal3-dimensional density
has peak at most C3 t^-3: its symbol is the full new symbol restricted
to that orthogonal complement. Summing the finitely many independent
triples bounds p_t uniformly as t approaches zero on the compact set.
The positive-time Fourier peak t^-6 handles large t. Hence G is finite
and continuous off S2 by dominated convergence. The same holds for K_T.

This rank-greedy proof preserves all dependent large-direction cases.
Explicit dependencies include a matching body pair and its support2
restriction, and the pair triangle
(e1+e2)-(e1+e3)-(e2-e3)=0. Such triples invalidate an any-three-independent
argument. Numerical g, short-potential minima, quantitative tube radii,
and kernel/source certificates for the new model remain unproved.
