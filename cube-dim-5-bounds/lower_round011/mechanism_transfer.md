# Corrected staggered-row witnesses and all-scale L1 transfer

Round 010's numerical improvement is withdrawn: the old y interval was the
union of the row center ranges rather than their intersection. The original
source, result and lock remain intact; `../lower_round010/withdrawal.json`
records the hashes. Only this round's corrected results may be interpreted
using this document.

Put u=2^(1/5), v=3^(1/5), a=1+u^3/2, s=(1+u)/2 and h=tau*a.
The periodic unit atoms are

    (j*a+h*(k mod 2), k*s, m3*s, m4*s, m5*s),
    j,k,m3,m4,m5 in Z.

The period lengths are T=(a,2s,s,s,s). Its half-open period cell has exactly
two atoms, represented by b0=(0,0,0,0,0), b1=(h,s,0,0,0). Because 0<h<a,
there is no collision: different k have different y, and for fixed k the
remaining integer indices have injective coordinates. A block indexed by
l=(l1,...,l5) has exactly N=product(l) distinct unit atoms.

Select l1 consecutive j values starting at j, l2 consecutive k values
starting at k, and consecutive tail integers centered about 0 for odd li,
or about 1/2 for even li. The first-coordinate extrema explicitly account
for the parity of every selected row. With the single common radius
r=N^(1/5)/2, the exact receiver range in every coordinate is

    [maximum selected coordinate - r, minimum selected coordinate + r].

In particular, y lies in [(k+l2-1)*s-r, k*s+r]. Intersect these ranges with
D=[0,a] x [0,2s] x [0,s/2]^3. The source directly enumerates every selected
atom of every retained box, verifies distinct cardinality N, and verifies
both receiver endpoints against every selected coordinate. Thus it guards
the geometric premise independently of the endpoint construction. Contained
boxes and duplicate boxes may be removed without altering the union.

Let V be the exact absolute volume of this finite union. Reflection
x_i -> s-x_i for i=3,4,5 preserves the periodic input and gives eight copies
inside one full period cell, overlapping only on null faces. The union has
volume 8V and the period mass is 2. Hence the finite periodic ratio is 4V.
No additional division by period-cell volume is appropriate.

For strict level 1, fix epsilon>0 smaller than half the width of every
retained box and smaller than every retained radius. Shrink each box by
epsilon on all ten faces. For any receiver in a shrunken box, every selected
atom lies at sup-norm distance at most r-epsilon. Choose the common cube
radius r'=r-epsilon/2. This is one radius for all five coordinates, and
N/(2r')^5>1. After replacing each atom by a nonnegative smooth bump of mass
one supported in a sup-norm ball of radius delta<epsilon/4, the same cube
contains every selected bump in its entirety. Its average remains strictly
larger than one. Reflections and translations preserve these margins.

Here is an explicit collar rather than an unspecified asymptotic boundary.
Use interior period indices q in {0,...,M-1}^5 and include the two atoms at
q*T+b0 and q*T+b1 for all q in {-10,...,M+9}^5. Before reflection every
registered witness has j,k in [-6,6], l1,l2<=4, so its first period index j'
is in [-6,9] and second period index floor(k'/2) is in [-3,4]. Tail indices
for li<=4 are in [-1,2]. Tail reflection replaces m by 1-m, which still
lies in [-1,2]. Therefore this collar of width c=10 contains every selected
atom used by all reflected witnesses in every interior period. These
period representatives are unique, so the bump input has L1 mass exactly
2(M+20)^5, including bumps whose support crosses the bounding rectangle.
There is no truncation of the bumps at that rectangle.

Let V_epsilon be the volume of the union of the shrunken boxes in D.
The strict level-1 superlevel set of this finite smooth L1 input has volume
at least 8 M^5 V_epsilon. Its weak ratio is at least

    4 V_epsilon * (M/(M+20))^5.

First take arbitrarily large M. Then let epsilon tend down to zero and
choose delta<epsilon/4 at each stage. A finite union of nondegenerate boxes
has V_epsilon -> V; all zero-width faces have measure zero. Thus the best
weak-(1,1) constant for the centered cube maximal operator on general L1
functions over R^5, with all positive scales, is at least the certified 4V.
This transfer proves validity of a correctly checked finite union; it does
not prove that this finite witness family exhausts a periodic superlevel
set and does not convert a ratio below the standing bound into improvement.
