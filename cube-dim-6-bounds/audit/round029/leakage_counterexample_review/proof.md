# Independent proof audit: a genuine Lebesgue counterexample to L6

The exact L6 statement is frozen in round025/theorem.md, Section4. It says
that for every f>=0 in L1 and lambda>0 one can choose countably many
pairwise disjoint positive-radius centered axis-parallel cubes Q_i such
that integral_Q_i f>=lambda|Q_i|, and, writing g_i=f1_Q_i,

    L_i={x:M_c g_i(x)>(1-theta)lambda},
    F={M_c f>lambda} minus union_i L_i,

one has lambda|F|<=(43/62500) integral f. Here theta=43/4000000, and M_c
uses one common positive radius for all six coordinates, at every scale.
We prove the asserted existential selection is impossible for one fixed f.

## Fixed input; no atoms

Set lambda=1, d=1/1000000, a=19/20, b=1-d, h=1-theta, and

    f(x)=(1/2)1_[-2,2]^6(x) + [32/(2d)^6]1_[-d,d]^6(x).

This is a nonnegative bounded compactly supported Lebesgue L1 function.
Its total mass is(1/2)4^6+32=2080. The spike has mass32 and positive
six-dimensional volume(2d)^6. No limit to a Dirac measure is taken.
Boundary conventions for cubes do not affect any integral.

## Every admissible selected family obeys the same mass restriction

Fix ANY countable family of positive-radius cubes Q_i=Q(c_i,R_i),
pairwise disjoint up to null sets, with integral_Q_i f>=|Q_i|. Put

    m_i=integral_Q_i [32/(2d)^6]1_[-d,d]^6.

Since the background density is at most1/2 everywhere,

    m_i >= |Q_i|-(1/2)|Q_i intersection [-2,2]^6|
         >= (1/2)|Q_i| =32 R_i^6.

This argument fully permits a cube to capture only part of the spike.
Countable additivity/Tonelli and disjointness give sum_i m_i<=32. Thus
R_i<=1 for every i. No cube need contain the origin, the entire spike,
or the entire background. If zero-radius cubes were permitted, their
components vanish a.e. and can be discarded.

Call a component big when m_i>112/5. There is at most one big component,
since2(112/5)>32. Every other component, including equality, is small.
The family may be empty, finite, infinite, or leave mass unselected;
the argument covers all of these possibilities.

## A slab of bad centers, chosen after the family

If there is a big component, choose the sign of x1 opposite to c_(big),1.
If that center coordinate is zero, choose either side. If no big component
exists, choose either side as well. By reflection symmetry it suffices to
write the negative slab

    S={x: -b<x1<-a, and |x_j|<b for j=2,...,6}.

In the big case we now have c_(big),1>=0. The slab has volume

    V=(b-a)(2b)^5.

The sign may depend on the selected family, but f and V do not. This is
consistent with the required quantifiers: for every family we will find
a subset of its own leakage set with the same lower volume.

For x in S let r=||x||infinity+d. Then0<r<b+d=1. The centered cube Q(x,r)
contains the entire spike because |x_j|+d<=r. It lies strictly inside
[-2,2]^6 since |x_j|+r<=2||x||infinity+d<2b+d<2. Its f-average is exactly

    1/2 +32/(2r)^6 =1/2 +1/(2r^6)>1.

Thus S is contained in E={M_c f>1}. The common radius and strict level-set
condition are both satisfied at every point of the open slab.

## Uniform exclusion from every component maximal set, at all radii

Put r0=a-d>0. For every x in S and0<r<r0, Q(x,r) misses the spike:
the first-coordinate distance to[-d,d] exceeds a-d. Therefore every
component g_i has average at most1/2<h at these radii.

For r>=r0 and a small component, its background average is at most1/2,
and its spike contribution is at most m_i/(2r)^6. Consequently

    average_Q(x,r) g_i
       <=1/2 +(112/5)/(64r^6)
       <=1/2 +(7/20)/r0^6=:S0<h.                    (1)

This bound is uniform in i and x, so countably many small cubes cannot
remove any point of S from the leakage. It includes partial-spike cubes.

It remains to treat the unique possible big component. Since c1>=0 and
R<=1, its first-coordinate interval has lower endpoint c1-R>=-1. The
first-coordinate upper endpoint of Q(x,r) is x1+r<r-a. Therefore the
intersection of these two first-coordinate intervals has length at most
r+1-a. This number is positive for r>=r0. If the intervals do not meet,
the intersection length is zero and the same upper bound still holds.
In the other five coordinates the intersection length is at most2r.
Restriction to the background cube can only reduce the volume. Hence

    background average of g_big
       <= (1/2)(r+1-a)(2r)^5/(2r)^6
       =1/4 +(1-a)/(4r).

The spike average is at most32/(2r)^6=1/(2r^6), regardless of how much
spike the big component actually captured. Thus

    average_Q(x,r) g_big
       <=1/(2r^6)+1/4+(1-a)/(4r)
       <=1/(2r0^6)+1/4+1/(80r0)=:B0<h.             (2)

The last inequalities in(1),(2) follow because their r-dependent positive
terms decrease for every r>0. The strict rational comparisons S0<h and
B0<h were independently verified by Fraction arithmetic, not by sampling
r. Together with the r<r0 case they control ALL positive radii. In
particular sup_r average g_i<h for each i at each x in S, because the
common uniform upper bound max(1/2,S0,B0) is strictly below h.

Thus S is disjoint from every L_i and S subset F for this arbitrary family.
No assumption on how the family was selected entered the proof.

## Strict violation and disposition

The independently checked exact values are

    V =49998750005499989500010249994950001
       /31250000000000000000000000000000000,

    (43/62500)*2080=4472/3125,

    V-(43/62500)*2080
      =5278750005499989500010249994950001
       /31250000000000000000000000000000000 >0.

Therefore EVERY admissible family has leakage exceeding the L6 allowance.
This proves L6 false for a genuine Lebesgue L1 input, with the prescribed
strict maximal level sets and all positive common radii.

Only L6 is refuted. The round025 support-conditioned corner lemma remains
valid, and can still be applied to each g_i individually. It cannot charge
the combined influence of the background and spike; this example exhibits
exactly that gap. The formal implication L6=>C6<=64-43/31250 remains a
correct conditional implication with a false premise, so it gives no bound.
The current general C6 upper bound remains64. Nothing here proves that
all other covering improvements, or all other positive-kernel routes, fail.
