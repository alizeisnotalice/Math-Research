# A bounded Lebesgue L1 counterexample to the proposed leakage lemma L6

This concerns only the UNPROVED selection target L6 of round025. It does
not contradict the proved single-support-cube result and gives no new C6
bound. All cubes below are centered, axis parallel, with one common radius;
all maximal operators include every positive radius.

Put d=1/1000000, a=19/20, b=1-d, theta=43/4000000 and h=1-theta. At
lambda=1 define the ordinary bounded, compactly supported L1 function

    f = (1/2) 1_[-2,2]^6 + [32/(2d)^6] 1_[-d,d]^6.

Write the second summand as s. Its mass is32, whereas the background mass
is2048. Thus M=2080. In particular this is not a point mass or a different
background measure. The containing support cube has volume4096>M, so the
support-compatible-height hypothesis of the round025 corollary fails here.

We prove that EVERY countable family of cubes Q_i, disjoint up to Lebesgue
null sets, with integral_(Q_i) f >= |Q_i| violates L6. Define

    g_i=f 1_(Q_i), L_i={M_c g_i>h}, E={M_c f>1},
    F=E minus union_i L_i.

## 1. Radius and mass restrictions hold for every allowed selection

Let R_i be the radius and c_i the center of Q_i, and put m_i=integral_(Q_i)s.
Since the background is everywhere at most1/2,

    |Q_i| <= integral_(Q_i) f <= |Q_i|/2+m_i.

Consequently m_i>=|Q_i|/2=32 R_i^6. Since m_i<=32, every R_i<=1.
Countable disjointness and s in L1 give sum_i m_i<=32. There is therefore
at most one index i with m_i>112/5=(7/10)*32. Call it big if present; all
others are small. This argument permits cubes that capture only part of
the spike, arbitrarily small cubes, and countably many cubes. There is no
atomic-boundary double counting.

## 2. A full slab of genuinely bad centers

If there is a big cube with first center coordinate c_(i,1)>=0, use

    S={-b<x1<-a, |xj|<b for j=2,...,6}.

If its first coordinate is negative use the reflected positive slab;
if no big cube exists choose either slab. In every case

    |S| = V=(b-a)(2b)^5.

For x in S choose the common radius r=||x||_infinity+d. It is strictly
less than1, its cube contains the entire spike support, and
||x||_infinity+r<2b+d=2-d<2, so its cube is inside the background cube.
The average of f on this witness is exactly

    1/2 + 32/(2r)^6 = 1/2+1/(2r^6)>1.

Thus S is contained in E. These are positive-volume bad centers with
strict levels; no finite-scale replacement of M_c is used.

## 3. No selected component is bad on that slab, at any radius

Put r0=a-d=949999/1000000. If a test cube about x in S meets the spike
in positive volume, its radius r is at least r0. At smaller radii the
component average is at most1/2<h.

For a small component and r>=r0, its background average is at most1/2
and its spike contribution is at most m_i/(2r)^6. Hence, for ALL radii,

    M_c g_i(x) <= max(1/2, 1/2+(7/20)/r0^6) =: U_small < h.       (1)

For the big component suppose first c_(i,1)>=0, so the negative slab
was chosen. Since R_i<=1, the lower endpoint of its first coordinate
interval is at least-1. A test cube about x has upper first endpoint
x1+r<r-a. The length of intersection in that coordinate is at most
r+1-a. In each of the other five coordinates it is at most2r. Thus
the background contribution to the normalized test average is at most

    (r+1-a)/(4r)=1/4+(1-a)/(4r).

The spike contribution is at most32/(2r)^6=1/(2r^6). Both terms are
decreasing for r>0. For r>=r0 their sum is therefore at most

    U_big = 1/(2r0^6)+1/4+1/(80r0) < h.                         (2)

The positive slab and c_(i,1)<0 case follows by reflection. For r<r0
the earlier spike-free bound applies. Consequently the ENTIRE slab S
is disjoint from EVERY L_i, including for an infinite selected family.
The bounds do not require the selected cube to contain the whole spike.

## 4. Strict quantitative violation

All inequalities in (1)-(2) and the following volume comparison are exact
rational inequalities, checked in verify.py:

    V=(1-10^-6-19/20)[2(1-10^-6)]^5
       > (43/62500)*2080 = 1.43104.

As S is contained in F, lambda|F|>=V>(43/62500)M. This contradicts L6,
which requires some allowed family with lambda|F|<=(43/62500)M.
The contradiction is uniform over all its allowed selections.

The missing foreign-background contribution is therefore not merely an
unpaid step that this particular L6 could resolve: L6 as stated is false.
The conditional prospective number64-43/31250 must not be registered as
a C6 upper bound. A revised localization statement would need different
assumptions or accounting and a new proof. The actual C6 bounds are unchanged.

The finite arithmetic checker supports the explicit rational margins and
symbolic monotonicity only; the selection quantifiers and maximal-function
implications are proved analytically above. Independent Astra review is
recorded separately, without changing this proof or its executed hash.
