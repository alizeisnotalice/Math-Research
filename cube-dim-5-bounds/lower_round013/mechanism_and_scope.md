# Round 013: body-centered five-axis periodic witnesses

The frozen reference throughout this round is the complete algebraic
C_cert in `../lower_round003/a2_equals_a3_results.json`. A rational display
threshold is not the reference. This round gives no improvement and makes no
change to the formal bound ledger.

## New geometry and finite formulas

For positive rational spacings a=(a1,...,a5), let

    P0={(a1*j1,...,a5*j5): j in Z^5},
    P1={(a1*(j1+1/2),...,a5*(j5+1/2)): j in Z^5}.

Put one unit atom at every point of P0 union P1. The two orbits are disjoint:
in half-spacing units the first has five even coordinates and the second
has five odd coordinates. Thus one half-open period cell product[0,ai)
contains exactly two distinct atoms and has mass two. The simultaneous
five-coordinate half shift gives a body-centered geometry, rather than the
previous row-staggered pattern or an alternating one-axis modulation.

Use D=product[0,ai/2]. For orbit 0, the selected l consecutive coordinate
indices are centered about 0 when l is odd and about 1/2 when l is even.
For orbit 1, physical half-integer coordinates are centered about 1/2 when
l is odd and about 0 when l is even. Explicit integer index ranges are in
`indices(orbit,l)`; none uses rounded geometry.

A single-orbit block has N=product(li) unit atoms. A cross-orbit block is the
disjoint union of two rectangular blocks with lengths l0,l1 and has

    N=product(l0_i)+product(l1_i),
    r_exact=N^(1/5)/2.

One centered radius r is used in all five coordinates. Its receiver box is

    product_q [max_selected p_q-r, min_selected p_q+r], intersected with D.

Every retained block is independently expanded into its actual distinct
half-integer coordinate tuples. The program checks cardinality exactly N,
then takes extrema from those enumerated tuples and checks all ten box
endpoint inequalities. These checks do not reuse the construction's extrema
table. Strict inner radius uses s_lo=m/2^80 with

    s_lo^5<N <= (s_lo+2^-80)^5,  r_lo=s_lo/2.

Thus each inner box is a true strict-level witness after the finite smoothing
step stated below. Empty boxes and boxes contained in another may be removed
without affecting the union volume. Two exact union engines, recursive
slicing and coordinate-mask intersection, agree in every executed candidate.

Reflection x_q -> a_q-x_q preserves both full periodic orbits separately.
The 32 reflected copies of D fill a period cell with only null-face overlaps.
If the inner union volume in D is V, its ratio is therefore 32V/2=16V.
There is no division by the geometric volume of the period cell.

## Registered scope and finite-family outer bound

Exactly six rational spacing vectors were registered before execution. Each
uses all 2*4^5=2048 single-orbit blocks, and all cross-orbit choices having
l0_i,l1_i in {1,...,4} and |l0_i-l1_i|<=1 in each coordinate. There are ten
per-axis pairs and 10^5 cross-orbit choices: M=102048 configurations total.
This is a finite family; no claim is made that it exhausts a periodic
superlevel set or the full BCC parameter space.

The rational inner radius differs from r_exact by at most

    delta=1/(2*2^80).

For each of the six domains every side length of D is at most 2. Replacing
r_lo by r_exact moves each unclipped box face outwards by at most delta.
For a nonempty original box, telescoping the five products of side lengths
bounds its added volume by 5*(2*delta)*2^4=160*delta, since every clipped
side remains at most 2. If the original box is empty but its enlargement is
nonempty, an originally nonpositive side becomes at most 2*delta, so the
whole enlarged volume is at most (2*delta)*2^4, also bounded by 160*delta.
Thus the finite union at the exact equality radii has the rigorous enclosure

    ratio_inner <= ratio_exact_family
        <= ratio_inner + 16*M*160*delta
         = ratio_inner + 1280*102048/2^80.

This accounts for all original configurations, including empty, duplicate
and dominated boxes. The exact upper endpoint was actually computed with
Fraction for each candidate and compared to the full C_cert lower isolating
endpoint. All six upper endpoints are strictly below it. Consequently these
six points and this finite family are excluded, rather than merely failing
because a dyadic radius was slightly too small. This argument gives no
exclusion of omitted blocks, other spacings or other periodic patterns.

## General L1 and all positive scales

The finite input can be transferred without changing the cube convention.
Every selected orbit index has magnitude at most two for lengths up to four.
After any of the five receiver reflections the orbit-0 index becomes 1-j
and the orbit-1 index becomes -j, which also stay in [-2,2]. For interior
period indices q in {0,...,K-1}^5, include both atoms from each input cell
q in {-2,...,K+1}^5. All selected atoms for every translated, reflected
witness are present. The resulting atomic mass is exactly 2(K+4)^5.

There are finitely many retained inner witnesses and each has (2r_lo)^5<N.
Choose a single eta>0 so that (2(r_lo+eta))^5<N for every retained witness.
Replace every input atom by a smooth nonnegative mass-one bump supported
in a sup-norm ball of radius less than eta/2. At every receiver in the
original inner box, the cube of common radius r_lo+eta contains all selected
bumps in their entirety, and its average is strictly larger than one.
No inward receiver shrink is needed because the dyadic radii already have
strict mass slack. The reflected receiving volume is at least 32 K^5 V.
Thus the ordinary finite L1 weak ratio is at least

    16 V * (K/(K+4))^5.

Let K increase without bound. This proves validity of the inner ratio for
the general L1 centered-cube operator over all positive scales. It gives a
lower bound only; the finite-family upper enclosure above is not a general
operator upper bound. No winning certificate was produced in this round.
