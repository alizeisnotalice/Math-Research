# D5 full-five-axis parity deletion: finite-family screen

This round is distinct from the sibling's two-coordinate parity weights and
from the previous BCC two-orbit point set. Here the unit atoms are

    P={ (a1*j1,...,a5*j5) : j in Z^5, sum(j) even }.

The rational spacings and six candidate names are frozen in registration.
The period lengths are 2*a_i. Each period has 16 residue vectors modulo two
of even parity, hence mass 16 and no repeated atom. Reflection of any one
coordinate j_i -> -j_i preserves sum parity. The positive box D=product
[0,a_i] therefore has 32 reflected copies filling one centered period cell,
with only null overlaps, and the finite ratio is 32V/16=2V.

For each coordinate choose a consecutive integer block with (length,start)
from ((1,0),(1,1),(2,0),(3,-1),(3,0),(4,-1)). The registered family has
6^5=7776 configurations. Only atoms whose complete five-index sum is even
are selected. If lengths are l_i and starts k_i, their exact number is

    N=(product(l_i)+(-1)^sum(k_i)*product(1_{l_i odd}))/2.

This follows by writing the parity indicator as (1+(-1)^sum(j))/2 and
factoring the rectangular sum. In one coordinate the alternating sum is
(-1)^k if the length is odd and zero if it is even. N=0 configurations are
discarded. Every nonzero configuration uses the common radius

    r_exact=N^(1/5)/2

in all five coordinates. Its receiver range is the product of
[max_selected p_q-r,min_selected p_q+r], intersected with D. If another
axis has length at least two, both parity choices on that other axis exist,
so both endpoints of the q-axis interval occur in selected points.
If all other lengths equal one, the q-axis interval is trimmed to its
first and last permitted parity positions. These rules compute the true
selected-point extrema rather than those of the unfiltered full rectangle.

For every retained box the source independently enumerates all even-sum
integer tuples, verifies their distinct cardinality N, takes extrema from
the tuples, and verifies all ten receiver endpoint inequalities. Strict
inner side s satisfies s^5<N<=(s+2^-80)^5. Direct guards and both exact
union engines pass for all six jobs.

To enclose the same finite family at exact equality radii, the radius error
is at most delta=1/(2*2^80). All six receiver domains have side lengths at
most two. Each configuration's new box volume is at most its old volume
plus 160*delta, including a box that was previously empty: if an old side
was nonpositive, its new positive width is at most 2*delta. For nonempty
boxes, telescope the five side-length products. Summing over all 7776
configurations and multiplying by the ratio factor two gives

    ratio_inner <= ratio_exact_family
      <= ratio_inner + 160*7776/2^80.

These exact rational upper bounds were actually computed and compared with
the complete frozen C_cert lower isolating endpoint. Every one is strictly
below C_cert. The conclusion excludes only these six rational spacing
vectors and these 7776 chosen configurations. It is not a global D5 bound
and does not exclude other blocks, other parameters or other point sets.

The strict inner model also has an explicit ordinary-L1 transfer. Selected
indices initially lie in [-1,2], and reflection puts them in [-2,2]. Use
centered period representatives j in {-1,0}^5 of even sum, at period offsets
2*q. For a receiver period q, each selected point has period offset at most
one in each direction because floor((j_i+1)/2) lies in [-1,1]. Thus K^5
interior periods need only the input cells q in {-1,...,K}^5. Their mass is
exactly 16(K+2)^5. For the finitely many inner witnesses choose a common
eta>0 such that (2(r_inner+eta))^5<N. Replace each atom by a smooth unit-mass
bump of sup-norm support radius smaller than eta/2. The same receiving
boxes still have centered cube averages strictly greater than one, using
one common radius r_inner+eta per witness. Their volume is at least
32 K^5 V. The L1 ratio is at least

    2 V * (K/(K+2))^5.

Let K tend to infinity. This validates only the constructed inner lower
ratio for the general-L1 all-positive-scale centered cube operator. Since
no candidate beats the frozen C_cert, no formal lower-bound promotion is
made and no new-breakthrough continuation counter is triggered here.
