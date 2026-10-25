# Full-displacement insertion kernel for the registered finite family

All coordinates below are physical coordinates. Let a be the exact standing
C lattice spacing vector, P=∏[0,a_j], and

    μ_(ε,z) = Σ_(k∈Z⁵) δ_(ak) + ε Σ_(k∈Z⁵) δ_(ak+z),   z∈P.

Coincident orbits, when z is a lattice point, are aggregated. Otherwise the
orbits are distinct and the exact period mass is1+ε. This is one common
periodic source measure. Overlap among different receivers never adds source
mass, and receiver volumes are always union volumes.

The finite family retains every primary consecutive coordinate block with
length vector in{1,2,3,4}⁵, plus mixed witnesses with satellite atoms captured
near its limiting radius, plus the isolated single-satellite witness.
For fixed sufficiently small ε, this is a finite family modulo periods:
all relevant atoms lie in a fixed collar and satellite coordinate counts
are finite. For a prescribed receiver, nearest consecutive blocks dominate
other consecutive blocks of the same count. The calculation does not claim
to include longer primary blocks or the complete periodic maximal halo.

Write U_ε(z) for its received volume per full period at threshold one, with
boundaries immaterial, and R_ε(z)=U_ε(z)/(1+ε). At ε=0 the original full
1024-block family gives U_0=C exactly. Hence, if B(z) is the first-order
mixed-boundary contribution and κ(z) the isolated-satellite contribution,

    R_ε(z) = C + ε D(z) + o(ε),
    D(z) = B(z) + κ(z) - C.

In particular, the denominator contributes−C. Dividing by the geometrical
period volume would be an erroneous second normalization.

## Common-radius boundary velocity

For a primary block W of N distinct unit atoms, r=N^(1/5)/2. At a limiting
receiver x let

    n_z(x,r)=∏_(j=1)^5 #{k∈Z : |a_j k+z_j-x_j|≤r}.

This counts distinct satellite atoms in the same closed five-dimensional
cube. A mixed witness has mass N+ε n_z, so its common radius satisfies

    r(ε)=r+ε n_z r/(5N)+O(ε²).

Every one of its ten faces moves with this same scalar radius variation;
there are no independently chosen coordinate radii. On a primary receiver
facet, the outward normal velocity is therefore

    v_W(x;z)=n_z(x,r) r/(5N).

At a facet where several primary witnesses have the same boundary plane,
the outer union envelope takes max_W v_W. Where an original positive-volume
box already covers the outward side, the contribution is zero. The code
implements these two facts as a weighted maximum integral with baseline
blockers. It starts from all51 feasible or degenerate raw primary boxes,
not only the42 positive-volume retained boxes. A box degenerate in exactly
one direction can open into a first-order strip and supplies both outward
faces. Degeneracy in at least two directions contributes only O(ε²).

Let F range over these oriented facets after the outward-covered parts are
removed. In one reflected half-period,

    B_H(z)=Σ_F ∫_F max_(W active on F) v_W(x;z) dH⁴(x).

Summing B_H over the32 coordinate reflections yields B(z). For the tested
one-tail slice only the sign of z_2 changes, so

    B(0,t,0,0,0)=16[B_H(0,t,0,0,0)+B_H(0,−t,0,0,0)].

Facet intersections and changes of satellite count within a facet have
codimension at least one inside the facet; their associated overlap of thin
normal strips is o(ε). Coincident entire facets are handled by the maximum,
not by summing their speeds. At aligned satellite boundary atoms, closed
cube counts are appropriate: the stationary cube endpoint stays attached
to the primary extreme and the opposite endpoint expands outward.

The usual finite-box symmetric-difference estimate justifies the expansion:
replacing each analytic radius by its value and linear term moves each face
by O(ε²); the volume error is bounded by the sum of these moves times bounded
cross-sectional areas. This also handles same-position/same-speed endpoint
ties whose first separation occurs only at higher order.

## Isolated satellite contribution

A satellite alone has a receiver cube of volume ε. If z lies in the interior
of the original finite halo, κ(z)=0; if outside its closure, κ(z)=1. On a
halo boundary, κ(z) is the fraction of the centered infinitesimal cube lying
outside the local union of coordinate orthants. It is a multiple of1/32.
The overlap with the O(ε)-thick mixed-boundary strips is o(ε), because the
isolated cube has side ε^(1/5). These statements concern the original finite
halo and include periodic copies near the period boundary.

For the registered slice interval,0<t<1/2 and all other coordinates vanish.
Thus z lies strictly inside the primary single-atom receiver cube and
κ(z)=0 throughout the interval. The pure satellite halo is entirely covered
for sufficiently small ε; no nonanalytic ε^(1/5) endpoint is omitted from
the derivative calculation.

## Displacement event cells

For an integration coordinate x_j, satellite counts change only at

    x_j=a_j k+z_j±r.

The original receiver facets and blockers have constant coordinate endpoints
e. The ordering changes only when a moving count cut meets such an endpoint:

    z_j=e±r−a_j k.

Within one reflection, all moving cuts in the same coordinate have the same
slope in z_j; their mutual ordering is constant. Reflection replaces z_j by
−z_j and is treated in a separate integral. Add the original halo endpoint
hyperplanes z_j=e to track κ. These coordinate hyperplanes partition the
full displacement cell into finitely many open rectangular cells. The
integer range−4,…,4 is safe here because r≤2 and every a_j>1.

Inside an open cell all count values, blocker masks, max-speed choices and
cut orderings are fixed. Each integration interval length is affine in its
own z_j. A facet integral uses only four integration coordinates, while its
normal-coordinate count is constant on the cell. Therefore B(z), and hence
D(z), is multiaffine in the displacement coordinates, with total degree at
most four; κ is constant there. This is a structural conclusion from the
endpoint arrangement, not an interpolation assumption. No sign is claimed
for the uncomputed cells or for boundary strata.

On the selected slice z=(0,t,0,0,0), fixed zero coordinates are retained as
aligned boundary strata, with closed count rules. Their combinatorics do
not vary with t. The same argument shows D(t) is affine on each open event
interval. The implementation includes all145 raw endpoint alternatives and
36 possible radii, including alternatives from infeasible blocks, so the
1196-cut arrangement is a safe refinement. It may split unnecessarily but
cannot omit a displacement event of the tested family.

## Exact selected interval

Let u⁵=2,v⁵=3. The interval containing the quarter-tail displacement a_2/4 is

    L=(-v+u v²)/2,
    U=(1+u+u²−u³v³)/2.

The code registers three physical positions L+θ(U−L), θ=1/4,1/2,3/4, then
computes their derivatives in Q(u,v). Their exact midpoint identity holds.
Since affine structure was already proved by the arrangement, the limiting
endpoint values are2D_(1/4)−D_(1/2) and2D_(3/4)−D_(1/2). Both have strictly
negative rational isolating intervals. Thus D(t)<0 for every L<t<U.

This proves a local obstruction for this finite witness functional at each
fixed t in the open interval. On any compact subinterval the expansion and
negative bound are uniform, giving an unquantified common ε_0>0 with
R_ε(t)<C for0<ε<ε_0. It does not include the two endpoint strata, all satellite
displacements, finite large ε, longer primary blocks or the full maximal
operator. A negative finite-family derivative is not an operator upper bound.

## Transfer to genuine L1 inputs

For any fixed small positive ε and finite witness union, first shrink every
positive-volume receiver box inward by a common τ>0. A corresponding radius
r_W−τ/2 still captures every selected atom strictly and satisfies
(2(r_W−τ/2))⁵<N_W+ε n_W. Replace each atomic mass by a smooth nonnegative bump
of the same mass with support radius less than τ/4; these same common radii
still capture the selected bumps with average strictly above one.

Use M⁵ receiving periods and both source orbits with indices in
{-5,…,M+5}⁵. This fixed collar safely contains the registered nearby witnesses.
The source mass is exactly(1+ε)(M+11)⁵, the shrunken received volume is
M⁵ U_(ε,τ), and the weak ratio is at least

    [U_(ε,τ)/(1+ε)] [M/(M+11)]⁵.

Let M→∞ and then τ↓0. This gives C₅≥R_ε(z) for the finite witness ratio,
using fixed Lebesgue measure, smooth compactly supported L¹ inputs and
admissible positive common radii. It does not turn the local negative result
into an upper bound. No new lower record is asserted in this round.
