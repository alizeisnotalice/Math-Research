# Free-threshold first variation at the selected halo-hole displacement

This concerns the finite primary-length family{1,2,3,4}⁵ with its actual
mixed satellite captures. It does not bound the full maximal operator.

## Scope correction and normalization

Round001 tested fixed geometry/threshold, h=0, and the common spatial
dilation s=(1+ε)^(1/5). Under a spatial dilation, weak ratios transform by

    s⁵ |{Mμ>s⁵}| / mass(μ).

Thus that second case corresponds exactly to λ=1+ε, or h=1 at first order.
It did not optimize h. More generally s=1+kε gives h=5k to first order,
with the normalized satellite position dilated together with the lattice.

Here the source measure has one primary unit mass and one satellite mass ε
per period; the threshold is λ=1+hε. Write U(ε,z,h) for the finite receiver
union volume per period. The weak ratio is

    R(ε,z,h) = λ U(ε,z,h)/(1+ε).

For a primary block of N distinct atoms with limiting common radius r,
let n_sat be the number of distinct satellite atoms in that same cube.
The equality radius is

    r(ε)=[(N+ε n_sat)/(1+hε)]^(1/5)/2,
    r'(0)=r(n_sat/N−h)/5.

The same radius drives all ten receiver faces. The denominator and threshold
together contribute(h−1)C to R'(0). There is no period-volume divisor.

## A genuine hole cell

The primary endpoint arrangement was searched only to locate one uncovered
open cell, with a deterministic symmetry-first/lexicographic rule. Its tail
indices are(0,0,8,8) and first-coordinate index is5. Each coordinate was
then refined to a satellite-count event interval. The exact endpoints and
centre coefficients are in `selected_cell.json`. The centre is approximately

    z=(0.5370890252561094,
       0.005487417743144651,0.005487417743144651,
       0.3126971573211914,0.3126971573211914).

The entire open fine cell is disjoint from the baseline finite halo. This
is checked with exact box intersections, not inferred from a decimal point.
At the centre, a positive distance separates the satellite from that halo.
The isolated satellite's receiver cube has volume ε/λ; its contribution to
λU is exactly ε for sufficiently small ε. Thus κ=1, independently of h.
This differentiates the present construction from round002's interior-halo
quarter-tail cell, where κ=0.

Only the centre's derivative was computed. No sign assertion is made for
every displacement in the selected open cell.

## Why the threshold response is concave here

Two finite geometric facts were checked exactly on all51 feasible primary
raw boxes:

1. No raw primary receiver box has a zero width.
2. All primary facets with the same oriented coordinate plane have the same
   limiting radius, hence the same N.

These are special properties of this finite baseline, not general facts
about unions of moving boxes. Their exhaustive checks are necessary: a
zero-width box could open a new strip, and competing different-radius faces
could change the outward envelope as h varies.

On an exposed oriented facet, all active blocks therefore have the same
velocity at a given receiver x:

    v(x,h)=r n_sat(x)/(5N)−(r/5)h.

Its integrated contribution is affine in h. Over the complete period put

    B0 = exposed contribution at h=0,
    K  = ∫_(exposed facets) r/5 dH⁴,
    S  = C−K.

Internal contacts must also be included when h>0; the h=0 exposed kernel
alone cannot account for gaps that appear as original boxes shrink. At a
contact plane, suppose the left and right boxes have limiting radii r_L,r_R
and counts N_L,N_R. If a baseline box crosses the plane in its interior,
it blocks any first-order gap. Otherwise let

    α=(r_L+r_R)/5>0,
    β=r_L n_L/(5N_L)+r_R n_R/(5N_R)≥0.

The two normal velocities are v_L and v_R. Their contact gap has first-order
length

    [−v_L−v_R]_+ = [α h−β]_+.

Satellite-count cuts partition the four-dimensional contact into pieces of
exact area A_j, with constant α_j,β_j. Different contact planes meet in lower
dimension and their thin-gap overlaps contribute only O(ε²). The resulting
normalized derivative is

    D(z,h)=B0+1−C+hS−Σ_j A_j[α_j h−β_j]_+.

Every hinge has a nonnegative coefficient. Consequently D(z,h) is concave
and piecewise affine for all fixed real h. This is a shape calculation, not
finite-difference fitting. For any fixed h, λ remains positive when ε is
sufficiently small.

The18 reflection classes use the two equal low-tail and two equal high-tail
coordinates; multiplicities sum to32. Five contact planes per halfcell and
333 nonzero count pieces yield20 distinct exact event ratios β_j/α_j.
Ratios are compared by cross multiplication in Q(u,v); field inversion and
displayed decimals are not used for decisions.

## Global maximizer of the first-order h profile

The exact values have the following decimal displays:

    K ≈1.903763626193357,
    S=C−K ≈0.0905528353751288,
    D(z,0)=B0+1−C ≈−0.1317150987176058.

Let W0=Σ_(β_j=0) A_j α_j. The exact computation gives

    S>0,      S−W0≈−0.263635885566582<0.

This gives an especially short all-h proof. For h≤0, every hinge vanishes,
so D(z,h)=D(z,0)+Sh≤D(z,0). For h≥0, retain only the zero-event hinges:

    D(z,h)≤D(z,0)+(S−W0)h≤D(z,0).

Thus h=0 is the unique global maximizer of this first-order profile and

    sup_(h∈R) D(z,h)=D(z,0)<0.

The adjacent event h≈0.089445370485576 has derivative approximately
−0.15529610817540165. Exactly two event candidates were evaluated; no h
grid was searched. The exact ratio representations, derivative numerators
and isolating intervals are in `free_threshold_results.json`.

This is an all-fixed-h first-order obstruction at one displacement. It is
not a bound for all finite ε, all thresholds λ, ε-dependent unbounded h,
all satellite displacements or the entire periodic maximal halo. For each
fixed h it implies a sufficiently small positive ε range with R<C; on
compact h ranges the finite-family expansion is uniform.

## Verification status

The exposed kernel is the previously validated finite-face implementation.
All contact pieces were integrated by recursive slicing. A separately coded
coordinate-mask intersection engine recomputed the333 nonzero areas
exactly. At the selected positive h event,333 literal one-dimensional
interval unions also reproduced the gap deficits, and the full normalized
derivative numerator agrees coefficient by coefficient. The baseline/contact
geometry is shared; this is explicitly a second integration-engine and gap
check, not independent geometric production.

The attempted direct five-dimensional mixed-receiver jet union reached its
registered120-second checkpoint in the first reflection, session5003,
exit status1. It produced no verification result. No other reflection was
started. `direct_jet_checkpoint.json` records this incomplete attempt and
`direct_jet.py` retains the source. No resumable in-memory union cache was
saved. A later registered chunked computation may resume this audit task;
the120-second checkpoint is not a permanent total-compute ceiling. The
failed attempt is not marked as independently passed.

## Finite mass and smooth L1 transfer

For fixed positive ε and λ=1+hε>0, any mixed witness uses the exact distinct
mass N+ε n_sat and common equality radius above. Shrink each positive-volume
receiver box by τ>0. Then choose a slightly smaller common radius retaining
all selected atoms, so the average is strictly greater than λ. Replace each
atom by a nonnegative smooth bump of the same mass with sufficiently small
support; the strict common-radius margin persists.

For M⁵ receiving periods retain both source orbits with indices
{-5,…,M+5}⁵. Their exact mass is(1+ε)(M+11)⁵, and the shrunken receiver
volume is M⁵ U_(ε,τ). The legal finite-L1 weak ratio is at least

    [λ U_(ε,τ)/(1+ε)] [M/(M+11)]⁵.

Let M→∞ and then τ↓0. All witness radii are among the allowed positive
common scales, and the background remains Lebesgue measure. This proves
lower-bound transfer for any successful finite candidate. No successful
candidate occurred here, and no new C5 record is declared.
