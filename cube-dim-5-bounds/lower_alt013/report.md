# lower_alt013: two-axis parity weights and holes

No improvement over the complete standing C is obtained. Four pre-registered
periodic candidates were actually evaluated with Python 3.15.0, using exact
rational inner and outer receiver boxes. Their *complete registered finite
block families*, including the ideal critical-radius version between those
boxes, are strictly below C. This is not an exclusion for the entire
periodic inputs at all scales, nor for general L1 inputs.

The clone is independent and only `lower_alt013/` contains new research files.
No commit or push was made. The initial local clone had an incomplete checkout
because unrelated historical blobs were unavailable; the project files used
here were readable, and the own-clone index was reconstructed with `read-tree`
without rewriting research files. No main-clone content was changed.

## Prior work and bounded screening

The lower004--009 retention index and their reports were consulted: common-tail
perturbations, length-five extension, first-axis alternating weights, local
common-tail event intervals, symmetric two-group tails, and alternating gaps
are already unfavorable. lower012's two-atom argument requires genuine
interaction beyond isolated pairs. This round changes dependence to two-axis
parity, retaining cross-group witnesses; it is not a finer first-axis scan.

Three mechanisms were registered, at most four actual inputs:

1. Two-axis parity weights, with the parity involving either coordinates
   (1,2) or (2,3), theta=1/8: two actual candidates.
2. Delete one parity class by theta=1, for the same two pairs: two actual
   candidates. Surviving atoms have mass two; absent sites have no mass.
3. Multiscale overlay: screen only, with no fifth input. If each selected
   witness uses mass from just one layer, its union is contained in the union
   of the layer witness sets. Hence its normalized certified volume is at
   most the mass-weighted maximum of the component certified ratios. A win
   from overlay requires mixed-layer capture or a component already winning;
   adding layer copies alone supplies no additive lower-bound gain. This
   statement concerns these selected layer witness sets, not the unknown full
   maximal level sets of arbitrary layers.

`registration.json` predates all four inner runs. `outer_registration.json`
records the subsequent, pre-run decision to bracket those *same four inputs*
with outer radii. No additional candidates were introduced.

## Exact geometry, mass and complete finite family

Use the exact rational spacings already constructed in lower005:

    a1 = 1932785622421 / 1099511627776,
    a2 = ... = a5 = 1181259412951 / 1099511627776.

At the distinct lattice point p_j=(a1*j1,...,a5*j5), assign

    w_j = 1 + theta*(-1)^(j_p+j_q).

The period is 2*a on the two active axes and a on each other axis, with
mass four per period cell. Reflection in each coordinate about zero leaves
the weighted input unchanged. Thus the reflection-reduced receiver domain
has upper coordinate a_k for the two active coordinates and a_k/2 otherwise.
Its full period cell consists of 32 copies, so its normalized level-one
witness ratio is **8*union_volume**, with no cell-volume denominator.

For every length vector l in {1,2,3,4}^5, enumerate both active starts in
{-3,...,3}. Inactive axes use their nearest canonical block: start
-(l-1)/2 for odd l and 1-l/2 for even l. This loses no box in the reduced
domain: the canonical interval minimizes the furthest endpoint distance
for any receiver x in [0,a/2], and the weights do not depend on those axes.
The active starts are complete for this finite family. The block mass is
at most 1024, so every critical radius is at most two; a positive
intersection with [0,a] requires starts strictly between -3 and 3, with
our inclusive range conservatively containing them. All a_k>1.

Let B denote the product block of distinct lattice indices. Its mass is

    m_B = product(l) + theta*(-1)^(j_p+j_q)*product(inactive l)

when both active lengths are odd, and product(l) otherwise. The code
independently counts even/odd active sites and multiplies by inactive counts,
verifying this distinct-site mass for every feasible raw box. Missing sites
in the deletion examples contribute zero, rather than a ghost unit atom.
The reported distinct_positive_atoms counts occupied points only.

There is exactly one common radius per block, not five independent radii:

    r* = (m_B/32)^(1/5).

Integer fifth-root bisection at denominator 2^40 constructs rlo<r*<=rhi,
with rhi-rlo=2^-40. At radius r the receiver box is

    product_k [max_{j in B} p_jk-r, min_{j in B} p_jk+r],

intersected with the reduced receiver domain. Inner boxes strictly witness
average >1; outer boxes contain every critical-radius box in the registered
family. At deletion theta=1 this formula conservatively retains zero sites
when taking endpoint extrema. This is harmless for inner validity. For
completeness of the finite-family enclosure, strip any all-zero outer
coordinate slices: the occupied support and mass stay unchanged, and its
physical bounding box is that of a smaller contiguous block, still having
lengths <=4 and already enumerated. If the entire block is zero it is omitted.
Thus deletion does not hide a larger receiver box within this finite family.

All raw boxes are used before exact deduplication and dominance removal.
Each result stores the complete retained union representation, including
starts, lengths, distinct mass, common radius and all ten rational endpoints.
The full raw family can be reconstructed from the code; omitted dominated
boxes cannot change the union. Recursive rational slicing and an independent
coordinate-mask union agree for all eight inner/outer runs.

## Results and exact comparisons

| Candidate | Inner ratio (display) | Outer ratio (display) | Gap C_lower - outer (display) | Retained inner/outer |
|---|---:|---:|---:|---:|
| weighted first-tail, theta=1/8 | 1.9472854746534198 | 1.9472854746662710 | 0.04703098690221468 | 96/96 |
| weighted tail-tail, theta=1/8 | 1.9563441893771545 | 1.9563441893898657 | 0.03797227217862009 | 95/95 |
| deleted first-tail, theta=1 | 1.6264996070010154 | 1.6264996070110578 | 0.3678168545574279 | 62/62 |
| deleted tail-tail, theta=1 | 1.6379763760647141 | 1.6379763760753059 | 0.3563400854931800 | 71/71 |

Decisions use full exact fractions, not these decimal displays. C's interval
is read directly from `lower_round003/a2_equals_a3_results.json` (SHA256
f385811e5460af63abc117804965461eed88d8ecb5b9d2cd30268387646f45d2).
Every outer ratio is strictly below its exact lower endpoint, which proves
strict inequality to the full algebraic C without relying on rounding its
display. The final code SHA256 is recorded in every result. Runs used at
most 2323 recursive nodes and less than one second each, below the
pre-registered per-run 120-second / 200000-node caps.

## Finite collar, smoothing, and all-scale scope

The following elementary transfer applies to any of the valid inner periodic
certificates; no winning transfer claim is needed for this failed round.
Let E be the inner witness union in a full period cell Q, whose period mass
is four. Retain n^5 translated receiver cells. Because there are finitely
many registered block offsets and radii, a fixed integer period collar k
contains every source atom used by every translated inner witness. Truncate
the periodic atomic input to those source cells. Its total mass is at most
4*(n+2*k)^5, and the level-one witness set has volume at least n^5*|E|.
Consequently the finite-input ratio tends from below to |E|/4=8*V.

For an actual L1 input, replace each distinct positive atom by a nonnegative
smooth bump with the same mass, supported in a sup-norm eta-neighborhood.
There are finitely many strict inequalities 32*rlo^5<m_B. Choose a common
eta>0 so that 32*(rlo+eta)^5<m_B for every retained witness. A receiver in its
inner box uses the centered cube with common radius rlo+eta, which captures
each selected bump completely. Its average remains >1. Thus the same
receiver witness volume transfers to smooth compactly supported L1 inputs
with unchanged mass. Duplicate selected atoms are never counted twice.

The maximal operator here is the full supremum over all positive common
cube radii. It dominates these selected witnessing averages, so restriction
to finitely many witness radii is legitimate for this **lower-bound transfer**.
The reverse inference is unavailable: failure of lengths <=4 cannot bound
the full all-scale level set, since longer capture blocks were not excluded.
Nor do these four input checks imply optimality for arbitrary theta,
arbitrary periodic deletions, or multiscale mixed captures.

## Reproduction and remaining question

From this directory run, for each NAME in weighted_first_tail,
weighted_tail_tail, deleted_first_tail, deleted_tail_tail:

    python3 parity.py NAME
    python3 parity.py NAME outer

Then run `python3 verify_results.py`. No external numerical library is
needed; exact integers and fractions supply the numerical evidence.

The finite screen disfavors these coarse theta choices. The next genuinely
new overlay question would require specifying a mixed-layer captured set,
counting distinct coincident atoms only once after combining their masses,
and showing a new receiver region beyond the component union. This round
has neither such a candidate nor a universal obstruction to it.
