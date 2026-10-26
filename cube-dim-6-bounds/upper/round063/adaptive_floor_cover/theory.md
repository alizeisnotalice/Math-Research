# C063: exact ordered-face cover and a resumable true-bound tree

This preparation uses true lower bounds on whole boxes. It does not use a
quasi-lower bound valid only on a box containing a minimizer, an incumbent
pruning argument, or a qBnB convergence theorem. The branch-and-bound skill
motivates the explicit coverage and boundary checks; the elementary statements
below are proved directly for this cover.

Let `D={y in [0,1]^5 : y1>=y2>=...>=y5}`, with all inequalities closed, and
identify `y` with `x=(1,y)`. The target is `g*=1/65000`. A node is the set
`D_B=B intersect D` for a closed rational rectangle
`B=product_i [l_i,u_i]`. A node's rectangle may contain unordered points.

## Exact order tightening and emptiness

Define `L_i=max_{j>=i} l_j` and `U_i=min_{j<=i} u_j`. For an ordered point in
`B`, `y_i>=y_j>=l_j` for every `j>=i`, and `y_i<=y_j<=u_j` for every `j<=i`.
Thus it lies in the tightened rectangle `B*=product_i[L_i,U_i]`. Conversely,
`L_i>=l_i` and `U_i<=u_i`, so `B* intersect D` lies in `B intersect D`. Therefore

    B intersect D = B* intersect D.

The lower vector `L` and upper vector `U` are descending. If all `L_i<=U_i`,
the point `L` is descending and lies in `B*`, proving nonemptiness. Consequently
the ordered node is empty if and only if some `L_i>U_i`. Equality must be kept:
it can describe a genuine lower-dimensional face, including all-equal ties or
coordinate zero. Tightening is idempotent. It does not remove S2 or S3.

The exact midpoint of `B*` is a legal six-dimensional face anchor after prefixing
1. In fact the tightened midpoint is descending as an average of descending
vectors; the analytic bound nevertheless requires no order check on an anchor.
The prepared root anchor is `(1,1/2,1/2,1/2,1/2,1/2)`, on S2, and is allowed.
For every fixed positive cutoff the band is smooth on all of R6. No full-G
singularity rule or off-S2 filter may reject this positive-band input.

## Closed splitting preserves every point

For a non-singleton tightened rectangle choose the largest positive width
`U_k-L_k`, breaking ties by the smallest coordinate index. Let
`m=(L_k+U_k)/2`, and form closed children by replacing the k-th interval with
`[L_k,m]` and `[m,U_k]`, respectively. Their union is exactly `B*`, and their
intersection includes the full seam `y_k=m`. Tighten each child independently.
By the preceding identity,

    D_B = D_lower union D_upper.

No point is lost when a child is certified truly empty. The closed seam may be
covered twice; uniqueness of geometric coverage is unnecessary, while omission
of a tie or split boundary would be an error. Singleton rectangles are never
split; an inconclusive singleton observation remains a blocked leaf.

Initially the only leaf is the root `[0,1]^5`. Replacing any live leaf by these
two children preserves the union. Marking a leaf passed retains it in the
passed-leaf ledger; it does not delete it from the partition. A halted or
inflight node also remains a leaf in the unresolved frontier. Induction proves
that passed, queued, inflight and blocked nonempty leaves jointly cover `D` at
every checkpoint, with empty leaves carrying exact emptiness witnesses. Queue
exhaustion alone is not a proof: inflight/blocked leaves and certificate validity
must also be accounted for.

## A valid positive-band bound on a node

Fix one `a` from the registered descending ladder
`(1/32,1/64,1/128,1/250,1/1000)` and `T=1/2`. Write
`F=B_(a,T)=integral_a^T p_t dt` for the unit-weight body32 model. The accepted
A057 fixed-right-inverse theorem gives, at an anchor `c=(1,c_face)`,

    F(c+(0,d)) >= F(c) exp((grad_face F(c)/F(c)).d - ||d||^2/(32a^2)).

This is a whole-space positive-time statement. It holds at S2/S3, ties, zeros
and all Fourier projection collisions. Because the exponent's spatial part is
concave, its minimum on the rectangle is attained at one of the 32 labelled
corners. Duplicate corners on zero-width coordinates are harmless.

Given the certified matching interval jet `0<L<=F(c)<=U` and
`g_i^-<=partial_(i+1)F(c)<=g_i^+`, the A057 robust corner expression is

    b_d = sum_i min(g_i^- d_i, g_i^+ d_i),
    f_d = min(U,max(L,b_d)),
    C_d = f_d exp(b_d/f_d - ||d||^2/(32a^2)).

The existing pinned A057 checker evaluates outward enclosures. A node passes
only when every corner lower endpoint is at least the exact target. The bound
is on the entire tightened rectangle, which is stronger than needed for its
intersection with `D`. Its value, five derivatives, model, normalization, time
window, anchor, complete box and external proof provenance must match. Full-G
artificial jets from 052, a different band's gradient or synthetic fixtures
cannot supply production evidence. The pure-Fraction adapter checks the exact
formula and geometry but does not authenticate an exponential or a native jet
merely because a dictionary contains a hash; main must verify the actual pinned
producer, original result and receipt files.

The main task has separately requested F061's independent genuine-gradient
crosscheck. Before it is accepted, 061's two local boxes remain conditional and
this new cover has no production release. No 061 leaf is imported into the
initial tree and no 061 seed is recomputed. The production release must bind the
accepted F061 evidence and the immutable source/audit inputs.

## Deterministic cutoff and branching policy

At a node's exact midpoint, use its current registered `a` and the fixed `T`.
If `F_lower>g*`, attempt the whole-node box certificate. If it passes, retain a
passed leaf. If it is inconclusive, split the node, and let both children inherit
this current cutoff. A small-radius certificate around the midpoint is not
accepted as coverage of a larger node. No arbitrary tiny box replaces the node.

Only `F_upper<g*` authorizes the next smaller `a` at the same node. Lowering `a`
increases the positive band, so it is an intelligible response to strict center
failure; it is not justified by curvature cost alone. The new window requires
its own matching value and gradient. At the last predeclared cutoff, strict
failure blocks the node and stops unresolved. It is not a full-G counterexample.
If the center interval overlaps or touches the target, or if the native/ABI/
input/hash/endpoint checks fail, stop without changing cutoff or precision.

Queue priority is `(depth,node_id)` (root `r`, lower suffix `0`, upper suffix `1`).
This gives deterministic breadth-first processing. One main scheduler owns the
queue; future parallel evaluation may only use uniquely reserved jobs issued by
that scheduler. The present prepared state machine has one inflight slot and no
worker/process launcher. It never recursively creates agents.

## What a completed cover would prove, and what is not proved

If a finite completed tree has only validly passed or truly empty leaves and no
queued/inflight/blocked leaves, its exact partition establishes `G(x)>=g*` on
the full ordered face: for each leaf use its own certified positive band and
the pointwise positivity inequality `G>=B_(a,T)`. A single common cutoff across
all leaves is unnecessary. Any subsequent symmetry/homogeneity transfer and
the relation to a C6 bound remain main's separate proof and integration steps.

No finite termination, complexity estimate or successful whole-face cover is
claimed here. A finite positive cutoff ladder can be inadequate. Even if full
G exceeds the target everywhere, the last fixed band's minimum may not do so.
Fixed precision may also be unable to separate a strict mathematical inequality.

A sufficient idealized termination condition would be a common strict uniform
margin for the final band's minimum on this compact face, uniformly reliable
value/gradient enclosures and exponential bounds, and fair longest-side
subdivision. Smoothness and compactness would then make the whole-box loss
uniformly small on sufficiently small boxes; each branch has only finitely many
cutoff decreases. These numerical and margin hypotheses have not been proved
or tested for this proposed cover. The finite first-batch cap is an execution
policy, not evidence of any such theorem.

The user's total research effort is unlimited. Each individual future program
still has a hard 180-second wall limit, with only its owned process reaped on
failure, original receipts/checkpoints retained and no automatic replay. The
proposed first batch reserves at most eight new geometric nodes and eight total
native attempts (including all cutoff attempts), followed by feedback and a new
main decision. No production kernel is authorized by this preparation.
