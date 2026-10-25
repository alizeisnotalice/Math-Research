# Exact acceptance contract for a finite band-Taylor simplex cover

This is the mathematical/serialization contract for the future global
checker, not a claim that a complete cover has been executed in this turn.
`calibrate.py` implements the local bound and two exact geometry checks on
frozen inputs. A kernel-proof reader and an actual root tree remain required
inputs for a future global acceptance.

## Required top-level records

All coordinates, times, beta, thresholds and final interval endpoints are
exact rational strings. NaN, infinity and binary64 claims are not accepted
as finite numeric inputs. Infinity is permitted only as a theorem-backed
S2 value, not a numerical endpoint.

    {
      version: "band-taylor-cover-v1",
      frame: {kind: "body32" or "body32_plus_axes", beta: rational-or-null,
              parameter_density: "standard Cauchy scale t",
              Fourier_normalization: "(2*pi)^(-6)", direction_order_hash: ...},
      T: positive rational,
      m_star: positive rational strictly greater than T/4096,
      root_id: node identifier,
      input_catalog: [...verified value-proof references...],
      transfer_proof_id: approved model-matched source/interface proof,
      refutation_catalog: [...approved fixed-model/time obstructions...],
      nodes: [...finite binary-tree nodes...],
      claimed_complete: boolean
    }

The checker computes the frame vectors and Gram itself; it does not trust
a supplied Gram scalar. The only allowed Gram scalars here are32 or
32+beta^2 with beta>0 and all six axis weights equal. Source records must
match that exact frame, parameter convention and normalization.

A hash binds a source file's identity; it does not prove a numerical value.
Each input_catalog reference must identify an independently verified kernel
certificate or a rigorous evaluator result approved by the main verifier,
including its validator identity and source/input hashes. Unapproved plain
Arb balls are not enough for a mathematical GLOBAL_PASS. This turn imports
only the already verified A007/B010 records.

Before attempting a cover, compare the exact frame and T with the approved
refutation catalog. A matching rigorous point obstruction returns
EXCLUDED_MODEL_TIME. B010 already excludes body32 T=.032, and B015 now
excludes body32 T=.048. The historical T48 tube example below is a local
component, not a claim that its full target can be certified.

## Geometry node format and exact tree checks

Every node stores id, parent_id, its ordered six rational six-vectors, and
exactly one of SPLIT, TUBE, BAND, UNRESOLVED. The six vectors must be affine
independent in x1=1 and satisfy1>=x2>=...>=x6>=0. The root must equal the
six prefix-one vertices w0,...,w5 specified in theorem.md, in that order.

A SPLIT record has edge_indices [i,j], exact midpoint, and two child IDs.
The checker recomputes every edge's Euclidean squared length exactly, chooses
the longest by lexicographic index tie-break, and requires that very edge.
The midpoint must be (vi+vj)/2. The first child replaces vj by the midpoint;
the second replaces vi. All other ordered vertices are unchanged. Every
split has exactly two present children and every non-root has exactly one
parent; cycles, missing nodes, duplicate IDs and orphan nodes reject the
certificate. The two children need not be asserted disjoint on their shared
face: the exact midpoint rule proves parent union and disjoint interiors.

Coverage is established only by these exact tree identities, not summed
floating volumes. A leaf cannot be deleted because its volume is small.

## TUBE leaf

Required fields: barrier_id, a single canonical_plane_k in1,...,5, delta>0,
and an exact rational lower claim.

For the present approved A009 barrier, the frame must be body32 and T>=4/125.
At every vertex v, the checker evaluates the squared distance to that same
fixed U_k=span(1_(1..k),1_(k+1..6)) and requires it <=delta^2. It must not
minimize k separately at each vertex. Convex projection norm then covers the
whole cell. The checker evaluates the exact barrier value

    L_tube=(1/5000000000)*max(1/delta-40000,0)

and requires L_tube>=m_star. At T=6/125, m_star=3/250000 and delta=1/100000
pass exactly. Old A009 proof/data cannot approve an augmented-frame TUBE
leaf; a separate model-matched barrier would be required.

## BAND leaf

Required fields: rational a with0<a<T; exact increasing partition
a=t0<...<tn=T; a source-data reference for every endpoint potential;
two source-gradient references at exactly a and T; center c; and a rational
barycentric witness c=sum lambda_j vj, lambda_j>=0, sum lambda_j=1.
Gradients use the fixed face-free coordinates x2,...,x6, with x1 fixed1.
The center in every value/gradient reference must equal c exactly.

Either source potential family is allowed, but not an unlabeled mixture:

* H data: B(c)=H_a-H_T, Delta_i=H_ti-H_ti+1,
  gradB=gradH_a-gradH_T.
* K data: B(c)=K_T-K_a, Delta_i=K_ti+1-K_ti,
  gradB=gradK_T-gradK_a.

The interval endpoints are outward rational enclosures from the verified
input catalog. Every band has exact nonnegative true mass. If subtraction
produces a negative upper endpoint, reject. Otherwise intersect its lower
endpoint with0 and keep its upper endpoint; a wide but consistent interval
is permitted, though it may fail the lower certificate.

The checker recomputes, without trusting cached claims:

    r=max_vertex ||v-c||inf,
    D=max_vertex ||v-c||2^2,
    L_B=lower(B(c)),
    ell=min_vertex lower(gradB(c).(v-c)),
    U2=sum_i ti^-2 Phi(r/ti)*upper(Delta_i).

Phi is evaluated by outward arithmetic at declared precision and rounded
up to a rational U2_upper. The certified rational leaf bound is

    L_band=L_B+ell-D*U2_upper/lambda,

where lambda=32 for body32 or32+beta^2 for the augmented frame. Require
L_band>=m_star; strict numerical separation is preferred for robustness.
There is no infinite-time-tail term and no epsilon^-2 H substitution.
The segments lie inside the convex cell by the barycentric center witness.

All source gradients and endpoint values must be for this same frame,
center, a,T and coordinate convention. A K_T gradient alone is insufficient
without the a-gradient. A Hessian upper bound cannot replace the negative
Hessian bound used here.

## Verdict and conditional termination

The checker visits the finite provided tree only; it does not start kernels,
refinement or a blind grid. Every interval calculation has a finite declared
precision cap. An unresolved sign, overflow, failed source proof, invalid
node or budget exhaustion yields NO_PASS, never a tolerant acceptance.

    verify model, T, m_star and approved source catalog
    reject an approved refuted model/time before accepting any global claim
    verify exact root and every split while traversing the finite tree
    for each terminal leaf:
        if TUBE, recompute and verify its same-plane bound
        if BAND, recompute and verify its rational Taylor bound
        otherwise append it to unresolved
    if any invalid record: INVALID
    elif unresolved nonempty: INCOMPLETE with all unresolved IDs
    else: GLOBAL_PASS, m_T>=m_star>T/4096

The last implication uses signed-permutation symmetry and the exact
K_T(r y)=r^-5 K_(T/r)(y)>=K_T(y) inward transfer. The separately verified
source theorem would then give cost<=T/(64m_star)<64. Until every leaf
passes, no global m_T, model cost or C6 improvement is reported. A builder
may checkpoint and stop at a finite authorized split/value budget, but must
retain every unresolved leaf and set claimed_complete=false.

The current eight overlapping local calibration cells are not a root
partition and cannot be submitted as a complete cover. T=.032 remains
globally refuted by B010 regardless of any passing local calibration.
T=.048 is also globally refuted by the later B015 point. The local Taylor
method remains a valid inequality for other cases, but supplies no numerical
certification for any new T or augmented frame in this turn.
