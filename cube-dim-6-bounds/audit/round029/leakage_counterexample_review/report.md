# A029 independent audit — PASS: L6 is false

The proposed counterexample is valid for the exact unproved L6 statement
in frozen round025. It is a bounded compactly supported Lebesgue L1
function, not an atomic surrogate. Every permitted countable disjoint
cube family fails L6. No defect was found in the selection quantifiers,
partial-spike treatment, boundary conventions, common-radius witnesses,
or all-positive-radius component estimates.

The independent complete proof is in proof.md. It was developed before
reading the main written theorem. Eight independent Python Fraction
checks ran once and passed; results.json retains their exact rational
margins. The main theorem, registration, checker source and saved results
were then read without executing the main checker. Their small/big bounds
and volume gap agree with the independently obtained rational expressions.
The source hashes are locked in source_receipt.json.

## Points that could have invalidated the argument, and their resolution

- A selected cube may capture only part of the spike. The bound
  m_i>=|Q_i|/2 uses only that the background density is at most1/2;
  it therefore still gives R_i<=1. Countable disjointness yields
  sum m_i<=32 by Tonelli, even with null boundary overlaps.
- No geometric assumption about containing the origin is needed. On the
  side opposite the unique big cube's first center coordinate, its
  first interval has the relevant endpoint no farther than radius1
  from its center. For c1>=0 this endpoint is at least-1, giving the
  intersection-length bound r+1-a and the stated background estimate.
- Small test radii cannot touch the spike because |x1|>a and its support
  lies in[-d,d]. Large test radii satisfy decreasing upper bounds below
  h. The bounds are uniform in x and in every small component, so taking
  all radii and a countable union of component exceptional sets is legal.
- The selected slab depends on the family only through one sign. Its
  volume is fixed. This proves the necessary forall-family statement;
  it does not improperly demand that one slab work for all families.
- Every point of the open slab has a genuine full-input bad cube: radius
  ||x||infinity+d contains the entire spike, stays inside the background,
  and is strictly less than1. Hence the average is strictly greater
  than1, while every component maximal value is strictly below h.

The leakage margin is the positive rational

    5278750005499989500010249994950001
    /31250000000000000000000000000000000.

This is the excess slab volume over (43/62500)*2080. It is not a numerical
estimate of a maximal-set union. No union engine or grid was used.

## Disposition

**Refute L6 as stated.** Its prospective conditional constant
64-43/31250 cannot be entered as a C6 bound. The round025 corner lemma and
its support-compatible-height corollary remain valid; this f at lambda1
has support-cube volume4096 but mass2080 and does not satisfy that
corollary's support-height hypothesis. Each selected component can satisfy
the local hypothesis without controlling the combined input's halo.

This does not exclude every revised leakage constant or all alternative
covering statements. No new C6 upper or lower bound results. Earlier
round025 files were kept unchanged; this audit records the refutation in
a new round rather than modifying the frozen conditional proposal.

Actual work: one8-check Fraction run, no repeated scalar calculation,
zero kernel/native/union/grid/MATLAB calls, no subagents, no stage/commit/push.
Only audit/round029/leakage_counterexample_review was written. The final
receipt records the runtime window and final deliverable hashes.
