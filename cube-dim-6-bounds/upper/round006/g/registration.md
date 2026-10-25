# A: finite first batch, 2026-10-25

This continues the existing body32 model. It seeks the **model** inequality
`g >= 1/60000`; it does not claim a C6 bound. Main-clone inputs are read only.
The frozen protocol, thresholds, role map, STATE, ledger, existing spatial
lemmas, and exact fan are hashed by `certify_local.py`.

Budget: at most 8 adaptive boxes, 24 new point/epsilon pairs (8 geometric
points), and 20 minutes of batch computation, whichever stops first. Precision
is Arb 384 bits; no binary64 value is used as a lower certificate. One Python
process is used, no MATLAB calls and no recursive agents. Stop after this
batch and report to the primary coordinator.

Points prioritize the previous minimum near `(1,1,.893,.807,.467,.097)`, an
adjacent interior point, tie/boundary points, axis and S2 strata, and one
generic control. For each point compare epsilon in `{1/200,1/100,1/50}`.
An exact rational radius is then chosen adaptively from the actual Arb result
and the proved spatial multiplicative comparison. Boxes are clipped to the
positive face. Their complement remains explicitly uncovered.

The existing right-inverse spatial estimates will also be inspected for a
tighter second-derivative bound. This seeks a new local estimate, not a repetition
of the excluded fixed coarse mesh. Full-face certification cannot be inferred
from any of the finite point data or the total volume of the boxes.

Ownership: only `upper/round006/g/` in the independent theory clone. No commits,
pushes, shared staging, or main-clone writes.
