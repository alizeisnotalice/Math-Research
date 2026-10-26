# C063 adaptive full-face floor cover: preparation ready for main review

The new directory prepares a cover of the entire closed ordered positive face,
not further isolated tiny-box pilots. It contains a direct proof of exact order
tightening, true emptiness and closed-seam subdivision; a deterministic pure
Fraction tree/adapter core; one recoverable authoritative queue format; and a
finite first-batch proposal. It contains no process or kernel launcher.

The initial authoritative snapshot contains only root `r`, representing
`[0,1]^5 intersect {y1>=...>=y5}`, queue `['r']`, no inflight task, no passed
leaves and environment `PREPARE_ONLY`. Its SHA256 is
`9cac6bc4fe23d9ae8f6b7c1a81130d1c6522769766b8142d2bc1f5506d797c3f`.
The exact root anchor `(1,1/2,1/2,1/2,1/2,1/2)` is legal even though it lies on
S2. Positive bands are smooth there; no off-S2 or full-G singularity rejection
is used. No 061 seed or leaf has been imported or recomputed, and 061's seal is
unchanged. Its two local boxes remain conditional pending main acceptance of
the separately requested F061 true-gradient crosscheck.

The single preparation invocation passed the 12 registered configurations.
It used only Python's local `fractions.Fraction` arithmetic, tree operations,
metadata adaptation and filesystem checkpoint fixtures. Checks included exact
finite-grid tightening equivalence (1215 box/grid cases), closed split-seam
coverage (3125 grid points), boundary/emptiness cases, strict cutoff transitions,
same-window/model/center/box binding, passed-leaf reuse rejection, snapshot
recovery, stale/duplicate/tampered state rejection and both finite budget gates.
The code and tests were pinned before that invocation. Program PID 66271 exited
0 under the 180-second alarm; its own receipt records 0.026919666008325294 seconds.
The original tool result is separately retained. There were 0 native calls,
0 scalar kernels, 0 production calls, 0 compiles, 0 retries and 0 061 recalculations.
No program from this preparation remained alive at the final process check.

The finite tests are software evidence. The exact geometry/partition arguments
are in `theory.md`; finite grids are not substituted for those proofs. All
test jets are explicitly synthetic constant-function fixtures and cannot enter
production state. Synthetic snapshots under `test_snapshots/` include a
deliberately malformed negative fixture; they are not an alternate queue and
are not in the final artifact whitelist. They remain locally retained as test
scratch. The only real initial queue is `state/snapshot_00000000.json`.

The main mathematical and scheduling decisions are:

- Tighten by `L_i=max_{j>=i}l_j`, `U_i=min_{j<=i}u_j`; discard only strict
  `L_i>U_i`. Ties and zero-width faces stay in the cover.
- Use one midpoint jet for the full tightened node rectangle. A center pass
  followed by an inconclusive whole-box bound causes a longest-positive-side,
  exact closed-midpoint split; children inherit the current cutoff. A small
  certified box around the center never substitutes for the full node.
- Only strict center failure `F_upper<1/65000` allows the next cutoff in
  `1/32,1/64,1/128,1/250,1/1000`, always with `T=1/2` and a matching new jet.
  Center overlap/equality, source/input error or timeout stops unresolved.
- Main owns the sole BFS queue and publishes an inflight reservation before
  dispatch. Immutable canonical snapshots form a verified SHA256 chain. Resume
  observes an original live handle or imports its original saved result; an
  unknown orphan never triggers a second launch. Passed nodes cannot change
  across snapshots.
- The proposed first batch starts at most eight new geometric nodes and makes
  at most eight total native reservations, including all cutoff attempts.
  This deliberately tighter total-call cap avoids interpreting eight nodes as
  forty hidden calls. Every result/stop returns feedback. No automatic new
  batch is implemented.

The main review should decide the proposed finite budget, authenticate the
pending derivative evidence, and approve the exact native producer/receipt
import before any release. The adapter authenticates neither native arithmetic
nor exponential evaluation merely from user-supplied hash strings: it verifies
exact geometry and rational corner algebra, while main must verify the actual
pinned C016 and A057 result/receipt artifacts. The initial implementation has
one inflight slot and full snapshots for this small audited batch. Parallel
inflight jobs or large-scale checkpoint compaction require an explicitly
reviewed continuation of the same unique scheduler, not recursive delegation.

No successful global cover, finite termination guarantee, global g certificate,
new C6 bound or closure of the main partial-I argument is claimed. A fixed
positive cutoff ladder may itself be insufficient, and fixed arithmetic may
stop on an inconclusive comparison. Total research effort remains unbounded;
the per-program 180-second wall limit and finite feedback batches remain firm.
No MATLAB, delegation, stage, commit or push was used.

Key files are `theory.md`, `state_format.md`, `first_batch_plan.json`,
`cover_core.py`, `check_fixed.py`, `preparation_results.json`,
`preparation_receipt.json`, and the initial authoritative snapshot. The input
manifest SHA256 is
`9db4af7ff6b77b94dec50619873a5035e5a822a4c6fdb48911fc5f8d61a0c299`.
The `math-c04-branch-and-bound` skill was read for true-bound versus quasi-bound
discipline; no source qBnB theorem or unverified literature conclusion is used
to claim termination or coverage here.
