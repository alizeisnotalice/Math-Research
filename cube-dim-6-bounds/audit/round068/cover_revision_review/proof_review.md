# F068: repaired state and mathematical-contract review

Frozen preparation seal: `a3b254ce4ba9b833ba369fa72c65b457e648d563b10c61bb969199a087e59db3`. The state/core and mathematical importer pass this independent review. Actual production authorization and launch/timeout handling are separate gates; the cleanup-inclusive runtime qualification in section 6 must not be omitted from a release claim.

## 1. Fixed origin closes R1

The genuine Context loads the pinned C063 source snapshot and derives the C068 genesis, including its execution-input-manifest hash. `verify_records` compares the entire first record to that exact derivation. The actual saved genesis hash is `10da5e4344f7876af1152b641bb9d02b913798c05542973ef6426b35a94b464a`; it is PREPARE_ONLY, queue [r], no release, no inflight attempt, zero reservations and zero native calls.

`validate_domain` fixes r's identity, parent, depth, raw box, tightened box and tightening witness, with root exactly [0,1]^5. It also fixes model, target, original input identity and both initial-batch caps. Thus the formerly accepted smaller-root state is not an admissible state. A PREPARE_ONLY state must equal the complete original unexecuted source state.

The previously proved closed ordered-face tightening and split identities can now be applied inductively: the fixed full root is the base case, and only the registered exact split transition can add child geometry. Lower-dimensional faces and closed split seams remain included. This proof concerns the region union; a passed status still additionally requires its actual mathematical certificate.

## 2. Typed whole-state replay closes R2 and protects topology

`append_event` accepts an event, not an arbitrary successor state. It locks, recovers the current history, checks the expected head hash, computes the successor through `apply_event`, then publishes that computed state. Event names and exact field sets are checked. There is no reset, refund, arbitrary requeue, batch replacement or subtree deletion event.

Recovery recomputes every successor from its preceding verified state and the same event, then compares the whole stored state to that result. Consequently a locally self-consistent queue cannot hide a missing sibling, an erased attempt or a counterfeit passed leaf. Hash consistency is necessary but no longer substitutes for transition legality.

Each reserve uses the existing deterministic BFS rule and consumes a new serial and a cumulative reservation. The same node at a new cutoff retains its earlier attempts. Existing node IDs cannot disappear; old token/input bindings, original handles and observations persist; a passed leaf remains exactly unchanged. An inflight slot is cleared only by a reconstructed original import or a retained explicit failure stop. A halted history permits no later transition.

The independent durable test retained `r.a0.n000001` and reservation count 1 after an unknown-orphan stop. Stale predecessor, reset, second inflight reserve, reinitialization and post-orphan reserve all failed. No native process existed in that fixture. The independent mixed-cutoff batch used nodes [r,r,r,r0,r1,r00,r01,r10], consuming eight total reservations for six distinct nodes. The ninth reservation failed. This checks cumulative attempts across cutoff changes, not merely eight different leaves.

## 3. Environment and release validation close R3

A TEST_ONLY Context admits only TEST_ONLY state with release=None. It cannot process AUTHORIZE_ONCE. A genuine production Context admits PREPARE_ONLY only as the original state; MAIN_RELEASED requires an actual release reconstructed by `verified_release` from the unique main_release.json and accepted source artifacts. The original illegal test-to-production successor is rejected by whole-state replay, and an injected arbitrary state field is rejected before transition processing.

The release contract binds the execution input manifest, current sealed preparation, finite 8/8 policy, the 180-second convention, accepted F061 manifest and main's revision acceptance artifact. It checks the referenced review artifacts and F061 sealed artifacts. The independent review itself creates no release or main-acceptance file. Current real-tree reserve is rejected because preparation does not authorize dispatch.

These guarantees apply to the actual pinned Context and production entry points. Passing an invented Context object or bypassing the entry points in arbitrary Python is not a substitute for main authorization. The production-branch unit fixture below is deliberately not an authorization test.

## 4. Original native artifacts and exact mathematical reconstruction

`math_contract.reconstruct` binds the current token, exact input, unique attempt folder, prelaunch intent, native bytes, stderr bytes and receipt. For genuine production imports it also binds the release and observed wrapper/native IDs to the original process record. Source/binary/fan/runtime hashes, exact command, model, normalization, time window and center are checked against the pinned context and reservation.

All six physical gradient intervals are read from canonical three-field dyadics. The importer requires positive band and endpoint lower bounds, supported counts and timing metadata, and native Euler containing zero. It independently reconstructs the interval

    x.grad(F)+5F+T*p_T-a*p_a

using exact rational interval arithmetic, including the first/normal derivative. Only after this check are the other five physical fields selected for the face certificate. The new mixed-sign interval fixture exercises all six fields. Corrupting only the normal derivative makes the independent Euler check fail even though the stored native Euler still contains zero.

These are consistency and binding checks under the accepted engine proof; Euler is not an independent six-component gradient reconstruction. Zero/repeat counts are checked for permitted ranges and parity, not independently re-enumerated. The code reports these limitations correctly. F061's two registered cases and twelve actual gradient fields do not mean future point/window inputs have already been independently rerun; every newly reconstructed jet explicitly denies that claim.

For a strict center pass, the actual pinned A057 checker runs at Arb384. Its complete 32 labelled corner output is then checked against exact geometry, gradient-dot endpoint minimum, shared-value clamp, fixed frame32 quadratic and exponent. The saved verification must equal the entire reconstructed result. Thus a rehashed PASS string or rehashed invented exponential endpoints does not become mathematical evidence.

The independent original-artifact attack first changed valid native inventory metadata and updated its receipt hash after a verification was saved. Reconstruction rejected the old verification. A separately rehashed false whole-box PASS was also rejected by the real Arb reconstruction. Both attacks changed only newly created review fixtures. No original author/native artifact was touched.

Configuration 9 independently exercises the production-only mathematical branch using fictional IDs and local unit artifacts. Correct branch fields reconstruct, while wrong release, source, time window and original PID are rejected. This does not authorize a real production Context, attest an actual process, or validate a real launch. The deliberate test marker adaptation is contained in that unit fixture and is never published to the actual tree.

## 5. Coverage and finite-run conclusions

The state induction preserves the full ordered-face geometry and the unique unresolved frontier. Whole-box evidence is tied to each leaf's own current cutoff. A split inherits that cutoff; only a strict center upper endpoint below target advances the ladder. Old-window commands and receipts were independently rejected after two successive cutoff imports at the same geometric node.

A completed structural tree still becomes a genuine body32 cover only after every original native/box result is authentic, all original attempts are settled, and main accepts the final partition and analytic inequalities. Synthetic fixtures do not prove such a cover. The real tree remains the unexecuted genesis. The all-scale/symmetry/homogeneity transfer and any C6 integration remain main's responsibility.

No automatic next batch, retry, precision escalation or native loop appears in the inspected entry points. `run_native_once.py` requires an already persisted matching reservation and creates an exclusive attempt directory before launch. Its existence prevents replay after an unobserved launch. `control.py` exposes separate bounded operations and does not spawn native processes. None of those actual launch/kill paths was executed in this review.

## 6. Operational runtime qualification, separate from core/math PASS

The source arms a 180-second real-time alarm and uses the remaining deadline for the native wait. It also targets only the owned child process group for SIGKILL. However, the exception path in `run_native_once.py` disables the alarm at line 95 and then calls `child.wait()` without a timeout at line 101. Receipt writing occurs after the alarm is disabled. Likewise `control.py` disables its alarm at line 111 before saving its failure receipt.

Therefore this source establishes a 180-second work deadline and a cleanup policy, not by itself a cleanup-inclusive absolute process-wall upper bound. No actual overrun was observed, because no native process was launched. This is a static limit of the implementation guarantee, not an invented failed execution or a refutation of event replay.

The user's literal per-program 180-second hard cap requires main to provide an independent outer enforcement mechanism that also retains ownership of cleanup, or a reviewed implementation with an earlier soft deadline and a bounded hard-cleanup deadline. A synthetic PASS cannot establish that operational guarantee. The current audit's core/math PASS is not an unconditional production release while this runtime convention/enforcement remains unsettled.

## 7. Validation scope

Nine independent configurations ran once; 31 actual A057 Arb384 checker calls evaluated 992 labelled corners. Native calls, author integration runs, old C063 writer runs, old formal replays, MATLAB calls and subagents were all zero. The verifier completed in approximately 0.384 seconds under its own hard 180-second alarm. All 43 input/source hashes were unchanged before and after the run.

The old R1/R2/R3 blockers are closed for the inspected actual API and fixed Context. No new state/core or mathematical-import blocker was found. Main still must handle the runtime qualification, actual release, genuine first launch and its original receipt, and final proof integration.
