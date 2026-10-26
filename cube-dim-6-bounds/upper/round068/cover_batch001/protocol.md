# C068 revision and finite-batch protocol

The independently reviewed C063 geometry remains valid, but its persistence API
was rejected for release. The final A/F063 audit seal is
`c4eb83e909789f52bdd9d908454466e5f96cf853b2405ac1da78b2eba0ac77a2`.
Its three actual public-API counterexamples were root shrinkage, rollback of an
inflight reservation followed by identical-token reissue, and TEST_ONLY
promotion to MAIN_RELEASED without a release. No author source is edited here.
The frozen C063 core is used only for exact geometry and deterministic ordinary
transitions; its old append/recover functions are not the C068 persistence path.

## Fixed origin and legal-history proof

The only C068 genesis is derived from the complete bytes of C063's sealed
`state/snapshot_00000000.json`, SHA256
`9cac6bc4fe23d9ae8f6b7c1a81130d1c6522769766b8142d2bc1f5506d797c3f`.
It carries that source hash and the C068 execution-input manifest hash. Its
state is exactly the old unexecuted PREPARE_ONLY root. A TEST_ONLY history used
by fixtures differs only in its explicit environment marker and cannot become
a production history. The real root, model, target, original C063 input identity
and 8-node/8-reservation caps are checked on every state. PREPARE_ONLY must equal
the unexecuted source state. MAIN_RELEASED must equal a release object derived
from the actual pinned main release and reviewed acceptance artifacts.

Every later immutable snapshot records one typed event as well as the complete
successor state. The only events are AUTHORIZE_ONCE, RESERVE_NEXT_ONCE,
ATTACH_ORIGINAL_HANDLE_ONCE, IMPORT_ORIGINAL_RESULT_ONCE and
HALT_ORIGINAL_ATTEMPT_NO_RETRY. Exact event fields are enforced; arbitrary
successor replacement, reset, delete, new batch and requeue are not operations.

Recovery starts at the fixed genesis and replays every event against the
previous verified state. It then compares the *entire* resulting state to the
stored successor. An intact hash chain without legal replay is insufficient.
Explicit prefix checks also preserve every old node, attempt token/input,
attached original handle, imported observation, passed leaf and nondecreasing
serial/reservation count. Existing topology can only grow by the deterministic
closed split already proved in C063. The queue is exactly the unprocessed leaf
frontier; an omitted child or a counterfeit passed leaf cannot be introduced by
a status edit. This induction preserves the full ordered-face root union.

Inflight can leave its active slot only through import of its original validated
result or a retained failure/orphan stop. The latter keeps the original attempt,
input and consumed reservation as a blocked leaf and halts the batch. It never
refunds a call or creates a replacement token. A halted history has no further
transition. A once-only origin marker and exclusive genesis directory prevent
the exposed initialization helper from resetting an existing history; the
production CLI does not expose initialization at all.

## Native and mathematical evidence

Each native attempt is reserved and its exact input durably published before
launch. The one-attempt runner creates an exclusive attempt directory and a
pre-launch input/intent; existing directories are never reused for execution.
Its exact command uses the unchanged C016 binary and pinned fan, model, source
and dylibs. Point numerators/denominator are exact, with the existing int64 /
int128 eligibility bounds. It computes the band, all six physical derivatives,
two endpoint densities and native Euler in one Arb512 invocation.

The native runner writes original stdout/stderr, immutable intent/process data
and an immutable final receipt. The receipt binds command parameters, all
source/binary/fan/runtime hashes, original wrapper/native PIDs, window/center,
release, exit and output hashes. It has a 180-second whole-program alarm, with
only its owned child group killed/reaped on timeout. Input failure consumes the
reservation, retains a failure record and does not launch or retry a native.

Before an import, `math_contract.reconstruct` reads those original files and
validates their relationships. It decodes full exact dyadic fields, requires all
six gradients, validates 967680 logical group terms and 483840 sign-pair terms,
records and checks the ranges/parity of all zero/repeated-node counts, and checks
finite timing metadata. It requires positive band/endpoint lower bounds, native
Euler containing zero, and an independent exact-rational interval reconstruction
of `x.grad(F)+5F+T*p_T-a*p_a` containing zero. It does not silently skip zero or
repeated projections. The zero/repeat inventory is bound to the original native
output and pinned complete-enumeration source; it is not independently
re-enumerated at every new point in this batch.

Only after these checks are the five physical face fields selected. For a strict
center pass, the pinned A057 checker recomputes the whole tightened node box at
Arb384, including all 32 labelled corners and actual exponential enclosures.
The frozen pure core independently checks exact corner displacement, gradient
dot minima, shared-value clamp, fixed frame32 curvature and exponent. A saved
verification is accepted only when it equals this entire native/Arb
reconstruction. Rehashing a false PASS or fake exponential endpoints does not
make them importable. Recovering a completed event can replay these *saved*
native/box checks; it never evaluates the old center kernel again.

The accepted F061 seal
`91b26762f71229fbaa4c011be12b5a678ce968afc298e23a2c5320403c51fbda`
checks the two registered C061 points and all twelve physical gradient fields.
It does not mean every future point/window has already been independently
evaluated. The release must explicitly accept future native interval jets under
the mathematical/implementation evidence for the pinned engine and acknowledge
that limited F061 scope. Every new jet marks
`independent_full_gradient_crosschecked_for_this_new_attempt=false`; Euler alone
is never described as a six-component independent proof.

## Reviewed release and bounded operations

Current authorization is preparation and eight integration fixtures only:
**zero native kernels**. A real release requires a new independent review of
this repaired revision, a main acceptance file that binds this preparation seal
and the review artifacts, and the unique `main_release.json`. A/F063's rejection
is not treated as a PASS. Neither synthetic test success nor a status string is
a production release or a mathematical certificate.

The future primary driver is a sequence of separate explicit short programs,
each under a 180-second whole-program limit; there is no encompassing loop:

1. `control.py authorize` appends the one reviewed release event.
2. `control.py reserve` persists exactly the next BFS reservation and prints its
   unique token/input identity. Exhausted 8/8 budgets retain the queue and stop.
3. `run_native_once.py TOKEN` launches that one original native only. The tool's
   original session is observed to completion; an interrupted call is polled
   through the same handle, never relaunched.
4. `control.py attach TOKEN TOOL_SESSION_ID` persists the observed original
   tool session and native/wrapper PIDs. An already-written matching handle
   artifact can finish its pending transaction without launching anything.
5. `control.py verify TOKEN` reconstructs and saves the real native/Arb result.
   Main can read it or repeat a read-only reconstruction; no kernel is called.
6. `control.py import TOKEN` verifies those original artifacts again and appends
   the unique mathematical transition. Strict center failure alone lowers a;
   center pass with whole-box failure splits at unchanged a; a whole-box pass
   retains the leaf; straddle/equality blocks and stops.
7. A native failure uses `control.py halt TOKEN` with its retained failure file.
   A missing or unknown original process is recorded as an orphan stop; it must
   not be replaced. Verification/input errors stop the operation with its own
   receipt and leave the existing frontier intact for main's retained-failure
   settlement. No next native is issued after such an error.

All controls record their own PID, operation, exit/error and before/after
snapshot hashes. They do not spawn processes. Each actual native command has
its own original process/receipt record. At most eight distinct geometric nodes
and eight total native reservations, including cutoff-ladder attempts, belong
to this batch. The ladder remains `1/32,1/64,1/128,1/250,1/1000`, with `T=1/2`.
No automatic retry, precision escalation, new batch or recursive delegation
exists. Full snapshots and one inflight slot are deliberately bounded here.

At batch end, freeze the exact queue, inflight/blocked status, passed/empty
frontier, source/release/receipt hashes and actual finite call count. Queue
exhaustion alone does not promote a global theorem. Main retains final
partition/proof acceptance, symmetry/homogeneity integration and any C6 claim.
No global g or C6 result is produced by this preparation or by synthetic tests.
