# C063 authoritative state and recovery contract

There is one authoritative format: canonical JSON records named
`state/snapshot_00000000.json`, `snapshot_00000001.json`, and so on. There is no
second CSV queue or mutable "latest state" file to diverge from it. Main alone
appends through the exclusive `.writer.lock`. An immutable snapshot is published
with an exclusive atomic hard link after its contents are flushed and fsynced;
the directory is then fsynced. Existing snapshots are never overwritten.

Each record is `{sequence,previous_sha256,state}`. The hash links the complete
bytes of the prior record. Canonical JSON uses sorted keys, compact separators,
UTF-8 and one final newline. Recovery verifies every sequence/hash link, state
invariant and unchanged input-manifest identity. A previously passed node must
be byte-for-byte structurally unchanged in every later state. Incomplete
`.pending` files are reported but never treated as published state. A missing
sequence, invalid JSON, noncanonical file, changed input or altered passed leaf
is an error and stops recovery.

The state contains:

- The immutable input-manifest SHA256, model, exact target, environment and main
  release identity. Environments are `PREPARE_ONLY`, `TEST_ONLY`, or
  `MAIN_RELEASED`. A test state cannot receive a production release.
- A map of nodes keyed by stable path ID. `r` is the root; suffixes `0` and `1`
  are the closed lower/upper children. Every node stores parent, depth, original
  rational rectangle, exact tightening data, tightened rectangle, current
  cutoff index, attempts, children and any pass proof.
- A unique BFS queue sorted by `(depth,node_id)`, one explicit inflight slot,
  persistent passed/empty leaf lists and a separate halted/blocked reason.
- The finite batch's distinct started-node IDs and its native attempt
  reservations. Reservations are consumed before dispatch, not after success;
  a failure or unknown launch cannot silently reclaim budget. The first format
  permits at most eight started nodes and eight total reservations.
- Observed native-call count, with per-attempt `native_launch_observed` and
  `launch_outcome_unknown`. An unknown interrupted launch is not silently
  recorded as an observed zero call; the reservation remains consumed. Actual
  receipts settle process counts before any final ledger is accepted.

An attempt token such as `r01.a2.n000007` binds a single node, cutoff and serial
number. Its canonical exact input includes body32 identity and normalization,
six-coordinate anchor, the positive band window, complete tightened face box,
target, fixed quadratic and immutable input-manifest SHA256. The token and
input digest accompany the native result, box result and external receipts.

The required dispatch/import transaction is:

1. Main verifies the actual source/audit/runtime files and its finite release,
   including accepted F061 gradient evidence. `reserve` removes only the next
   queued node, records the token/input, consumes the budget and sets inflight.
   Main publishes and fsyncs this snapshot **before any native launch**.
2. Main or a uniquely assigned future worker launches only that reserved input
   under a 180-second whole-program limit. It writes native stdout/stderr and
   its original receipt. Main appends the observed original tool session/PID/
   deadline handle; that handle is immutable. No launcher exists in this
   preparation module.
3. Main authenticates the saved native jet and the pinned A057 box checker
   execution. The adapter binds model/window/center/box and independently
   recomputes rational corner data. It cannot establish provenance or exp
   correctness from hash-shaped strings alone. Main must verify actual files
   and compare their hashes to the release. `ingest` accepts only the current
   token and exact input digest and imports the result once.
4. Main appends the resulting queued, split, passed or blocked state. If the
   process crashes after output is saved but before import, recovery imports
   the original completed result. It never runs that attempt again.

If recovery finds inflight, it must inspect only the recorded original live
handle or its original completed receipt. If the crash occurred after dispatch
but before the handle was saved, the launch outcome is unknown: stop and seek
that original receipt/process identity, without reissuing the token. Timeout,
hash/ABI/input failure and target overlap preserve a blocked leaf, receipt and
frontier. Only the owning process group may be stopped/reaped. No automatic
precision increase, retry, job stealing or replay is permitted.

`check_saved_box` checks 32 labels including duplicate geometric corners,
physical displacement, gradient-dot interval minimum, the shared-value clamp,
the fixed `||d||^2/32` curvature, the exact exponent and the lower-endpoint
minimum. Its certificate SHA256 is over canonical parsed certificate content;
the external artifact's byte SHA256 remains a separate receipt/manifest field.
The pinned checker source SHA256 is also required. A real jet requires
`C016_NATIVE_SAME_WINDOW_POSITIVE_BAND` provenance, six-coordinate matching
metadata, five face derivative intervals, Arb512/full-group counts and actual
native/receipt/accepted-gradient evidence hashes. Synthetic fixtures are allowed
only in `TEST_ONLY` state; full-G/artificial source kinds are rejected.

Budget exhaustion does not delete queued work or imply mathematical failure.
It returns control to main. A later batch requires a new main decision and a
reviewed continuation; this initial version deliberately does not include an
unbounded driver or automatic batch-reset function. Full snapshots make this
first finite batch easy to audit. Large-scale checkpoint compaction or a
multi-inflight extension would need an explicit reviewed format change; neither
is silently enabled here.

The global theorem flags remain false even when the structural tree verifier
reports that all nonempty leaves have passed. Main separately authenticates all
proofs, settles unknown attempts, verifies the final partition and performs any
symmetry/homogeneity and C6 integration. Empty queue, synthetic test completion
or two local 061 boxes alone cannot promote a global statement.
