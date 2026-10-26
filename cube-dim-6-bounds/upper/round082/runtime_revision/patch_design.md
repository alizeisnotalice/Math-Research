# C082 minimal runtime revision

Scope: replace only the failed cleanup/receipt boundary of the C069 runtime
supplement (archive071). C069 sources, fixtures and seal and all C068 core/math
files remain frozen. No production operation is released. No native/body/Arb
call or old normal/soft/math fixture is part of this batch.

The fixed A071 manifest is
7cf470bc51ac1f79771c65fe4abac79eaeb4eff145c1734674f1586d42392ae7.
Its retained-reader full-pipe case actually blocked, then observed an uncaught
PermissionError at the supervisor's redundant second killpg. Before-receipt,
receipt and completion were absent. This was a receipt-loss finding, with no
observed180-second overrun. The OS-internal cause of EPERM is not inferred.

## Exact changes

1. Accept group acknowledgements only from the original watchdog pipe and PID.
   When that acknowledgement already records group kill/not-found, omit the
   redundant supervisor killpg. The original unreaped worker is still the identity
   anchor; no new PID/name lookup is used for kill selection.
2. Without acknowledgement, catch OSError including EPERM, persist the error,
   mark the task stopped and termination unconfirmed. Do not send a fabricated
   kill-ack request to the watchdog. Bound further original-event reading and
   observe only the original PID with a bounded ps call. A ps observation alone
   never promotes the group to confirmed/reaped. Without ACK retain the original
   leader until the original watchdog exits or its bounded reap remains unresolved.
3. Receipt exceptions preserve before_receipt and a separate receipt_failure;
   they are never labelled final_receipt_written=true. Receipt blocking remains
   under the external supervisor-kill fuse and default total alarm. Watchdog
   remains armed until either the final receipt or a receipt-failure record is
   written. When group termination is unconfirmed, no disarm message is sent;
   closing the control pipe leaves cleanup to the original watchdog/worker fuse.
4. On a future released native stop, create failure.json only if absent, bound to
   the original inflight token/input and runtime/original artifacts. Preserve an
   original worker failure rather than replacing it. The existing reservation is
   already consumed; this helper never refunds, resets, relaunches or edits the
   replay state. Main must perform the original-token HALT after observing the
   actual handles. If file I/O itself blocks, the already persisted reservation
   and runtime start/handles remain orphan evidence; no success is claimed.

Production constants remain soft170/group178/parent179/watchdog179.5/
supervisor179.8 seconds; native start_new_session=False, worker own PG, and fixed
bounded reap remain. The revised supervisor uses a single cleanup budget at most
1 second (fixture0.12 seconds), including its ACK/reap sequence. There is no
timeout override. Default total alarms are not cancelled during errors/writes.
The watchdog, worker and native/control launcher sources themselves are exact
copies of the frozen runtime versions; only supervisor/support and isolated test
code change. Real production still needs this seal, independent acceptance and
main runtime/C068 release. No release file is created.

## Four once-only registered configurations

1. retained_reader_cleanup: soft deadline enters a real full anonymous pipe with
   its read end retained and never drained. Require independent group kill event,
   no second supervisor kill, preserved receipts and direct-child reaps.
2. early_return_descendants: a task returns after a real child/grandchild are
   ready in its PG. Require their termination while a harness-owned unrelated PG
   sentinel remains alive; then separately stop/reap only that sentinel.
3. receipt_retained_reader_hang: complete the task, write before_receipt, then
   block the supervisor in a real retained-reader full pipe. Require independent
   supervisor SIGKILL and partial checkpoint with no final receipt. Do not claim
   the killed supervisor reaped its watchdog.
4. receipt_error_gate_unacked_eperm: dynamically reject the absent production
   release and invalid timeout/reset arguments before mathematical imports.
   Fixture-only injection raises EPERM before any supervisor group kill/ACK;
   require an unresolved error plus original-PID observation. Make receipt.json
   a directory to force a real write exception and require receipt_failure plus
   final_receipt_written=false. The injected EPERM is a controlled branch test,
   distinct from A071's actual OS EPERM.

Each configured supervisor deadline is1.6 seconds; auxiliaries self-limit1.7,
harness wait1.95, and the entire harness has179.8-second default alarm. At most
four configurations run once, first failure stops with no automatic continuation.
The harness writes each observation even on an assertion failure. This batch is
runtime evidence only; independent audit and main release remain required.
