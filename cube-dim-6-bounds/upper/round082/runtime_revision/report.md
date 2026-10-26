# C082 runtime revision — four author configurations passed

Status: **AUTHOR PASS; independent review and main release still required**.
Production remains closed. This round contains zero native/body/Arb calls, zero
old normal/soft/math fixture reruns, zero retries and no global g or new C6 claim.
Only upper/round082/runtime_revision was written. C069/071 and C068 frozen inputs
remain hash-identical. No stage, commit, push, MATLAB or delegation was used.

## Problem and minimal revision

A071 manifest7cf470bc51ac1f79771c65fe4abac79eaeb4eff145c1734674f1586d42392ae7
recorded a real retained-reader full-pipe hang followed by an uncaught EPERM at
the supervisor's redundant killpg, before any final receipts. It observed no
180-second overrun. The revision now accepts termination acknowledgements only
from the original watchdog pipe/PID and omits a second kill after such an ACK.
Without ACK, a kill OSError is retained as unconfirmed termination, with bounded
original-PID observation and no false disarm. The original worker is not reaped
while it remains the identity anchor for a possible watchdog group action.

Receipt-write errors retain before_receipt and a separate receipt_failure,
with final_receipt_written=false. A blocked write is still covered by the
independent supervisor kill and total alarm. A future released native failure
can create an original-token failure record without refund/reset/relaunch;
main remains responsible for the immutable queue's HALT transition. This helper
was not production-executed in this round.

Production timing remains soft170/group178/parent179/watchdog179.5/
supervisor179.8 seconds. Native inherits the worker's PG; bounded cleanup remains
at most1 second, and no total alarm is cancelled during cleanup or receipt work.
worker.py, watchdog.py, native_task.py and control_task.py are byte-identical
copies of the prior frozen launch layer. C068 core/math are neither changed nor
rerun. The changes and exact finite tests were fixed in patch_design.md before
execution; source SHA pins did not change after registration.

## Single actual execution

Registration SHA:
d719404a3de868d8ea5248b251f3af55adac066aeeaf74a608478fd43f1ec204.
Input manifest SHA:
4a127182772c7b5281d911ef6c10e334ca199d330fb35c652e875ac165be6542.
Patch design SHA:
ac82c7117a65b55da8a33baea2f9f03c3bed7919466135ed8894507b64210a84.

Exactly one check_revision.py invocation ran through the registered local venv
Python: original tool session75593, harness PID89838, exit0. Harness monotonic
elapsed2.453384375s. Four configurations started and all four passed their
registered assertions. Only the same tool handle was polled, never reissued.

| New configuration | Supervisor / worker / watchdog | Expected supervisor exit | Configuration wall seconds |
|---|---|---|---:|
| Retained-reader cleanup hang | 89846 / 89847 / 89848 | 1: stopped task | 0.884832083 |
| Early return with descendants | 89850 / 89851 / 89852 | 0 | 0.148521042 |
| Retained-reader receipt hang | 89857 / 89858 / 89859 | -9: independent SIGKILL | 1.271716542 |
| Unacked EPERM, receipt error, gate | 89864 / 89865 / 89866 | 1: explicit failure | 0.129980042 |

The first case filled65536 bytes, retained its read end and entered blocking
cleanup at0.353446542s. Its persisted events include the original watchdog's
GROUP_SIGKILL_SENT; group_ack=true and group_sigkill_sent_by_supervisor=false.
Both direct children were reaped (worker-9, watchdog0), and before/receipt/
completion exist. This closes the observed redundant-kill receipt-loss path for
this actual scenario without requiring the OS to repeat its earlier EPERM.

The second case started actual child89853 and grandchild89854 inside PG89851,
then returned while they were alive. After owned-group cleanup, both known PID
observations were non-running. Unrelated sentinel PID/PG89849 remained alive
across that cleanup; the harness separately terminated/reaped only this original
sentinel, returncode-15. Worker/watchdog reaps were-9/0. This tests the selected
same-PG descendant and unrelated-own-PG boundary; it does not claim coverage of a
malicious descendant that deliberately creates another session.

The third case completed work, saved before_receipt, filled65536 bytes in the
supervisor's own retained-reader pipe at0.040844750s and blocked. The supervisor
exited-9 at about1.272s, before its default1.6-second SIGALRM. The harness sent no
emergency kill. There is intentionally no final receipt or completion. The
known worker had already been reaped; no claim is made that the killed supervisor
reaped its watchdog. This expected partial result is a test PASS, not a claim
that a final receipt can be guaranteed while its filesystem/write path blocks.

The fourth case dynamically rejected the absent runtime_release.json and
timeout/reset argument overrides before importing mathematical modules. A
fixture-only injected EPERM occurred before any supervisor group kill or ACK.
The receipt records errno1, identity_anchor_unreaped=true,
termination_confirmed=false, group_ack=false and the bounded original-PID ps
observation. Creating a directory at receipt.json caused the actual exception
FileExistsError(17,'File exists') under exclusive write mode; receipt_failure
records it and final_receipt_written=false. This differs from the design's
informal EISDIR label but is the intended real write-exception path. The original
watchdog subsequently handled pipe closure; completion observed watchdog0 and
worker-9 reaped. The earlier unconfirmed group fields are preserved rather than
retroactively claiming successful initial cleanup. Injected EPERM is explicitly
distinct from A071's real OS return.

A later read-only ps query of the16 known original PIDs returned no rows (exit1);
it killed no process. This is absence at observation time, not proof that every
descendant was reaped by the supervisor. Every record retains
all_descendants_reaped_claim=false. No actual180-second overrun was observed.

## Release, evidence and remaining scope

No runtime_release.json exists and the template is unauthorized. The eight-node /
eight-native cap and the original C068 queue remain unchanged. Main must obtain
an independent review of this new seal and explicitly release it before any
production operation. The tests exercised fixed shortened deadlines and stub
processes; actual body execution, every OS fault, deliberate session escape,
watchdog death and filesystem failure of all checkpoint paths were not tested.
The timers are OS enforcement mechanisms, not a mathematical guarantee for a
suspended machine or unscheduled kernel.

The whitelist contains source/configuration, patch design, original run marker,
results, execution receipt, this report and merged_fixture_evidence.json. That
merged evidence stores every original regular fixture file and each observation
with exact bytes and SHA, including empty stdout/stderr; the raw fixture tree
remains local as a redundant convenience. No special files need to be archived.
The original marker prevents rerunning this batch in place. Further tests require
a separately registered independent review or new revision; this author batch
is complete and does not automatically launch another one.
