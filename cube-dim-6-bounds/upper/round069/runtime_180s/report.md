# Runtime supplement — partial author run, sealed after first failure

Directory: upper/round069/runtime_180s. Main schedules this object and its
independent review for archive round071. Status is **PARTIAL; runtime acceptance
not established; production remains closed**. This is not a mathematical result
and yields neither a global g certificate nor a new C6 bound.

## Scope and preserved inputs

The task was to supplement C068's accepted replay/math implementation with a
180-second whole-program envelope covering native work, cleanup and receipts.
C068's core_v2.py, math_contract.py, original launch/control sources, preparation
seal and genesis snapshot remain byte-for-byte pinned and unchanged. All writes
in this task are under this supplement directory. There was no compile, body
oracle call, native call, old mathematical fixture rerun, state transition,
commit or push. Original first-batch policy remains eight new nodes/eight native
reservations and requires main's explicit later release.

The independent A068 manifest is
517abb6583d09006b56e22f4fe16efbd679551756a62c2c997b4e25a4fe1d5da.
It accepts core/math while requesting this runtime supplement. The current
supplement is not that supplement's independent acceptance.

## Implementation prepared, not production released

The envelope creates an own worker process group, with native children inheriting
that PG. An external watchdog acknowledges readiness before work starts. Fixed
production offsets from supervisor entry are soft170, worker-group kill178,
supervisor kill179, watchdog default alarm179.5 and supervisor default alarm179.8
seconds. Worker also has a group fallback at178. Total alarms are never cancelled
before cleanup or receipt writes. Native cleanup uses a direct original-child
kill and wait bounded by at most1 second. Supervisor reap budgets total at most1
second. A logical task result is distinct from the SIGKILL used to close the held
group leader; direct-child reaping is distinct from descendant-group termination.

The copied native/control layers preserve the frozen C068 artifact and mathematical
contract. New runtime release validation also binds a future independent runtime
acceptance, this supplement seal, and the original C068 main release. No actual
runtime_release.json exists. The template says production_authorized=false.
These are source-level descriptions, not a claim that the incomplete test batch
has validated every production path. OS scheduling/suspension is outside a
mathematical wall-time guarantee; the source uses explicit wall-clock termination
mechanisms without disabling them during recovery.

## Exact registered execution and outcome

Registration SHA256:
dd266049c572de60f832328e0ca82dc85c393db3b136b8337e47bb5a5d7fc8e1.
Input manifest SHA256:
686ebbd1ebe1d73ee3090b00692504d8f375f3aceb7b33900001bfa79c154322.

Only one invocation of check_runtime.py occurred: tool session87047, harness
PID79440, observed final exit1. The original handle was polled; it was not
reissued. Six configurations were registered, three started, two passed the
registered assertions, and the third stopped the harness. Configurations4–6 were
never started. No source, registration, marker or fixture artifact was replaced.

| Configuration | Observed result | Supervisor PID / exit | Envelope completion elapsed |
|---|---|---|---|
| normal_child | PASS | 79448 / 0 | see original completion.json; harness wall0.107589s |
| soft_timeout | PASS | 79453 / 1, expected stop | 0.568s approximately; harness wall0.592914s |
| cleanup_hang | FAIL: required independent group-kill event absent | 79459 / 1 | 0.572226s |
| descendant_group | NOT RUN | none | none |
| receipt_write_hang | NOT RUN | none | none |
| gate_and_deadline_isolation | NOT RUN | none | none |

The soft case observed a fixed work-cutoff signal and a child wait timeout,
followed by bounded native-like cleanup. It does not prove signal handling is
instantaneous. The cleanup case recorded soft cutoff at0.352453s but its task
returned BrokenPipeError(32,'Broken pipe'). The supervisor then issued owned-group
SIGKILL, received watchdog acknowledgement and reaped worker79460 and watchdog79461.
The independent watchdog's own178-equivalent group-kill event did not occur.
The fixture's intended prolonged FIFO blocking scenario was therefore not covered.
This is neither a runtime PASS nor evidence of an observed real guard overrun.
The cause of the FIFO's broken-pipe behavior is not established by these receipts.

All three completed envelopes record worker and watchdog reaped=true with
worker returncode-9 and watchdog returncode0. These are observed direct-child
results; all_descendants_reaped_claim remains false. A later read-only ps query
of only the twelve known original PID receipts returned no rows (exit1):
79448,79449,79450,79451,79453,79454,79455,79456,79459,79460,79461,79465.
That is evidence that they were absent at observation time, not an assertion that
the supervisor reaped grandchildren. The read-only check killed no process.

## Unexecuted intentions and handoff

Configuration4 intended to start a real child and grandchild in the worker PG,
block their parent cleanup and observe the independent watchdog terminate that
group. Configuration5 intended to finish work, preserve before_receipt.json,
block the supervisor on receipt FIFO open, and observe the independent watchdog
kill the supervisor with no completed receipt. Configuration6 intended to reject
timeout/reset argument overrides and an actual absent-release production gate
before any mathematical import. None of these intentions is claimed as tested.

Main directs an independent auditor to use new fault injection, including a full
pipe with a retained reader that never drains, real descendants and unrelated-PG
boundaries. That avoids falsely repeating the author's unfulfilled FIFO scenario.
Author will not rerun normal/soft, repair this frozen source in place, or start
the remaining configurations. A substantive production-source finding requires
a separate revision directory and finite registered batch. Production remains
closed pending the independent review and main's explicit release.

## Reproduction and evidence whitelist

The original command was the registered local venv Python with -B check_runtime.py.
The one-shot marker deliberately prevents rerunning it in this directory. Inspect
registration.json, input_manifest.json, fixture_run_started.json,
execution_receipt.json, results_partial.json and each original fixture receipt /
completion / stub handle for the actual evidence. There is intentionally no
results.json that would claim all six passed. FIFO special files are not read,
hashed or included in the archive whitelist. The manifest lists regular source,
configuration, report and receipt files only. Parent owns any eventual archival
commit; this agent performed no git mutation.
