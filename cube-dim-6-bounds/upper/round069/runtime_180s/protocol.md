# Runtime supplement (directory round069; scheduled archive round071)

This supplement changes only the process launch/control envelope. C068 core_v2,
math_contract, immutable genesis and full-chain replay remain frozen. No source
or native binary is compiled. No body oracle is called during preparation.

Production entry point is one explicit `envelope.py production native TOKEN` or
`envelope.py production control OP ...`; there is no batch loop or automatic
next operation. Main retains the eight-new-node/eight-native-reservation limit.
The existing C068 queue, reservation, unique inflight and mathematical evidence
rules still decide every transition. A killed/partial operation remains an
original orphan until main records and reviews its actual evidence; it cannot
reset, refund, or re-launch its token.

The supervisor starts a new-session worker and an independent watchdog outside
that worker's group. Worker/native stay in the worker's original process group.
Native Popen explicitly has start_new_session=False. Work starts only after the
watchdog sends an armed acknowledgement. The worker remains alive after its
logical result, preserving the original group-leader identity until group kill
and watchdog acknowledgement. No process is found or killed by name.

Fixed production times, relative to supervisor monotonic entry: soft work cutoff
170 seconds; watchdog group SIGKILL at 178; watchdog supervisor SIGKILL at 179;
watchdog default SIGALRM at 179.5; supervisor default SIGALRM at 179.8. Worker has
an additional own-group fallback at 178. The supervisor never cancels its total
alarm, and the independent parent fuse remains armed through the final receipt
write. Kill and watchdog messages precede any potentially blocking disk writes.
Direct-child reap calls use explicit remaining-time bounds, at most 1 second
for native cleanup; supervisor direct-child reap budgets total at most 1 second.
Partial, unobserved, or unreaped outcomes are represented explicitly. No claim
that all descendants were reaped is made. Absolute wall-time timers are OS
enforcement requests, not a theorem about a suspended machine or kernel.

There are exactly six fixed short runtime configurations, each with supervisor
timer 1.6 seconds and no mathematical imports/calls. The standalone child stubs
also have a 1.55-second default termination timer. The single fixture harness
refuses an existing started marker or fixture directory and stops on first
failure. Normal and soft-cleanup child exits are distinguished from forced group
termination. A FIFO with no reader deliberately blocks cleanup or receipt open;
the latter must leave before_receipt.json but no final receipt.json. A real
child/grandchild verifies PG inheritance. Known stub PIDs are observed via ps;
missing or zombie is reported as non-running, never as proof of reaping.

Fixture constants are selected solely by the separate fixture entry mode and
six names. There is no CLI timeout override. Fixtures are disabled if a runtime
release file exists. Production requires a main-owned runtime_release.json
binding this preparation seal, the accepted independent runtime review and the
actual C068 main_release.json; C068's own independent math/replay release is
validated again. The template is explicitly unauthorized. No production release
exists during this preparation.

Evidence limitations: six stub configurations exercise runtime control paths,
not the native body engine or all possible OS faults. No claim of global g or
new C6 bound follows from these tests. The mathematical fixtures already passed
in C068/F068 are not rerun. Main must review the runtime supplement before any
real node reservation/launch through this entry point.
