# C160 v1 prepare failure and minimal v2 delta

This is additive. The original seven source files, PACKAGE_MANIFEST.json,
Terra inputs.json and quantization_certificate.json are immutable. Sol
has not executed the v2 resume, fixture, mathematical consumer checks,
optimizer, roots or numerical preflight.

## Evidence actually read

The first terminal directory read found the two complete partial JSON
outputs, the original source bytes and an empty phaseA directory. There
was no registration, release, started marker or numerical result. The
original package hashes matched the source. The later actual operator
file `operator_prepare_failure_observation.json` is pinned in the
read-only partial observation at SHA-256
`dbf3de868713d0252298eaa0be13829f14e647a4f6fa87a96f960a74de78f1c6`.

That operator file labels itself a post hoc observation, not a complete
raw capture. Its available excerpt is: `precheck KeyError: path;
prepare UnboundLocalError: cannot access local variable importlib where
it is not associated with a value`. It reports tool exit1, no rerun,
and the original command/source hash; no full stdout/stderr was captured
because the original tool output was truncated. We do not reconstruct a
traceback, PID, wall/CPU receipt or process activity from that excerpt.
The separate agent-capacity error does not establish an active program.

Static source inspection confirms the coding defect: an `import
importlib.util` statement inside phaseA makes importlib local throughout
that function. The earlier `importlib.metadata.version` access therefore
raises UnboundLocalError. This defect is in the prepared code I supplied.
The operator's wrong `path` lookup for a package row using `name` is a
separate wrapper defect; its shell lacked a failure barrier and continued
to preparation. Neither engineering failure is a mathematical-family
exclusion. The existing strict source identity and quantization proof
are retained, not weakened to obtain a smaller experiment.

The complete reused bytes are:

* inputs.json SHA-256
  `b74e867ab92449f8604a5f0dc07a2ce9896163410f2b62444b943ed5a25250aa`;
* quantization_certificate.json SHA-256
  `737479ae78c4a841a77c9af23ec9b9f85a70d19d976ba0c6c65b6029d26159ad`;
* original PACKAGE_MANIFEST.json SHA-256
  `7719ad8534d9f3bb4048bb7bf2c0f4960ef76e7f6be6f8677816a39412e378ab`.

PARTIAL_FAILURE_OBSERVATION.json preserves the earlier read-only snapshot.
Its `original_tool_trace_or_operator_receipt_read=false` records the
state before opening the newly arrived operator observation; this later
paragraph records that observation read without changing the snapshot.
The complete raw trace remains unavailable.

## Changes limited to the missing preparation suffix and its gate

`repair_prepare_v2.py` uses explicit global aliases imported from
importlib.metadata and importlib.util. It never calls the original
phaseA prepare again. It preserves the old partial byte hashes and old
source hashes, then validates the already saved N/D witness as a new
consumer: for B=n/d, existing D=2^96 and existing N it checks exactly
`0<nD`, `nD<=Nd`, `Nd<nD+d`, plus exact identity of the already saved cap
in both JSON outputs. These three integer cross-product inequalities
establish positivity, upward inclusion and less than one grid step.
It does not recompute ceil, select a cap, compare an old M/H objective,
or run a scalar baseline. This is strict validation of existing data,
not a new quantization experiment.

It adds only the previously unexecuted version metadata, phaseA
registration and versioned release requirements, plus prepare_v2
start/exit records. It writes no new inputs or quantization certificate.
All writes are exclusive; a genuine failure retains the new partial
suffix and stops, without retry. The original failure record is retained.

`run_phaseA_v2.py` is an additive gate entrypoint. It checks the recovery
receipt, fixture receipt and new v2 exact-hash contract, then calls the
unchanged original phaseA runner. **The numerical input and output paths
do not change**: inputs.json and phaseA/* remain the original locations.
The original run_phaseA.py and common.py bytes are unchanged. Original
source156 saved-audit acceptance and all mathematical/runtime checks
remain in that runner. Versioned fields augment the primary and
independent release; old approval alone is insufficient.

The scope remains49 features, 735471 prefix plus16 true tail rows, one
SLSQP start, <=150 total objective/gradient calls, one CPU, nice10, no
wall/CPU cap, no probability renormalization. No science is authorized
by preparation. PhaseB remains separately selected and gated.

## Engineering regression and independent wrapper barrier

`scope_regression_fixture_v2.py` is prepared for Terra, not executed by
Sol. It checks actual symbol-table binding of the original failing
importlib and the corrected global aliases. It then invokes the new
metadata and module-loader helpers with engineering stubs, checking
actual dispatch rather than syntax parsing alone. It also invokes the
actual operator main on a deliberately wrong package schema in an
isolated temporary directory, with a recording subprocess stub, and
requires rejection before any child. It calls no mathematical function,
quantization check, model or numerical program.

`operator_prepare_v2.py` accepts only package `files[].name` rows; external
source `path` rows are a distinct format. It validates both original and
v2 source packages before any child, then runs the engineering fixture
and the metadata recovery as separate argv-based subprocesses with
`check=True`. Any exception or child failure exits immediately. There
is no shell concatenation and no automatic science command. Its wrong
schema failure barrier is exercised by the fixture described above.

An independent review should check these deltas, the actual fixture
receipt after Terra executes it, the source/partial identity, witness
consumer inequalities, complete unchanged scope and new release fields.
It should not convert the incomplete old failure capture into a claim
that a full original receipt or mathematical run exists.
