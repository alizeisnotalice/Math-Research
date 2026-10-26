# Round 102 source scope

Registered UTC: 2026-10-26T15:28:33.951672+00:00

This finite read-only source scope was registered before any source contents were read. Source cap: 30; registered: 30. No registered source may be executed.

| ID | Root | Relative path | Role |
|---|---|---|---|
| S01 | main | `upper/round068/cover_batch001/control.py` | 068 wrapper controller |
| S02 | main | `upper/round068/cover_batch001/run_native_once.py` | 068 production entry |
| S03 | main | `upper/round068/cover_batch001/integration_receipt.json` | 068 observed receipt |
| S04 | main | `audit/round068/cover_revision_review/report.md` | 068 independent review |
| S05 | main | `upper/round082/runtime_revision/worker.py` | 082 worker |
| S06 | main | `upper/round082/runtime_revision/control_task.py` | 082 control task |
| S07 | main | `upper/round082/runtime_revision/watchdog.py` | 082 watchdog |
| S08 | main | `upper/round082/runtime_revision/execution_receipt.json` | 082 execution receipt |
| S09 | main | `audit/round082/runtime_revision_review/report.md` | 082 independent review |
| S10 | main | `lower/round089/ordered_a1_derivative/ordered_derivative.py` | 089 final checker |
| S11 | endpointfix | `audit/round089/ordered_derivative_review/report.md` | 089 independent review |
| S12 | main | `upper/round090/cover_batch001_actual/release_once.py` | 090 actual entry; run_native_once.py absent |
| S13 | coveraudit | `audit/round090/actual_cover_review/report.md` | 090 independent review |
| S14 | coveraudit | `audit/round090/actual_cover_review/runtime_run_receipt.json` | 090 observed receipt |
| S15 | main | `upper/round091/cost_interval/check_exact.py` | 091 final checker |
| S16 | curvature | `audit/round091/cost_interval_review/report.md` | 091 independent review |
| S17 | main | `upper/round093/full_gamma_account/run_once.py` | 093 final entry |
| S18 | greenaudit | `audit/round093/full_gamma_review/report.md` | 093 independent review |
| S19 | main | `upper/round095/smooth_kernel_correction/check_constants.py` | 095 final checker |
| S20 | theory | `audit/round095/smooth_correction_review/report.md` | 095 independent review |
| S21 | main | `upper/round096/joint_middle_budget/run_guard.py` | 096 final entry |
| S22 | localaudit | `audit/round096/joint_middle_review/report.md` | 096 independent review |
| S23 | main | `lower/round097/full_spacing_local_max/check_saved.py` | 097 final checker |
| S24 | localaudit | `audit/round097/full_spacing_review/report.md` | 097 independent review |
| S25 | upper098 | `upper/round098/direct_green_geometry/run_once.py` | 098 production entry |
| S26 | upper098 | `upper/round098/direct_green_geometry/run_receipt.json` | 098 observed receipt |
| S27 | localaudit | `audit/round098/direct_green_review/report.md` | 098 independent review |
| S28 | coveraudit | `audit/round099/parity_mass_review/check_exact_saved.py` | 099 auditor checker |
| S29 | coveraudit | `audit/round099/parity_mass_review/timing_binding_run_receipt.json` | 099 timing binding receipt |
| S30 | coveraudit | `audit/round099/parity_mass_review/report.md` | 099 independent review |

## Limits

- At most two bounded source-listing/search batches; no recursive project-wide scan.
- Static text/JSON reads and AST parsing only; zero scientific recomputation and zero execution of legacy scripts.
- No writes outside this directory; no stage, commit, push, MATLAB, or subagents.
- Each helper process is guarded with `signal.alarm(179)`.
