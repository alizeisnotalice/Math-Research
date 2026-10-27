# R038 Phase A delivery review

Date: 2026-10-28. Read-only cross-check of the compute report and summary, the softclone joint-cover SPEC, and the main general-cover ledger. No production file was changed, no pipeline or LP was rerun, and no Git operation was performed.

## Sources and identity

- Compute source: `/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round038/compute/R038_REPORT.md` and `phasea_summary.json`.
- Soft-agent SPEC: `/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round038/joint_cover/SPEC.md`. The path under `cube-r030-luna` is absent; the SPEC is in `cube-r030-soft`.
- Main ledger: `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/results/general_cover_ledger.json`.

`R038_REPORT.md:21` says the Phase A production registration/code hashes were frozen before execution and the theoretical SPEC was completed after launch. The pre-run `phasea_registration.json` hashes the actual production source files and 15 input JSONL files, but contains no SPEC hash. The current softclone SPEC SHA-256 is `69ef33a9c6ea5689390a81f311f8ca453c1d9c699c4525a7a32daa6563c8319c`; this identity is not bound by the pre-run registration. File timestamps also show the SPEC later than the registration (01:37:37 vs. 01:36:32 +0800). The report claims a subsequent read-only check found the implemented pair/plateau rule matched the SPEC; this is a post-run comparison, not pre-run SPEC identity. Do not describe the executed run as cryptographically frozen against this later SPEC.

The three-LP stage has a separate post-Phase-A registration. `phasea_cert3_registration.json` binds the phase summary hash, its own source/core hashes, and three selected case hashes before the three LP runs (`R038_REPORT.md:21`). This identity is distinct from the original Phase A registration.

## Phase A counts

The report (`R038_REPORT.md:3-9`) and summary agree on 79 complete cached R037 H_FAIL patterns. Parsing the 79 summary rows and comparing their stored rational `L_after` and `L_before` values exactly gives 62 strict load reductions and 17 equal loads; none increased. Terminal statuses are 50 `stop_pair_coordinate` and 29 `iteration_cap`, matching the report.

The 50 coordinate stops mean a full nonedge-pair scan found no strict improvement for the specified lexicographic pair move. The 29 capped cases reached 128 accepted moves and were not fully rescanned at the end. Neither status proves global payment optimality (SPEC §5, especially lines 84-88; report lines 5-7).

## Three LP certifications and scope

Only three selected cases received fresh full mass-LP primal/dual certificates: `C/2/5 full_disjoint`, `C/1/6 full_disjoint`, and `A/3/15 partial` (`R038_REPORT.md:9-17`; `phasea_cert3_registration.json`; `phasea_cert3_results.json`). The outputs contain exact rational primal/dual vectors and feasibility/objective intervals, and all three are labeled `H_fail`; the report says they satisfy original gap eligibility but fail the half-gap condition at the attained final loads. The independent exact LP audit now reports PASS for all three (`/Users/zhengzhihao/Desktop/T/cube-r030-num-audit/cube-dim-7-bounds/audit/round038/lp/README.md:3-13`, `audit_lp.json`): exact primal/dual nonnegativity, covering/packing constraints, objectives, bound directions, eligibility margin, and strict H-failure margin were checked. This is a bounded certificate audit, not a new LP solve or a global claim.

The 23/56 split based on old R037 scalar LP intervals is marked unverified historical diagnostic data in `R038_REPORT.md:9`; it is not 79 new LP certifications. The three selected cases do not support a conclusion about the remaining 76 patterns.

## Claims that do not follow

- R038 does not prove an improved operator bound (C_n). Its pair objective reduces finite payment-certificate load on a subset of old R037 patterns; the report states it asserts no operator bound (`R038_REPORT.md:3,5`).
- It does not show that pair balancing reaches the global LP optimum or that the cover family is globally optimal (SPEC §5, lines 84-86; report lines 5,7).
- The main ledger remains at `current_round: 36`, with `global optimality of cover family` unaccepted (`results/general_cover_ledger.json`). No R038 conclusion is entered in that ledger.
- The run is explicitly not a new geometry batch: it reuses 79 cached R037 patterns and does not generate 240 new geometries or rerun the original atoms/selection/baseline LP (`R038_REPORT.md:3,21`; SPEC lines 7-11). Describe this as a post-processing run on outputs from the earlier R037 pipeline; do not say the original pipeline was never run.
- No seven-dimensional lower-bound work is authorized or reported. SPEC line 11 excludes dimension-seven cases/lower-bound work; main `STATE.md:3,5,7` records `LOWER_PAUSED` and `N7_SPECIFIC_STOPPED`.

## Independent case audit and compact replay

The num-audit Phase A report (`/Users/zhengzhihao/Desktop/T/cube-r030-num-audit/cube-dim-7-bounds/audit/round038/phasea/README.md`; `summary.json`) records **79/79 PASS**, 3,827 accepted steps, 50 coordinate stops, 29 iteration caps, and 129,033 independently checked final nonedges across the 50 coordinate-stop cases. The independent checker reconstructed each witness graph and complete maximal-clique list, checked each recorded payment/load/energy step, and verified final vectors and stop tags. The intermediate pair choice was not checked against all possible pair rankings, so the audit does not certify that each chosen pair was globally lex-best. It also does not reconstruct old atoms, original histories, or old LP certificates. The initial three integer-string digit-limit failures were retried; final audit input hashes show no changed source cases.

The producer-authored `replay_compact.py` archival receipt reports PASS for all 79 compacted cases and 3,827 recorded steps, and validates exact recorded-step transitions, payment conservation, lexicographic decrease, final-q hashes, and the three LP certificates. The compact bundle maps the 79 original raw hashes and certificate hashes; it is lossless for the replay fields described in `R038_REPORT.md:27-29`, but is not byte-identical compression of the 71 MB raw outputs. This replay is standalone in the operational sense that it consumes the compact archive, but its script is producer-written, so it is not a second independent author. The independent original-case mathematical audit is the separate num-audit Phase A work described above. Replay did not request the optional final terminal scan; the independent audit did check all final nonedges for the 50 coordinate stops. Neither replay nor audit establishes pair-choice optimality along the path or a search-wide result.

## Current integration status

The main repository is at R037 commit `a3cd1967` as reported by the coordinator; the accepted upper bound remains R036, and no new `C_n` result is recorded. R038's improved delivery evidence strengthens the audit trail for the finite cached-pattern post-processing and three exact LP certificates, while leaving those mathematical scope limits intact. The earlier post-launch SPEC provenance caveat remains: R038 production inputs/code were registered before execution, but the subsequently completed SPEC itself was not hashed into that pre-run registration.

## Resource wording

`phasea_summary.json` records 79 tasks, 328.3319 s wall time, 1691.0232 s summed task CPU, 15 distinct worker PIDs, and maximum per-worker RSS 92,225,536 bytes (about 92 MB decimal). The report (`R038_REPORT.md:19`) describes a dynamic queue, an early sample of 15 active workers at 87.7–95.4% CPU, and later drainage to four busy workers. These are task/resource receipts; they do not establish 15-core saturation throughout the batch or justify claiming a sustained full-machine load. The report also says work was not extended to fill the machine.

## Delivery gap

The main ledger has no R038 entry; the three LP certificates now have independent exact-feasibility PASS, but no R038 conclusion has been accepted into the ledger. The review substantiates the finite-run counts, exact audit scope, replay scope, statuses, and stated resource receipts only; it does not elevate the mathematical conclusion beyond those records.
