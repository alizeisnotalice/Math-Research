# C160 additive v2 engineering-recovery static review

Verdict: PASS for Terra's one engineering fixture plus missing-metadata recovery command only. This is not a rerun of v1 preparation, a mathematical fixture, or a PhaseA science release. Actual fixture/recovery receipts must be inspected before a new concrete v2 science gate. The original numerical model and run_phaseA.py bytes remain unchanged.

The earlier code-only review missed the function-local importlib binding defect. Static inspection now confirms it: import importlib.util inside v1 phaseA makes importlib local throughout that function, so the earlier importlib.metadata.version access raises UnboundLocalError. The additive v2 aliases remove that lexical ambiguity from the missing metadata suffix. This correction does not weaken the original quantization/model proof, and the original failed source/evidence remains preserved.

## Exact delta and preserved evidence

Reviewed V2_PACKAGE_MANIFEST SHA256: a70ff961fdffa282d1126820e865caae09e3b6c1c16153726e46f84aa83f7310, with exactly seven new members. Their names/hashes/bytes and the original seven-file package hashes match. All11 partial snapshot members match the real Terra root. In particular:

- original PACKAGE_MANIFEST:7719ad8534d9f3bb4048bb7bf2c0f4960ef76e7f6be6f8677816a39412e378ab;
- existing quantization certificate:737479ae78c4a841a77c9af23ec9b9f85a70d19d976ba0c6c65b6029d26159ad;
- existing inputs:b74e867ab92449f8604a5f0dc07a2ce9896163410f2b62444b943ed5a25250aa;
- real operator failure observation:dbf3de868713d0252298eaa0be13829f14e647a4f6fa87a96f960a74de78f1c6.

The real phaseA directory is still empty at this review; no prepare_v2 directory or fixture receipt exists yet. The observed failure file is explicitly a post-hoc record containing a truncated excerpt, not complete raw stdout/stderr. It reports the separate package name/path precheck error and missing shell failure barrier. No full traceback, original PID, duration or process-state claim is reconstructed from that excerpt. The additive contract correctly requires acknowledgement that the original capture is incomplete.

## Metadata recovery scope

repair_prepare_v2 uses direct global importlib.metadata/util aliases. It never calls original prepare_for_terra.phaseA. It restricts the target to the exact frozen partial root, checks the original source/input/quantization/failure-observation hashes, requires the observed empty phaseA directory and a fresh prepare_v2 marker, and copies only new versioned file names. Existing original files are not overwritten. A genuine recovery failure preserves the additive suffix and stops without retry.

The only arithmetic consumer checks during recovery use the ALREADY SAVED N,D and exact full source B=n/d:0<nD, nD<=Nd, Nd<nD+d, plus D=2^96 and exact saved cap identities. These three cross-product inequalities are the positive/upward/less-than-one-grid-step witness, equivalent to the original ceiling property. The code does not compute N again, perform ceil, choose another cap or compare an old M/H objective. Source B is matched by the complete exact string/selector and frozen full156 result hash. Existing inputs/quantization bytes are checked again at completion.

Only package-version metadata, phaseA registration/release requirements and additive prepare_v2 start/exit records are newly written. The numerical scope, paths, h values/caps, warm start,49 features,150-call limit and all original mathematical checks remain unchanged. A partial metadata failure excludes no mathematical family.

## Engineering fixture and fail-fast wrapper

The prepared engineering fixture checks symbol-table binding, then calls the corrected metadata/util helper functions with local engineering stubs. Its util loader stub sets a namespace flag; it does not load/execute common's model preflight or any mathematical routine. The fixture also calls the actual operator main against a deliberately invalid package row using path instead of name, with a recording subprocess stub, and requires rejection with zero child calls. Temporary fixture data is isolated, sys.argv restored and the stub module assignment is not a production-source edit.

operator_prepare_v2 distinguishes package files[].name from external-source path rows. It validates BOTH complete source packages before any child, then executes the engineering fixture and recovery with explicit argv and check=True in sequence. A schema/validation exception starts no child; a fixture failure prevents the recovery child. There is no shell concatenation, fallback or automatic science command. The new repair requires the actual fixture PASS receipt and binds its repair/operator source hashes; source V2 hashes bind the fixture itself. This fixes both the lexical defect and the independent wrapper failure-barrier error without replaying the original preparation prefix.

Sol has inspected this control flow only. Neither symtable fixture/stub execution nor the three saved-witness arithmetic checks has been run by this reviewer. Terra's actual receipts remain required evidence.

## New PhaseA gate and actual156 acceptance

run_phaseA_v2 adds only recovery/fixture/output identity checks and the new exact-hash contract, then invokes the unchanged original run_phaseA.run function. The old mathematical runner still enforces its source156 saved acceptance, exact model/cap checks, exclusive science marker and one<=150-call selector. New contract fields bind V2 manifest/entrypoint, partial input/quantization bytes and the failed-observation acknowledgement, so old approval alone is insufficient for the authorized v2 entry.

The sole authorized subsequent science entry must be run_phaseA_v2.py; direct invocation of the preserved unversioned runner/prepare is not authorized. The original runner remains on disk for mathematical byte preservation, not as an alternative execution route. PhaseB remains separately selected/prepared/released later.

C156 is no longer merely a producer candidate pending a preflight. This review read the actual main acceptance file SHA b4743d20f78e1f4320aeeb602e84d803e6a2fac16c899633219012b9327f9381 and actual independent saved seal SHA18eca9a3a6b33c987952e42244221286efb143c038c54335a5bb02125e9713e4. Their structured status/acceptance flags bind the same frozen full156 producer result f9a0689da14fb34d48c3ffb65b6abaef06193906ff71cc3cab77420074bee21f. Main acceptance evidence-file hashes match. No156 mathematics was rerun here. These actual saved records, not this code-only review, must populate the unchanged science acceptance gate.

No science is authorized by the fixture/recovery approval. The future concrete PhaseA release must include its actual new registration/requirements/fixture/recovery evidence plus the genuine156 acceptance records. Candidate floats still do not announce M/global-g/C6 bounds; the earlier mathematical scope limits remain unchanged.

## Reviewed new file pins

| File | SHA256 |
|---|---|
| repair_prepare_v2.py |20f3771e62d18a2ac4e0c82ca88874889a64e0123ca7fa83a1c69b4f9204804f|
| run_phaseA_v2.py |8636a5f55dccda427fc17d2ed76e0b062eb731264e1af872a69fb810174c02d2|
| operator_prepare_v2.py |61b1fa47b79ad34697106ffb46f8c47f0bc23685f5334b7a7aa49833f9278cbd|
| scope_regression_fixture_v2.py |7ab8f037648111ebc1e4967ca20798e45f5b004dfe901c591bdd9712b3358ff6|
| PARTIAL_FAILURE_OBSERVATION.json |b1c9a59d967c617efc3beb216a306291b4ca3c34eae0a79208278961944339cd|
| V2_FAILURE_AND_DELTA.md |6ad30500dd9427fedd237dbf3d4e1e767a4ea1422146d84ce475e198b9831a33|
| HANDOFF_V2.md |57d50d4200e82dec790234c6db753af261c6b873e0a3d946455e5f18fab6f73a|

Only source/proof reads, AST syntax parsing, JSON identity/schema and hash/byte checks were performed. No engineering fixture, symtable/stub run, recovery, existing-witness arithmetic, ceil, optimization, model/H/root/native evaluation or old functional was executed. This reviewer changed only its own audit/upper160_preflight append-only v2 report files. Old audit/source/partial evidence bytes remain unchanged; no stage, commit, push or recursive agent occurred.
