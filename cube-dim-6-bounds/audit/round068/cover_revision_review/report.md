# F068: core/math PASS; actual release and runtime enforcement remain pending

Final scope decisions:

- **Core/math: ACCEPT.** The repaired state history and mathematical importer satisfy the reviewed contract.
- **Runtime: CONDITIONAL / NEEDS_WRAPPER.** Main keeps the real production gate closed until a separately sealed hard-wall-time supplement is reviewed. The original C068 seal is unchanged; no passed mathematical configuration is rerun for that future runtime-only review.

The repaired event-replayed core and mathematical importer pass independent review against preparation manifest `a3b254ce4ba9b833ba369fa72c65b457e648d563b10c61bb969199a087e59db3`. The original R1/R2/R3 counterexamples no longer pass the actual API: root/environment substitution is rejected, stale/orphan history cannot reset reservations, and test history cannot acquire a production release.

Nine independent configurations completed once in 0.3835530420183204 seconds, under a 180-second hard alarm. They included a mixed-cutoff eight-attempt/six-node history, omitted sibling and counterfeit passed leaf, mixed-sign six-gradient intervals, corruption of the normal derivative, replacement of original artifacts after verification, rehashed false exponential endpoints, and wrong release/source/window/PID bindings in a production-contract unit branch. All expected rejections occurred.

The actual pinned A057 Arb384 checker ran 31 times, covering 992 labelled corners. No native kernel ran. The author's eight tests, old C063 writer and old twelve formal tests were not executed. All 43 pinned files were unchanged. No author or previous-audit file was modified; nothing was staged, committed or pushed.

There is one explicit operational qualification before production: `run_native_once.py` disables its alarm before an unbounded cleanup `child.wait()`, and `control.py` disables its alarm before failure-receipt writing. The code supplies a 180-second work deadline but does not alone guarantee that cleanup-inclusive process wall time is at most 180 seconds. Main must provide independent outer hard enforcement or a reviewed bounded cleanup arrangement. No real timeout/kill overrun was observed or tested, and no hypothetical overrun is being reported as an actual failure.

The real state is still PREPARE_ONLY genesis: no release, attempt, inflight job or native call. The unit fixtures—including the production-specific branch test—are not actual native receipts or a body32 cover. F061's accepted finite gradient checks are bound as source evidence, with no universal future-input rerun claim. Genuine launch/kill handling, actual main release, authentic new native results and final global proof remain outside this review's completed execution.

See `proof_review.md` for the state induction, mathematical binding analysis, exact runtime source locations and release boundary. `results.json` records all nine configurations and rejection messages. The independent reproducer is `verify_revision.py`; it uses the frozen new modules and actual A057 checker, with newly created temporary fixtures only under this directory.

Replay in a clean copy without results.json:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B verify_revision.py
```

The final receipt and manifest bind the execution and artifacts. Core/math PASS does not authorize a production batch, certify global g or establish a new C6 bound.
