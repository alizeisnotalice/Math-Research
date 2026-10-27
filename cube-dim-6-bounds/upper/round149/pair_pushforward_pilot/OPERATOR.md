# C149 operator handoff — no science release yet

Designated execution directory:
`/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round149/pair_pushforward_pilot`

Check the target first. Do not overwrite existing files or outputs. Copy only the exact source files listed in `execution_whitelist.json`, preserving bytes. Their P-relative writes land in the target; every old scientific reference is an absolute read-only path. No source directory move or C147 edit is needed.

Phase A, only after main authorizes the code-only whitelist: run `prepare.py` once with the existing finite-dimension Python. It does Fraction/hash/canonical source bookkeeping, not density algebra. Preserve its terminal receipt. It writes case_packet.json, source_manifest.json, registration.json and prepare_receipt.json and refuses a preexisting packet/start/release. Do not automatically rerun on error. This phase is explicitly not a science release.

Return the frozen prepared packet and hashes to main and the independent reviewer. They may run the saved-only preflight and inspect the exact selection. Do not create parent_release.json yourself. Do not call the bare run_prepared.py science entry until main provides its concrete dual-approved exact-hash release, including workers between1 and15.

Phase B is queued behind round148 native work and any main resource allocation. Main owns when and with how many workers it starts. Execute one released batch, one-thread environment, nice>=10, no wall/CPU timeout, no automatic restart. Capture original controller PID, terminal exit/reap, complete child CPU, stdout/stderr hashes and release hash. Monitor original processes and checkpoints without reissuing failed/unknown cases.

The allowed science entry compiles only the selected NEW density laws. There are no k_h values, special-function integrations, old6D grids, beta2/C6 computations, source edits, stage, push, recursive agents, or expansion of the selected set. A stop/failure preserves all per-case checkpoints and receipts for main. Keep operator outputs in the designated directory.
