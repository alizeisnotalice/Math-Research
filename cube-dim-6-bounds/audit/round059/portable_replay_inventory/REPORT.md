# Round 059 — portable replay inventory (static review)

## Scope and evidence

This is a static source inventory of 20 final verifier/review files associated with rounds 050–056. The finite search commands, UTC receipt/exit times, source paths, SHA-256 values, roles, outputs, path dependencies and wrapper proposals are recorded in [MANIFEST.md](MANIFEST.md). The inventory uses Python standard-library AST parsing, targeted text reads and SHA-256. No source script, checker, writer, native binary, kernel, matrix recomputation, or MATLAB session was run.

The prior independent clone’s HEAD and untracked `cube-dim-6-bounds/` were retained. Only this round’s `portable_replay_inventory/` directory was created and written.

## Static findings

- Five absolute-path literals occur in four source files (one in round 051; three in round 054 across its two scripts; one in round 056). Each literal is used to construct a `Path` that is later read. They are runtime dependencies, not provenance-only strings.
- Many other checkers read filesystem paths stored in JSON registrations, manifests, receipts, release records, and build records. These are runtime path inputs, but their values were not inspected, so this review does not claim each currently contains an absolute path.
- The round-050 checker and independent fan review both traverse the active ledgers and recompute determinants. The source describes up to 783,994 active systems. The round-050 review also repeats fan and cell structure work. Their “no kernel/no solve” labels do not mean they are cheap saved-data/hash checks.
- The round-055 structural review imports a helper that recomputes two exact Gram determinants for 11,476 pairs (22,952 matrix operations), then checks frame/group/fan structure. It does not call the full G/H kernel. The round-055 closure verifier checks the saved operation count and hashes without repeating that calculation.
- Round 052’s 12-case formal scalar check constructs and inverts 6×6 SymPy matrices. This is a bounded symbolic/matrix replay, but not a full production kernel. Round 053’s formal checker also performs symbolic algebra and small exact elimination. Rounds 054 and 056 perform bounded exact scalar/root checks, not full matrix/kernel passes.
- Several “check/replay” scripts write `checked_results.json`, `results.json`, receipts, or failure records. Guards vary: some refuse existing outputs; others only guard one receipt or no output at all. The manifest records the visible guard and write target. None was run, so no output was created or overwritten by this audit.

## Suggested portable wrapper shape

Create new wrappers outside the frozen source files. Each wrapper should accept explicit `--bundle-root`/`--author-root`/`--predecessor-root` mappings, resolve manifest-listed paths under those roots, and write any new receipt/result only into a separately selected `--output-dir`. Give compute-heavy modes explicit names and budgets: round 050 ledger/fan replay and round 055 structural replay should be opt-in and report their operation caps before work. Keep pure saved-evidence verifiers read-only. Do not modify the old source scripts to achieve this.

## Limits

These are static findings from the listed source files only. No manifest/registration JSON values were checked, no path was resolved by executing the scripts, no writer behavior was exercised, and no timing was measured. “Portable wrapper suggested” does not mean a wrapper exists. Nothing here verifies certificate contents or establishes any mathematical claim; portability is separate from mathematical correctness.
