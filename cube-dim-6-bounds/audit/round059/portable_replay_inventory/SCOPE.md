# Round 059 portable replay inventory — scope registration

- Registered: 2026-10-26 09:03:10 UTC (shell `date -u`).
- Writer clone: `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics` (HEAD `d70bd06`; retained as-is; pre-existing untracked `cube-dim-6-bounds/` retained).
- Sole write area: `cube-dim-6-bounds/audit/round059/portable_replay_inventory/`.
- Read-only sources: `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds` and sibling author trees `Math-Research-cube6-upper-astra`, `Math-Research-cube6-audit`, `Math-Research-cube6-curvature-astra`, `Math-Research-cube6-literature` under `/Users/zhengzhihao/Desktop/T`.
- Registered batch: final minimal replay/check scripts associated with rounds 050–056 only. Exclude construction/generation pipelines, broad archival crawls, and any unrelated round files. Hard cap: 40 source files total.
- Search budget: at most two batch source-enumeration/search passes after this registration. Record exact commands, receipt and exit times, and counts in the manifest.
- Method: source code is read only as text/AST; compute SHA-256. Classify runtime absolute-path access, provenance-only absolute-path strings, writer behavior (including overwrite/create targets), and replay operations that visibly repeat expensive matrix/kernel work. Static indicators are leads, not execution evidence.
- Prohibited: run any old writer/checker/replay; execute kernels or scalar/matrix recomputation; MATLAB; modify/checkout/update source clones; alter old frozen scripts; stage/commit/push; recursive delegation.
- Deliverables: manifest and report in this directory. Findings must label static discovery vs unverified behavior and must not imply that portability establishes mathematical validity.
