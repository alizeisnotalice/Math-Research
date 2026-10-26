# Round 091 n5 comparison erratum scope

This additive correction does not alter the sealed 40-source `SCOPE.md`, source manifest, inventory, or receipts. It uses exactly four additional read-only originals (maximum four):

1. C6 main-clone current guidance: `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/coordination/guidance_20261026/next_steps_original.md`.
2. Current n5 on-disk ledger: `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/results/bounds_ledger.json`.
3. Its complete certified lower-construction interval source: `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/lower_round003/a2_equals_a3_results.json`.
4. Main C6 saved comparison: `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/lower/round066/primary/n5_dynamic_check.json`.

One Python standard-library static hash/field-extraction pass, max 180 seconds. It may read these four sources only, hash their bytes, compare the n5 certificate hash to the existing C6 dynamic-check binding, and extract existing JSON values/text. No certificate checker, source code, root/kernel/scalar/matrix recomputation, network, decompression, MATLAB, stage, or push. The C6 075 reliable lower endpoint and source hash are reused from the sealed 40-source manifest; no reread. Record n5 git HEAD and tracked/untracked dirtiness from git metadata separately; do not interpret in-progress files as accepted certificates.
