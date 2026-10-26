# Round 091 cost-source inventory scope

- Writer clone: `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics`; sole write area is `cube-dim-6-bounds/audit/round091/cost_source_inventory/`.
- Main six-dimensional source root: `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds`, observed HEAD `c03fed2f1590b71ea036f8e5da94765065418aaa`. Its worktree was observed dirty; the manifest hashes the on-disk bytes at scan time, not an assertion that each file is committed at HEAD.
- Five-dimensional source root: `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds`, observed HEAD `b6acac941cb0dc16c3fa87460de711121779cfdb`. Its working tree is dirty; hashes below bind current on-disk STATE/ledger, not a claim that those edits are committed.
- Scope ceiling: exactly the 40 read-only paths enumerated below. Includes C6 README/STATE/ledger and existing round084 registry/manifest; selected existing C6 reports/data for rounds 032, 052, 059, 061, 074 and the complete 075 review; plus n5 README omitted by design, with only current n5 STATE and ledger included. No other text/JSON source may be opened for this inventory.
- One Python standard-library static read/hash/extraction program, max 180s including cleanup. It may parse registered JSON and extract verbatim fields only; it must not import project/science code, execute an old writer/checker, evaluate roots/kernels/scalars/matrices, decompress bundles, or call MATLAB/network.
- Extract exact stored intervals and source JSON pointers. Do not reconstruct intervals from display decimals. Separate global g lower, pointwise g upper, full-face I upper/lower, and model-cost diagnostic. Model identity must match the fixed α=1, 32-direction Cauchy body32 (q=6/5) frame.
- n5 is read-only context: record its HEAD and file hashes, compare its formal strongest lower against the C6 complete upper 64. No n5 upper inheritance and no mathematical recalculation. The 078 literature conclusion is already embedded in the registered C6 state/084 context; no network or new literature search.
- No staging, commit, push, MATLAB, recursion, or writes outside this directory.

## Exact source allowlist (40)

C6 main clone (except the two following numerics-clone files):
1. `README.md`
2. `STATE.md`
3. `results/bounds_ledger.json`
4. `upper/round032/green_target_diagnostic/report.md`
5. `upper/round032/green_target_diagnostic/checked_results.json`
6. `upper/round032/green_target_diagnostic/results.json`
7. `upper/round032/green_target_diagnostic/candidate_01.json`
8. `upper/round032/green_target_diagnostic/manifest.json`
9. `upper/round052/zero_node_green_finite_part/theorem.md`
10. `upper/round052/zero_node_green_finite_part/report.md`
11. `audit/round052/zero_node_finite_part_review/report.md`
12. `audit/round052/zero_node_finite_part_review/results.json`
13. `upper/round059/global_region_account/theorem.md`
14. `upper/round059/global_region_account/report.md`
15. `upper/round059/global_region_account/source_index.json`
16. `upper/round059/global_region_account/actual_run/results.json`
17. `audit/round059/global_region_account_review/report.md`
18. `audit/round059/global_region_account_review/results.json`
19. `upper/round061/cutoff_gradient_pilot/pilot_report.md`
20. `upper/round061/cutoff_gradient_pilot/production_release.json`
21. `upper/round061/cutoff_gradient_pilot/production/low-primary/result.json`
22. `upper/round061/cutoff_gradient_pilot/production/low-primary/box_1.json`
23. `upper/round061/cutoff_gradient_pilot/production/low-primary/box_2.json`
24. `upper/round061/cutoff_gradient_pilot/production/low-primary/jet.json`
25. `upper/round061/cutoff_gradient_pilot/production/S3-primary/result.json`
26. `upper/round061/cutoff_gradient_pilot/production/S3-primary/box_1.json`
27. `upper/round061/cutoff_gradient_pilot/production/S3-primary/box_2.json`
28. `upper/round061/cutoff_gradient_pilot/production/S3-primary/jet.json`
29. `audit/round061/cutoff_gradient_pilot_review/report.md`
30. `upper/round074/middle_long_budget/theorem.md`
31. `upper/round074/middle_long_budget/report.md`
32. `upper/round074/middle_long_budget/results.json`
33. `audit/round074/middle_long_budget_review/report.md`
34. `audit/round074/middle_long_budget_review/results.json`
35. `audit/round075/post_s8_overlap_review/report.md`
36. `audit/round075/post_s8_overlap_review/review_results.json`

Numerics clone, round084 registry (counted within the same 40-source ceiling):
37. `cube-dim-6-bounds/audit/round084/upper_mechanism_registry/registry.md`
38. `cube-dim-6-bounds/audit/round084/upper_mechanism_registry/source_manifest.json`

n5 main clone:
39. `STATE.md`
40. `results/bounds_ledger.json`
