# B011: one justified midpoint cutoff, at most two H passes

Owner B Astra; only upper/round011/cutoff48/ writable. B008–B010 remain frozen. No MATLAB, agents, staging/commit/push, gradient, new point, new G, other T, grid or Monte Carlo.

Exactly one new fixed cutoff T=6/125=.048, the single midpoint between the already refuted fixed times .032 and .064. This is a necessary-condition test, not a time grid. At most two full H_.048 passes; maximum elapsed batch time1200seconds.

Point order and stopping, fixed before evaluation:
1. The exact second B010 point from point2.json/batch2_results.json. Reconstruct its G from the saved exact rational/pi^5 certificate; never evaluate G again. Compute H_.048 with the audited real Arb384 full-fan evaluator. If the upper endpoint of true K_T=G-H_T is <=T/4096=3/256000, stop immediately and refute this fixed T's strict-64 target.
2. Only if point1 does not refute the strict-64 target, compute one H_.048 at the old A007 reference center, likewise reusing its saved exact G. Stop after this second necessary condition regardless of result.

The sufficient target m0=T/4000=3/250000=1.2e-5 is separate: a true K_T<m0 refutes the62.5target, but if K_T>T/4096 it does not refute all<64 improvements and the registered second point is still allowed. Any surviving point tests are not a global minimum or full-cube certificate.

Before kernels: freeze source/input SHA, exact coordinates, existing G arithmetic, and inspect available matching point/time records for duplicates. Every full H pass must account for all967680 group/cone terms, including real confluent/zero handling; exact sign pairing may count483840 evaluated terms with multiplicity2. No tolerance merging. Python bytecode writes disabled before imports. Primary separately controls any cross-implementation audit; B does not launch it.
