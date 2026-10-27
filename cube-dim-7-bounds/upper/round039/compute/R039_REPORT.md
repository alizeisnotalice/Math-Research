# R039: one-step pair scans with (L, Σq, E)

This run reads exactly three frozen R038 final payment vectors and their complete cliques. It performs six scans: each vector once over nonedges and once over all unordered pairs. Every scan selects at most one best step and stops. No old geometry, atomization, selection or LP was repeated, and no seven-dimensional calculation or new operator bound is claimed.

For each pair, the complete cliques are partitioned into the four membership classes 00, 10, 01, 11. On the rational payment interval, each nonempty class yields one affine load branch. Evaluating endpoints and all in-domain crossings gives the exact minimum load and its full optimal interval. The scan then minimizes S=Σq over that interval, breaking a remaining tie by E=Σh q². It compares the exact global merit (L,S,E) and only accepts a strict decrease; nonnegative payments and Σhq=U are checked after the chosen step. A zero-payment pair has no possible changed state and is counted but need not be evaluated.

| R038 final state | Domain | Pairs scanned | Best pair | Outcome |
|---|---|---:|---|---|
| C/2/5 full | nonedge | 3060 | (0,7) | L unchanged, S decreases, E increases |
| C/2/5 full | all | 3160 | (0,7) | L unchanged, S decreases, E increases |
| C/1/6 full | nonedge | 1840 | (24,55) | L unchanged, S decreases, E increases |
| C/1/6 full | all | 2016 | (24,55) | L unchanged, S decreases, E increases |
| A/3/15 partial | nonedge | 3969 | (66,79) | L unchanged, S decreases, E increases |
| A/3/15 partial | all | 4560 | (3,11) | L unchanged, S decreases, E increases |

Thus all three fixed states admit a strict one-step improvement of the revised merit even over nonedges, but no scan found a strict immediate load decrease. This does not establish that a later sequence lowers L, that the state is globally optimal, or that a multi-coordinate move is necessary. Every claim in the table is an exact arithmetic result from the frozen complete clique list, subject to independent audit. The best-pair record for each scan includes its rational branch coefficients, candidate abscissae, entire load-optimal interval, new coefficients, exact S/E changes, final q and payment check.

The seed inputs and code were registered before execution in `registration.json` (SHA-256 `cc24240275ff240cad9fe18665fb817554de321a3c019cf6270ee803a910a639`), matching the three cached input hashes in the soft proof. The single-pair code SHA is `372a1e11eb2969399e5b75ec73bc67b048c76cf33ff6f66ca0a7f9f07112560e`. The three-node equal-h P3 example was checked separately as a synthetic algorithm control; it is not an additional production geometry and was not used in the six scan statistics.

Actual parent PID was 70895, six worker PIDs 70904–70909, overall wall time 12.42 s. At one five-second sample, the two A-case workers were near 99.7% CPU and four shorter C scans had completed. Individual CPU, wall, RSS and exact output hashes are in `summary.json`. No iteration, grid extension or extra load was run. The theoretical all-pair formula and scope are in the independent soft-agent R039 proof; an independent numerical audit of the six records is pending.
