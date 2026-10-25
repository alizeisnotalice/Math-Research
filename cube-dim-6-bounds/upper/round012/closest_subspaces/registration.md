# A012: exact closest-subspace geometry, finite independent test

Only this directory is writable. Keep every completed source directory
unchanged. No kernel, density, MC, grid, MATLAB, agents, stage, commit or push.

Prove or refute the five contiguous two-cluster formula for body32 S2 and
the proposed eight-term formula for the augmented frame V union {beta e_i}.
Handle full sign/permutation symmetry, zeros, ties and nonunique nearest
planes. Do not assume plane tubes are disjoint or that three augmented
directions are independent. Geometry only; no G/K finite/singular theorem
or C6 certificate is inferred.

Finite batch: exactly twelve rational vectors, including zero, equal values,
signed values, coordinate axes, a two-axis vector, repeated levels, two
generic controls, a body-axis vector, and the two read-only exact B010 inputs.
For every point compare the formulas with all496 body-body pairs and all703
augmented38-direction pairs, using exact Fraction Gram projection arithmetic.
The base496 comparison may reuse the corresponding subset of the703 pairs.
Use one fixed rational beta=3/7 in the complete38-direction check; prove
scaling invariance for every beta>0 analytically. No extra test point.

Record exact formulas, all input/script hashes, the twelve vector definitions,
the B010 distance-square comparisons and actual runtime. Stop on a mismatch
and preserve its exact counterexample; otherwise stop after the twelve checks.
Do not extend to another beta, another grid or kernel evaluation.
