# R011 deterministic truncated-kernel calibration and pilot

The fixed calibration and pilot completed at 256-bit precision, with no promotion and no new Monte Carlo sample. The following are complete symmetry-orbit enclosures for the original truncated kernel Hε; rounded values are shown here, and exact outward dyadic endpoints are in results.json.

| x | ε=1/64 | ε=1/256 |
|---|---:|---:|
| (1,0,0,0,0,0,0) | 3.3873083502829714e-5 | 1.3413989812390446e-4 |
| (1,1,1,1,1,1,1) | 1.8565697319391483e-6 | 3.7543925813411231e-4 |

All four raw intervals are positive, satisfy the expected cutoff monotonicity and the origin-peak checks, and have absolute width below 1e-12; their actual widths range from about 4.51e-72 to 1.72e-69. They are point values, not a continuous face lower bound g.

At the generic R003 exact-dyadic point of index 2710 and ε=1/1024, only half-phase IDs 0–255 were evaluated, as four lexicographic permutation chunks of 64. Their signed partial contribution is approximately 5.7177628981218954e-6. **This is not Hε**: the remaining 322304 half-phases were not evaluated and may cancel or add. The actual partial interval width is about 1.06e-71, below the planning threshold 7.94e-16, but this does not certify the width of unseen phases. No full generic evaluation was started.

The 46 exact primitive rays and 212 rational α weights were compiled from the accepted R010 fan without changing scale; q=ψ(ξ) was recomputed exactly from all 64 prefix forms, and q's least common multiple is 2514318975148780800. Exact A_sum is stored in results.json, with the rational sanity check 5040 A_sum≤1/512. The sharp ψ≥32√2|ξ| input is supported by the separately hashed primary 21-orbit functional receipt, which was read and verified rather than re-enumerated.

Four scalar divided-difference cases were checked using common exact integer knots: all zero, all one, partial confluence and seven distinct nearby rational knots. The real-Arb path overlaps an independently computed complex-DD path; all-zero and all-one formulas pass, and the near-distinct case also overlaps a complete-homogeneous expansion with 32 extra terms and its explicit rational tail. Equal nodes are detected only by exact integer equality. Three origin checks at ε=1/64,1/256,1/1024 agree with 840 A_sum/(π⁷ε⁶); the origin receives the full group factor once, with no opposite-pair factor.

The complete symmetry calculations use seven positive-axis phases with fibre 46080 for e1, and 64 sign phases with fibre 5040 for the all-ones point. The generic point has fibre 1; the partial sum includes the opposite partner factor 2 but no extrapolation to the full orbit. Per-phase top/sub-DD cache, log, Hermite, zero-knot and repeat counts are retained. Across 142 complete-symmetry phase evaluations there were 30,104 top requests, 26,778 unique tops, 228,724 unique sub-DDs, 4,022 logs and 1,690 confluent Hermite evaluations. Across 256 generic pilot phases there were 54,272 top requests, 54,272 unique tops, 517,440 unique sub-DDs, 11,776 logs and no Hermite or zero-knot events. These counts include both epsilons for each symmetry case. counts.json separately records the scalar/origin counts and fixed-input checker work.

Wall time was 0.670 s; summed worker-task CPU was 0.962 s and parent CPU 0.058 s. The largest task-process RSS high-water mark was 40,206,336 bytes. The three-second sampler caught only the parent because the run finished sooner; it does not establish a peak concurrency or system-wide memory figure. Sampled system swap usage was 6256.38 MiB. Fifteen workers were configured for eight tasks, but no saturation claim is made.

The checker verifies all source hashes, 212 exact weights and normalization, independently repeats the four fixed scalar complex checks and three origin constants, and reassembles the four symmetry outputs and generic partial from their saved per-phase balls without evaluating another phase. It also checks orbit fibres, phase IDs, width targets and the primary sharp-floor receipt. The completed calibration establishes no global g, I or C7 result, and a complete generic run requires a separate release.
