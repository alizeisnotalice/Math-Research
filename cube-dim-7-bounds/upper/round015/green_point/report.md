# R015 true Green potential at one exact off-S3 point

For the exact R003 index-2710 point recorded in testcases.json, the original potential G—not a positive-cutoff approximation—satisfies

    0.000003782255829017294943880908 < G(x)
                                     < 0.000003782255829017294943880909.

The actual outward dyadic interval in result.json has width about 3.3054e-68. Its upper endpoint is strictly below each of 1/219000, 1/194000 and 1/179000 by exact rational comparison. The first coordinate is 1, all six remaining coordinates lie strictly between 0 and 1, and all seven absolute coordinates are distinct exact dyadic numbers. A span of at most three body directions has at most four distinct absolute coordinate values, so x is off S3. R001 proves G is continuous locally off S3; hence the strict upper inequality holds on a relative-open face neighborhood of positive face measure. Each of the three proposed global or essential face lower targets is therefore false. This is not an operator C7 counterexample or lower bound.

The pointwise computation is genuinely at zero time cutoff. The completed R012 whole half-orbit had zero exact integer ray knots in all 322560 phases × 46 rays. The same source geometry and exact point were verified by hash, and every R015 phase tested all 46 knots nonzero again. Consequently every cone's time integral is absolutely integrable near zero; the finite full cone/group sum permits the termwise ε↓0 limit. Monotone convergence of positive p_t gives the original G. The implementation evaluates real divided differences of log|a| at exact common-denominator integer knots; a repeated nonzero knot uses (-1)^(k−1)/(k a^k), and the B^6 scaling is applied exactly once. No zero node was deleted, clipped, or replaced by a small ε. As an independent consistency check, the R012 H_(1/1024) upper endpoint is strictly below the present G lower endpoint.

All 5040 new G chunks and 322560 half-phases completed once, including the first four G chunks as pilot; no R012 H interval was reused as a G value. The pilot passed its planning-width gate. The 256-bit full sum is positive and meets the 1e-12 width target. It comprises 68,382,720 top DD requests, 659,175,346 unique sub-DD computations and 14,837,760 log evaluations; zero knots, repeated top nodes and Hermite evaluations were all zero. The signed normalized chunk intervals are summed once without an extra group factor. Source hashes, exact per-chunk intervals and counts are retained.

For the particular fixed complete-I budget I₀=7723627498344/10^18, the exact rational test

    7^7 I₀^6 / G_upper^7 > 128^6

holds; its sixth root is approximately 158.38947. Every legally certified face floor g is at most this point's G value, so the specific bound-expression using this fixed numerator cannot yield a result below 128. This statement does not rule out a sharper I budget or a different proof, and **does not assert C7≥128 or C7≥158.38947**.

Run wall time was 126.925 s, with 869.716 s summed chunk CPU. The four pilot chunks ran in the parent and account for 0.670 s of that sum; excluding this overlap, worker plus parent CPU was approximately 873.074 s, or 6.88 core-equivalents over wall. Twenty-four worker snapshots showed 15 workers totaling 603–720% CPU and summed RSS 354272–494832 KiB, which can double-count shared pages. Individual peak RSS was 39,550,976 bytes; sampled system swap use ranged about 7173–7205 MiB. These are observed shared-machine figures, not full-occupancy or system-wide peak claims.

The checker verifies source hashes, the preceding R012 all-phase nonzero-node receipt, the exact point and off-S3 distinctness, all 5040 new G chunk terminals, the exact dyadic outward sum, the H<G consistency, and the three rational thresholds and fixed-I₀ power comparison. It does not rerun the 322560 phases or independently prove the accepted fan and limit identities. No new physical point, Green grid, whole-face cover, MC, old shape cache, MATLAB, commit or push was used.
