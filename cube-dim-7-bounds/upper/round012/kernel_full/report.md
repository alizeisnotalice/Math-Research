# R012 complete generic truncated-kernel point

For the exact R003 index-2710 dyadic point specified in testcases.json, at ε=1/1024, the full signed-permutation calculation gives the strictly positive enclosure

    0.000003778291519696606186963321 < Hε(x)
                                  < 0.000003778291519696606186963322.

The actual outward dyadic interval, substantially narrower than the displayed decimal brackets, has width approximately 3.3053e-68, below the prescribed 1e-12 target. This is a **single point** of the true truncated kernel; it is not a lower bound on the whole face, nor a g, I or C7 certificate.

The accepted R011 half-phase chunks 0–3 were copied verbatim and verified by hash; their 256 phases were not recomputed. The remaining 5036 chunks were evaluated once, each with 64 lexicographic sign phases, giving all 5040 chunks and 322560 opposite-pair representatives. Every chunk interval was already normalized by 2/(2π)^7. The final assembly adds those signed intervals once and applies no second group or normalization factor. The finite complete-fan and phase-orbit identities are supplied by the independently checked R010/R011 proof sources; all 212 simplex contributions are retained at each phase.

The full accounting is 68,382,720 top DD requests, all unique within their exact phase contexts; 659,175,346 unique sub-DDs, 697,882,546 cache hits, and 14,837,760 real log evaluations. This point has no zero ray knots, repeated top nodes or Hermite evaluations in the full half-orbit, as recorded by the executed finite checks. Counts include the four reused pilot chunks. Knot sorting and differences use exact common-denominator integers, never a near-equality threshold. The arithmetic kernel remains the pre-hashed 256-bit R011 real-Arb implementation; no automatic 512-bit full rerun took place.

The run used 133.679 s wall time, 898.730 s summed new-chunk CPU and 3.333 s parent CPU, or about 6.75 core-equivalents combined over wall. Twenty-six worker-phase snapshots each showed 15 child workers, with summed process CPU 549–714% and summed RSS 354992–491968 KiB; summed RSS can count shared pages repeatedly. Maximum individual worker RSS was 40,550,400 bytes. The sampled system swap use ranged approximately 6615–6693 MiB. The machine was shared with other work, so these figures do not imply full-core occupancy or a system-wide peak.

The checker verifies source hashes and the 5040 unique chunk terminals, checks that reused intervals are exactly the R011 records, validates phase ranges/counts, and independently sums all chunk lower and upper dyadic endpoints as exact rational numbers. The final reported interval encloses that exact sum. It does not re-evaluate any phase or by itself prove the fan or Fourier identities; those depend on the accepted structural proof and pilot. No ε=0 limit, new physical point, full-face cover, Monte Carlo, old shape cache, MATLAB, commit or push was used.
