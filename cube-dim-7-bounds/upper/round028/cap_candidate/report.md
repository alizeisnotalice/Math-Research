# R028 local cap candidate

The bounded SciPy exploration on [4,10] completed after 15 new double evaluations and zero cache hits; raw IEEE inputs and results are recorded in trace.jsonl. The five initial cached values are double conversions of the exact midpoints of the complete R027 intervals. This exploratory result does not certify a continuous minimum.

The returned t is exactly the IEEE float 0x1.9e0829abe1420p+2. Rounding its exact binary rational times 2^20 gives integer 6783498, hence the exact rational t=3391749/524288. One new 512-bit Arb row gives F(t) approximately 24.55496052596470030163. Its exact outward upper endpoint is strictly below the full R027 winner lower endpoint; upward rounding yields the certified C7 < 24.554960525965, subject to the R027 cube-cap theorem. There is no global optimality or C7 lower-bound claim.

The primary run used 0.008756 s wall and 0.004700 s CPU, with peak RSS 84,377,600 bytes including loaded libraries. The separately reported checker repeated the one exact candidate, not the exploratory points, and passed. Source and code were locked before running; no neighbors, extra precision, phase, MC, MATLAB, or lower search were used.
