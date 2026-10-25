# A038 bounded Green floor diagnostic

Started UTC 2026-10-25 20:50:39; absolute deadline 21:00:39, 600 seconds including preparation. Own only this directory. All previous sources frozen.

Reuse B032 complete saved G and six gradient dyadic intervals at w; never recompute that seed. Use frozen B030 proposal: exact midpoint free gradient normalized by maximum absolute free component, negative step 1/20, ordered face gap 1e-8 clipping, nearest-even denominator 1e14, then only the predeclared inward gap repair k=1,...,6 if an exact Fourier node is zero. No optimizer, grid, changed step or seed, H/K/p, MATLAB, recompilation, delegation, stage, commit or push.

At most three new full C013 G+six-gradient calls, each at most120 seconds and bounded by the total absolute deadline. Before every call save the exact input, full preceding endpoint source/hash, preflight and command. Require all9308 primitive zero forms nonzero, int64-safe inputs, all967680 logical terms, Arb512 precise endpoints and Euler checks. Stop at the first reliable G upper endpoint <1/65000, any unsupported candidate, duplicate, hash mismatch, failure, timeout, or after call3. No retry.

This batch diagnoses the new conditional g>=1/65000 target only. A surviving finite list is not a global lower bound. A strict counterexample would not exclude every I/g configuration or change C6.
