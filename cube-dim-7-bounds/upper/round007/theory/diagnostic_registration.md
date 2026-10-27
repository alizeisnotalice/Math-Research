# R007 matrix_slack fixed diagnostic contract

2026-10-27. Main has explicitly released this finite diagnostic to Terra. Theory source: `hierarchy_proof.md`; this Sol owner executes no numerics. Numeric output belongs to Terra's separate `upper/round007/matrix_slack/`.

Fixed original64-direction process and same8 Hadamard blocks/pairing as R004-R006. 15 processes, each2^20 fresh rawGamma(1/2,1)^64 rows, P=Y/sumY; total15*2^20. No old optimizer, Sobol substitution, alternate model or adaptive extra samples. The numerical owner locks15 distinct original-generator seeds before the run; recommended new SeedSequence identifiers[202610277,7007,j],j0..14, if no earlier main-authorised seed scheme is already frozen. Record exact scheme in the actual numeric registration and never switch/reseed mid-run.

Per row compute all6 observables Z0,Z1,Z2,Z3,Zpairmat,ZpairSum from proof(1). Assert only the diamond(2); Zpairmat and Z1 are unordered. Record the8 nonnegative paired gaps listed in proof and the signed difference Zpairmat−Z1. Exact envelope B=D/64 covers all6levels/gaps; signed range is[−B,B]. q0/no additional R variable is involved nowhere here; D analytically conditions out R.

Streaming: fixed2048-row chunks,512 chunks per worker; preserve every row count and RNG/checkpoint handle. Summaries retain per-worker/chunk counts, compensated means/sumsquares/cross-gap statistics, extrema and flags. Full15-million-row shape cache need not be archived; retain code/seeds/checkpoint identities and all fallback exact input hex. Do not silently rerun completed terminal chunks.

Arithmetic: log-determinants and log-sum-exp,15 workers, BLAS/OpenMP1. Flag same-input fallback if globalYmax/Ymin>2^40, Cholesky/nonfinite failure, or a proved range/hierarchy violation exceeds tau_gate=(D/64)*2^(-32). Retain tiny signed floating violations and count them; never clip. A nonpositive/nonfinite original Y stops as a generation error, not a redraw. Fixed fallback recommendation:128 decimal digits on the exact original binary Y inputs, normalise there, recompute all6observables. Verify invariants to the documented fallback tolerance; if not resolved report and stop.

Maximum4096 fallback rows globally, unique(worker,row) keys; shared locked allocation prevents exceeding it. Required fallback4097 stops the batch, no row deletion or new seed. Unflagged floating arithmetic is not a rigorous integration proof. All statistical results remain diagnostic; no strict I/g/C7 update.

Statistical list fixed:6levels+8nonnegativegaps+1signedgap. Record block means, paired gap distributions, bounded-sampling confidence expressions with declared alpha (suggest1/10^6) and rangeB or2B as appropriate. State their ideal-iid/arithmetical scope; no claim that empirical confidence endpoints are deterministic expectation bounds.

Stop after exactly15 terminal worker streams/complete assembly, source/seed/hash checks and main review. Report actual CPU/wall/RSS/concurrency; requesting15 does not prove saturation. Any failure/cap exhaustion is explicit and does not permit automatically increasing samples/fallback cap or changing arithmetic gates.
