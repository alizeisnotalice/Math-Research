# R009 fixed exact ray-catalog registration for Terra

2026-10-27. Source `fan_proof.md` and frozen literal `fan_normals.json`. Main explicitly authorises only the finite ray catalog after source/hash lock; no automatic chambers, triangulations, H/G evaluations, MC or new C7 claim. Sol theory owner runs no math script; Terra owns its separate numerical directory. No MATLAB/commit/push.

Input:36 integer w normals, row0..6=e1..e7; row7..35=29listed cut-prefix rows in the exact JSON order. Full V is(1,tail), tail lex in{-1,+1}^6 with−1 first. Verify regeneration by first-negative-hit templates;64directions,35nonnegativeprefix,29internalcuts, no global duplicate normals/sections. Coordinate w=xi_i−xi_(i+1), xi=Mw with detM1, lambda=sum_i i w_i. Do not enumerate using rounded u_i/k normals.

Finite jobs: every lexicographic6-subset of36 row indices, exactly1,947,792. Integer6x7 kernel/cofactors. Rank<6/allzero cofactors is not a unique ray; mixed-sign kernels rejected; otherwise orient nonnegative including zeros, gcd-reduce to primitivew. Storexi=Mw,lambda,normalised section rational point, all36incidences, all64body signs, one6row/rank6minor witness. Deduplicate only by primitivew. Walls/degenerate rays retained; no positivity-strict filter. Optional bodyzero-rank label must use all64body normals, not only29cuts.

Source/code/API probe before scientific run. Pilot exactly first4096 and last4096 lex combination IDs (8192 unique); compare optimised kernel to signed6minor reference on these, checkMkernel0,rank witness,orientation/gcd,7axis rays and coordinate transform. Pilot records are reused.

After pilot acceptance,15 single-thread processes compute remaining1,939,600 combinations, no wall cap and no adaptive additional normals/model. Optimised fraction-free kernel or wall-stratum dimensional reduction is allowed if identical to reference. Audit exact coverage via lex IDs/counts and the wall-stratum counts475020,831285,498771,127890,14210,609,7 for c0..6, respectively. Do not silently rerun terminal IDs.

Output minimal structural certificate: source normals/generation hashes; complete finite coverage/terminal counts;rank-deficient/mixedsign/acceptedhit counts; sorted unique primitive rays with incidences/signs/witnesses; merge/dedup receipt; versions/CPU/wall/RSS/actual concurrency. No need archive all rejected cofactor arrays. A checker validates every stored ray/witness and coverage manifests; no sampled-only ray validity assertion. Source hashes rechecked on completion.

Stop on any exact mismatch, coverage loss, inconsistent primitive ray or source mutation. A failed implementation may be repaired only with transparent new receipt; no midpoint/rank-tolerance substitute. Stop after the complete ray catalog/main review. This is a new arrangement batch, not a rerun of the earlier spectral-unit-ball catalog, and not a complete fan triangulation or H/g certificate.
