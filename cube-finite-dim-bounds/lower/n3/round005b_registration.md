# Round005 stage b, after independent checkpoint audit PASS

2026-10-25. Parent and independent auditor explicitly confirmed first round005 checkpoint PASS against full prior T3 upper; original checkpoint input/result remain immutable. This authorizes the additional finite stage below. Use distinct round005b_* files and hashes.

Fixed start is round005_certificate_input.json, target its full exact upper endpoint from round005_certificate_result.json. Same primary+small-satellite construction, all three periods, all three phases and epsilon free. Preserve every original 170 key and update all box radii/masses/geometries at each evaluation; full joint-union volume/(1+epsilon) only.

At most 1000 Nelder-Mead evaluations, split into up to two 500-evaluation chunks so selected best points can be independently recaptured and their added keys included before the second chunk. Initial recapture: 32768 centers, seed 20261025511. Between chunks/final recapture: 65536 centers, seeds 20261025512 and 20261025513. Entry radii<=3.2, offsets[-4,4]; no all-scale capture exhaustion claim. Explicit simplex steps in log periods (.0003,.00005,.00005), phases (.005,.005,.005), log epsilon .2, halved for the second chunk.

Stop optimizing immediately on full-union gain >1e-6 above this checkpoint's complete exact upper endpoint; then recapture and certify. Smaller gains are retained and certified at the batch end. Exact comparison must beat the fixed checkpoint's full exact upper, not its shorter display value. End after this certification; no third stage without further authorization. Actual Python+MATLAB; explicit b/finite mass/strict-L1 transfer. No old evidence changes, commits/pushes/ledger edits.
