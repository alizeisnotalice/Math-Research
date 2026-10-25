# Round 007: local common-tail event intervals

Set `a1=1+u^3/2` and `a2=a3=a4=a5=(1+u)/2+delta`, retaining all
`l in {1,2,3,4}^5`.  This round enumerated all 175 distinct *raw* affine
endpoint alternatives before min/max choices or dominance removal.  Their
pairwise crossings in `delta in [-1/16,1/16]`, together with the window and
`delta=0`, gave 369 exact event values.  The audit selected exactly the four
nearest event intervals on each side of zero.

On an open selected event interval all raw endpoint orderings, min/max choices
and grid-cell coverage masks are fixed.  The first-coordinate lengths are
constant in `delta`; the other four are affine.  Thus the exact ratio is a
polynomial of degree at most four.  This structural conclusion follows from
the complete raw-endpoint event list, **not** from interpolation.  Six exact
dual-union evaluations at interior rational nodes per interval were used to
construct and redundantly validate the polynomial implementation.

The comparison target is the full 25-coefficient algebraic C certificate in
`lower_round003/a2_equals_a3_results.json`, not a displayed decimal and not
the upper endpoint of its isolating interval.  On each of the eight selected
intervals, the Bernstein coefficients of `P(delta)-C` are all nonpositive
without subdivision.  Consequently `P(delta) <= C` on every one of those
closed event intervals, so none can strictly improve C.  The four left
intervals have a nonnegative derivative Bernstein certificate and the four
right intervals have a nonpositive one, giving the expected one-sided local
peak behavior at C.

This is an accurately delimited local exclusion, not a global common-tail
optimization or a new weak `(1,1)` claim.  No finite-collar or smoothing
transfer was needed because there is no new winning finite model.

## Status and audit correction

An initial selector accidentally joined the `-1/16` window endpoint to the
four nearest left events, crossing unselected events.  Its sixth interpolation
check failed and no exclusion was made; those six samples are preserved as
diagnostic evidence in `event_interval_0_sample_*.json`.  The corrected
selector uses exactly four intervals per side.  A preliminary corrected
result for interval zero compared to C's isolating *upper endpoint* and is
withdrawn; the `exactC_event_interval_*_results.json` files below are the
only valid finalizations.

## Reproduction and hashes

Run six commands `python3 common_tail_event_intervals.py --sample I K` for
`K=0,...,5`, then `python3 common_tail_event_intervals.py --finalize I`, for
each `I=0,...,7`.  Each sample is an independently exact dual-union job with
the registered 120-second and 200,000-node limits.

* registration: `1b8550a2ba88bfeb79d95b3f7b879301c533e289c4d241e6628c0bc0673157d7`
* event/interpolation/Bernstein source: `7b9b688982af07cfc784c5ab8d16795697e84be4382913a17faac19ce17ac383`
* exact C coefficient source: `f385811e5460af63abc117804965461eed88d8ecb5b9d2cd30268387646f45d2`
* finalized interval results, indices 0--7 respectively:
  `8228d160a1c0b54c5469952124a7afbaef1d204581f9bd3a90680ca4f0f8a901`,
  `0c811894d788e611410df45f76e5ba8ad694b0d9d0256023aa861e6676832ce2`,
  `6e4f8d82ae5c732b75cee059aead1d3dc053c6b24dace82f870f920085bdb192`,
  `708643ef54f5808e039fedd20081debc64da72e5094ca3fc36e210631c3a4aa8`,
  `beec6b02d36c4877c3fcf7d162485112eb40c1fb8aa20fc0d0e507a5bc9e48b5`,
  `17f16aabc9676225760965e5d859d818df30163e3dc685fd8f6c6e5c2e904c82`,
  `c9a8f22208d910598e58f4b07dcae081ec5c286b70601bea459ddc53d91d48d4`,
  `0cabb7c5f94b1bd618c41bfa7557dfd668c70beef66f15e7473f0a6336fcddc8`.

## Next finite candidate

This batch excludes only the nearest eight common-tail order intervals.  The
next distinct lower-bound mechanism should vary two tail spacing groups (for
example `a2=a3` and `a4=a5`) and first screen for a new non-dominated witness
box before attempting a two-parameter event-cell certificate.  It should not
repeat either the excluded one-parameter intervals or the failed alternating
mass model.
