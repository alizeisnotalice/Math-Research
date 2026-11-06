# r054 final independent review

## Result

Accepted, with one scope-language cleanup recommended for the interpolated
chain.  The scalar and finite-cone computations were checked separately; no
M512 receiver computation was rerun.

## Scalar staircase

The staged `code/scalar54_cert.py` was actually run against
`rounds/r054/scalar54.json`, writing only `/tmp/r054_scalar54_recheck.json`.
It completed with the same strict directed margins:

\[
 E H\in[1.0003040150239913507\ldots,\;1.0003040150239913508\ldots],
\]
\[
 J_{L=10^8}\le6.5837639397339196\ldots\times10^{-32}
 <e^{-71.79}=6.6374176241009982\ldots\times10^{-32}.
\]

The exact-Decimal audit is sound: `O` is explicitly `Decimal(1)`, so
`gpoly` initializes both `t` and `s` with the `j=0` term and then adds
`j=1,...,63`; it is the full Gamma(64) survival polynomial.  All finite
input decimals are protected by the 500-digit `Inexact` trap, while only
`exp`/`ln` are widened one ulp outward.  The finite transfer uses
`alpha=1-63/10^8` and
\((1+t/L)^{63}\le e^{63t/L}\) on the positive staircase support, exactly
matching the computed `Jalpha`.

This certifies the explicitly stated finite scalar construction; it is not
an optimality or general endpoint assertion.

## Cone and polyline score

`code/cone_score54.py` was run against the immutable existing D1 replay
record, without rerunning M512.  Startup/final hashes agreed.  It reproduced
the finite-anchor score upper bound `71.432003... < 71.5` and the straight
old-anchor chain score `71.668105... < 71.79`.  The directed interval formula
for `F(z)=(1-e^{-z})/z` is sign-safe for the observed range crossing zero,
and the Pareto-volume prefix-minimum formula is correct.  The script checks
the complete anchor list before using the diagonal/chain ordering facts.

Recommended wording change only: the `linear_pareto_chain.scope` JSON text
should say the certified witness is **exactly the particular interpolated
chain containing the old anchors**.  A score computed from that chain does
not automatically upper-bound a larger halo or a receiver family with extra
points, even if it retains the old constraints.  The updated prose proof
already makes this distinction; align the output text with it.

## Limits

Nothing here proves a full halo result, an arbitrary receiver-family result,
or a statement about `C_N`.  The D1 record remains a separate frozen,
directed lower certificate; the cone calculation correctly validates its
input/source hashes.
