# r040: hard localized right clip audit

Date: 2026-11-06.  This records the valid reduction supplied by the old
inner rate; it is not a proof of the remaining intermediate one-step bound.

Let `alpha=1+eps` and let

`J_old=I_mu-(1+eps)D`, `Delta=J_old-I_*=F-eps D`.

For sufficiently small fixed `delta>0` with `2 delta<A3/K4`, the expansions
at zero give a uniform strict buffer

`Delta(eta)<=-c eps eta^3<0`, `delta eps<=eta<=2delta eps`.        (1)

The lower endpoint is a fixed positive number once `eps` is fixed.  Thus
the buffer is bounded away from zero on this band by a constant depending
on `eps`, and saddlepoint/Fourier prefactors (`O(log t)`) are eventually
dominated by `ct`.  Rate continuity in `a_N->mu` is enough uniformly up to
`t=N`; one must not require `t|a_N-mu|->0`.  The finite-call perturbation
is `O(t/N)`, hence bounded, and is also harmless against the linear-in-t
buffer.

Define the localized candidate as follows:

* retain `min(O^-,O^+,I_old)` for `eta<=eta_b=2delta eps`;
* set it equal to `O=O^- wedge O^+` for `eta>eta_b`.

By (1), throughout the buffer `[delta eps,eta_b]` the old inner is larger
than `O^+` for all `t>=t0(eps)` and large matched N.  Therefore the minimum
already equals `O` in a full neighbourhood of the artificial interface;
the forced removal causes no jump or positive clipping residual.

At a current state with `eta>=delta eps`, the localized value equals `O`,
so positivity directly gives `T W_prev<=T O_prev<=O=W`; this resolves all
fixed-positive-tilt re-entry.  If a current core state has `eta<delta eps`,
then for every accessible predecessor `r-z` (`z>=0`) its saddle coordinate
is at most `delta eps+O(1/t)<eta_b` once t is large.  Hence `W_prev<=I_old,prev`
on every such predecessor, so any independently proved local inequality
`T I_old<=I_old` transfers without a crossing term.

**Critical scope:** the last sentence applies only on whatever subset the
inner one-step inequality is already known (the r036 natural tube).  It
does not prove that inequality for the remaining intermediate strip
`A sqrt(log t/t)<eta<delta eps`.  The clipping lemma removes the large
positive-tail/re-entry issue, but the uniform one-step estimate in that
strip and the left clip remain open.
