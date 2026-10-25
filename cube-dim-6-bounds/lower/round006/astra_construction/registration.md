# D Astra: density-preserving dimerization, finite first batch

This is the handoff from the failed fixed-spacing parity and staggered
screens, not a restarted project. Only this directory is owned by D Astra.
No MATLAB, other agents, commits, or pushes are invoked here.

Let `a` be the exact R6-003 spacing vector. Choose one active axis `i`.
For `0<=e<1` put `s=(1+e)^(1/6)`, `b_i=a_i*s/(1+e)`, and
`b_j=a_j*s` for `j!=i`. In the active coordinate place the consecutive
points `p(k)=b_i*(k+e*(k mod 2))`; assign the full product atom the mass
`1+d*(-1)^k`, `0<=d<1`. Other axes use `b_j Z`.

The two alternating gaps are `b_i*(1+e)` and `b_i*(1-e)`.
The period cell has mass 2 and volume `2*prod(a)`, independent of e,d.
This is a genuine spacing/mass joint family: even-length block spans change
by `+/-e*b_i` while odd-length spans do not; odd-length block masses change
by `+/-d` while even-length masses do not. It is not any previous
fixed-spacing parity or cross-coordinate staggered family.

At e=0,d=0 it is exactly the baseline. At e tending to 1, the nearby
opposite-parity atoms merge to mass 2, and all effective lattice spacings
tend to `2^(1/6)*a`; global dilation/mass normalization recovers the baseline
ratio. Thus the chosen density-preserving path compares nontrivial
two-scale interiors against two equivalent endpoints, rather than losing
only through an arbitrary density drift. This is motivation, not a proof
that an interior improves the ratio.

For active start q and length l, the selected active mass is
`l+d*(-1)^q` when l is odd and l when l is even. Multiply by the other
five lengths, set the one common radius `r=M^(1/6)/2`, and use center
intervals `[source_max-r,source_min+r]`. Intersect with
`H=[0,2b_i] x prod_{j!=i}[0,b_j/2]`. For inactive axes take the same
central starts as the baseline. Every retained box is a valid common-radius
block witness. Reflect only the five inactive coordinates; these 32 copies
fill the full period cell, so the ratio is `32|U_H|/2=16|U_H|`.

The first registered batch has 24 cases: active axes 1,2,3 (3 represents
the four equal tails); e=1/8,1/3,2/3,7/8; d=0,1/5. The active lengths are
1..8, inactive lengths 1..4. Active length 8 is needed to include four
merged pairs near the second endpoint. Enumerate every active start whose
center interval meets H using a radius-based finite outer bound, then
filter by actual intersection. A separate e=d=0 calibration is required;
it must agree with the complete current certified endpoint to 2e-9.

Hard limits: 24 new candidates, 5000 raw boxes per candidate, 1,000,000
recursive union states globally, and 1200 seconds globally. Stop at the
first ratio more than 1e-9 above the current complete exact upper endpoint,
and export all raw boxes for E; no blind follow-up scan. Stop on mismatch
or any resource cap. Fraction arithmetic parses and retains the complete
ledger and frozen upper endpoints; float values are only discovery triggers.
No current C6 claim results from this screen.

`dimer_screen.py` hashes the protocol, thresholds, ledger, old D reports,
current source, and this registration before calculation. Its checkpoint
records PID, elapsed time, candidate counts and node counts after each case.
The actual execution session is recorded in the final report. Local Python
is the existing `cube-finite-dim-scratch/venv/bin/python`.
