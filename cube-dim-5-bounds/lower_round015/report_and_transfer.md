# Geometry 015: one rational tail-plane rotation

The single pre-registered angle cos(theta)=4/5, sin(theta)=3/5 completed
successfully. The model uses the exact algebraic C_cert spacings, not rational
approximations substituted as a different input. The exact finite-family
ratio is enclosed by

    [1.47332918344428616026493779757201420421195251226882888270847,
     1.47332918346304105892696890440200012879207403292133703539896]

(decimal displays; exact rational endpoints are in the result). Its upper
endpoint lies strictly below the complete C_cert lower isolating endpoint.
This is no new lower bound. It excludes only this angle and the explicitly
registered rotated consecutive-block family, not the rotated lattice's
complete superlevel set and not a neighborhood of zero angle.

The run actually processed 230400 configurations. It retained 779 inner and
779 outer boxes. Exact recursive slicing and coordinate-mask intersection
agree separately for both unions. Direct integer-index cardinality checks
and all ten geometric endpoint inequalities passed for every inner retained
box, for a total of 7790 checked inequalities. The run used 55554 nodes and
15.930261792 seconds against its 200000-node/120-second limit. Actual local
Python Fraction arithmetic, integer fifth roots and outward interval
propagation were used; no MATLAB process was used.

`supercell_theory.md` was saved and hash-locked before execution. It proves
the 25-point supercell, integer representative uniqueness, complete planar
translation cover, the nine neighboring translation sufficiency for the
registered block shapes, and the 8/25 normalization. The calculation does
not assume individual reflections of the rotated coordinates are symmetries.
The pre-execution theory and source files remain unchanged.

## Explicit finite L1 transfer for the strict inner witnesses

Use the period lengths (a,s,s,5s,5s) in the computational coordinate order.
The unrotated selected indices belong to [-1,2]; reflection across the middle
of an unrotated period sends j to 1-j and preserves this range. In the
rotated plane, each representative start is in [0,5s)^2 and neighboring-cell
translations add at most 5s in either direction. The physical displacement
of a block atom from its start belongs to [-9s/5,12s/5] in the first rotated
coordinate and [0,21s/5] in the second. Hence every selected planar atom
lies between -2 and 2 neighboring period indices. A uniform collar width
three in all five period directions is therefore sufficient.

For K^5 interior receiver periods include all 25 atoms from input period
indices in {-3,...,K+2}^5. The mass is exactly 25(K+6)^5, with no duplicate
atoms because representatives and half-open period translations are unique.
Every atom used by every reflected and translated inner witness is present.
For the finitely many retained radii r_lo, choose a common eta>0 such that
(2(r_lo+eta))^5<N for each selected mass N. Replace every input atom by a
smooth nonnegative unit-mass bump supported in a sup-norm ball of radius
less than eta/2. Each inner receiving box still admits a single centered
cube radius r_lo+eta in all coordinates containing all selected bumps, so
its average is strictly larger than one. The eight unrotated-axis reflected
copies have total receiving volume at least 8 K^5 V_inner. Thus ordinary
finite smooth L1 inputs attain ratios at least

    (8/25) V_inner * (K/(K+6))^5.

Let K tend to infinity. This validates the strict inner ratio for the general
L1 centered cube operator with all positive scales; it creates no upper
bound for that operator and does not promote this failing candidate.

## Why the next question concerns witness shape

After rotation, a physical axis-aligned cube cuts the original integer
(m,n) plane in a tilted parallelogram. The captured lattice points need not
form a consecutive rectangular (m,n) block. The low result above can therefore
reflect an overly restrictive witness family. No general one-sided angular
derivative or local rotation obstruction has been proved, and no blind
angle grid is launched.

The next finite geometric question is to fix the unrotated three-axis length
product K, choose a planar atom count n, use the common radius
r=(K*n)^(1/5)/2, and compute the planar centers lying in at least n of the
atom-window squares [p_x-r,p_x+r] x [p_y-r,p_y+r]. Multiplying that exact
planar overlap region by the three unrotated receiver intervals gives valid
nonrectangular-index witnesses. This describes an unresolved next mechanism;
it has not been registered or executed in this round.

Reproduce from the repository root:

    python3 cube-dim-5-bounds/lower_round015/rotated_supercell.py

All rational endpoints, hashes, source limits, reference hashes and completion
checks are in `rotation_4_3_5_results.json` and `completion_receipt.json`.
No commit or formal ledger edit was made. The existing bound does not count
as a newly discovered breakthrough in this round.
