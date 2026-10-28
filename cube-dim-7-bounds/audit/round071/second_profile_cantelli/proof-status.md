# R071 independent preflight — full theorem PENDING

Stage note: the preparation status below is preserved as historical scope.
It is superseded by the full-chain ACCEPT in review.md, locked to author
source 0e36488654bfa1413c19335a3eb9cf7b900e16b25cf9b73b2b8f36efc7381a10.

This is an adversarial preparation note, not acceptance of a new bound for C_n.
The author has communicated a tuple but has not yet supplied a frozen full proof.
No earlier R068 manuscript or accepted source was changed.

## Interfaces checked independently

Write L=log n >=8192, x=.91L, h=1/(256L), trim=2^-19,
rlo=1+trim, rhi=1.365-trim. The tuple uses eN=eB=eT=2^-19,
eJ=eRad=eS=eD=2^-20. At the base trim=4h. Thus trim>2h on the
whole tail. The complete matching interval [(r-h)/s,(r+h)/s] stays
inside the source interval, its baseline q_b=(r+h)/s is available,
and an active matched source has R=sq-h>rlo-2h>1. These facts would
fail if the R066 trim 2^-30 were retained with the new h.

The source-count endpoint calibration has drift at most 2h+1/n.
The new eN gives 2h=eN/2 at the base, leaving room below 3eN/4.
The old eN=2^-20 would not suffice. The exact negative derivative of
ell_pin still has the correct sign for the endpoint-ratio comparison;
this does not depend on h being of order 2^-80/L.

The rigorous relative-coordinate bound |partial_q log J|<3/rho,
not the old coarse 256/rho bound, is needed. Across matching width
2h/s, summing n coordinates gives at most 6xh <=6Lh=3/128.
Consequently the proposed P0=1-eB-3/128 pays the product loss via
exp(-u)>=1-u. The old nearly-unit P0 is unavailable. Individual J
requires only the pointwise bound and an n lower bound; it must not
be replaced by the n-fold product estimate.

Radial cancellation is exact. The remaining width factor is
atanh(h/r)/(h/r), between 1 and 1/(1-(h/r)^2), and the new h still
fits eRad. The finite-a center guard is algebraically unchanged,
provided the eventual manuscript verifies t<n, beta<2, the actual
source radius >4 and target item radius <3 in the new domain.

## Profile and probability status

The geometry author confirms that Cantelli uses only the per-source-face
variance upper bound sum_(i != j) J2_i/J_i <= sum_i J2_i/(1-eJ).
Nonnegativity permits adding the omitted source pin back in this upper
bound. No unproved variance lower bound or subtraction of a target-law
mean is needed. The source pin in the mean is paid separately by .183.
Cantelli must be applied within each fixed source-face probability law
and only the resulting tail probabilities may then be mixed.

The roots auditor has communicated a local H2 result (formal report
pending): exact J2=R/(3rho)(Qg^3+Qv^3), a.e./AC derivative <.1/rho,
and matching transport <=xh/5. The coarser 2xh budget already fits
10^-6 x at the proposed base. Even the old 16xh first-profile transport
fits 10^-5 x here. These checks do not establish target-good probability:
the frozen proof must show six bad-event budgets, with the two-sided J2
condition counted as the event whose two tails are jointly paid.

An independent Fraction run, PID 99742, checked twelve compressed
critical margins in preflight.py, all PASS. It used n>=L^32/32!, not the
author's 64-term Taylor envelope. It checks the active trim, endpoint
drift, product loss, both profile bridges, pointwise J, radial width,
and Cantelli directly as v/(v+x gap^2)<14/25, as well as the conservative
gap 1/1330. All adverse terms decrease for larger L: h decreases,
1/(L^32/32!) decreases, x increases, and xh is constant. The conditional
energy-budget rational inequality is also positive. This arithmetic
does not certify the profile laws or the operator transfer.

## Still required before full ACCEPT

1. Lock and read the complete author source, including the exact H2 local
   dependency and its actual domain; no acceptance from a parameter receipt.
2. Reconstruct both target profile expectations, pin handling, target-face
   independence and six-event probability union for every L>=8192.
3. Check source-face J denominators, pin mean, H2 upper bound and Cantelli
   separately before mixture, then clipping and original corner-cap loss.
4. Prove the source-cap and tau exponential budgets for the entire tail,
   including Csrc*rhi/h, finite-n exp correction and Taylor-prefix monotonicity.
5. Reconstruct the actual same-history weighted row bound, finite-a ca^2/ca
   energy transfer, general L1 witnesses and all positive common radii.
6. Verify the complete F(A,B-gamma n) scalar cost, including positive A/cross
   buffers, rather than compare only gamma with a displayed log price.

Current disposition: PREPARATION ONLY / full C_n theorem PENDING.
