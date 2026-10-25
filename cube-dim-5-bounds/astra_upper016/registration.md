# Astra016: two ridge-weighted Cauchy Young necessary screens

Parent authorization 2026-10-26 after015. Only this directory may change. Old015 backend and evidence are read-only. No commits, no other agent directory changes.

Model alpha=1, beta=1/2, fixed u=(W+2)/12 in[0,1]. Exactly two new shapes Z=a*u^m with m in {2,3}; no m1 rerun, no eta grid, no other m.

Stage1: necessary corner x=(1,1,1,1,1) only. New moments E[u^(2m)], E[d u^(2m)], E[B u^(3m)] are produced under beta=1/2. Cauchy–Schwarz bounds U=E[sqrt(d)u^(2m)] by sqrt(AE), so all-a family ceiling is c1(AE)^(3/2)/V². Compare with014's current I threshold using outward Arb balls and rational endpoint comparisons. This bounds certificate families, not Green. Compute the degree<=13 moment arrays in a single finite batch, each completed W power checkpointed. Reuse only the matching beta=1/2 cache/header from015, SHA bfdf502cac067d2a1d9465c96ac15076c2e897cee89f8121bfed40977fc3fbf0. No degree guard change.

If both shapes are excluded, stop and lock. If one or both remain unexcluded, select at most one m: the one with the larger certified necessary ceiling lower endpoint; ties choose smaller m. This selection is only a finite pilot priority, not an assertion of achievability.

Stage2 for only that m: use actual new-beta E[d^k u^(2m)], k1..4, highest degree<=26. A two-node Gauss pair for d*u^(2m) dR selects one rationalized pair; Hermite-interpolate d^(-1/2), then multiply by d for a lower bound on sqrt(d). Direct cubic Hermite of sqrt(d) has the wrong direction. Round each coefficient downward with verified intervals, choose one rational a near U_lower/V. Check this actual corner candidate before any spatial production. Negative polynomial coefficients are permitted and carried exactly.

If warranted afterStage2, a subsequent recorded checkpoint may build the selected full spatial certificate, with either the one fixed a or adaptive a(x)=U(x)/V(x), requiring U>0,V>0 and cL U³-gV²>=0 everywhere. Neither a corner value nor floating minimum is a certificate. No other parameter/shape search is included.

Each process has120-second checkpoint, new cache states persist; a partial spatial/corner accumulation is never double counted. Existing exact source and graph fingerprints are recorded. Key format remains ten5-bit edge counts, total degree<=31; all moment arithmetic is arbitrary-precision GMP. General alpha1 L1/all-scale transfer and S1/S2 issues remain separate proof obligations for any C5 claim.
