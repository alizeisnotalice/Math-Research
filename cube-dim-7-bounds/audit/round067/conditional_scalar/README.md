# R067 independent conditional scalar audit

The source hashes and finite gates were registered before running the independent stdlib Fraction checker. The checker reconstructs the fixed parameters, radius trim, error constants, P0, U bounds, row budget and M. It uses 40 terms of the atanh series with the rational remainder bound 2y^81/[81(1−y²)], independently of the producer's 48-term implementation. Because μ0=1/2−(5/2)log(6/5), its **upper** bound uses the lower logarithm endpoint.

The log-mass lower bound is approximately 0.318812435>637/2000, and gamma is approximately 0.318419204>317/1000. The μ0 upper bound is approximately 0.044196108<221/5000; the separate profile ceiling is approximately 0.044220042<223/5000, and the variance bound is approximately 0.006061426<1/150. Exact identities are η−χ=47/183400 and τ=47/200000, and 317/1000−83/917=207689/917000 exceeds 9/40 by 341/229250. Exact rational results, hashes and resource receipt are in `audit.json`.

This is an audit of one conditional scalar candidate, not a target-profile, geometry, or operator proof and not an optimization of the complete Phi. Existing R066 results have their own proof and audit; they are not independently re-established here.
