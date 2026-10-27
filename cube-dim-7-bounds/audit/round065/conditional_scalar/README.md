# R065 independent conditional scalar audit

The fixed finite gates and producer source hashes were registered before the check. `audit.py` uses independent stdlib Fraction arithmetic and 40 terms of the positive atanh series with the rigorous tail 2y^81/[81(1−y²)], without importing the producer helper. The stated radius endpoints, trim, P0, U bounds, M, row budget and all scalar margins were reconstructed. The 40-term log interval contains the saved 48-term producer endpoint, while all rational fields agree exactly.

The certified mass lower bound is approximately 0.30975213>3097/10000, and gamma is approximately 0.30638875>61/200. Also 18233/100000>log(6/5), μ<999/20000<χ, and v_rate<1/150. Exact calculations give η−χ=1/18180, σ/2−χ(1−σ)=1/20000, 61/200−91/909=37249/181800, and subtracting 1/5 gives 889/181800. Exact fractions, hashes and runtime are in `audit.json`.

The audit concerns these conditional scalar values only; geometry and the operator theorem are outside its scope. Earlier artifacts involving the incorrect sign in the tau expression are not used as evidence of a mathematical exclusion. Different valid conservative constants in R064 need not be equal to these fixed values.
