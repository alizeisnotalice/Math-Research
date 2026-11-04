#!/usr/bin/env python3
"""Exact rational checks for the cube example in angular_grouping.tex."""
from fractions import Fraction
from decimal import Decimal, getcontext

getcontext().prec = 50

def dec(q: Fraction) -> Decimal:
    return Decimal(q.numerator) / Decimal(q.denominator)

for n in range(1, 4):
    q = Fraction(n + 2, 3 * n)  # E_sigma theta_1^2
    for delta in (Fraction(1, 2), Fraction(3, 4)):
        individual_sq = 1 + delta * delta * q
        balanced_sq = Fraction(1, 1)
        # The normalized L2 Minkowski gap is 2(sqrt(individual_sq)-1).
        # Positivity is certified exactly by individual_sq > 1.
        assert individual_sq > balanced_sq
        for cplus, cminus in ((Fraction(1), Fraction(1)),
                              (Fraction(2), Fraction(1)),
                              (Fraction(3), Fraction(5))):
            deff = delta * (cplus - cminus) / (cplus + cminus)
            grouped_sq = 1 + deff * deff * q
            assert grouped_sq < individual_sq
            if cplus != cminus:
                assert grouped_sq > 1
        print(
            f"n={n}, delta={delta}: q={q}, "
            f"S_individual^2={individual_sq}≈{dec(individual_sq)}, "
            f"S_balanced^2=1"
        )

print("PASS: exact Fraction identities, strict balanced Minkowski gaps, and coefficient grouping checks")
