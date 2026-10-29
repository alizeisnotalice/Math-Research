#!/usr/bin/env python3
"""Exact finite checks of the truncation factor; the general proof is in the note."""
from fractions import Fraction

count = 0
for dimension in range(1, 21):
    for collar in (1, 2, 5):
        for size in (2 * collar + 1, 10 * collar, 100 * collar):
            loss = 1 - (1 - Fraction(2 * collar, size)) ** dimension
            bound = Fraction(2 * dimension * collar, size)
            assert 0 <= loss <= bound
            count += 1
print(f"PASS: {count} exact rational checks of (T2); general proof uses Bernoulli.")
