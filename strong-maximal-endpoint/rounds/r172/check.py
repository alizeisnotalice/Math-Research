"""Exact volume and two binomial tails for the ordered-cover obstruction."""
from fractions import Fraction
from math import comb, log

for n in (12, 120, 1200):
    threshold = 7 * n // 12
    source = sum(comb(n, i) * 2**i for i in range(n // 2 + 1))
    outside = sum(comb(n, i) * 2**i for i in range(threshold))
    inside = sum(comb(n, i) for i in range(threshold, n + 1))
    bound = Fraction(outside, 3**n) + n * Fraction(inside, 2**n)
    if n == 1200:
        assert bound < Fraction(5122, 10**9)
    print(n, n * log(3) - log(source), float(min(1, bound)))
