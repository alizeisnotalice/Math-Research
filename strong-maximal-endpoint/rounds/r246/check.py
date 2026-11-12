from itertools import combinations
from math import comb
from fractions import Fraction

# All 2-subsets of two groups of 3 hit each 3-subset of six axes.
blocks = [c for group in ((0, 1, 2), (3, 4, 5)) for c in combinations(group, 2)]
assert all(any(set(b) <= set(t) for b in blocks) for t in combinations(range(6), 3))
for n in range(100, 1501):
    k, d = n // 100, n // 4
    assert Fraction(comb(n, k), comb(d + 1, k)) <= Fraction(25, 6)**k
print('six-axis covering and 1401 exact binomial ratios passed')
