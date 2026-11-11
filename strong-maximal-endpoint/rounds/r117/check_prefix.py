"""Exact arithmetic for the tagged strict-cover example."""
from fractions import Fraction as F
from itertools import combinations

p, eta, delta, q = F(1, 2), F(1, 4), F(1, 4), 16
lengths = [q**j for j in range(1, 9)]
source = sum((1 - eta) * size * F(1, size) * F(4, 5) for size in lengths)
covered = sum(size * F(1, size) for size in lengths)
density = (1 - eta) * F(4, 5)
pair_bound = max((1 - eta) * a / (eta * b - a) for a, b in combinations(lengths, 2))
assert density - p == F(1, 10)
assert covered / source == F(5, 3)
assert pair_bound == F(1, 4) < p / (1 + delta)
print('source, union, strict margin, cross-label upper:', source, covered, density-p, pair_bound)
