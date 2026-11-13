from fractions import Fraction as F
from math import comb

def knots(n):
    mass = F(0)
    cost = 0
    points = [(mass, cost)]
    for j in range(n, -1, -1):
        mass += F(comb(n, j), 2**n)
        cost += comb(n, j) * 2**(n-j)
        points.append((mass, cost))
    return points

def value(points, p):
    for (x, a), (y, b) in zip(points, points[1:]):
        if p <= y:
            return a + (b-a)*(p-x)/(y-x)

for n in range(1, 17):
    low = knots(n-1)
    for p in (F(i, 16) for i in range(17)):
        left, right = max(0, 2*p-1), min(1, 2*p)
        candidates = {F(left), F(right)}
        for x, _ in low:
            candidates.update(v for v in (x, 2*p-x) if left <= v <= right)
        assert value(knots(n), p) == min(
            2*value(low, a) + value(low, 2*p-a) for a in candidates
        )
