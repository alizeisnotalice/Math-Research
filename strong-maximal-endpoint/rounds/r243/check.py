from fractions import Fraction
from itertools import combinations, product

intervals = list(combinations(range(-4, 5), 2))
count = 0
for (a, b), (c, d) in product(intervals, repeat=2):
    h = min(b, d) - max(a, c)
    if h <= 0:
        continue
    t = max(Fraction(1), Fraction(a + b - 2*c, b - a),
            Fraction(2*d - a - b, b - a))
    assert t <= 1 + Fraction(2*(d - c - h), b - a)
    assert t <= Fraction(d - c, h)**2
    count += 1
print(f'{count} ordered interval pairs passed')
