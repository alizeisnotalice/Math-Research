from fractions import Fraction
from itertools import permutations, product
from collections import Counter

n = 4
bits = list(product(range(2), repeat=n))
cases = 0
for x in bits:
    if sum(x) % 2:
        continue
    for order in permutations(range(n)):
        counts = Counter()
        for y in bits:
            hybrid = list(y)
            pattern = []
            for i in order[:-1]:
                hybrid[i] = x[i]
                pattern.append(sum(hybrid) % 2)
            counts[tuple(pattern)] += 1
        assert len(counts) == 2 ** (n - 1) and set(counts.values()) == {2}
        cases += 1
for arm, density in [(Fraction(1, 4), Fraction(1, 5)), (Fraction(1, 3), Fraction(1, 4))]:
    assert arm / (1 + arm) == density > Fraction(1, 8)
assert Fraction(1, 4) + Fraction(1, 3) < 1
print(f'{cases} source-order pairs passed')
