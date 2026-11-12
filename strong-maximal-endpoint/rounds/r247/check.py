from fractions import Fraction as Q
from itertools import product
from math import prod

records = list(product(range(3), repeat=2))  # short, outer long, inner long
axes = range(2)
for q in [(Q(2, 3), Q(1, 4)), (Q(1, 3), Q(1, 5))]:
    h = [1 - v for v in q]
    for code in range(1, 1 << len(records)):
        family = [w for j, w in enumerate(records) if (code >> j) & 1]
        shells = {}
        for w in family:
            outer = tuple(i for i in axes if w[i] == 1)
            weight = prod(q[i] for i in axes if w[i] == 2)
            shells[outer] = max(shells.get(outer, 0), weight)
        envelope = sum(prod(h[i] for i in outer) * v for outer, v in shells.items())
        mean = Q(0)
        for bits in product(range(2), repeat=2):
            t = [i for i in axes if bits[i]]
            projections = {tuple(w[i] for i in t) for w in family
                           if all(w[i] != 2 for i in t)}
            mean += len(projections) * prod(h[i] if bits[i] else q[i] for i in axes)
        assert envelope <= mean
print('1022 exact ternary-family comparisons passed')
