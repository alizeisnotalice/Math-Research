from fractions import Fraction as Q
from itertools import product
from math import factorial, prod

h = [Q(1, 3), Q(1, 2), Q(3, 4)]
masks = range(8)
weights = {f: prod(q if f & (1 << i) else 1-q for i, q in enumerate(h)) for f in masks}
for code_mask in range(1, 256):
    code = [l for l in masks if code_mask & (1 << l)]
    shell = sum(max(prod((h[i] if b & (1 << i) else 1-h[i])
                         for i in range(3) if l & (1 << i)) if b & l == b else 0
                    for l in code) for b in masks)
    proj = sum(weights[f] * len({l & f for l in code}) for f in masks)
    assert shell <= proj
    d = int(len({l & 4 for l in code}) == 2)
    assert proj <= Q(3, 2)**3 * 2**d
# log(3) + log(2)/4 < 1.2724 by taking the fourth power.
x = 4 * Q(12724, 10000)
assert sum(x**i / factorial(i) for i in range(20)) > 3**4 * 2
print('255 codes and the rational exponent comparison passed')
