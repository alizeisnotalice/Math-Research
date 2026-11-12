from fractions import Fraction as Q
from itertools import product
from math import comb

for j in range(101):
    t = Q(j, 100)
    outside = max(0, t-Q(1, 2)) + max(0, Q(1, 2)-t)
    assert outside == Q(1, 2)-min(t, 1-t)
for n in range(1, 7):
    cells = list(product(range(4), repeat=n))
    for k in range(n+1):
        assert sum(sum(a == 0 for a in cell) == k for cell in cells) == comb(n,k)*3**(n-k)
        assert sum(sum(a < 2 for a in cell) == k for cell in cells) == comb(n,k)*2**n
print('101 endpoint positions and six exact product-cell counts passed')
