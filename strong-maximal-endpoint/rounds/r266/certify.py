from fractions import Fraction as Q
from math import prod, factorial
from itertools import product

count = 0
for m in range(2,18):
    N = 4*m
    for k in range(N+m):
        owners = [j for j in range(N+1) if j<=k<j+m]
        cells = {t for j in owners for t in range(j,j+(1 if j<N else m))}
        assert len(cells) <= 2*m
        if m<=k<N:
            assert len(cells)==m
        count += 1
for shift in product((Q(0),Q(1,4),Q(1,2),Q(1),Q(3,2)), repeat=3):
    assert 1-prod(max(0,1-v) for v in shift)<=sum(shift)
x = Q(1271,1000)
assert sum(x**k/factorial(k) for k in range(12))>Q(356,100)
print(count, 'reverse-volume cells, 125 translation bounds and exponent verified')
