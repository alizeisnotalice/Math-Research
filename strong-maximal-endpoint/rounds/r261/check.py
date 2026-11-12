from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import factorial

points = [F(j,8) for j in range(5)]
for depths in combinations_with_replacement(points,3):
    for j in range(8):
        t = F(2*j+1,16)
        cdf = sum(x<t for x in depths)
        assert sum(x-1<t<x for x in depths) == 3-cdf
        assert sum(x-1<t-1<x for x in depths) == cdf
    for x in depths:
        assert min(x,0)-max(x-1,-F(1,2)) == F(1,2)
assert sum(F(6932,10000)**j/factorial(j) for j in range(9)) > 2
assert F(12784,10000)-F(3,2)*F(6932,10000) > F(2385,10000)
print('35 depth profiles, signed bins, half-segments and exponent margin verified')
