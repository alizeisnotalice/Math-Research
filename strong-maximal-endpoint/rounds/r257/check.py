from fractions import Fraction as F
from math import factorial

for m in range(1, 33):
    starts = [F(j,m) for j in range(-2*m, 2*m+1)]
    for j in range(m):
        x = F(2*j+1,2*m)
        cover = [s for s in starts if s < x < s+1]
        left, right = min(cover), max(cover)
        assert (right,left+1) == (F(j,m),F(j+1,m))
        core = left+1-right
        union = right+1-left
        assert (core+union)/2 == 1
x = F(6932,10000)
assert sum(x**i/factorial(i) for i in range(9)) > 2
assert F(184,100)*x < F(12784,10000)
assert F(11,6)*x < F(12784,10000)
print('528 grid cells and both exponent comparisons verified')
