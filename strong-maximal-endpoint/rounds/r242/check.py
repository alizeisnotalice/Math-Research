from fractions import Fraction
from math import comb, isqrt

for n in (4, 64, 256, 1024):
    k = isqrt(n) // 2
    h = Fraction(1, 2 * k)
    assert (1 - h)**k >= Fraction(1, 2)
    assert comb(n, k) * h**k >= 2**k
    print(n, k, 2**k)

b = Fraction(1, 10000)
assert b + b / (1 - b) < Fraction(201, 1000000)
