from fractions import Fraction as F
from math import factorial

values = [F(i, j) for i in range(1, 8) for j in range(1, 5)]
for length in values:
    for prototype in values:
        u = length  # the kernel is smallest at the largest allowed displacement
        kernel = 1 / (prototype * max(F(1), u / prototype)**4)
        factor = max(length / prototype, prototype / length)**3
        assert factor * kernel >= 1 / length
x = F(981, 1000)
assert sum(x**i / factorial(i) for i in range(9)) > F(8, 3)
assert x + F(6, 25) + F(1, 25) + F(1, 100) + F(7, 1000) < F(12784, 10000)
print('784 scalar comparisons and both rational exponent bounds passed')
