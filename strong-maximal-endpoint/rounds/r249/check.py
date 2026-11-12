from fractions import Fraction as F
from math import comb, log

for j in range(101):
    a = F(j, 100)
    tails = sum(t - t*t/2 for t in (a, 1-a))
    assert 1 + tails == F(3, 2) + a*(1-a)
# Simpson integrates the quadratic anchor term exactly.
assert (F(3, 2) + 4*F(7, 4) + F(3, 2))/6 == F(5, 3)
for m in range(2, 19):
    # Total gap s has density m(m-1)s(1-s)^(m-2); conditional error is s/4.
    mean = sum(F((-1)**j * comb(m-2,j)*m*(m-1), 4*(3+j)) for j in range(m-1))
    assert mean == F(1, 2*(m+1))
print('101 anchors, exact quadratic integral and 17 midrange means passed')
