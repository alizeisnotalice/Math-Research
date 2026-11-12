from fractions import Fraction as F
from math import factorial

for lam in (F(1),F(3,2),F(2),F(3)):
    a = 2*lam/(1+lam)**2
    heights = (a,lam*a,a,a/lam)
    assert sum(heights)/2 == 1
    assert all(h <= lam*j for h,j in zip(heights,reversed(heights)))
    for j in range(101):
        b = F(j,200)
        assert min(0,b)-max(-F(1,2),b-1) == F(1,2)
# 2*H(1/3) = (1/3)*log(729/16).
def exp_lower(x,n):
    return sum(x**i/factorial(i) for i in range(n+1))
assert exp_lower(3*F(12731,10000),18) > F(729,16)
assert exp_lower(F(7,10),6) > 2
assert F(12731,10000)+F(7,10000)+F(4,1000) < F(12784,10000)
print('Four direction densities, 404 interval tests and exponent margins verified')
