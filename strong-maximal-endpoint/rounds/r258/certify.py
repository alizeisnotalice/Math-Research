from fractions import Fraction as F
from math import factorial

eps = F(1,400)
ln2lo = 2*sum(F(1,3)**(2*j+1)/(2*j+1) for j in range(6))
x = (F(1,6)-2*eps)*ln2lo
expneg_upper = sum((-x)**j/factorial(j) for j in range(5))
zupper = 1-(1-12*eps)*(1-expneg_upper)
assert 4*zupper < sum(F(1278,1000)**j/factorial(j) for j in range(12))
for m in range(1,7):
    h = F(1,m)
    starts = [(j,j*h+shift) for j in range(-8,9) for shift in (-eps,0,eps)]
    endpoints = {e for _,s in starts for e in (s,s+1)}
    for j,s in starts:
        for k in range(j,j+m):
            lo, hi = k*h+eps,(k+1)*h-eps
            assert s <= lo < hi <= s+1
            assert all(not lo < e < hi for e in endpoints)
print('Strict exponent and all six perturbed grids verified')
