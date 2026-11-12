from fractions import Fraction as F
from math import comb
from itertools import product
for m in (2,4,8,16,32,64):
    loss=sum(F(comb(2*m,i)*(i-m),4**m*m) for i in range(m+1,2*m+1))
    assert loss==F(comb(2*m,m),2*4**m)
# Three unit cells; intervals consisting of cells 0,1 and 1,2.
# At p=3/4 the common-source infimum is 2 (entire middle, half of each end).
for m in (1,2,4,8):
    val=max((F(3,4)*2*(a+b)-max(0,a+b-m))/m
            for a,b in product(range(m+1),repeat=2))
    assert val==2
    print(m,val)
