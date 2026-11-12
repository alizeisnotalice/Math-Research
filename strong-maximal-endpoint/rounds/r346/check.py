from itertools import combinations_with_replacement,product
from fractions import Fraction as F
from math import comb
for n,m in [(3,1),(4,2),(6,2),(6,3)]:
    def value(a):
        sums=[0]
        for v in a: sums+= [x+v for x in sums]
        return F(sum(a),4*m)-F(sum(max(0,x-m) for x in sums),m*2**n)
    best=max(value(a) for a in combinations_with_replacement(range(m+1),n))
    assert best==F(1,2)-F(comb(2*m,m),2*4**m)
    print(n,m,best)
assert F(3,8)-sum(max(0,F(sum(bits),2)-1) for bits in product(range(2),repeat=3))/8==F(5,16)
