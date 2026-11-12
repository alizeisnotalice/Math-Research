from fractions import Fraction as F
from math import factorial, comb, prod
from itertools import product

for k in range(1,10):
    for s in (3,8,15):
        a=F(s-1)
        L=F(7,3)
        vals=[sum(comb(j,i)*L**(j-i)*factorial(i)/a**(i+1) for i in range(j+1))
              for j in range(k)]
        for j in range(1,k):
            assert a*vals[j]==L**j+j*vals[j-1]
for k in range(1,7):
    b,q=F(7,16),F(1,2)
    H=(1-b)/(q-b)
    for side in product((F(1),F(3,2),F(2)),repeat=k):
        P=prod(side)
        assert (b+(1-b)/P>q)==(P<H)
print('Layer integral recurrences and six dimensions of strict hull averages passed')
