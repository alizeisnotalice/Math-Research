from fractions import Fraction as Q
from math import comb
for s in range(2,16):
    for p in (Q(1,4),Q(1,2),Q(3,4)):
        # The halo of the interval has length 2/q-1 for 0<q<1.
        via_layers = 2*Q(s,s-1)/p - 1
        assert via_layers == 2/p-1+2/(p*(s-1))
        for j in range(15):
            f=Q(j,10)
            h=min(f,p)
            primitive=s*f*h**(s-1)-(s-1)*h**s
            assert primitive==h**s+s*p**(s-1)*max(0,f-p)
print('42 interval prices and 630 truncated layer primitives passed')
