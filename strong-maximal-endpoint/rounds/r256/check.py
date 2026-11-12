from fractions import Fraction as Q
from itertools import product
from math import prod

a, lam = Q(1,10), Q(1,100)
mgf = (1-a)*2/(2-lam) + a*6/((2-lam)*(3-lam))
# log(mgf) <= mgf-1
assert Q(54,100)*lam - (mgf-1) > Q(35,10**6)
assert Q(6932,10000)+Q(1,100)+Q(140,10000)+Q(2,140)+Q(54,100) < Q(1272,1000)
u = [Q(0), Q(1,9), Q(2,3), Q(1), Q(7,6)]
by_mask = sum(prod(a*max(0,1-v) if bit else 1-a for bit,v in zip(mask,u))
              for mask in product((0,1), repeat=len(u)))
assert by_mask == prod(1-a*min(1,v) for v in u)
print('Both exponent margins and all 32 masks verified')
