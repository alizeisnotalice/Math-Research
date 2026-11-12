from fractions import Fraction as Q
from itertools import product
from math import factorial

n, radius = 4, 1
options = range(2**n)
sets = ((0, 15), (3, 7))
near = [{x for x in options if any((x ^ a).bit_count() <= radius for a in s)} for s in sets]
actual = sum(any(x in s for x,s in zip(pair,near)) for pair in product(options, repeat=2))
assert Q(actual,2**(2*n)) == 1-(1-Q(len(near[0]),2**n))*(1-Q(len(near[1]),2**n))
# 6*(ln2+2H(1/3)-H(1/6)) = ln(3^6*5^5/2^8).
def exp_lower(x):
    return sum(x**k/factorial(k) for k in range(40))
assert exp_lower(6*(Q(12784,10000)+Q(24,100))) > Q(3**6*5**5,2**8)
assert exp_lower(3*Q(12784,10000)) > Q(729,16)
print('All 256 template pairs and both strict exponent gaps verified')
