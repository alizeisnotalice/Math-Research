from fractions import Fraction as Q
from math import factorial
lam, u, b = Q(1,2000), Q(6,1000), Q(19,25)
assert lam*u - 4/b*lam**2/(1-lam) > Q(16, 10**7)
def exp_lower(x, n):
    return sum(x**i / factorial(i) for i in range(n+1))
assert exp_lower(Q(6932,10000), 8) > 2
assert Q(6932,10000) + 2*u + Q(572,1000) < Q(12784,10000)
assert exp_lower(Q(8), 14) > 2000
assert Q(1+8,2000) < u  # H(t) <= t(1+log(1/t))
print('Tail, volume exponent and Hamming bound verified')
