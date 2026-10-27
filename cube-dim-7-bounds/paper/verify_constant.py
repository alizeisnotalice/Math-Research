from fractions import Fraction as Q
from math import comb
import json

t = Q(3391749, 524288)
U = Q(4910992105193, 200000000000)
N = 128
x = (t - 1) / (t + 1)
s = 2 * sum((x**(2*j+1) / (2*j+1)
             for j in range(N)), Q())
r = 2*x**(2*N+1) / ((2*N+1)*(1-x*x))
z = (s+r)*10**18
ell = Q((z.numerator+z.denominator-1)//z.denominator,
        10**18)
q = 1 + 1/(2*t)
h = (1+1/t)/2
beta = 7 - Q(7, 2)*q**6
terms = sum((Q(comb(6, j), j)*q**(6-j)
             * Q(-1, 2)**j*(1-t**(-j))
             for j in range(1, 7)), Q())
gamma = 1 - h**7 - Q(7, 2)*terms
B_low = gamma + beta*s
B_up = gamma + beta*ell
A = (1+1/t)**7
D = U - A - B_up
assert beta > 0 and B_low > 0
assert D > 0 and D*D > 4*A*B_up
print(json.dumps({
    'status': 'PASS exact rational certificate',
    't': str(t), 'U': str(U), 'N': N,
    'log_upper': str(ell),
    'B_lower': str(B_low), 'B_upper': str(B_up),
    'D': str(D), 'D_squared_minus_4AB': str(D*D-4*A*B_up)
}, indent=2))
