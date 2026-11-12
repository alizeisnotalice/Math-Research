from fractions import Fraction as F
from math import factorial, comb

def exp_lower(x, degree):
    return sum(x**j/factorial(j) for j in range(degree+1))

assert exp_lower(F(6932,10000), 8) > 2
assert F(6932,10000)+F(2,100)+F(24,1000)+F(25,1000)+F(51,100) <= F(12722,10000)
x = F(1,50)
neglog_upper = x + x*x/(2*(1-x))
assert F(4,100)*F(51,100)-neglog_upper > F(19,100000)
for n in range(1, 25):
    s = F(3*n,10000)
    exact = sum(comb(n,j)*(2*s)**j/factorial(j) for j in range(n+1))
    assert exact <= exp_lower(80*s,n)*(1+F(2,80))**n
print('Volume, gamma tail, and 24 exact dilation comparisons passed')
