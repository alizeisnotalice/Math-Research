from fractions import Fraction as F
from math import factorial

def exp_lower(x):
    return sum(x**j/factorial(j) for j in range(13))

assert exp_lower(F(811,1000)) > F(9,4)
assert exp_lower(F(6932,10000)) > 2
assert 2*sum(F(1,3)**(2*j+1)/(2*j+1) for j in range(5)) > F(6931,10000)
assert exp_lower(F(954,10000)) > F(11,10)
assert F(49,100)+F(1,10)+F(954,10000)+F(811,1000)-F(6931,30000) < F(1266,1000)
assert F(811,1000)+F(3,5)*F(6932,10000)+F(45,2000)+F(1,45) < F(1272,1000)
tests = 0
for scale in [F(j,12) for j in range(1,49)]:
    for radius in (F(1,2), F(1)):
        for power in (1,3,10,45):
            for delta in (F(0), radius/2, radius*scale, radius*scale*2):
                actual = F(delta <= radius*scale)/scale
                tail = max(F(1),delta/radius)**(-1-power)
                price = scale**power if scale>1 else 1/scale
                assert actual <= tail*price
                tests += 1
print(tests, 'component comparisons and both exponent bounds passed')
