from fractions import Fraction as Q
from math import factorial

def length(a, b, c, d):
    return max(Q(0), min(b,d)-max(a,c))

cases = 0
for m in range(1,9):
    h = Q(1,m)
    eta = Q(1,1000)
    assert h*(1-eta)+(2-2*h)*(1-eta)/2+h*(eta/h) == 1
    for j in range(20):
        a = h*j/20
        for k in range(m+1):
            cell = k*h
            same = length(a,a+1,cell,cell+h)
            if not same:
                continue
            far = sum(length(a,a+1,(k+v)*h,(k+v+1)*h) for v in (-m,m))
            assert same + far == h
            cases += 1
def exp_lower(x, n):
    return sum(x**j/factorial(j) for j in range(n+1))
assert exp_lower(Q(6932,10000),8) > 2
assert exp_lower(Q(512,100),20) > Q(999,6)
assert Q(11,6)*Q(6932,10000)+Q(1,999)+Q(512,100*800) < Q(12783,10000)
print(cases, 'receiver-cell cases and both logarithm certificates passed')
