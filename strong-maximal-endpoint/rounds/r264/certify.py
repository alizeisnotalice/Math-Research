from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial

intervals = list(combinations(range(6), 2))
count = 0
for (a,b),(c,d) in product(intervals, repeat=2):
    overlap = max(0, min(b,d)-max(a,c))
    for y in (Q(c),Q(c+d,2),Q(d)):
        dist = max(0,a-y,y-b)
        for t in (1,2,6):
            assert Q(overlap,d-c)**(t+1) <= (1+Q(dist,b-a))**(-1-t)
            count += 1
# For t=3 the symmetric weight h^4 and both asymmetric weights are rational.
weights = [Q(1,7),Q(4,7),Q(2,7)]
lengths = (2,5,8)
sym = asym = 0
for i,j in product(range(3),repeat=2):
    v,w = lengths[i],lengths[j]
    h = min(v,w)
    sym += weights[i]*weights[j]*Q(h*h,v*w)**2
    asym += weights[i]*weights[j]*Q(h,v)*Q(h,w)**3
assert sym <= asym
def exp_low(x):
    return sum(x**k/factorial(k) for k in range(14))
assert exp_low(Q(6932,10000)) > 2
assert Q(12784,10000)-Q(6932,10000)-Q(57,100) > Q(15,1000)
assert exp_low(Q(981,1000)) > Q(8,3)
print(count,'interval comparisons, posterior symmetry and exponent bounds verified')
