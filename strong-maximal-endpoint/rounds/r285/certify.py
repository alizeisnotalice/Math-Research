from fractions import Fraction as Q
from itertools import combinations_with_replacement
amounts=[Q(0),Q(1,3),Q(2,3),Q(1)]
sides=list(combinations_with_replacement(amounts,3))
cases=0
for left in sides:
    for right in sides:
        for a in range(4):
            for b in range(4):
                soft=sum(left[3-a:])+sum(right[3-b:])
                hard=min(a,sum(left))+min(b,sum(right))
                assert soft<=hard
                cases+=1
# For rates 1,2, g(t) is proportional to v(1-v), v=exp(-t).
for j in range(2,10):
    d=Q(j)  # exp(D)
    ratios=[]
    for m in range(2,80):
        v=Q(1,m)
        ratios.append(d*(1-d*v)/(1-v) if d*v<1 else 0)
    assert ratios==sorted(ratios)
print(cases,'exact interval comparisons and eight shifted two-rate laws passed')
