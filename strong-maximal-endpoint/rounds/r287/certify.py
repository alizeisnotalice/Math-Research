from fractions import Fraction as F
from itertools import combinations
intervals=list(combinations(range(6),2))
for mask in range(1,32):
    A={j for j in range(5) if mask>>j&1}
    h=[F(i+1,9) for i in range(5)]
    q=max(sum(h[a:b])/(b-a) for a,b in intervals if A.intersection(range(a,b)))
    for s in (F(1,5),F(1,3),F(1,2),F(3,4)):
        cover=set().union(*(set(range(a,b)) for a,b in intervals if len(A.intersection(range(a,b)))>s*(b-a)))
        assert sum(h[i] for i in cover)<=2*q*len(A)/s
d=F(8)
for a,b in intervals:
    for c,e in intervals:
        length=e-c
        mid=F(c+e,2)
        if max(a,mid-d*length/2)>min(b,mid+d*length/2):continue
        if length<=b-a:
            assert a-d*(b-a)<=c<=e<=b+d*(b-a)
        else:
            assert mid-(d+2)*length/2<=a<=b<=mid+(d+2)*length/2
print('124 weighted covers and interval dilation comparisons passed')
