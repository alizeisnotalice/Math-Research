from fractions import Fraction as F
from itertools import product

for k in range(1,9):
    for shift in (F(0),F(1,3),F(2)):
        ts=[F(i,12) for i in range(1,97)]
        values=[(1-shift/t)**(k-1) if t>shift else F(0) for t in ts]
        assert values==sorted(values)
weights=[F(1,6),F(1,3),F(1,2)]
for h in product([F(0),F(1,3),F(2,3),F(1)],repeat=3):
    total=sum(a*b for a,b in zip(h,weights))
    left=total
    used=[]
    for w in reversed(weights):
        used.append(min(left,w))
        left-=used[-1]
    used.reverse()
    assert sum(used)==total
    for slope in ((1,2,4),(0,0,1),(1,1,1)):
        assert sum(a*w*v for a,w,v in zip(h,weights,slope))<=sum(u*v for u,v in zip(used,slope))
print('24 likelihood grids and 192 exact mass-preserving comparisons passed')
