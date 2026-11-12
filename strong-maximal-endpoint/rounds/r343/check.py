from fractions import Fraction as F
from itertools import combinations
tags=list(map(set,combinations(range(8),4)))
for A in tags:
    for B in tags:
        d=len(A-B)
        ratios=[F(3,4) if (i in A)!=(i in B) else F(1) for i in range(8)]
        hit=F(1)
        for r in ratios: hit*=r
        assert hit==F(9,16)**d
        increment=sum((1 if r<1 else F(3,4))-F(3,4) for r in ratios)
        assert increment==F(d,2)
print(len(tags)**2,"exact overlap and marginal checks passed")
