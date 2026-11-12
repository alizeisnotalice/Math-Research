from fractions import Fraction as F
from itertools import combinations, product
intervals=[(a,b) for a,b in combinations(range(4),2)]
def cells(rs):
    return {i for a,b in rs for i in range(a,b)}
count=0
for mask in range(1,1<<len(intervals)):
    menu=[r for j,r in enumerate(intervals) if mask>>j & 1]
    for f in product((0,1,3),repeat=3):
        mean={r:F(sum(f[r[0]:r[1]]),r[1]-r[0]) for r in menu}
        g=[max((mean[r] for r in menu if r[0]<=i<r[1]),default=0) for i in range(3)]
        levels=sorted(set(mean.values())|set(f)|{0})
        charge=F(0)
        for lo,hi in zip(levels,levels[1:]):
            t=(lo+hi)/2
            cover=[r for r in menu if mean[r]>t]
            for r in cover.copy():
                rest=[j for j in cover if j!=r]
                if cells(rest)==cells(cover):
                    cover=rest
            mult=[sum(a<=i<b for a,b in cover) for i in range(3)]
            assert max(mult)<=2
            # For s=2 the layer integrand is affine between all listed breakpoints.
            charge+=2*(hi-lo)*sum(max(0,f[i]-t) for i in range(3) if mult[i]==2)
        assert 0<=charge<=sum(v*v for v in f)
        assert sum(v*v for v in g)<=2*sum(u*v for u,v in zip(f,g))+charge
        count+=1
print(count,'exact interval-menu and input pairs passed')
