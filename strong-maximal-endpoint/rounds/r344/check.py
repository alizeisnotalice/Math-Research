from itertools import product,combinations
from fractions import Fraction as F
for h,r in [(2,3),(3,2),(4,2)]:
    n=h*r
    code=[v for v in product((0,1),repeat=n) if all(sum(v[j:j+h])%2==0 for j in range(0,n,h))]
    f=F(len(code),2**n)
    assert f==F(1,2**r)
    for I in combinations(range(n),h-1):
        counts={}
        for v in code:
            z=tuple(v[i] for i in I)
            counts[z]=counts.get(z,0)+1
        assert len(counts)==2**len(I) and len(set(counts.values()))==1
    for eta in (f/3,f,2*f):
        kept=max(F(0),1-eta/f)
        assert sum((kept/F(2**n) for _ in code),F(0))==max(F(0),f-eta)
    print(n,len(code),f,"passed")
