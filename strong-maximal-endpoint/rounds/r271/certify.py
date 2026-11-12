from fractions import Fraction as F
from itertools import product
from math import comb

atoms = list(product(range(3),repeat=2))
intervals = ((0,2),(1,2),(1,3))
labels = {(0,1):(2,1),(1,0):(1,2)}
label = lambda x: labels.get(x,(2,2))
boxes = [(tuple(b-a for a,b in r),{x for x in atoms if all(a<=u<b for u,(a,b) in zip(x,r))})
         for r in product(intervals,repeat=2)]
coef = [F(0) for _ in range(3)]
for y in atoms:
    for s in {label(x) for x in atoms}:
        rows = [row for shape,row in boxes if shape==s and y in row]
        if not rows:
            continue
        def score(row):
            hits = [x for x in row if label(x)==s]
            return F(len(hits),len(row)), -F(sum(sum(u!=1 for u in x) for x in hits),len(row))
        row = max(rows,key=score)
        for x in row:
            if label(x)==s:
                coef[sum(u!=1 for u in x)]+=F(1,len(row))
assert coef==[F(9,4),5,F(9,4)]
for k in (12,24,48,96):
    t=7*k//12
    proportions=[F(sum(comb(i,j) for j in range(min(i,t)+1)),2**i) for i in range(k+1)]
    assert min(proportions)==proportions[-1]>F(1,2)
    poly=[F(1)]
    for _ in range(k//2):
        out=[F(0)]*(len(poly)+2)
        for i,v in enumerate(poly):
            for j,c in enumerate(coef):
                out[i+j]+=c*v
        poly=out
    assert sum(poly[:t+1])/3**k>=F(1,2)*F(19,18)**(k//2)
print('Exact block polynomial, four threshold checks and four retained gains passed')
