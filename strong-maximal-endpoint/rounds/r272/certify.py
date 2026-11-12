from fractions import Fraction as F
from itertools import product
from math import comb, factorial

cells = list(product(range(3),repeat=2))
ivs = ((0,2),(1,2),(1,3))
boxes = [(tuple(b-a for a,b in r), {i for i,x in enumerate(cells) if
          all(a<=z<b for z,(a,b) in zip(x,r))}) for r in product(ivs,repeat=2)]
choices = [sorted({s for s,h in boxes if i in h}) for i in range(9)]
best=0
for lab in product(*choices):
    score=F(0)
    for s in set(lab):
        chosen={i for i,v in enumerate(lab) if v==s}
        for y in range(9):
            score+=max((F(len(chosen & h),len(h)) for sh,h in boxes if sh==s and y in h),default=0)
    best=max(best,score/9)
assert best==F(19,18)<F(3,2)
for k in range(1,41):
    for j in range(k+1):
        assert F(comb(k,j),2**k)*F(2,3)**k*2**j==F(comb(k,j)*2**j,3**k)
# log(9/(4 sqrt(2))) < .465
assert sum(F(930,1000)**j/factorial(j) for j in range(9))>F(81,32)
print('64 joint assignments, 40 binomial laws and the exponent bound passed')
