from fractions import Fraction as F
from itertools import product
from math import factorial

segments = [(0,2),(0,1),(2,4),(2,3),(0,4)]
for f in product((0,1,3),repeat=4):
    mx=[max(F(sum(f[a:b]),b-a) for a,b in segments if a<=i<b) for i in range(4)]
    for t in (F(1,2),F(1),F(3,2),F(2)):
        points=[i for i in range(4) if mx[i]>t]
        assert len(points)*t <= sum(f[i] for i in points)
boxes=[({0,1},{0}),({1,2},{1}),({0},{0,1})]
for mask in range(512):
    f={(i,j):F(mask>>(3*i+j)&1) for i,j in product(range(3),repeat=2)}
    g={(i,j):max((sum(f[a,j] for a in I)/len(I) for I,J in boxes if i in I and j in J),default=0)
       for i,j in f}
    for I,J in boxes:
        avg=sum(f[i,j] for i in I for j in J)/(len(I)*len(J))
        for i,j in product(I,J):
            assert sum(g[i,b] for b in J)/len(J)>=avg
assert sum(F(6932,10000)**j/factorial(j) for j in range(9))>2
assert 1+F(2,5)*F(6932,10000)<F(12784,10000)
print('81 laminar inputs, 512 variable-fiber masks and the exponent margin passed')
