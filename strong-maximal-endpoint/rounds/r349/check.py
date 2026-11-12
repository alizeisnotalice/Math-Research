from fractions import Fraction as F
from itertools import combinations
# R1=[0,3], R2=[2,5], E=[0,1] union [2,4].
# Assign receivers [0,3] to R1, [3,5] to R2. Source atoms have equal unit length.
a=(F(3,5),F(2,5))
rows=((F(1,3),F(1,3),0),(0,F(1,3),F(1,3)))
delta=F(1,2)
ratios=[]
for size in range(4):
    for A in combinations(range(3),size):
        left=sum(weight*min(1,sum(row[j] for j in range(3) if j not in A)/delta)
                 for weight,row in zip(a,rows))
        if size: ratios.append((1-left)/F(size,3))
        else: assert left==1
k=max(ratios)
assert k==F(11,10)
flow=((F(11,30),F(7,30),0),(0,F(2,15),F(4,15)))
assert all(sum(row)==a[i] and all(v<=a[i]*rows[i][j]/delta for j,v in enumerate(row))
           for i,row in enumerate(flow))
assert max(sum(row[j] for row in flow)*3 for j in range(3))==k
print("8 cuts checked; optimal density factor", k)
