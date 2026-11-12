from fractions import Fraction as F
from itertools import combinations
from math import comb
n,a=9,4
tags=[set(t) for i,t in enumerate(combinations(range(n),a)) if i%4==0]
for j in (4,5,6,7,9):
    measured=0
    for size in range(j+1):
        for t in combinations(range(n),size):
            w=sum(label<=set(t) for label in tags)
            measured+=w*(w-1)
    pairs=sum(sum(comb(n-b,k) for k in range(j-b+1))
              for u in tags for v in tags if u!=v for b in [len(u|v)])
    assert measured==pairs
    print(j,F(pairs,2**n))
