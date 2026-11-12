from fractions import Fraction
from math import comb
def half_tail(n, j):
    return Fraction(sum(comb(n,i) for i in range(max(j,0),n+1)),2**n)
for n in (8,16,32):
    k=n//2
    j=3*n//4
    hit=half_tail(k,j-k)
    assert hit>Fraction(1,2)
    print(n,"union",half_tail(n,k),"source",half_tail(n,j),"hit",hit)
