from fractions import Fraction as F
from math import comb, log
def tail(n,t,j):
    return sum(comb(n,i)*t**i*(1-t)**(n-i) for i in range(j,n+1))
f=F(3,8)
q=2*f/(1+f)
n,j=8,4
mass,hit=tail(n,f,j),tail(n,q,j)
assert hit>F(1,2) and q==F(6,11)
print("source volume",mass,"each box average",hit,"union volume",1)
print("sharp rate",2*log(2)-1-log(2*log(2)))
