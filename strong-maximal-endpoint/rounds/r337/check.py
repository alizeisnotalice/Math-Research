from fractions import Fraction
from math import comb, log
D=.75*log(3)-log(2)
for n in (32,64,128,256,512):
    k=n//2
    j=3*n//4
    # Symmetric Bin(k, 1/2) upper half, including half of the central atom.
    mass=Fraction(2*sum(comb(n,i) for i in range(j+1,n+1))+comb(n,j),2**(n+1))
    hit=Fraction(2*sum(comb(k,i) for i in range(k//2+1,k+1))+comb(k,k//2),2**(k+1))
    assert hit==Fraction(1,2)
    print(n, "boundary hit",hit,"log residual",round(-log(mass)-n*D-.5*log(n),6))
