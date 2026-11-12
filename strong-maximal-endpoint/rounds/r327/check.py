from fractions import Fraction
from itertools import product
for n,k in [(2,3),(3,4),(4,5)]:
    hits=[label for label in product(range(n), repeat=k) if len(set(label))==1]
    assert len(hits)==n and len(set(hits))==n
    assert Fraction(len(hits), n**k)==Fraction(1,n**(k-1))
    print(n,k,"union",len(hits),"incidence multiplicity",1,"success fraction",Fraction(len(hits),n**k))
