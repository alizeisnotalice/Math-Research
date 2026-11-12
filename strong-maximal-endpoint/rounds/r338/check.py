from fractions import Fraction
from math import comb, log
rate=sum(float(q)*log(float(3*q)) for q in (Fraction(1,2),Fraction(1,3),Fraction(1,6)))
for n in (24,48,96,192,384):
    count=sum(comb(n,a)*comb(n-a,c) for a in range(n//2,n+1)
              for c in range(min(n-a,n//6)+1))
    union=Fraction(count,3**n)
    assert 0<union<1
    print(n, "log union + nD + log n =", round(log(union)+n*rate+log(n),6))
