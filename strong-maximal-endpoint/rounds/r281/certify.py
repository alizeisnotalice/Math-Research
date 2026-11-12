from fractions import Fraction as F
from math import comb, factorial

for b in (F(1,4),F(1),F(3)):
    poly=[F(1)]
    for j in range(2,17):
        new=[sum(c*comb(i+1,h)*b**(i+1-h)/(i+1)
                 for i,c in enumerate(poly) if i+1>=h) for h in range(j)]
        for x in (F(0),F(2,3),F(5)):
            assert sum(c*x**i for i,c in enumerate(new)) == (x+b)*(x+j*b)**(j-2)/factorial(j-1)
        poly=new
for s in (5,15,50):
    for l in (F(4),F(8)):
        for k in range(1,13):
            ratio=sum(F(comb(k-1,i)*factorial(i),1)/((s-1)*l)**i for i in range(k))
            a=(k-1)/((s-1)*l)
            assert 1<=ratio<=1/(1-a)
print('Shifted polynomial recurrences and 72 exact integration ratios passed')
