from fractions import Fraction as F
from math import factorial

#  log(4) - phi(4/5)/6 < 1.264 after multiplying by 30 and exponentiating.
ratio = F(4)**30 / (F(9,5)**9 * F(1,5))
x = 30*F(1264,1000)
assert sum(x**j/factorial(j) for j in range(90)) > ratio
for m in range(1,13):
    h = F(1,m)
    for j in range(101):
        t = F(j,100)
        cuts = sorted({F(0),F(1),t,max(0,t-(1-h)),min(1,t+(1-h))})
        avg = sum((v-u)*max(0,abs((v+u)/2-t)-1+h)/h for u,v in zip(cuts,cuts[1:]))
        b = min(t,1-t)
        assert avg == (max(0,h-b)**2+max(0,h-1+b)**2)/(2*h)
print('Strict rate and 1212 phase averages verified')
