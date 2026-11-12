from fractions import Fraction as Q
for j in range(1,100):
    e=Q(j,1000)
    mass=8*(1-e+e/20)
    assert 8-mass==Q(38,5)*e
    assert Q(4,5)+Q(1,5)*(1-Q(33,10)*e)==1-Q(66,100)*e
    assert mass/(8-e)<=1-Q(825,1000)*e
    assert (1+Q(2,5)*e)*(1-e/3)>1
    assert (1+Q(55,100)*e)*(1-Q(66,100)*e)<=1-Q(11,100)*e
print('99 exact moving-threshold comparisons passed')
