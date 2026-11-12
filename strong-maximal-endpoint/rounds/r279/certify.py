from fractions import Fraction as Q
for p in (Q(1,4),Q(1,2),Q(3,4)):
    for gap in (Q(1,10),Q(1,20),Q(1,50)):
        b=p*(1-gap)
        for q in (p,p*(1+gap/2),p*(1-gap/2)):
            D=(1-q)/(q-b)
            assert (1+b*D)/(1+D)==q
            shorter=D*Q(99,100)
            assert (1+b*shorter)/(1+shorter)>q
            for s in (2,5,10):
                assert 1+2*shorter*(b/p)**s < 1+2*shorter
                # Integrate (2/(epsilon*t))*s*t**(s-1)/q**s exactly.
                eps=(q-b)/q
                assert 2*s*q**(s-1)/((s-1)*eps*q**s)==2*s/((s-1)*(q-b))
print('27 two-rectangle averages and 81 layer charges verified')
