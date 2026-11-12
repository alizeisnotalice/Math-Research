from fractions import Fraction
from itertools import combinations
from math import comb,log
def H(x):
    return -x*log(x)-(1-x)*log(1-x)
for d in (1,2,3):
    for k in range(5):
        sets=list(combinations(range(5),k))
        actual=Fraction(sum(set(range(d))<=set(v) for v in sets),len(sets))
        assert actual==Fraction(comb(k,d),comb(5,d))
delta=1/16
old=delta*log(2)
new=.5*H(2*delta)-.25*H(4*delta)
packing=log(2)-H(2*delta)
assert old<new<packing
print("15 exact intersection tests passed",old,new,packing)
