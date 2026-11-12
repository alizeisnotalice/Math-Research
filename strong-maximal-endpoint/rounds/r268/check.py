from fractions import Fraction as F
from itertools import product

intervals = ((0,4),(1,5),(1,7),(3,6),(4,8))
source = {0,1,3,4,5,7}
receiver = range(8)
p = F(1,3)
sizes = sorted({b-a for a,b in intervals})
def candidates(z):
    return {b-a for a,b in intervals if a<=z<b}
for weights in product((1,2,5), repeat=len(sizes)):
    w = dict(zip(sizes,weights))
    def mass(z):
        return sum(w[l] for l in candidates(z))
    values = {x:F(0) for x in source}
    for y in receiver:
        for l in candidates(y):
            a,b = next((a,b) for a,b in intervals if b-a==l and a<=y<b)
            xs = source.intersection(range(a,b))
            assert len(xs)>p*l
            for x in xs:
                values[x] += F(w[l],len(receiver)*mass(y)*len(xs))*F(mass(y),mass(x))
    assert all(v <= F(2,len(receiver))/p for v in values.values())
    if len(set(weights))==1:
        assert sum(values.values()) <= F(2*len(source),len(receiver))/p
print('27 positive weight assignments checked on every source atom')
