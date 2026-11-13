from fractions import Fraction as F
from itertools import product

def intervals(x, prefix):
    return [(a,b) for a in range(x+1) for b in range(x,4)
            if a<b and (not prefix or x in (a,b))]

def values(g,x,y,prefix):
    columns = [[F(sum(g[3*t+i] for i in range(a,b)), b-a)
                for t in range(3)] for a,b in intervals(x,prefix)]
    maxima = list(map(max, zip(*columns)))
    def vertical(h):
        return max(sum(h[a:b])/(b-a) for a,b in intervals(y,prefix))
    return vertical(maxima), max(map(vertical,columns))

for g in product(range(2), repeat=9):
    for x,y in product((1,2),repeat=2):
        assert values(g,x,y,True) == values(g,x,y,False)
