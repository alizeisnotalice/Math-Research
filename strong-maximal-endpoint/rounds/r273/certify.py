from fractions import Fraction as F
from itertools import product
J = ({0,1},{1},{1,2})
cols = [sum(max((F(1,len(S)) for S in J if y in S and z in S),default=0)
            for y in range(3)) for z in range(3)]
assert cols == [1,2,1]
def eq(b):
    return 8*b*b+3*b-1
assert eq(F(2126,10000))<0<eq(F(2128,10000))
assert 1+2*F(2128,10000)<F(3,2)
for mask in range(1,2**9):
    E = {(x,z) for x,z in product(range(3),repeat=2) if mask & (1<<(3*x+z))}
    total = 0
    for y in range(3):
        g = [max(F(sum((x,z) in E for z in S),len(S)) for S in J if y in S)
             for x in range(3)]
        total += sum(g)
        # On each unit atom use its left segment of length g[x].
        for R in J:
            for S in J:
                if y in S:
                    assert sum(g[x] for x in R)>=F(len(E & set(product(R,S))),len(S))
    assert total<=2*len(E)
print('Kernel columns, algebraic root bracket and 511 arbitrary source masks verified')
