from fractions import Fraction as F
from itertools import permutations,product

boxes = ((0,2),(1,3),(1,4),(2,6),(0,4))
source = {0,1,3,4,5}
shapes = {b-a for a,b in boxes}
def cost(l,S):
    old = {x for a,b in boxes if b-a in S for x in range(a,b)}
    val = F(0)
    for y in range(6):
        atoms = [source.intersection(range(a,b)) for a,b in boxes if b-a==l and a<=y<b]
        val += max((F(len(c-old),len(c)) for c in atoms),default=0)
    return val
G = {frozenset():F(0)}
for size in range(1,4):
    for path in permutations(shapes,size):
        S=frozenset(path)
        G[S]=max(G[S-{l}]+cost(l,S-{l}) for l in S)
direct=max(sum(cost(l,set(path[:j])) for j,l in enumerate(path)) for path in permutations(shapes))
assert G[frozenset(shapes)] == direct
starts = (0,1,2)
owners = ({0,1},{2},{3})
small = [sum(F(len(own),2) for a,own in zip(starts,owners) if a<=x<a+2) for x in range(4)]
assert small == [1,F(3,2),1,F(1,2)]
for n in range(1,5):
    total = F(0)
    for x in product(range(4),repeat=n):
        value = F(1)
        for u in x:
            value *= max(1,small[u])
        total += value
    assert total/4**n == F(9,8)**n
print('All six priority orders, the subset recurrence and four product dimensions verified')
