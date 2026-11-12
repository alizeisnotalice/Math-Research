from fractions import Fraction as F
from itertools import product

cells=list(product(range(3),repeat=2))
rects=[({0,1},{0,1}),({1,2},{1,2}),({2},{0,1,2}),({0,1,2},{2})]
def step(f,axis):
    g={}
    for x in cells:
        vals=[]
        for r in rects:
            if all(x[i] in r[i] for i in range(2)):
                vals.append(sum(f[(v,x[1]) if axis==0 else (x[0],v)] for v in r[axis])/len(r[axis]))
        g[x]=max(vals,default=F(0))
    return g
t=F(1,4)
for mask in range(1,512):
    f={x:F((mask>>j)&1) for j,x in enumerate(cells)}
    mass=sum(f.values())
    for i in range(2):
        prev={x:F(v>1-t**i) for x,v in f.items()}
        f=step(f,i)
        aux=step(prev,i)
        level=[x for x in cells if f[x]>1-t**(i+1)]
        assert all(aux[x]>1-t for x in level)
        assert len(level)<=F(5,3)**(i+1)*mass
# s=4, P and Q proportional to fourth powers; clear the normalization factors.
for a,b in product(product(range(1,5),repeat=2),repeat=2):
    lhs=sum(x*y**3 for x,y in zip(a,b))
    mid=sum(x*x*y*y for x,y in zip(a,b))
    assert lhs*lhs<=mid*sum(y**4 for y in b)
print('511 indicator sources and 256 exact Holder comparisons passed')
