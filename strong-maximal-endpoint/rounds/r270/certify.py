from fractions import Fraction as F
from itertools import product

intervals = [(0,2),(1,2),(1,3)]
for n, optimum in [(1,F(3)),(2,F(19,2))]:
    points = list(product(range(3),repeat=n))
    boxes = []
    for rect in product(intervals, repeat=n):
        shape = tuple(b-a for a,b in rect)
        atoms = {i for i,x in enumerate(points) if all(a<=v<b for (a,b),v in zip(rect,x))}
        boxes.append((shape,atoms))
    labels = [sorted({s for s,atoms in boxes if i in atoms}) for i in range(len(points))]
    best = 0
    cases = 0
    for lab in product(*labels):
        reward = F(0)
        for s in set(lab):
            chosen = {i for i,v in enumerate(lab) if s==v}
            for y in range(len(points)):
                reward += max((F(len(atoms & chosen),len(atoms))
                               for shape,atoms in boxes if shape==s and y in atoms),default=0)
        if n==2 and lab==tuple({(0,1):(2,1),(1,0):(1,2)}.get(x,(2,2)) for x in points):
            assert reward==optimum
        best=max(best,reward)
        cases+=1
    assert best==optimum
    print(n,cases,best/len(points))
short = [{0,1},{1,2},{2,3}]
owners = [{0},{1,2},{3}]
for x in range(4):
    envelope = F(len(set().union(*(r for r in short if x in r))),2)
    score = sum(F(len(o),2) for r,o in zip(short,owners) if x in r)
    assert max(1,score)==max(1,envelope)
assert sum(max(1,F(len(set().union(*(r for r in short if x in r))),2))
           for x in range(4))/4==F(5,4)
print('Four-interval envelope attained on every source atom')
