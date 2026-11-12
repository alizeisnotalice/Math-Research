from fractions import Fraction as F
from itertools import product

boxes = [set(range(i, i+6)) for i in (0, 2, 4, 6)]
source = {0, 2, 3, 5, 6, 7, 9, 10, 11}
owners = [{x for x in box if all(x not in prev for prev in boxes[:j])}
          for j, box in enumerate(boxes)]
src = [r & source for r in boxes]
vol = [len(r) for r in owners]
v = [[sum((x >= min(box)+3) == b for x in own) for b in (0,1)]
     for box, own in zip(boxes, owners)]
for r, ar in enumerate(src):
    lam = sum(F(vol[s],vol[r])*F(len(ar & other),len(other))**2
              for s,other in enumerate(src) if s != r)
    for o in range(2):
        if not v[r][o]:
            continue
        threshold = F(v[r][o],len(ar)**2)
        candidates = [s for s in range(len(src)) if F(v[s][1-o],len(src[s])**2)>=threshold]
        union = set().union(*(set(product(ar & src[s],repeat=2)) for s in candidates))
        J = F(len(union),len(ar)**2)
        lost = 0
        for x,z in product(ar,repeat=2):
            scores = [max(F(v[s][b],len(a)**2) for s,a in enumerate(src) if x in a and z in a)
                      for b in range(2)]
            lost += scores.index(max(scores)) != o
        assert F(lost,len(ar)**2) <= J <= (1+lam)*F(vol[r],v[r][o])
print('All positive-volume box-direction states passed')
