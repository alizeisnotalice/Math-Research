from fractions import Fraction as F
from itertools import product

checked=0
for masses in product((0,1,2), repeat=6):
    if sum(masses)<8:continue
    j0=next(i for i in range(6) if sum(masses[:i+1])>=2)
    assert sum(masses[j0+1:])>=4
    checked+=1
for seq in product(range(3),repeat=6):
    for gap in range(1,6):
        assert 2*sum(x*y for x,y in zip(seq,seq[gap:]))<=2*sum(x*x for x in seq)
for j in range(-6,7):
    scale=F(2)**j
    for ratio in (F(1),F(3,2),F(19,10)):
        for D in (F(6),F(10),F(20)):
            assert (D+1)*scale*ratio/2 <= (D+1)*scale
print(checked,'scale splits and all finite correlation comparisons passed')
