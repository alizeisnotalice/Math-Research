from fractions import Fraction as F
from random import Random

rng = Random(533)
for trial in range(100):
    rows = [[rng.randrange(2) for _ in range(16)] for _ in range(5)]
    ends = [i for i in range(1, 17) if any(
        row[i-1] and (i == 16 or not row[i]) for row in rows)]
    if not ends:
        continue
    values = [[F(sum(row[:i]), i) for i in range(1, 17)] for row in rows]
    loss = [sum(max(v)-v[i] for v in values)/5 for i in range(16)]
    D = min(loss)
    assert D == min(loss[i-1] for i in ends)
    gaps = [[max(v)-v[i-1] for v in values] for i in ends]
    levels = sorted({F(0), F(1)}.union(*map(set, gaps)))
    def bad(h):
        return min(F(sum(x > h for x in col), 5) for col in gaps)
    assert sum((b-a)*bad(a) for a,b in zip(levels,levels[1:])) <= D
    assert D <= min(h+(1-h)*bad(h) for h in levels)
