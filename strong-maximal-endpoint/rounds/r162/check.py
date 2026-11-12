"""Exhaust all subsets of a binary tree of seven vertices with unit active masses."""
from itertools import product

def ancestors(i):
    path = {i}
    while i:
        i = (i-1)//2
        path.add(i)
    return path

paths = [ancestors(i) for i in range(7)]
for weights in product(range(3), repeat=7):
    masses = list(weights)
    for i in range(6, 0, -1):
        masses[(i-1)//2] += masses[i]
    end = 0
    while 2*end+1 < 7:
        end = max((2*end+1, 2*end+2), key=masses.__getitem__)
    loss = sum(weights[i] for i in range(7) if i not in paths[end])
    conflict = sum(weights[i]*weights[j] for i in range(7) for j in range(7)
                   if i not in paths[j] and j not in paths[i])
    assert loss**2 <= 2*conflict
print(3**7, 'integer tree profiles passed')
