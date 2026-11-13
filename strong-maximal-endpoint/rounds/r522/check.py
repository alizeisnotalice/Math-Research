from fractions import Fraction as F

for n in range(3, 25):
    rows = []
    for i in range(1, n+1):
        cells = []
        for j in range(1, n+1):
            size = 2 + int(i == j-1) - int(i == j)
            cells.extend([0]*(4-size) + [1]*size)
        mass = 0
        row = []
        for end, occupied in enumerate(cells, 1):
            mass += occupied
            row.append(F(mass, end))
        assert max(row) == F(1, 2)
        rows.append(row)
    for j in range(n):
        assert all(row[4*j+3] == F(1, 2) for i, row in enumerate(rows) if i != j)
    best = max(sum(column)/n for column in zip(*rows))
    assert F(1, 2)-best == F(1, 4*n*n)
