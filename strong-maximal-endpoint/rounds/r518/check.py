from random import Random

rng = Random(518)
grid = {(x, y) for x in range(12) for y in range(12)}
for trial in range(128):
    source = {point for point in sorted(grid) if rng.randrange(2)}
    rows = []
    for _ in range(60):
        x, y = rng.randrange(12), rng.randrange(12)
        w, h = rng.randint(1, 12-x), rng.randint(1, 12-y)
        box = {(i, j) for i in range(x, x+w) for j in range(y, y+h)}
        if 2*len(box & source) > len(box):
            rows.append((w, box))
    selected = set()
    union = set()
    for width, box in sorted(rows, key=lambda row: -row[0]):
        union.update(box)
        if 4*len(box & selected) <= len(box):
            assert 4*len((box-selected) & source) > len(box)
            selected.update(box)
    assert len(union) <= 24*len(selected)
