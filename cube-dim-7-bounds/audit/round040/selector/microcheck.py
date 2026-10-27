#!/usr/bin/env python3
"""Independent, new exact case for the frozen selector proof; no halo scan."""

from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json
import platform

SOURCE = Path('/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round040/selector/PROOF.md')
EXPECTED = '19d057f510e0f8b529f1090014cfc5a56f85efb9cd72fe6bdc28bfffe6b11b40'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED


def volume(box):
    out = F(1)
    for lo, hi in box:
        assert lo < hi
        out *= hi - lo
    return out


def overlap(left, right):
    out = F(1)
    for (lo, hi), (bot, top) in zip(left, right):
        out *= max(F(0), min(hi, top) - max(lo, bot))
    return out


n, N = 4, 3  # Different from all three frozen baseline diagnostics.
theta, level, a = F(5, 8), F(7, 5), F(9, 4)
target = tuple((F(0), F(2)) for _ in range(n))
vstar = volume(target)
A = (1 - theta / 2) * vstar
B = 3 * theta * vstar / 4
sources, supports, radii = [], [], []
previous_delta = F(0)

for j in range(1, N + 1):
    R = 3 - F(j, 2 * N)
    delta = F(j, 2 * (N + 1))
    cube = ((delta - 2 * R, delta),) + ((-R, R),) * (n - 1)
    P = ((previous_delta, delta),) + ((F(0), F(1)),) * (n - 1)
    bulk = ((previous_delta, delta), (F(-1), F(0))) + ((F(0), F(1)),) * (n - 2)
    sources.append(cube)
    radii.append(R)
    pmass = level * A / N
    bmass = level * (2 * volume(cube) - A / N)
    assert pmass > 0 and bmass > 0
    supports.extend([(f'P{j}', P, pmass), (f'B{j}', bulk, bmass)])
    previous_delta = delta

Pstar = ((F(1), F(2)),) * n
supports.append(('Pstar', Pstar, level * B))
assert all(radii[j] > radii[j + 1] for j in range(N - 1))
assert radii[-1] > 1
for (_, C, _), (_, D, _) in combinations(supports, 2):
    assert overlap(C, D) == 0

# For this example every positive support piece is wholly in or wholly out of
# each witness. This permits exact independent union/fresh calculations.
for _, C, _ in supports:
    for V in sources + [target]:
        assert overlap(C, V) in (0, volume(C))

selected, ownership, records = [], [], []
for index, V in enumerate(sources + [target]):
    witness_mass = sum(m * overlap(C, V) / volume(C) for _, C, m in supports)
    fresh_pieces = [(name, C, m) for name, C, m in supports
                    if overlap(C, V) == volume(C)
                    and not any(overlap(C, old) == volume(C) for old in selected)]
    fresh_mass = sum(m for _, _, m in fresh_pieces)
    assert witness_mass > level * volume(V)
    decision = fresh_mass >= theta * level * volume(V)
    assert decision == (index < N)
    if decision:
        assert {name for name, _, _ in fresh_pieces} == {f'P{index + 1}', f'B{index + 1}'}
        assert fresh_mass == 2 * level * volume(V)
        ownership.append(fresh_pieces)
        selected.append(V)
    else:
        assert {name for name, _, _ in fresh_pieces} == {'Pstar'}
        assert fresh_mass == level * B
    records.append({'index': index + 1, 'selected': decision,
                    'witness_mass': str(witness_mass), 'fresh_mass': str(fresh_mass),
                    'threshold': str(theta * level * volume(V))})

contributions = [sum(m * overlap(C, target) / volume(C) for _, C, m in pieces)
                 for pieces in ownership]
assert contributions == [level * A / N] * N
assert sum(contributions) == level * A
assert level * A > (1 - theta) * level * vstar
assert level * B < theta * level * vstar
assert level * (A + B) > level * vstar

# The halo estimates use the same normalization as the frozen statement.
ca = 1 + 1 / a
halo_threshold = (1 - theta) * level / ca ** n
assert level * A / (ca ** n * vstar) > halo_threshold
assert theta * level / ca ** n > halo_threshold
source_mass = sum(m for pieces in ownership for _, _, m in pieces)
input_mass = sum(m for _, _, m in supports)
assert input_mass == source_mass + level * B

result = {
    'status': 'all independent exact checks passed', 'python': platform.python_version(),
    'source': str(SOURCE), 'source_sha256': EXPECTED,
    'parameters': {'n': n, 'N': N, 'theta': str(theta), 'lambda': str(level), 'a': str(a)},
    'normalized_target_witness_mass': str((A + B) / vstar),
    'normalized_target_fresh_mass': str(B / vstar),
    'normalized_target_captured_mass': str(A / vstar),
    'each_captured_share': str(F(1, N)),
    'each_share_of_rejection_threshold': str(A / (N * (1 - theta) * vstar)),
    'halo_threshold': str(halo_threshold), 'source_mass': str(source_mass),
    'input_mass': str(input_mass), 'decisions': records,
    'scope': 'one new rational finite case; no production solve, halo-volume approximation, or uniform-bound inference',
}
out = Path(__file__).with_name('microcheck.json')
out.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ('status', 'parameters',
    'normalized_target_witness_mass', 'normalized_target_fresh_mass',
    'normalized_target_captured_mass', 'each_captured_share')}, indent=2))
