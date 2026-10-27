#!/usr/bin/env python3
"""Exact rational diagnostics for the fresh-witness-mass selector. No scan."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform

ROOT = Path(__file__).resolve().parent


def volume(box):
    v = F(1)
    for lo, hi in box:
        v *= hi - lo
    return v


def intersection(a, b):
    out = tuple((max(x[0], y[0]), min(x[1], y[1])) for x, y in zip(a, b))
    return out if all(lo < hi for lo, hi in out) else None


def subtract(a, b):
    c = intersection(a, b)
    if c is None:
        return [a]
    core = list(a)
    out = []
    for k, (lo, hi) in enumerate(c):
        if core[k][0] < lo:
            piece = core.copy()
            piece[k] = (core[k][0], lo)
            out.append(tuple(piece))
        if hi < core[k][1]:
            piece = core.copy()
            piece[k] = (hi, core[k][1])
            out.append(tuple(piece))
        core[k] = (lo, hi)
    assert sum(map(volume, out), F(0)) + volume(c) == volume(a)
    return out


def integral(pieces, domain, previous=()):
    ans = F(0)
    for box, mass in pieces:
        cut = intersection(box, domain)
        remaining = [cut] if cut else []
        for old in previous:
            remaining = [r for p in remaining for r in subtract(p, old)]
        ans += mass * sum(map(volume, remaining), F(0)) / volume(box)
    return ans


def test_fragmentation(n, N):
    theta = F(3, 4)
    target = tuple((F(0), F(2)) for _ in range(n))
    vt = volume(target)
    captured = (1 - theta / 2) * vt
    fresh = 3 * theta * vt / 4
    witnesses = []
    pieces = []
    previous_delta = F(0)
    for j in range(1, N + 1):
        radius = 3 - F(j, 2 * N)
        delta = F(j, 2 * (N + 1))
        witness = ((delta - 2 * radius, delta),) + ((-radius, radius),) * (n - 1)
        witnesses.append(witness)
        hit = ((previous_delta, delta),) + ((F(0), F(1)),) * (n - 1)
        bulk = ((previous_delta, delta), (F(-1), F(0))) + ((F(0), F(1)),) * (n - 2)
        pieces += [(hit, captured / N), (bulk, 2 * volume(witness) - captured / N)]
        previous_delta = delta
    pieces.append((((F(1), F(2)),) * n, fresh))
    selected = []
    residual_masses = []
    for witness in witnesses:
        total = integral(pieces, witness)
        residual = integral(pieces, witness, selected)
        assert total > volume(witness)
        assert residual == 2 * volume(witness)
        assert residual >= theta * volume(witness)
        selected.append(witness)
        residual_masses.append(residual)
    target_total = integral(pieces, target)
    target_fresh = integral(pieces, target, selected)
    assert target_total == captured + fresh > vt
    assert target_fresh == fresh < theta * vt
    contributions = []
    for j, witness in enumerate(selected):
        cut = intersection(target, witness)
        contribution = integral(pieces, cut, selected[:j]) if cut else F(0)
        assert contribution == captured / N
        contributions.append(contribution)
    assert sum(contributions, F(0)) == captured
    assert captured > (1 - theta) * vt
    return {
        'dimension': n, 'sources': N, 'theta': str(theta),
        'target_mass_over_target_volume': str(target_total / vt),
        'target_fresh_over_target_volume': str(target_fresh / vt),
        'captured_over_target_volume': str(captured / vt),
        'each_source_over_captured': str(F(1, N)),
        'each_source_over_rejection_threshold': str((captured / N) / ((1 - theta) * vt)),
        'all_source_witnesses_strict': True,
        'all_selection_decisions_exact': True,
        'mass_conservation_exact': sum(residual_masses, F(0)) + fresh == sum(m for _, m in pieces),
    }


sources = {
    'R032': ROOT.parents[1] / 'round032/joint_cover_theory/proof.md',
    'R035': ROOT.parents[1] / 'round035/joint_cover/proof.md',
    'R037': Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round037/general_cover/PROOF.md'),
    'R037_clarifications': Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round037/general_cover/CLARIFICATIONS.md'),
}
result = {
    'python': platform.python_version(),
    'arithmetic': 'Python standard-library fractions.Fraction, exact rational arithmetic',
    'checks': [test_fragmentation(2, 8), test_fragmentation(2, 32), test_fragmentation(3, 8)],
    'source_sha256': {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in sources.items()},
    'scope': 'Finite selector interface checks only; no uniform upper-bound improvement.',
}
(ROOT / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
