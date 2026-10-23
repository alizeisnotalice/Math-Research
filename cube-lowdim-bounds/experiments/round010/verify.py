#!/usr/bin/env python3
"""Verify the saved rational leaves, without running search or bisection."""
import importlib.util
import json
import sys
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def verify():
    data = json.loads((ROOT / 'certificate.json').read_text())
    source = (ROOT / data['core_dependency']).resolve()
    assert sha256(source.read_bytes()).hexdigest() == data['core_sha256']
    spec = importlib.util.spec_from_file_location('verified_round008_core', source)
    core = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = core
    spec.loader.exec_module(core)
    p = data['parameters']
    core.R, core.E2, core.E3 = map(Q, (p['R'], p['e2'], p['e3']))
    assert core.R > 1 and core.E2 >= 0 and core.E3 >= 0
    core.PREVIOUS_UPPER = Q(21773155252261, 5615797500000)
    leaves = data['accepted_leaves']
    assert leaves
    cursor = Q(1, 10)
    lower_bounds = []
    for row in leaves:
        lo, hi, saved = map(Q, (row['lo'], row['hi'], row['lower']))
        assert lo == cursor and lo < hi <= 1
        assert hi <= Q(3, 4) or lo >= Q(3, 4)
        actual = core.lower_on_box(core.I(lo, hi)).lo
        assert actual == saved and saved > 0
        lower_bounds.append(saved)
        cursor = hi
    assert cursor == 1
    digest = sha256('\n'.join(f"{r['lo']}/{r['hi']}" for r in leaves).encode()).hexdigest()
    assert digest == data['inner']['accepted_leaf_sha256']
    assert len(leaves) == data['inner']['accepted_boxes']
    assert min(lower_bounds) == Q(data['inner']['minimum_accepted_lower'])
    small, outer = core.small_r_certificate(), core.outer_checks()
    assert small['strictly_negative'] and outer['strict']
    assert outer == data['outer_and_cost']
    assert Q(outer['cost']) == Q(6938554331, 1790000000)
    print(json.dumps({'status': 'all_saved_leaves_and_scalar_conditions_verified',
                      'leaf_count': len(leaves), 'cost': outer['cost'],
                      'improvement': outer['improvement'],
                      'analytic_dependency': '../round008/theorem.md plus parameter_review.md'},
                     indent=2))


if __name__ == '__main__':
    verify()
