#!/usr/bin/env python3
"""Primary, read-only exact spectral enumeration and saved lower interfaces.

This checker does not recalculate a previously accepted lower union or turn
the numerical H diagnostic into an expectation certificate.
"""
import hashlib
import itertools
import json
import math
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import sympy as sp


def main():
    tic = time.perf_counter()
    root = Path(__file__).resolve().parents[1]
    folders = ['lower/round003', 'audit/round003/lower',
               'upper/round003/spectral_catalog',
               'upper/round003/truncated_green',
               'upper/round003/truncated_numeric']
    bindings = {}
    for folder in folders:
        base = root / folder
        data = (base / 'manifest.json').read_bytes()
        bindings[folder] = hashlib.sha256(data).hexdigest()
        entries = json.loads(data)['files']
        if isinstance(entries, dict):
            entries = [dict(path=k, sha256=v) if isinstance(v, str)
                       else dict(path=k, **v) for k, v in entries.items()]
        for item in entries:
            name = item.get('path', item.get('file', item.get('name')))
            path = (base/name).resolve()
            assert path.is_relative_to(base.resolve())
            assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256']
    prior = json.loads((root/'lower/round002/results.json').read_text())
    new = json.loads((root/'lower/round003/results.json').read_text())
    l, u = map(F, new['full_construction_comparison_interval'])
    old = F(prior['full_construction_comparison_interval'][1])
    finite = F(new['finite_L1']['exact_ratio_lower'])
    assert u >= l > finite > old > F(2)
    assert finite > F(2136429741978, 10**12)
    assert F(new['complete_round002_upper']) == old
    state = json.loads((root/'lower/round003/continuation_state.json').read_text())
    assert state['this_round_evaluated_indices'] == [51, 39]
    assert set(state['this_round_evaluated_indices']).isdisjoint(
        state['all_prior_R002_evaluated_event_indices'])
    assert state['remaining_new_events'] == 59 and not state['family_exhausted']
    v = [(1,)+s for s in itertools.product((-1, 1), repeat=6)]
    assert sp.Matrix(v).T*sp.Matrix(v) == 64*sp.eye(7)
    found, histogram = {}, Counter()
    for ordered in itertools.combinations_with_replacement(range(7), 7):
        x = tuple(reversed(ordered))
        if math.gcd(*x) != 1:
            continue
        dots = [sum(a*b for a, b in zip(row, x)) for row in v]
        zero = [row for row, dot in zip(v, dots) if dot == 0]
        rank = sp.Matrix(zero).rank() if zero else 0
        histogram[rank] += 1
        if rank == 6:
            mult = Counter(x)
            orbit = (2**sum(c != 0 for c in x)*math.factorial(7)
                     // math.prod(math.factorial(k) for k in mult.values()))
            found[x] = (sum(map(abs, dots)), orbit)
    catalog = json.loads((root/'upper/round003/spectral_catalog/catalog.json').read_text())
    saved = {tuple(r['primitive_xi']): (r['psi_rate'], r['orbit_cardinality'])
             for r in catalog['orbits']}
    assert found == saved and len(found) == 21
    assert sum(t[1] for t in found.values()) == 504868
    assert dict(histogram) == {0:954, 1:13, 2:30, 3:109, 4:242, 5:192, 6:21}
    for row in catalog['orbits']:
        x = sp.Matrix(row['primitive_xi'])
        b = sp.Matrix([v[i] for i in row['rank6_witness_indices']])
        assert b.rank() == 6 and b*x == sp.zeros(6, 1)
        cof = [(-1)**j*b[:, [k for k in range(7) if k != j]].det() for j in range(7)]
        assert cof == row['witness_signed_cofactors']
        assert all(c % 32 == 0 and abs(c) <= 216 for c in cof)
        assert [F(c, row['psi_rate']) for c in x] == list(map(F, row['unit_ball_vertex']))
    print(json.dumps(dict(status='PASS', manifest_bindings=bindings,
        new_lower_simple='2136429741978/1000000000000',
        finite_margin_vs_complete_R002=str(finite-old),
        spectral_orbits=21, spectral_vertices=504868,
        independent_rank_backend='SymPy exact rational rank',
        original_green_new_upper=False, numerical_H_expectation_certified=False,
        elapsed_seconds=time.perf_counter()-tic), indent=2))


if __name__ == '__main__':
    main()
