#!/usr/bin/env python3
"""Exact saved certificates and S3 witnesses; stdout only, no kernel calls."""
import hashlib
import json
import time
from fractions import Fraction as Q
from pathlib import Path

import sympy as s


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parents[1]
    old = json.loads((root/'lower/round001/results.json').read_text())
    new = json.loads((root/'lower/round002/results.json').read_text())
    lo, hi = map(Q, new['full_construction_comparison_interval'])
    finite = Q(new['finite_L1']['exact_ratio_lower'])
    assert lo <= hi and finite > Q(old['full_construction_comparison_interval'][1])
    assert finite > Q(2136208981884, 10**12)
    moment = json.loads((root/'upper/round002/hadamard_moment/results.json').read_text())
    iu = Q(moment['rigorous_dyadic_I_upper'])
    old_i = json.loads((root/'audit/round001/exact_budget_results.json').read_text())
    assert iu < Q(old_i['I_upper_expression_outward_interval'][0])
    target = Q(1, 194000)
    cost6 = 7**7*iu**6/target**7
    assert cost6 < 128**6
    directory = root/'upper/round002/g_probe'
    diag = json.loads((directory/'diagnostics.json').read_text())
    body = [s.Matrix([1]+[1 if (b>>k)&1 else -1 for k in range(6)])
            for b in range(64)]
    witnesses = []
    for candidate in diag['candidate_bundles']:
        rounded = candidate['rounded_rational_point']
        x = s.Matrix([s.Rational(a) for a in rounded['rational']])
        proof = rounded['S3']
        matrix = s.Matrix.hstack(*(body[j] for j in proof['body_indices']))
        coeff = s.Matrix([s.Rational(a) for a in proof['coefficients']])
        assert matrix*coeff == x and matrix.rank() <= 3
        original = s.Matrix([s.Rational(Q(float.fromhex(a))) for a in
                            candidate['original_optimized_point']['point']['floathex']])
        assert matrix.row_join(original).rank() == matrix.rank()
        near = [Q(a) for a in candidate['fixed_nearby_offS3_probe']['rational']]
        assert len({abs(a) for a in near}) == 7
        witnesses.append({'body_indices': proof['body_indices'],
                          'rounded_and_original_in_S3': True,
                          'nearby_off_S3_by_distinct_absolute_values': True})
    assert len(witnesses) == 8
    assert diag['optimization_calls'] == 868
    assert diag['diagnostic_full_point_calls'] == 216
    assert diag['high_precision_same_point_subset_verification_calls'] == 24
    assert diag['total_budget_charged_calls'] == 1108
    handles = []
    for path in sorted((directory/'starts').glob('*.json')):
        item = json.loads(path.read_text())
        assert item['status'] == 'finished'
        assert item['full_point_calls'] == len(item['trace']) <= 256
        handles.append(item['handle'])
    assert len(handles) == len(set(handles)) == 64
    print(json.dumps({
        'status': 'PASS primary exact round002 saved checks',
        'lower_complete_interval': [str(lo), str(hi)],
        'finite_L1_exact_lower': str(finite),
        'strict_margin_above_complete_round001_upper': str(finite-Q(old['full_construction_comparison_interval'][1])),
        'new_complete_I_upper': str(iu),
        'conditional_next_g_target': str(target),
        'conditional_cost_sixth_power': str(cost6),
        'conditional_gap_below_128_sixth_power': str(128**6-cost6),
        'conditional_global_g_proved': False,
        'singular_witnesses': witnesses,
        'all_64_handles_terminal': True,
        'new_C7_lower': True, 'new_C7_upper': False,
        'new_Green_or_shape_evaluations': 0,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'seconds': time.perf_counter()-started,
    }, indent=2))


if __name__ == '__main__':
    main()
