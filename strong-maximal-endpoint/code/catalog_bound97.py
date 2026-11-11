#!/usr/bin/env python3
"""Exact finite catalog-density lower bounds for the r097 grid cases.

This enumerates rational cell volumes and the elementary minimal-subcover
bound. It does not invoke or modify any MILP solver or solver output.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction as F
from pathlib import Path

DYADIC6 = (F(0), F(1, 16), F(1, 8), F(1, 4), F(1, 2), F(3, 4), F(1))
CASES = (
    {'dimension': 2, 'cells_per_axis': 6, 'grid': 'dyadic6', 'p': '19/20', 'theta': '9/10'},
    {'dimension': 2, 'cells_per_axis': 6, 'grid': 'dyadic6', 'p': '39/40', 'theta': '9/10'},
    {'dimension': 2, 'cells_per_axis': 6, 'grid': 'dyadic6', 'p': '99/100', 'theta': '99/100'},
    {'dimension': 2, 'cells_per_axis': 8, 'grid': 'uniform', 'p': '19/20', 'theta': '9/10'},
    {'dimension': 2, 'cells_per_axis': 8, 'grid': 'uniform', 'p': '39/40', 'theta': '99/100'},
    {'dimension': 2, 'cells_per_axis': 8, 'grid': 'uniform', 'p': '99/100', 'theta': '3/4'},
    {'dimension': 3, 'cells_per_axis': 3, 'grid': 'uniform', 'p': '19/20', 'theta': '99/100'},
    {'dimension': 3, 'cells_per_axis': 3, 'grid': 'uniform', 'p': '39/40', 'theta': '3/4'},
    {'dimension': 3, 'cells_per_axis': 3, 'grid': 'uniform', 'p': '99/100', 'theta': '9/10'},
    {'dimension': 3, 'cells_per_axis': 4, 'grid': 'uniform', 'p': '19/20', 'theta': '3/4'},
    {'dimension': 3, 'cells_per_axis': 4, 'grid': 'uniform', 'p': '39/40', 'theta': '9/10'},
    {'dimension': 3, 'cells_per_axis': 4, 'grid': 'uniform', 'p': '99/100', 'theta': '99/100'},
)


def cell_volumes(spec):
    m, dim = spec['cells_per_axis'], spec['dimension']
    if spec['grid'] == 'dyadic6':
        assert m == 6
        breaks = DYADIC6
    elif spec['grid'] == 'uniform':
        breaks = tuple(F(i, m) for i in range(m + 1))
    else:
        raise ValueError(spec['grid'])
    widths = tuple(breaks[i + 1] - breaks[i] for i in range(m))
    vols = [F(1)]
    for _ in range(dim):
        vols = [a * w for a in vols for w in widths]
    assert sum(vols, F(0)) == 1
    return sorted(vols)


def one_case(spec):
    p, theta = F(spec['p']), F(spec['theta'])
    volumes = cell_volumes(spec)
    prefix = [F(0)]
    for v in volumes:
        prefix.append(prefix[-1] + v)

    h_rows = []
    for h in range(1, len(volumes) + 1):
        s_h = prefix[h]
        c_lo = max(theta, s_h)
        c_hi = F(1)
        a = 1 - (1 - p) * h
        b = (1 - p) * (h - 1) * s_h
        # g1(C)-g2(C) = (a-p/h) C+b. A convex max of two
        # affine functions reaches its minimum at an endpoint or a crossing.
        candidates = [c_lo, c_hi]
        slope = a - p / h
        if slope:
            cross = -b / slope
            if c_lo <= cross <= c_hi:
                candidates.append(cross)
        evaluated = []
        for c in candidates:
            g1 = a * c + b
            g2 = p * c / h
            evaluated.append((max(g1, g2), c, g1, g2))
        best = min(evaluated, key=lambda q: q[0])
        h_rows.append({
            'h': h,
            's_h': str(s_h),
            'C_interval': [str(c_lo), str(c_hi)],
            'bound': str(best[0]),
            'C_minimizer': str(best[1]),
            'affine_loss_bound_at_minimizer': str(best[2]),
            'max_box_bound_at_minimizer': str(best[3]),
        })

    best_row = min(h_rows, key=lambda r: F(r['bound']))
    return {
        **spec,
        'p_fraction': str(p), 'theta_fraction': str(theta),
        'cell_count': len(volumes),
        'cell_volume_min': str(volumes[0]),
        'cell_volume_max': str(volumes[-1]),
        'bound_fraction': best_row['bound'],
        'bound_decimal_12': f'{float(F(best_row["bound"])):.12f}',
        'minimizing_h': best_row['h'],
        'minimizing_C': best_row['C_minimizer'],
        'minimizing_s_h': best_row['s_h'],
        'all_h_bounds': h_rows,
    }


def main():
    rows = [one_case(case) for case in CASES]
    checkpoint_path = Path(__file__).resolve().parent / 'periodic-search' / 'periodic_search97.checkpoint.json'
    candidate_upper_bounds = {}
    if checkpoint_path.exists():
        checkpoint = json.loads(checkpoint_path.read_text())
        checkpoint_keys = {
            tuple(row.get(k) for k in ('dimension', 'cells_per_axis', 'grid', 'p', 'theta'))
            for row in checkpoint.get('completed_rows', [])
        }
        expected_keys = {
            tuple(row.get(k) for k in ('dimension', 'cells_per_axis', 'grid', 'p', 'theta'))
            for row in rows
        }
        if (checkpoint.get('completed_cases') != len(CASES)
                or checkpoint_keys != expected_keys
                or len(checkpoint_keys) != len(CASES)):
            raise AssertionError('r097 checkpoint case keys do not exactly match the catalog cases')
        for row in checkpoint.get('completed_rows', []):
            candidate = row.get('candidate') or {}
            if candidate.get('certificate_pass'):
                key = tuple(row.get(k) for k in ('dimension', 'cells_per_axis', 'grid', 'p', 'theta'))
                candidate_upper_bounds[key] = candidate['source_density_fraction']
    for row in rows:
        key = tuple(row.get(k) for k in ('dimension', 'cells_per_axis', 'grid', 'p', 'theta'))
        upper = candidate_upper_bounds.get(key)
        row['certified_candidate_density_upper_bound_fraction'] = upper
        if upper is not None and F(row['bound_fraction']) > F(upper):
            raise AssertionError(f"catalog lower bound exceeds witness upper bound: {key}")
    result = {
        'scope': 'Exact rational lower bound for the periodic finite grid catalog only; no continuous-operator or all-dimension extrapolation.',
        'method': 'For each h, minimize max((1-(1-p)h)C+(1-p)(h-1)s_h, pC/h) over C in [max(theta,s_h),1], then minimize over h.',
        'case_count': len(rows),
        'cases': [{k: v for k, v in row.items() if k != 'all_h_bounds'} for row in rows],
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
