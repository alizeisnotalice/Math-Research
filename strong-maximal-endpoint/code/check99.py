"""Independent exact replay and comparison for the compact r099 grid search."""
import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def qstr(x):
    return f'{x.numerator}/{x.denominator}'


def log_upper(x, terms=32):
    """Directed rational upper bound for log(x), using log x=2 atanh z."""
    x = F(x)
    if x < 1:
        raise ValueError('x must be at least 1')
    z = (x-1)/(x+1)
    partial = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    tail = (2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))) if z else F(0)
    return partial+tail


def cell_geometry(row):
    n, m = int(row['dimension']), int(row['cells_per_axis'])
    breaks = tuple(F(x) for x in row['axis_breaks'])
    assert n > 0 and m > 0 and len(breaks) == m+1
    assert breaks[0] == 0 and breaks[-1] == 1
    assert all(a < b for a, b in zip(breaks, breaks[1:]))
    coords = tuple(itertools.product(range(m), repeat=n))
    index = {x: i for i, x in enumerate(coords)}
    widths = tuple(breaks[i+1]-breaks[i] for i in range(m))
    volumes = tuple(math.prod(widths[j] for j in c) for c in coords)
    assert sum(volumes, F(0)) == 1
    return n, m, breaks, coords, index, volumes


def arc_id(origin, length, m, compact):
    if compact:
        if length == m:
            assert origin == 0
            return m*(m-1)
        assert 1 <= length < m
        return (length-1)*m+origin
    return origin*m+(length-1)


def replay_case(row, compact):
    n, m, breaks, coords, index, volumes = cell_geometry(row)
    p, theta = F(row['p']), F(row['theta'])
    assert 0 < p < 1 and 0 < theta <= 1
    cand = row['candidate']
    den = int(cand['source_rounding_denominator'])
    nums = cand['source_cell_fraction_numerators']
    assert den > 0 and len(nums) == len(coords)
    assert all(type(x) is int and 0 <= x <= den for x in nums)
    assert cand['source_cell_coordinates_lexicographic'] == [list(c) for c in coords]
    assert cand['certificate_pass'] is True
    density = sum((volumes[i]*F(x, den) for i, x in enumerate(nums)), F(0))
    assert 0 < density <= 1
    assert density == F(cand['source_density_fraction'])
    assert (density.numerator, density.denominator) == (
        int(cand['source_density_numerator']), int(cand['source_density_denominator']))

    base = m*(m-1)+1 if compact else m*m
    union = set()
    margins, rect_records = [], []
    for rect in cand['certified_selected_rectangles']:
        axes = [tuple(map(int, a)) for a in rect['axis_origin_length_cells']]
        assert len(axes) == n
        axis_cells = []
        rid = 0
        for origin, length in axes:
            assert 0 <= origin < m and 1 <= length <= m
            if compact and length == m:
                assert origin == 0
            axis_cells.append(tuple((origin+j) % m for j in range(length)))
            rid = rid*base+arc_id(origin, length, m, compact)
        assert int(rect['id']) == rid
        ids = {index[c] for c in itertools.product(*axis_cells)}
        listed = rect['receiver_cell_ids']
        assert len(listed) == len(set(listed)) and set(listed) == ids
        volume = sum((volumes[i] for i in ids), F(0))
        mass = sum((volumes[i]*F(nums[i], den) for i in ids), F(0))
        mean, margin = mass/volume, mass/volume-p
        assert mean == F(rect['mean_fraction'])
        assert margin == F(rect['mean_minus_p_fraction']) and margin > 0
        margins.append(margin)
        union.update(ids)
        rect_records.append({'id': rid, 'axis_origin_length_cells': [list(a) for a in axes],
                             'receiver_cell_ids': sorted(ids),
                             'volume_fraction': qstr(volume),
                             'mean_fraction': qstr(mean),
                             'mean_minus_p_fraction': qstr(margin)})
    assert rect_records and cand['all_final_rectangles_strictly_above_p'] is True
    assert union == set(cand['receiver_cell_ids_exact_union'])
    covered = sum((volumes[i] for i in union), F(0))
    assert covered == F(cand['theta_actual_fraction']) and covered >= theta
    assert covered == F(int(cand['theta_actual_numerator']), int(cand['theta_actual_denominator']))
    ratio = covered/density
    assert ratio == F(cand['theta_actual_over_density_fraction'])
    assert ratio <= 1/density

    K = 1000*n
    boundary = F(K-2, K)**n
    finite_halo_lower = ratio*boundary
    log_hi = log_upper(1/p)
    phi_hi = (1/p)*(1+log_hi)**(n-1)
    finite_C_lower = finite_halo_lower/phi_hi

    return {
        'case': {'dimension': n, 'cells_per_axis': m, 'grid': row['grid'],
                 'p': qstr(p), 'theta': qstr(theta)},
        'axis_breaks': [qstr(x) for x in breaks],
        'source_rounding_denominator': den,
        'source_cell_fraction_numerators': nums,
        'source_density': qstr(density),
        'theta_actual': qstr(covered),
        'selected_rectangles': rect_records,
        'selected_rectangle_count': len(rect_records),
        'minimum_strict_mean_margin': qstr(min(margins)),
        'periodic_halo_ratio_lower': qstr(ratio),
        'periodic_halo_ratio_upper_trivial': qstr(1/density),
        'finite_box_periods_per_axis': K,
        'finite_box_boundary_factor': qstr(boundary),
        'finite_box_indicator_ratio_lower': qstr(finite_halo_lower),
        'log_1_over_p_upper': qstr(log_hi),
        'orlicz_denominator_upper': qstr(phi_hi),
        'finite_endpoint_C_lower': qstr(finite_C_lower),
        'solver_status': row.get('status'),
        'solver_optimality_certified': False,
        'solver_mip_gap_numeric_only': row.get('mip_gap_numeric'),
        'solver_dual_bound_numeric_only': row.get('dual_bound_numeric_only'),
    }


def decimal(x, digits=12):
    return format(float(F(x)), f'.{digits}g')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project-root', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    root = Path(args.project_root)
    rnd = root/'rounds/r099'
    raw_new = json.loads((rnd/'periodic99.json').read_text())
    raw_old = json.loads((rnd/'periodic96_reference.json').read_text())
    assert raw_new['pass_all_candidates'] and raw_old['pass_all_candidates']
    assert raw_new['case_count'] == raw_old['case_count'] == 12
    assert len(raw_new['cases']) == len(raw_old['cases']) == 12

    new_map, old_map = {}, {}
    for row in raw_new['cases']:
        cert = replay_case(row, compact=True)
        key = tuple(cert['case'][x] for x in ('dimension','cells_per_axis','grid','p','theta'))
        assert key not in new_map
        new_map[key] = cert
    for row in raw_old['cases']:
        cert = replay_case(row, compact=False)
        key = tuple(cert['case'][x] for x in ('dimension','cells_per_axis','grid','p','theta'))
        assert key not in old_map
        old_map[key] = cert
    assert set(new_map) == set(old_map)

    comparisons = []
    better = tied = worse = 0
    total_rectangles = 0
    for key in sorted(new_map):
        new, old = new_map[key], old_map[key]
        dn, do = F(new['source_density']), F(old['source_density'])
        delta = do-dn
        if delta > 0:
            better += 1
            state = 'improved'
        elif delta == 0:
            tied += 1
            state = 'tied'
        else:
            worse += 1
            state = 'worse'
        assert new['case'] == old['case']
        total_rectangles += new['selected_rectangle_count']
        comparisons.append({
            'case': new['case'],
            'density_r096': old['source_density'],
            'density_r099': new['source_density'],
            'r096_minus_r099_density': qstr(delta),
            'cost_comparison': state,
            'selected_rectangles_r096': old['selected_rectangle_count'],
            'selected_rectangles_r099': new['selected_rectangle_count'],
            'r099_minimum_strict_mean_margin': new['minimum_strict_mean_margin'],
            'r099_theta_actual': new['theta_actual'],
            'r099_periodic_halo_ratio_lower': new['periodic_halo_ratio_lower'],
            'r099_finite_box_boundary_factor': new['finite_box_boundary_factor'],
            'r099_finite_box_indicator_ratio_lower': new['finite_box_indicator_ratio_lower'],
            'r099_orlicz_denominator_upper': new['orlicz_denominator_upper'],
            'r099_finite_endpoint_C_lower': new['finite_endpoint_C_lower'],
            'r099_finite_endpoint_C_lower_decimal': decimal(new['finite_endpoint_C_lower']),
            'r099_source_rounding_denominator': new['source_rounding_denominator'],
            'r099_source_cell_fraction_numerators': new['source_cell_fraction_numerators'],
            'r099_axis_breaks': new['axis_breaks'],
            'r099_selected_rectangles': new['selected_rectangles'],
            'solver_status': new['solver_status'],
            'solver_optimality_certified': False,
        })

    assert better+ tied+worse == 12
    sources = ['periodic99.py', 'check99.py']
    inputs = ['periodic99.json', 'periodic96_reference.json']
    documents = ['report.md', 'audit99.md']
    validation = {
        'schema': 'r099-package-validation-v1',
        'pass_all': True,
        'case_count': 12,
        'certified_compact_rectangles_replayed': total_rectangles,
        'density_cost_comparison': {'improved': better, 'tied': tied, 'worse': worse,
                                    'meaning': 'lower exact source density is a better candidate'},
        'solver_optimality_certified': False,
        'best_state': 'best certified finite-grid candidate within each parameter case; no solver or global optimality proof',
        'scope': ('Exact finite periodic candidate replay, comparison with the r096 witnesses, and '
                  'finite-box lower / periodic-volume upper ledger. No general maximal-operator upper '
                  'bound, onset dimension, or continuous d_N claim.'),
        'cases': comparisons,
        'source_sha256': {name: sha256(root/'code'/name) for name in sources},
        'input_sha256': {name: sha256(rnd/name) for name in inputs},
        'document_sha256': {name: sha256(rnd/name) for name in documents},
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(validation, indent=2)+'\n')
    print(json.dumps({k:v for k,v in validation.items() if k != 'cases'} |
                     {'cases': [{'N': c['case']['dimension'], 'm': c['case']['cells_per_axis'],
                                 'p': c['case']['p'], 'theta': c['case']['theta'],
                                 'comparison': c['cost_comparison'],
                                 'd96': decimal(c['density_r096']),
                                 'd99': decimal(c['density_r099']),
                                 'C99_lower': c['r099_finite_endpoint_C_lower_decimal']}
                                for c in comparisons]}, indent=2))


if __name__ == '__main__':
    main()
