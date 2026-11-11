"""Exact portable replayer for finite periodic-grid candidate witnesses."""
import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction as Q
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frac_text(value):
    return f'{value.numerator}/{value.denominator}'


def log_upper_atanh(x, terms=32):
    """Exact rational upper bound for log(x), x >= 1."""
    x = Q(x)
    if x < 1:
        raise ValueError('log upper routine requires x >= 1')
    z = (x-1)/(x+1)
    partial = 2*sum((z**(2*j+1)/Q(2*j+1) for j in range(terms)), Q(0))
    tail = 2*z**(2*terms+1)/(Q(2*terms+1)*(1-z*z)) if z else Q(0)
    return partial+tail


def axis_arcs(m):
    """Original full torus catalog, with axis enumeration matching periodic97."""
    return [(origin, length,
             tuple((origin+offset) % m for offset in range(length)))
            for origin in range(m) for length in range(1, m+1)]


def replay_case(row):
    dim, m = int(row['dimension']), int(row['cells_per_axis'])
    p, theta = Q(row['p']), Q(row['theta'])
    breaks = tuple(Q(x) for x in row['axis_breaks'])
    if not (dim >= 1 and m >= 1 and len(breaks) == m+1):
        raise AssertionError('invalid dimension, grid, or axis break count')
    if breaks[0] != 0 or breaks[-1] != 1 or any(a >= b for a, b in zip(breaks, breaks[1:])):
        raise AssertionError('axis breaks must strictly partition [0,1]')
    if not (0 < p < 1 and 0 < theta <= 1):
        raise AssertionError('p and theta outside the supported range')

    coords = tuple(itertools.product(range(m), repeat=dim))
    cell_index = {c: i for i, c in enumerate(coords)}
    widths = tuple(breaks[i+1]-breaks[i] for i in range(m))
    volumes = tuple(math.prod(widths[i] for i in c) for c in coords)
    if sum(volumes, Q(0)) != 1:
        raise AssertionError('grid cell volumes do not sum to one')

    cand = row['candidate']
    den = int(cand['source_rounding_denominator'])
    nums = cand['source_cell_fraction_numerators']
    if den <= 0 or len(nums) != len(coords) or any(type(x) is not int or not 0 <= x <= den for x in nums):
        raise AssertionError('invalid rational source cell data')
    if cand['source_cell_coordinates_lexicographic'] != [list(c) for c in coords]:
        raise AssertionError('source-cell coordinate order mismatch')
    if not cand['certificate_pass']:
        raise AssertionError('source candidate is not marked as certified')

    density = sum((v*Q(x, den) for v, x in zip(volumes, nums)), Q(0))
    if density <= 0 or density > 1 or density != Q(cand['source_density_fraction']):
        raise AssertionError('exact source density mismatch')
    if (density.numerator != int(cand['source_density_numerator']) or
            density.denominator != int(cand['source_density_denominator'])):
        raise AssertionError('source-density numerator/denominator mismatch')

    arc_index = {(o, length): i for i, (o, length, _) in enumerate(axis_arcs(m))}
    union = set()
    strict_margins = []
    rects_out = []
    for rect in cand['certified_selected_rectangles']:
        axes = rect['axis_origin_length_cells']
        if len(axes) != dim:
            raise AssertionError('rectangle dimension mismatch')
        axis_cells = []
        for pair in axes:
            origin, length = map(int, pair)
            if not (0 <= origin < m and 1 <= length <= m):
                raise AssertionError('rectangle endpoint is outside the torus catalog')
            axis_cells.append(tuple((origin+j) % m for j in range(length)))
        rect_id = 0
        for pair in axes:
            rect_id = rect_id*m*m + arc_index[tuple(map(int, pair))]
        if int(rect['id']) != rect_id:
            raise AssertionError('rectangle id does not match the original catalog enumeration')

        cell_ids = {cell_index[c] for c in itertools.product(*axis_cells)}
        listed_ids = rect['receiver_cell_ids']
        if len(listed_ids) != len(set(listed_ids)) or set(listed_ids) != cell_ids:
            raise AssertionError('rectangle receiver-cell set mismatch')
        volume = sum((volumes[i] for i in cell_ids), Q(0))
        mass = sum((volumes[i]*Q(nums[i], den) for i in cell_ids), Q(0))
        mean = mass/volume
        margin = mean-p
        if mean != Q(rect['mean_fraction']) or margin != Q(rect['mean_minus_p_fraction']) or margin <= 0:
            raise AssertionError('rectangle mean or strict margin mismatch')
        strict_margins.append(margin)
        union.update(cell_ids)
        rects_out.append({'id': int(rect['id']), 'axes': [list(map(int, a)) for a in axes],
                          'receiver_cell_ids': sorted(cell_ids),
                          'volume_fraction': frac_text(volume),
                          'mean_fraction': frac_text(mean),
                          'mean_minus_p_fraction': frac_text(margin)})

    if not rects_out:
        raise AssertionError('candidate has no selected certified rectangles')
    if union != set(cand['receiver_cell_ids_exact_union']):
        raise AssertionError('exact receiver union mismatch')
    covered = sum((volumes[i] for i in union), Q(0))
    if covered != Q(cand['theta_actual_fraction']) or covered < theta:
        raise AssertionError('receiver coverage fraction mismatch or below theta')
    if covered != Q(int(cand['theta_actual_numerator']), int(cand['theta_actual_denominator'])):
        raise AssertionError('receiver coverage numerator/denominator mismatch')
    periodic_ratio = covered/density
    if periodic_ratio != Q(cand['theta_actual_over_density_fraction']):
        raise AssertionError('periodic halo density ratio mismatch')
    if periodic_ratio > 1/density:
        raise AssertionError('periodic halo ratio exceeds the unit torus volume bound')
    if cand.get('all_final_rectangles_strictly_above_p') is not True:
        raise AssertionError('raw strictness summary flag is not true')

    # Tiling a K^N box by K periods and discarding one period at each face
    # leaves every selected box available to interior receiver points.
    K = 1000*dim
    boundary_factor = Q(K-2, K)**dim
    finite_halo_lower = periodic_ratio*boundary_factor
    log_upper = log_upper_atanh(1/p)
    phi_upper = (1/p)*(1+log_upper)**(dim-1)
    finite_C_lower = finite_halo_lower/phi_upper

    return {
        'dimension': dim,
        'cells_per_axis': m,
        'grid': row['grid'],
        'p': frac_text(p),
        'theta_required': frac_text(theta),
        'axis_breaks': [frac_text(x) for x in breaks],
        'source_rounding_denominator': den,
        'source_cell_fraction_numerators': nums,
        'source_density': frac_text(density),
        'source_realization': 'In each grid cell, occupy a half-open coordinate-1 slab of the listed relative length.',
        'selected_rectangles': rects_out,
        'selected_rectangle_count': len(rects_out),
        'minimum_strict_mean_margin': frac_text(min(strict_margins)),
        'receiver_cell_ids_union': sorted(union),
        'receiver_fraction_actual': frac_text(covered),
        'periodic_halo_ratio_lower': frac_text(periodic_ratio),
        'periodic_halo_ratio_upper_trivial': frac_text(1/density),
        'finite_box_periods_per_axis': K,
        'finite_box_boundary_factor': frac_text(boundary_factor),
        'finite_box_indicator_ratio_lower': frac_text(finite_halo_lower),
        'log_1_over_p_upper_rational': frac_text(log_upper),
        'orlicz_denominator_upper': frac_text(phi_upper),
        'finite_endpoint_C_lower': frac_text(finite_C_lower),
        'solver_diagnostics_noncertifying': {
            'status': row.get('status'),
            'has_incumbent': row.get('has_incumbent'),
            'seconds': row.get('solver_seconds'),
            'mip_gap_numeric_only': row.get('mip_gap_numeric'),
            'dual_bound_numeric_only': row.get('dual_bound_numeric_only')},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_json')
    parser.add_argument('summary_json')
    args = parser.parse_args()
    raw_path, summary_path = Path(args.raw_json), Path(args.summary_json)
    raw = json.loads(raw_path.read_text())
    rows = raw.get('completed_rows', raw.get('cases'))
    if not isinstance(rows, list) or not rows:
        raise AssertionError('raw file has no completed candidate rows')
    completed = int(raw.get('completed_cases', len(rows)))
    expected = int(raw.get('expected_cases', raw.get('case_count', len(rows))))
    if completed != len(rows) or not 0 < completed <= expected:
        raise AssertionError('checkpoint row count metadata mismatch')
    summaries = [replay_case(row) for row in rows]
    keys = [(r['dimension'], r['cells_per_axis'], r['grid'], r['p'], r['theta_required']) for r in summaries]
    if len(set(keys)) != len(keys):
        raise AssertionError('duplicate completed parameter case')
    out = {
        'schema': 'periodic97-exact-finite-candidate-summary-v1',
        'pass': True,
        'complete': completed == expected,
        'completed_cases': completed,
        'expected_cases': expected,
        'raw_sha256': sha256(raw_path),
        'scope': ('Exact finite-grid candidate replay, a finite-box lower ledger, and a periodic-volume upper check. '
                  'Solver diagnostics are non-certifying; this does not prove a general H_N upper bound, '
                  'an onset dimension, or optimality of the finite search.'),
        'cases': summaries,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'cases'} |
                     {'case_summaries': [{'N': r['dimension'], 'm': r['cells_per_axis'],
                                          'p': r['p'], 'theta': r['theta_required'],
                                          'source_density': r['source_density'],
                                          'finite_C_lower_approx': float(Q(r['finite_endpoint_C_lower'])),
                                          'rectangles': r['selected_rectangle_count']}
                                         for r in summaries]}, indent=2))


if __name__ == '__main__':
    main()
