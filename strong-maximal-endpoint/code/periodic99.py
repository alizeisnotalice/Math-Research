"""Finite periodic torus-grid MILP search with exact candidate certificates.

This is a finite-grid search for candidates in the periodic source-density
proxy. HiGHS bounds are numerical only. Every reported witness candidate is
rechecked with rational arithmetic after rounding source cell fractions upward
to a 2**24 denominator.
"""
import argparse
import hashlib
import itertools
import json
import math
import multiprocessing as mp
import os
import time
import warnings
from fractions import Fraction as Q
from pathlib import Path

for _thread_var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                    'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_thread_var] = '1'

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

EPS = Q(1, 10**6)
SOURCE_DENOMINATOR = 2**24
SOLVER_THREADS = 1
TIME_LIMIT = 120
MAX_MODEL_BYTES = 100*2**20
WORKERS = 12
DYADIC6_BREAKS = (Q(0), Q(1, 16), Q(1, 8), Q(1, 4), Q(1, 2), Q(3, 4), Q(1))

# Twelve distinct cases spanning each requested grid family and all threshold
# values, without creating the excluded 4D-by-4-cells 10,000-rectangle model.
CASES = (
    {'dimension': 2, 'cells_per_axis': 6, 'grid': 'dyadic6', 'p': '4/5', 'theta': '1/2'},
    {'dimension': 2, 'cells_per_axis': 6, 'grid': 'dyadic6', 'p': '4/5', 'theta': '3/4'},
    {'dimension': 2, 'cells_per_axis': 6, 'grid': 'dyadic6', 'p': '9/10', 'theta': '1/2'},
    {'dimension': 2, 'cells_per_axis': 6, 'grid': 'dyadic6', 'p': '9/10', 'theta': '3/4'},
    {'dimension': 2, 'cells_per_axis': 6, 'grid': 'dyadic6', 'p': '19/20', 'theta': '1/2'},
    {'dimension': 2, 'cells_per_axis': 6, 'grid': 'dyadic6', 'p': '19/20', 'theta': '3/4'},
    {'dimension': 2, 'cells_per_axis': 8, 'grid': 'uniform', 'p': '4/5', 'theta': '1/2'},
    {'dimension': 2, 'cells_per_axis': 8, 'grid': 'uniform', 'p': '9/10', 'theta': '1/2'},
    {'dimension': 3, 'cells_per_axis': 3, 'grid': 'uniform', 'p': '4/5', 'theta': '1/2'},
    {'dimension': 3, 'cells_per_axis': 4, 'grid': 'uniform', 'p': '9/10', 'theta': '3/4'},
    {'dimension': 4, 'cells_per_axis': 3, 'grid': 'uniform', 'p': '4/5', 'theta': '1/2'},
    {'dimension': 4, 'cells_per_axis': 3, 'grid': 'uniform', 'p': '9/10', 'theta': '1/2'},
)


def _axis_breaks(m, kind):
    if kind == 'dyadic6':
        assert m == 6
        return DYADIC6_BREAKS
    assert kind == 'uniform'
    return tuple(Q(i, m) for i in range(m+1))


def build_geometry(spec):
    dim, m, kind = spec['dimension'], spec['cells_per_axis'], spec['grid']
    breaks = _axis_breaks(m, kind)
    widths = tuple(breaks[i+1]-breaks[i] for i in range(m))
    coordinates = tuple(itertools.product(range(m), repeat=dim))
    cell_index = {coord: i for i, coord in enumerate(coordinates)}
    cell_volumes = tuple(math.prod(widths[index] for index in coord)
                         for coord in coordinates)
    assert sum(cell_volumes, Q(0)) == 1

    # Proper torus arcs have a unique (origin, length), with length 1..m-1.
    # The full circle has m equivalent origins, so retain only origin 0.
    # This removes only identical cell sets; it does not change the selectable
    # geometric rectangles or any source/receiver constraint.
    axis_arcs = []
    for length in range(1, m):
        for origin in range(m):
            cells = tuple((origin+offset) % m for offset in range(length))
            measure = sum((widths[i] for i in cells), Q(0))
            axis_arcs.append({'origin': origin, 'length': length,
                              'cells': cells, 'measure': measure})
    cells = tuple(range(m))
    axis_arcs.append({'origin': 0, 'length': m, 'cells': cells,
                      'measure': sum(widths, Q(0))})

    rectangles = []
    for arc_tuple in itertools.product(axis_arcs, repeat=dim):
        axis_cells = [arc['cells'] for arc in arc_tuple]
        coords = tuple(itertools.product(*axis_cells))
        ids = tuple(cell_index[coord] for coord in coords)
        volume = math.prod(arc['measure'] for arc in arc_tuple)
        assert sum((cell_volumes[i] for i in ids), Q(0)) == volume
        rectangles.append({
            'axes': tuple((arc['origin'], arc['length']) for arc in arc_tuple),
            'cell_ids': ids,
            'volume': volume,
        })

    return {
        'breaks': breaks,
        'widths': widths,
        'coordinates': coordinates,
        'cell_volumes': cell_volumes,
        'rectangles': rectangles,
    }


def _full_axis_arc_geometry(spec):
    """Build the old catalog as a reference for exact compact-catalog replay."""
    dim, m, kind = spec['dimension'], spec['cells_per_axis'], spec['grid']
    breaks = _axis_breaks(m, kind)
    widths = tuple(breaks[i+1]-breaks[i] for i in range(m))
    coordinates = tuple(itertools.product(range(m), repeat=dim))
    cell_index = {coord: i for i, coord in enumerate(coordinates)}
    cell_volumes = tuple(math.prod(widths[index] for index in coord)
                         for coord in coordinates)
    axis_arcs = []
    for origin in range(m):
        for length in range(1, m+1):
            cells = tuple((origin+offset) % m for offset in range(length))
            axis_arcs.append({'origin': origin, 'length': length,
                              'cells': cells,
                              'measure': sum((widths[i] for i in cells), Q(0))})
    rectangles = []
    for arcs in itertools.product(axis_arcs, repeat=dim):
        ids = tuple(cell_index[x] for x in itertools.product(*(a['cells'] for a in arcs)))
        rectangles.append({'cell_ids': ids,
                           'volume': math.prod(a['measure'] for a in arcs)})
    return {'coordinates': coordinates, 'cell_volumes': cell_volumes,
            'rectangles': rectangles}


def _fraction_record(value):
    return [value.numerator, value.denominator]


def catalog_equivalence(input_path):
    """Compare exact set, volume, and candidate means against the old catalog."""
    old = json.loads(Path(input_path).read_text())
    groups = {}
    for row in old['cases']:
        if not row.get('candidate'):
            raise AssertionError('reference case has no exact candidate source')
        key = (row['dimension'], row['cells_per_axis'], row['grid'])
        groups.setdefault(key, []).append(row)
    expected = {(2, 6): 961, (2, 8): 3249, (3, 3): 343,
                (3, 4): 2197, (4, 3): 2401}
    summaries = []
    for (dim, m, kind), rows in sorted(groups.items()):
        spec = {'dimension': dim, 'cells_per_axis': m, 'grid': kind}
        compact = build_geometry(spec)
        full = _full_axis_arc_geometry(spec)

        def exact_map(geometry, candidate_rows, allow_duplicate_sets):
            out = {}
            duplicate_count = 0
            for rect in geometry['rectangles']:
                ids = tuple(sorted(set(rect['cell_ids'])))
                volume = sum((geometry['cell_volumes'][i] for i in ids), Q(0))
                means = []
                for row in candidate_rows:
                    cand = row['candidate']
                    den = cand['source_rounding_denominator']
                    nums = cand['source_cell_fraction_numerators']
                    mass = sum((geometry['cell_volumes'][i]*Q(nums[i], den)
                                for i in ids), Q(0))
                    means.append(mass/volume)
                value = (volume, tuple(means))
                if ids in out:
                    duplicate_count += 1
                    if not allow_duplicate_sets or out[ids] != value:
                        raise AssertionError('duplicate cell set has inconsistent exact data')
                else:
                    out[ids] = value
            return out, duplicate_count

        old_map, old_duplicates = exact_map(full, rows, True)
        new_map, new_duplicates = exact_map(compact, rows, False)
        assert len(new_map) == expected[(dim, m)]
        assert new_duplicates == 0
        assert old_map == new_map
        canonical = [
            [list(ids), _fraction_record(volume),
             [_fraction_record(mean) for mean in means]]
            for ids, (volume, means) in sorted(new_map.items())
        ]
        digest = hashlib.sha256(json.dumps(canonical, separators=(',', ':')).encode()).hexdigest()
        summaries.append({'dimension': dim, 'cells_per_axis': m, 'grid': kind,
                          'reference_cases': len(rows),
                          'full_catalog_rows': len(full['rectangles']),
                          'unique_catalog_rows': len(new_map),
                          'removed_duplicate_rows': old_duplicates,
                          'exact_set_volume_mean_hash': digest})
    assert len(summaries) == 5
    return {'pass': True, 'families': summaries}


def build_model(spec, geometry):
    from scipy.optimize import Bounds

    ncell = len(geometry['coordinates'])
    rectangles = geometry['rectangles']
    nrect = len(rectangles)
    x0, w0, z0 = 0, ncell, ncell+nrect
    nvar = ncell+nrect+ncell
    threshold = Q(spec['p'])+EPS

    rr, cc, dd = [], [], []
    lower, upper = [], []

    # One strict-target source-density constraint per selectable torus box:
    # sum_j vol_j x_j - (p+eps)|R| w_R >= 0.
    for rid, rectangle in enumerate(rectangles):
        row = len(lower)
        for cell_id in rectangle['cell_ids']:
            rr.append(row)
            cc.append(x0+cell_id)
            dd.append(float(geometry['cell_volumes'][cell_id]))
        rr.append(row)
        cc.append(w0+rid)
        dd.append(-float(threshold*rectangle['volume']))
        lower.append(0.0)
        upper.append(np.inf)

    rectangles_by_cell = [[] for _ in range(ncell)]
    for rid, rectangle in enumerate(rectangles):
        for cell_id in rectangle['cell_ids']:
            rectangles_by_cell[cell_id].append(rid)

    # z_i <= sum_{R containing the whole receiver cell i} w_R.
    for cell_id, covering in enumerate(rectangles_by_cell):
        row = len(lower)
        rr.append(row)
        cc.append(z0+cell_id)
        dd.append(1.0)
        for rid in covering:
            rr.append(row)
            cc.append(w0+rid)
            dd.append(-1.0)
        lower.append(-np.inf)
        upper.append(0.0)

    # At least theta of receiver cells, measured by exact cell volumes.
    row = len(lower)
    for cell_id, volume in enumerate(geometry['cell_volumes']):
        rr.append(row)
        cc.append(z0+cell_id)
        dd.append(float(volume))
    lower.append(float(Q(spec['theta'])))
    upper.append(np.inf)

    matrix = coo_matrix((np.asarray(dd, dtype=np.float64),
                         (np.asarray(rr, dtype=np.int32), np.asarray(cc, dtype=np.int32))),
                        shape=(len(lower), nvar)).tocsr()
    objective = np.zeros(nvar, dtype=np.float64)
    objective[x0:w0] = np.asarray([float(v) for v in geometry['cell_volumes']])
    integrality = np.zeros(nvar, dtype=np.uint8)
    integrality[w0:z0] = 1
    lb = np.zeros(nvar, dtype=np.float64)
    ub = np.ones(nvar, dtype=np.float64)
    constraint = LinearConstraint(matrix, np.asarray(lower), np.asarray(upper))
    model_bytes = (matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes+
                   objective.nbytes+integrality.nbytes+lb.nbytes+ub.nbytes+
                   np.asarray(lower).nbytes+np.asarray(upper).nbytes)
    assert model_bytes < MAX_MODEL_BYTES
    return {
        'constraint': constraint,
        'matrix': matrix,
        'objective': objective,
        'integrality': integrality,
        'bounds': Bounds(lb, ub),
        'var_slices': (x0, w0, z0, nvar),
        'n_source': ncell,
        'n_rectangles': nrect,
        'n_receiver': ncell,
        'n_constraints': len(lower),
        'nnz': matrix.nnz,
        'sparse_model_bytes': model_bytes,
    }


def _ceil_source_fraction(value):
    if not math.isfinite(float(value)) or value < -1e-7 or value > 1+1e-7:
        raise ArithmeticError(f'solver source variable outside [0,1]: {value}')
    clipped = min(1.0, max(0.0, float(value)))
    numerator = min(SOURCE_DENOMINATOR,
                    max(0, math.ceil(clipped*SOURCE_DENOMINATOR)))
    return Q(numerator, SOURCE_DENOMINATOR), numerator


def _candidate_certificate(spec, geometry, model, solution):
    ncell = model['n_source']
    nrect = model['n_rectangles']
    x0, w0, z0, _ = model['var_slices']
    vector = solution.x
    source_fracs, source_numerators = [], []
    for value in vector[x0:w0]:
        frac, num = _ceil_source_fraction(value)
        source_fracs.append(frac)
        source_numerators.append(num)

    raw_w = vector[w0:z0]
    binary_errors = [abs(float(value)-round(float(value))) for value in raw_w]
    binary_ok = max(binary_errors, default=0.0) <= 1e-5
    raw_selected = [rid for rid, value in enumerate(raw_w) if value >= 0.5]

    p = Q(spec['p'])
    certified_rectangles = []
    failed_rectangles = []
    covered = set()
    for rid in raw_selected:
        rect = geometry['rectangles'][rid]
        numerator = sum((geometry['cell_volumes'][j]*source_fracs[j]
                         for j in rect['cell_ids']), Q(0))
        mean = numerator/rect['volume']
        if mean > p:
            certified_rectangles.append({
                'id': rid,
                'axis_origin_length_cells': [list(pair) for pair in rect['axes']],
                'mean_fraction': str(mean),
                'mean_minus_p_fraction': str(mean-p),
                'receiver_cell_ids': list(rect['cell_ids']),
            })
            covered.update(rect['cell_ids'])
        else:
            failed_rectangles.append(rid)

    density = sum((volume*source_fracs[i]
                   for i, volume in enumerate(geometry['cell_volumes'])), Q(0))
    receiver_volume = sum((geometry['cell_volumes'][i] for i in covered), Q(0))
    theta = Q(spec['theta'])
    cert_pass = bool(binary_ok and certified_rectangles and
                     receiver_volume >= theta and density > 0)
    solver_z_volume = sum(float(v)*float(vector[z0+i])
                          for i, v in enumerate(geometry['cell_volumes']))

    candidate = {
        'certificate_pass': cert_pass,
        'source_rounding_denominator': SOURCE_DENOMINATOR,
        'source_cell_coordinates_lexicographic': [list(c) for c in geometry['coordinates']],
        'source_cell_fraction_numerators': source_numerators,
        'source_realization': ('In each cell occupy a half-open slab in coordinate 1 '
                               'of relative length numerator/2^24; this realizes its exact '
                               'cell-average density.'),
        'source_density_fraction': str(density),
        'source_density_numerator': density.numerator,
        'source_density_denominator': density.denominator,
        'solver_receiver_volume_numeric': solver_z_volume,
        'solver_positive_w_count': len(raw_selected),
        'solver_binary_integrality_max_error': max(binary_errors, default=0.0),
        'failed_solver_selected_rectangle_ids': failed_rectangles,
        'certified_selected_rectangles': certified_rectangles,
        'receiver_cell_ids_exact_union': sorted(covered),
        'theta_required_fraction': str(theta),
        'theta_actual_fraction': str(receiver_volume),
        'theta_actual_numerator': receiver_volume.numerator,
        'theta_actual_denominator': receiver_volume.denominator,
        'theta_actual_over_density_fraction': str(receiver_volume/density) if density else None,
        'theta_actual_over_density_numerator': (receiver_volume/density).numerator if density else None,
        'theta_actual_over_density_denominator': (receiver_volume/density).denominator if density else None,
        'all_final_rectangles_strictly_above_p': all(
            Q(row['mean_fraction']) > p for row in certified_rectangles),
    }
    return candidate


def solve_case(spec):
    started = time.monotonic()
    try:
        geometry = build_geometry(spec)
        model = build_model(spec, geometry)
        options = {'presolve': True, 'time_limit': TIME_LIMIT,
                   'mip_rel_gap': 0.0, 'threads': SOLVER_THREADS}
        # SciPy exposes time_limit/presolve/mip_rel_gap; its local wrapper
        # forwards the HiGHS-specific `threads` option verbatim.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='Unrecognized options detected:.*threads.*',
                                    category=RuntimeWarning)
            result = milp(
                c=model['objective'], integrality=model['integrality'],
                bounds=model['bounds'], constraints=model['constraint'],
                options=options)

        has_incumbent = getattr(result, 'x', None) is not None
        candidate = (_candidate_certificate(spec, geometry, model, result)
                     if has_incumbent else None)
        axes_breaks = [str(x) for x in geometry['breaks']]
        return {
            **spec,
            'epsilon_generation_fraction': str(EPS),
            'rectangle_count': model['n_rectangles'],
            'source_cell_count': model['n_source'],
            'receiver_cell_count': model['n_receiver'],
            'variable_count': len(model['objective']),
            'constraint_count': model['n_constraints'],
            'constraint_nnz': model['nnz'],
            'sparse_model_bytes': model['sparse_model_bytes'],
            'axis_breaks': axes_breaks,
            'solver': 'scipy.optimize.milp / HiGHS',
            'scipy_version': scipy.__version__,
            'solver_options': options,
            'solver_threads_option_forwarded': True,
            'time_limit_seconds': TIME_LIMIT,
            'status': int(result.status),
            'success': bool(result.success),
            'message': str(result.message),
            'has_incumbent': has_incumbent,
            'objective_numeric': (float(result.fun) if getattr(result, 'fun', None) is not None else None),
            'mip_gap_numeric': (float(result.mip_gap) if getattr(result, 'mip_gap', None) is not None else None),
            'dual_bound_numeric_only': (float(result.mip_dual_bound)
                                        if getattr(result, 'mip_dual_bound', None) is not None else None),
            'node_count': (int(result.mip_node_count)
                           if getattr(result, 'mip_node_count', None) is not None else None),
            'solver_seconds': time.monotonic()-started,
            'candidate': candidate,
            'pid': os.getpid(),
        }
    except Exception as exc:
        return {
            **spec,
            'status': None,
            'success': False,
            'message': f'{type(exc).__name__}: {exc}',
            'has_incumbent': False,
            'objective_numeric': None,
            'mip_gap_numeric': None,
            'dual_bound_numeric_only': None,
            'node_count': None,
            'candidate': None,
            'solver_seconds': time.monotonic()-started,
            'pid': os.getpid(),
        }


def run(output_dir, pilot=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = CASES[:1] if pilot else CASES
    start = time.monotonic()
    with mp.get_context('spawn').Pool(WORKERS) as pool:
        rows = pool.map(solve_case, cases)
    pids = sorted({row.pop('pid') for row in rows})
    result = {
        'schema': 'periodic-torus-grid-milp-candidate-v2',
        'pass_all_candidates': all(row['candidate'] is not None and
                                   row['candidate']['certificate_pass'] for row in rows),
        'configured_workers': WORKERS,
        'actual_worker_pids': pids,
        'actual_workers': len(pids),
        'solver_threads_per_worker': SOLVER_THREADS,
        'case_count': len(rows),
        'pilot_mode': pilot,
        'source_rounding_denominator': SOURCE_DENOMINATOR,
        'threshold_epsilon_generation_only': str(EPS),
        'all_candidates_certified_at_p_not_p_plus_epsilon': True,
        'scipy_version': scipy.__version__,
        'cases': rows,
        'seconds': time.monotonic()-start,
        'scope': ('Finite torus-grid candidate search for the core95 periodic source-density proxy. '
                  'Each certified theta_actual/d is a lower bound witnessed by the listed periodic '
                  'indicator candidate. It is not a lower certificate for continuous d_N or a general '
                  'H_N bound; finite-box transference would require a separate boundary-tolerance argument.'),
        'solver_dual_bound_is_numerical_only': True,
        'general_H_N_bound_claimed': False,
        'continuous_d_N_lower_bound_claimed': False,
    }
    (output_dir/'periodic99.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key != 'cases'}, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_directory')
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--catalog-audit', metavar='PERIODIC96_JSON',
                        help='exactly compare compact and full catalogs without solving MILPs')
    args = parser.parse_args()
    if args.catalog_audit:
        result = catalog_equivalence(args.catalog_audit)
        print(json.dumps(result, indent=2))
        return
    run(args.output_directory, pilot=args.pilot)


if __name__ == '__main__':
    main()
