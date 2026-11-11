"""Finite periodic torus-grid MILP search with exact candidate certificates.

This is a finite-grid search for candidates in the periodic source-density
proxy. HiGHS bounds are numerical only. Every reported witness candidate is
rechecked with rational arithmetic after rounding source cell fractions upward
to a 2**24 denominator.
"""
import argparse
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
TIME_LIMIT = 900
MAX_MODEL_BYTES = 100*2**20
WORKERS = 12
DYADIC6_BREAKS = (Q(0), Q(1, 16), Q(1, 8), Q(1, 4), Q(1, 2), Q(3, 4), Q(1))

# Twelve distinct cases spanning each requested grid family and all threshold
# values, without creating the excluded 4D-by-4-cells 10,000-rectangle model.
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

    # Each torus interval is indexed by its start cell and number of wrapped
    # cells, 1..m. Its set consists of whole partition cells.
    axis_arcs = []
    for origin in range(m):
        for length in range(1, m+1):
            cells = tuple((origin+offset) % m for offset in range(length))
            measure = sum((widths[i] for i in cells), Q(0))
            axis_arcs.append({'origin': origin, 'length': length,
                              'cells': cells, 'measure': measure})

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


def _atomic_json(path, payload):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def _case_result_path(output_dir, row):
    safe = lambda value: str(value).replace('/', '-').replace('.', '-')
    return Path(output_dir) / (
        f"case_d{row['dimension']}_m{row['cells_per_axis']}_{row['grid']}"
        f"_p{safe(row['p'])}_theta{safe(row['theta'])}.json")


def _solve_case_persist(args):
    spec, output_dir = args
    result = solve_case(spec)
    cert = result.get('candidate')
    _atomic_json(_case_result_path(output_dir, result), {
        'schema': 'periodic-torus-grid-milp-case-v1',
        'case': result,
        'candidate_exact_certificate_present': cert is not None,
        'candidate_exact_certificate_pass': bool(cert and cert.get('certificate_pass')),
        'solver_dual_bound_is_numerical_only': True,
    })
    return result


def _record_completed_case(output_dir, rows, result, expected_cases):
    """Persist one solver result before exposing it to the aggregate checkpoint."""
    case_path = _case_result_path(output_dir, result)
    cert = result.get('candidate')
    _atomic_json(case_path, {
        'schema': 'periodic-torus-grid-milp-case-v1',
        'case': result,
        'candidate_exact_certificate_present': cert is not None,
        'candidate_exact_certificate_pass': bool(cert and cert.get('certificate_pass')),
        'solver_dual_bound_is_numerical_only': True,
    })
    checkpoint_path = Path(output_dir) / 'periodic_search97.checkpoint.json'
    checkpoint = {
        'schema': 'periodic-torus-grid-milp-candidate-v1-checkpoint',
        'completed_cases': len(rows), 'expected_cases': expected_cases,
        'case_count': expected_cases, 'configured_workers': WORKERS,
        'solver_threads_per_worker': SOLVER_THREADS,
        'time_limit_seconds': TIME_LIMIT,
        'completed_rows': rows,
        'case_files': [str(_case_result_path(output_dir, row)) for row in rows],
        'completed_at_unix': time.time(),
        'resource_note': 'Per-case solver/certificate results are persisted before aggregate checkpoint; sampled run resource measurements are in math_guard report.',
    }
    _atomic_json(checkpoint_path, checkpoint)
    log_record = {
        'checkpoint': len(rows), 'total': expected_cases,
        'case_file': str(case_path),
        'spec': {k: result[k] for k in ('dimension','cells_per_axis','grid','p','theta')},
        'status': result.get('status'), 'has_incumbent': result.get('has_incumbent'),
        'certificate_pass': bool(cert and cert.get('certificate_pass')),
        'solver_seconds': result.get('solver_seconds'),
    }
    print(json.dumps(log_record), flush=True)


def _aggregate_result(cases, rows, pilot, elapsed):
    keys = [tuple(row.get(k) for k in ('dimension','cells_per_axis','grid','p','theta'))
            for row in rows]
    if len(rows) != len(cases) or len(set(keys)) != len(cases):
        raise RuntimeError(f'case aggregate mismatch: rows={len(rows)} cases={len(cases)} unique={len(set(keys))}')
    pids = sorted({row.get('pid') for row in rows if row.get('pid') is not None})
    return {
        'schema': 'periodic-torus-grid-milp-candidate-v1',
        'pass_all_candidates': all(row['candidate'] is not None and
                                   row['candidate']['certificate_pass'] for row in rows),
        'configured_workers': WORKERS,
        'actual_worker_pids': pids,
        'actual_workers': len(pids),
        'solver_threads_per_worker': SOLVER_THREADS,
        'case_count': len(rows),
        'pilot_mode': pilot,
        'time_limit_seconds': TIME_LIMIT,
        'source_rounding_denominator': SOURCE_DENOMINATOR,
        'threshold_epsilon_generation_only': str(EPS),
        'all_candidates_certified_at_p_not_p_plus_epsilon': True,
        'scipy_version': scipy.__version__,
        'cases': rows,
        'seconds': elapsed,
        'scope': ('Finite torus-grid candidate search for the periodic source-density proxy. '
                  'Each certified theta_actual/d is a lower bound witnessed by the listed periodic '
                  'indicator candidate. It is not a lower certificate for continuous d_N or a general '
                  'H_N bound; finite-box transference would require a separate boundary-tolerance argument.'),
        'solver_dual_bound_is_numerical_only': True,
        'general_H_N_bound_claimed': False,
        'continuous_d_N_lower_bound_claimed': False,
    }


def run(output_dir, pilot=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = CASES[:1] if pilot else CASES
    start = time.monotonic()
    checkpoint_path = output_dir/'periodic_search97.checkpoint.json'
    rows_by_key = {}
    if checkpoint_path.exists():
        old = json.loads(checkpoint_path.read_text())
        for row in old.get('completed_rows', []):
            key = tuple(row.get(k) for k in ('dimension','cells_per_axis','grid','p','theta'))
            rows_by_key[key] = row
    # A case file is written before the aggregate checkpoint. Recover a file
    # that was committed just before an interrupted parent could update it.
    for case_path in output_dir.glob('case_d*_m*_*.json'):
        payload = json.loads(case_path.read_text())
        row = payload.get('case', {})
        if not row:
            continue
        key = tuple(row.get(k) for k in ('dimension','cells_per_axis','grid','p','theta'))
        if key in {tuple(spec.get(k) for k in ('dimension','cells_per_axis','grid','p','theta')) for spec in cases}:
            rows_by_key[key] = row
    rows = list(rows_by_key.values())
    completed_keys = set(rows_by_key)
    pending_cases = [spec for spec in cases
                     if tuple(spec.get(k) for k in ('dimension','cells_per_axis','grid','p','theta'))
                     not in completed_keys]
    if pending_cases:
        with mp.get_context('spawn').Pool(WORKERS) as pool:
            for result in pool.imap_unordered(_solve_case_persist, [(x,str(output_dir)) for x in pending_cases], chunksize=1):
                key = tuple(result.get(k) for k in ('dimension','cells_per_axis','grid','p','theta'))
                if key in rows_by_key:
                    raise RuntimeError(f'duplicate case result: {key}')
                rows_by_key[key] = result
                rows.append(result)
                _record_completed_case(output_dir, rows, result, len(cases))
    rows.sort(key=lambda row: next(i for i, spec in enumerate(cases)
                                   if all(row.get(key) == value for key, value in spec.items())))
    result = _aggregate_result(cases, rows, pilot, time.monotonic()-start)
    pids = result['actual_worker_pids']
    _atomic_json(output_dir/'periodic_search97.json', result)
    _atomic_json(checkpoint_path, {
        'schema': 'periodic-torus-grid-milp-candidate-v1-checkpoint',
        'completed_cases': len(rows), 'expected_cases': len(cases),
        'case_count': len(rows), 'configured_workers': WORKERS,
        'actual_worker_pids': pids, 'actual_workers': len(pids),
        'solver_threads_per_worker': SOLVER_THREADS,
        'time_limit_seconds': TIME_LIMIT, 'final_result_path': str(output_dir/'periodic_search97.json'),
        'completed_rows': rows, 'completed_at_unix': time.time(),
    })
    print(json.dumps({key: value for key, value in result.items()
                      if key != 'cases'}, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_directory')
    parser.add_argument('--pilot', action='store_true')
    args = parser.parse_args()
    run(args.output_directory, pilot=args.pilot)


if __name__ == '__main__':
    main()
