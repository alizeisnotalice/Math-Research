#!/usr/bin/env python3
"""Round03: one 26-direction LP. Importing performs no numerical experiment."""
from __future__ import annotations
from time import monotonic
_IMPORT_START = monotonic()
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
import scipy
import scipy.sparse as sp
import highspy
_IMPORT_SECONDS = monotonic()-_IMPORT_START
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
ALPHA_Q, R_Q, H_Q = Q(47, 40), Q(33, 20), Q(1, 32)
ALPHA, RADIUS = float(ALPHA_Q), float(R_Q)
TOTAL_SECONDS, LP_SOFT_SECONDS = 120., 12.
K_TOL, A_TOL, GAIN_TOL = 1e-9, 1e-8, 5e-7
LOCAL_CLASSES = [(i, j) for i in (48, 49, 50, 51) for j in (3, 4, 5, 6)]+[
    (20, 12), (24, 8), (28, 4), (31, 1), (20, 13), (24, 9), (28, 5), (32, 1)]
LOWER = np.asarray([-1/20000, -1/1000]+[-1/4096]*24)
UPPER = -LOWER


class StopBudget(Exception):
    pass


class ModelRejected(Exception):
    pass


class BaselineRejected(Exception):
    pass


def dump(path, value):
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)


def finite(name, values):
    values = np.asarray(values)
    if not np.all(np.isfinite(values)):
        raise ModelRejected(name+'_nonfinite')
    return values


def read_helper():
    for path in (ROOT/'dependency_round02_search.py', ROOT.parent/'round02/search.py'):
        if path.is_file():
            spec = importlib.util.spec_from_file_location('round03_round02_definitions', path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            return module, path
    raise FileNotFoundError('Reviewed round02 definition dependency unavailable.')


def input_file(local, historical):
    a, b = ROOT/local, ROOT.parent/historical
    return a if a.is_file() else b


def orbit(i, j):
    return sorted({(x, y) for x in (i, -i) for y in (j, -j)}|
                  {(y, x) for x in (i, -i) for y in (j, -j)})


def full_coefficient_matrix(coefficients, N):
    out = np.zeros((2*N+1, 2*N+1))
    for key, value in coefficients.items():
        i, j = map(int, key.split(','))
        for x, y in orbit(i, j):
            out[x+N, y+N] = float(value)
    return out


def main():
    result_path = ROOT/'result.json'
    if result_path.exists():
        raise SystemExit('Existing result: the one finite round03 run must not be repeated.')
    start = monotonic()
    counts = {'LP_started': 0, 'LP_completed': 0,
              'baseline_K_started': 0, 'baseline_K_completed': 0,
              'baseline_A_started': 0, 'baseline_A_completed': 0,
              'candidate_K_started': 0, 'candidate_K_completed': 0,
              'candidate_A_started': 0, 'candidate_A_completed': 0,
              'K_direction_rows_built': 0, 'A_direction_rows_built': 0}
    result = {'status': 'started', 'phase': 'prepare', 'counts': counts, 'candidate_exported': False,
              'limits': {'total_seconds_internal': TOTAL_SECONDS,
                         'LP_time_limit_soft_seconds': LP_SOFT_SECONDS, 'LP_calls': 1,
                         'point_candidate_evaluations': 62414, 'rationalized_candidates': 1},
              'outer_deadline': 'The human-authorized external wrapper must SIGKILL the process group at 120s.'}
    dump(result_path, result)

    def save():
        result['elapsed_seconds'] = monotonic()-start
        dump(result_path, result)

    def guard():
        if monotonic()-start >= TOTAL_SECONDS:
            raise StopBudget('internal_120s_deadline')

    try:
        if len(LOCAL_CLASSES) != 24 or len(set(LOCAL_CLASSES)) != 24:
            raise ValueError('Fixed 24-class list changed.')
        helper, helper_path = read_helper()
        source = input_file('input_baseline.json', 'round01/scaled_certificate.json')
        arithmetic_source = input_file('input_scaling_arithmetic.json', 'round01/scaling_arithmetic.json')
        basis = helper.Basis(source)
        if Q(basis.spec['alpha']) != ALPHA_Q or Q(basis.spec['R']) != R_Q or basis.hq != H_Q:
            raise ValueError('Only the reviewed fixed alpha/R/h baseline is approved.')
        baseline_c = Q(basis.spec['c'])
        baseline_coefficients = {k: Q(v) for k, v in basis.spec['a'].items()}
        baseline_matrix = full_coefficient_matrix(baseline_coefficients, basis.N)
        arithmetic = json.loads(arithmetic_source.read_text())
        if arithmetic['new_certificate_sha256'] != sha256(source.read_bytes()).hexdigest():
            raise ValueError('Baseline cost/certificate mismatch.')
        baseline_cost_upper, baseline_cost_lower = Q(arithmetic['new_cost_upper']), Q(arithmetic['new_cost_lower'])
        bspl, singular = helper.load_formulas()
        formula_dir = helper.formula_directory()
        pools, pool_manifest = helper.make_pools(str(R_Q))
        Kpoints = np.unique(np.vstack((pools['initial_K'], pools['validation_K'])), axis=0)
        Apoints = np.unique(np.vstack((pools['initial_A'], pools['validation_A'])), axis=0)
        if len(Kpoints) != 15816 or len(Apoints) != 14693:
            raise ModelRejected('Reviewed deterministic-union counts do not match.')
        if len(Kpoints) > 15849 or len(Apoints) > 15358 or len(Kpoints)+len(Apoints) > 31207:
            raise ModelRejected('Fixed pool caps exceeded; no truncation/refill permitted.')
        if np.any(Kpoints.sum(axis=1) <= 0):
            raise ModelRejected('Origin in K pool.')
        if np.any(Apoints[:, 1] <= 0) or np.any(abs(Apoints.sum(axis=1)-RADIUS) < 1e-12):
            raise ModelRejected('Axis or truncation kink in A pool.')
        target_K = np.where(Kpoints.sum(axis=1) <= 1., 1.+1e-7, 0.)
        target_A = np.zeros(len(Apoints))
        np.savez(ROOT/'point_pools.npz', K=Kpoints, A=Apoints)
        result.update(variable_order=['t0_all_K1', 't1_spline_only']+[f'd_{i}_{j}' for i, j in LOCAL_CLASSES],
                      variable_bounds=list(zip(LOWER.tolist(), UPPER.tolist())),
                      pool_counts={'K': len(Kpoints), 'A': len(Apoints), 'round02_source_pools': pool_manifest},
                      dependencies={'round02_definition_sha256': sha256(helper_path.read_bytes()).hexdigest(),
                                    'baseline_certificate_sha256': sha256(source.read_bytes()).hexdigest(),
                                    'baseline_arithmetic_sha256': sha256(arithmetic_source.read_bytes()).hexdigest(),
                                    'formula_hashes': {name: sha256((formula_dir/name).read_bytes()).hexdigest()
                                     for name in ('bspl_closed.py', 'sing2.py', 'sing3.py')}},
                      baseline_cost_interval={'lower': str(baseline_cost_lower), 'upper': str(baseline_cost_upper)},
                      fixed_parameters={'alpha': str(ALPHA_Q), 'R': str(R_Q), 'h': str(H_Q), 'N': basis.N},
                      point_pools_sha256=sha256((ROOT/'point_pools.npz').read_bytes()).hexdigest())
        guard(); save()
        local_orbits = [orbit(i, j) for i, j in LOCAL_CLASSES]
        spline_mass_q = sum((H_Q**2*len(orbit(*map(int, k.split(','))))*value/2
                             for k, value in baseline_coefficients.items()), Q(0))
        # All objective columns are HALF integrals, with the affine B1 constant.
        numeric_baseline_cost = float(baseline_c)*ALPHA*RADIUS**(2-ALPHA)/(2-ALPHA)+float(spline_mass_q)
        objective = finite('objective', [numeric_baseline_cost, float(spline_mass_q)]+
                           [float(H_Q**2*len(o)/2) for o in local_orbits])
        result['mass_model'] = {'constant_B1_display': numeric_baseline_cost,
                                'spline_half_mass_exact': str(spline_mass_q),
                                'direction_half_integrals': objective.tolist()}
        all_indices = np.arange(-basis.N, basis.N+1)
        cached_Apsi = None

        def evaluate(kind, coefficient_matrix, cvalue, directions=False, Apsi_cache=None):
            nonlocal cached_Apsi
            points = Kpoints if kind.endswith('_K') else Apoints
            is_A = kind.endswith('_A')
            counts[kind+'_started'] += len(points)
            result['phase'] = kind; save(); guard()
            values, columns, singular_values = [], [], []
            for first in range(0, len(points), 256):
                guard()
                p = points[first:first+256]; r = p.sum(axis=1)
                Bu = finite('Bu', bspl.Bv(p[:, 0, None]/basis.h-all_indices[None, :]))
                Bv = finite('Bv', bspl.Bv(p[:, 1, None]/basis.h-all_indices[None, :]))
                BA = Bu@coefficient_matrix
                spline_K = np.sum(BA*Bv, axis=1)
                if is_A:
                    Lu = finite('Lu', bspl.LBv(p[:, 0, None]/basis.h-all_indices[None, :], ALPHA)*basis.h**(-ALPHA))
                    Lv = finite('Lv', bspl.LBv(p[:, 1, None]/basis.h-all_indices[None, :], ALPHA)*basis.h**(-ALPHA))
                    spline_value = np.sum((Lu@coefficient_matrix)*Bv, axis=1)+np.sum(BA*Lv, axis=1)
                    aps = finite('true_Apsi', singular.A_psi_vec(ALPHA, RADIUS, p[:, 0], p[:, 1])
                                 if Apsi_cache is None else Apsi_cache[first:first+len(p)])
                    singular_values.append(aps)
                    full_value = finite('true_AK', spline_value+float(cvalue)*aps)
                else:
                    spline_value = spline_K
                    psi = finite('Kpsi', np.maximum(r**(-ALPHA)-RADIUS**(-ALPHA), 0.))
                    full_value = finite('true_K', spline_value+float(cvalue)*psi)
                values.append(full_value)
                if directions:
                    cols = np.empty((len(p), 26))
                    cols[:, 0], cols[:, 1] = full_value, spline_value
                    for j, orb in enumerate(local_orbits):
                        column = np.zeros(len(p))
                        for i, k in orb:
                            if is_A:
                                column += Lu[:, i+basis.N]*Bv[:, k+basis.N]+Bu[:, i+basis.N]*Lv[:, k+basis.N]
                            else:
                                column += Bu[:, i+basis.N]*Bv[:, k+basis.N]
                        cols[:, j+2] = column
                    columns.append(finite('direction_columns', cols))
                    counts['A_direction_rows_built' if is_A else 'K_direction_rows_built'] += len(p)
                guard(); counts[kind+'_completed'] += len(p); save()
            if is_A and Apsi_cache is None:
                cached_Apsi = np.concatenate(singular_values)
            return np.concatenate(values), np.vstack(columns) if directions else None

        base_K, Kcolumns = evaluate('baseline_K', baseline_matrix, baseline_c, True)
        base_A, Acolumns = evaluate('baseline_A', baseline_matrix, baseline_c, True)
        np.savez(ROOT/'baseline_rows.npz', K_values=base_K, A_values=base_A, K_columns=Kcolumns,
                 A_columns=Acolumns, Apsi=cached_Apsi, K_targets=target_K, A_targets=target_A)

        def residual_record(Kvalues, Avalues):
            vk, va = target_K-Kvalues, target_A-Avalues
            bad_K, bad_A = np.flatnonzero(vk > K_TOL), np.flatnonzero(va > A_TOL)
            return {'K_bad_indices': bad_K.tolist(), 'A_bad_indices': bad_A.tolist(),
                    'K_max_violation': float(vk.max()), 'A_max_violation': float(va.max()),
                    'minimum_K': float(Kvalues.min()),
                    'minimum_K_on_D': float(Kvalues[Kpoints.sum(axis=1) <= 1].min()),
                    'minimum_true_A': float(Avalues.min()), 'K_tolerance': K_TOL, 'A_tolerance': A_TOL,
                    'passed': len(bad_K) == 0 and len(bad_A) == 0}

        result['baseline_residual'] = residual_record(base_K, base_A); save(); guard()
        if not result['baseline_residual']['passed']:
            result['baseline_failure_points'] = {
                'K': [{'index': i, 'point': Kpoints[i].tolist(), 'value': float(base_K[i]), 'target': float(target_K[i])}
                      for i in result['baseline_residual']['K_bad_indices']],
                'A': [{'index': i, 'point': Apoints[i].tolist(), 'value': float(base_A[i])}
                      for i in result['baseline_residual']['A_bad_indices']]}
            raise BaselineRejected('Unmodified baseline failed fixed numerical residual check; do not relax or repair.')
        # True AK only: no near-origin reduction, cap, or positive A margin.
        matrix = sp.csr_matrix(np.vstack((Kcolumns, Acolumns)))
        lower = np.r_[target_K-base_K, target_A-base_A]
        scales = 1./np.maximum(1., abs(matrix).max(axis=1).toarray().ravel())
        matrix = (sp.diags(scales)@matrix).tocsr(); lower = lower*scales
        finite('scaled_LP_rows', matrix.data); finite('scaled_LP_rhs', lower)
        solver = highspy.Highs()
        for key, value in (('output_flag', False), ('threads', 1), ('solver', 'simplex'),
                           ('simplex_strategy', 1), ('primal_feasibility_tolerance', 1e-10),
                           ('dual_feasibility_tolerance', 1e-10)):
            solver.setOptionValue(key, value)
        solver.addVars(26, LOWER, UPPER)
        solver.changeColsCost(26, np.arange(26, dtype=np.int32), objective)
        solver.addRows(matrix.shape[0], lower, np.full(matrix.shape[0], highspy.kHighsInf), matrix.nnz,
                       matrix.indptr[:-1].astype(np.int32), matrix.indices.astype(np.int32), matrix.data)
        status = solver.setSolution(26, np.arange(26, dtype=np.int32), np.zeros(26))
        result['zero_primal_start_accepted'] = status == highspy.HighsStatus.kOk
        if not result['zero_primal_start_accepted']:
            raise ModelRejected('Checked zero primal start was not accepted; no fallback solve.')
        guard()
        solver.setOptionValue('time_limit', min(LP_SOFT_SECONDS, TOTAL_SECONDS-(monotonic()-start)))
        result['phase'] = 'single_LP'; counts['LP_started'] += 1; save()
        solve_start = monotonic()
        solver.run()  # The only HiGHS.run, never retried.
        solve_wall = monotonic()-solve_start
        counts['LP_completed'] += 1
        model_status = solver.getModelStatus(); sol = solver.getSolution()
        infeasibility = float(solver.getInfo().max_primal_infeasibility)
        result['LP'] = {'status': solver.modelStatusToString(model_status),
                        'optimal': model_status == highspy.HighsModelStatus.kOptimal,
                        'valid_primal': bool(sol.value_valid), 'actual_call_wall_seconds': solve_wall,
                        'time_limit_soft_seconds': LP_SOFT_SECONDS,
                        'max_scaled_primal_infeasibility': infeasibility if math.isfinite(infeasibility) else None}
        save(); guard()
        if not result['LP']['optimal'] or not result['LP']['valid_primal']:
            raise ModelRejected('Only optimal with valid primal may enter the unique candidate exit.')
        t = finite('26_direction_solution', sol.col_value)
        result['floating_direction_variables'] = t.tolist()
        bound_error = max(float(np.max(LOWER-t)), float(np.max(t-UPPER)), 0.)
        result['variable_bound_violation'] = bound_error
        if bound_error > 1e-10:
            raise ModelRejected('Direction variable bounds violated.')
        candidate_c = (1.+t[0])*float(baseline_c)
        candidate_coefficients = {k: (1.+t[0]+t[1])*float(value) for k, value in baseline_coefficients.items()}
        for delta, (i, j) in zip(t[2:], LOCAL_CLASSES):
            key = f'{i},{j}'; candidate_coefficients[key] = candidate_coefficients.get(key, 0.)+float(delta)
        candidate_matrix = full_coefficient_matrix(candidate_coefficients, basis.N)
        Kcandidate, _ = evaluate('candidate_K', candidate_matrix, candidate_c)
        Acandidate, _ = evaluate('candidate_A', candidate_matrix, candidate_c, Apsi_cache=cached_Apsi)
        result['candidate_residual'] = residual_record(Kcandidate, Acandidate)
        result['affine_vs_full_evaluation'] = {
            'K_max_absolute_difference': float(np.max(abs(Kcandidate-(base_K+Kcolumns@t)))),
            'A_max_absolute_difference': float(np.max(abs(Acandidate-(base_A+Acolumns@t))))}
        np.savez(ROOT/'candidate_evaluations.npz', K_values=Kcandidate, A_values=Acandidate)
        dump(ROOT/'floating_coefficients.json', {'c': candidate_c, 'a': candidate_coefficients, 't': t.tolist()})
        # Recompute the cost from the expanded full coefficients, not only LP objective.
        candidate_spline_mass = sum(len(orbit(*map(int, k.split(','))))*value*float(H_Q)**2/2
                                    for k, value in candidate_coefficients.items())
        numeric_cost = candidate_c*ALPHA*RADIUS**(2-ALPHA)/(2-ALPHA)+candidate_spline_mass
        affine_cost = numeric_baseline_cost+float(objective@t)
        finite('expanded_candidate_cost', [numeric_cost, affine_cost])
        result.update(candidate_float_half_cost=numeric_cost, affine_float_half_cost=affine_cost,
                      cost_affine_vs_expanded_difference=abs(numeric_cost-affine_cost),
                      numerical_gain_vs_baseline_display=float(baseline_cost_upper)-numeric_cost)
        save(); guard()
        if not result['candidate_residual']['passed']:
            raise ModelRejected('Full raw K/true-AK residual failed; no repair or replacement.')
        if numeric_cost < float(baseline_cost_upper)-GAIN_TOL:
            denominator = 2**50
            c_q = Q(int(round(candidate_c*denominator)), denominator)
            a_q = {key: Q(int(round(value*denominator)), denominator) for key, value in candidate_coefficients.items()}
            nonzero = {key: str(value) for key, value in a_q.items() if value}
            if c_q <= 0 or not nonzero:
                raise ModelRejected('Invalid rationalized c/support.')
            actual_max = max(max(map(int, key.split(','))) for key in nonzero)
            support = (actual_max+2)*H_Q
            if R_Q > support:
                raise ModelRejected('R exceeds actual rationalized support scan bound.')
            certificate = {'alpha': str(ALPHA_Q), 'R': str(R_Q), 'h': str(H_Q), 'N': basis.N,
                           'c': str(c_q), 'a': nonzero}
            target = ROOT/'candidate.json'
            if target.exists(): raise ModelRejected('Refusing a second rationalized candidate.')
            guard(); dump(target, certificate)
            result.update(candidate_exported=True, candidate_sha256=sha256(target.read_bytes()).hexdigest(),
                          rationalization={'denominator': str(denominator), 'clip_bump_feedback': False,
                                           'zero_classes_omitted': len(a_q)-len(nonzero), 'support_U': str(support)})
        result['status'] = 'finite_single_LP_discovery_closed'
    except BaselineRejected as exc:
        result.update(status='baseline_numerical_check_rejected', reason=str(exc))
    except StopBudget as exc:
        result.update(status='interrupted_saved_checkpoint', reason=str(exc))
    except Exception as exc:
        result.update(status='failed_or_candidate_rejected', reason=f'{type(exc).__name__}: {exc}')
    assert counts['LP_started'] <= 1
    total_points = sum(counts[k] for k in ('baseline_K_started', 'baseline_A_started', 'candidate_K_started', 'candidate_A_started'))
    assert total_points <= 62414
    result.update(point_candidate_evaluations_started=total_points,
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  module_import_seconds_separate=_IMPORT_SECONDS,
                  environment={'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__,
                               'highspy': highspy.Highs().version()},
                  scope='Finite numerical discovery only. No automatic checker, cost certificate, feedback, or repeat run.')
    save()
    print(json.dumps({'status': result['status'], 'counts': counts, 'seconds': result['elapsed_seconds'],
                      'candidate_exported': result['candidate_exported'],
                      'float_half_cost': result.get('candidate_float_half_cost'),
                      'baseline_residual': result.get('baseline_residual')}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
