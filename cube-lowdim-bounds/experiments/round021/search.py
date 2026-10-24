#!/usr/bin/env python3
"""One finite round02 experiment; importing performs no rows, solves or scans."""
from __future__ import annotations
from time import monotonic
_IMPORT_START = monotonic()
from fractions import Fraction as Q
from hashlib import sha256
import argparse
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys

# Set before NumPy/SciPy/HiGHS are imported in both parent and pair children.
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import numpy as np
import scipy
import scipy.sparse as sp
import highspy
_IMPORT_SECONDS = monotonic()-_IMPORT_START
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
CACHE = ROOT/'cache'
U = Q(9, 4)
ALPHAS = ('23/20', '47/40', '6/5')
RADII = ('13/8', '33/20', '67/40')
PAIRS = [('47/40', '33/20')]+[(a, r) for a in ALPHAS for r in RADII if (a, r) != ('47/40', '33/20')]
TOTAL_SECONDS, PREP_SECONDS, PAIR_SECONDS, LP_SECONDS = 300., 45., 25., 8.
PUBLIC_FIRST_SCALE_SHA = '21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
K_TOL, A_TOL = 1e-8, 1e-7


class StopBudget(Exception):
    pass


class PairRejected(Exception):
    pass


def dump(path, value):
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    tmp.replace(path)


def finite(label, values):
    values = np.asarray(values)
    if not np.all(np.isfinite(values)):
        raise PairRejected(label+'_nonfinite')
    return values


def formula_directory():
    for path in (ROOT/'formula_modules', ROOT.parent/'upstream-v4/search'):
        if all((path/n).is_file() for n in ('bspl_closed.py', 'sing2.py', 'sing3.py')):
            return path
    raise FileNotFoundError('Unmodified formula modules not found.')


def load_formulas():
    directory = formula_directory()
    sys.path.insert(0, str(directory))
    result = []
    for name in ('bspl_closed', 'sing3'):
        spec = importlib.util.spec_from_file_location('round02_'+name, directory/(name+'.py'))
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        result.append(module)
    return result


class Basis:
    def __init__(self, path):
        self.bytes = path.read_bytes()
        if sha256(self.bytes).hexdigest() != PUBLIC_FIRST_SCALE_SHA:
            raise ValueError('First-round input certificate hash mismatch.')
        self.spec = d = json.loads(self.bytes)
        if d.get('extra'):
            raise ValueError('No extra component is approved.')
        self.hq, self.N = Q(d['h']), int(d['N'])
        if self.hq != Q(1, 32) or self.N != 80:
            raise ValueError('Fixed h/N changed.')
        self.h = float(self.hq)
        self.classes = sorted(tuple(map(int, s.split(','))) for s in d['a'])
        if len(self.classes) != 1181 or any(not 0 <= j <= i <= self.N for i, j in self.classes):
            raise ValueError('Fixed 1181-class basis changed.')
        self.n, self.nvar = 1181, 1182
        self.I = np.arange(-self.N, self.N+1)
        self.mapping = np.full((2*self.N+1, 2*self.N+1), -1, dtype=np.int64)
        self.mult = np.zeros(self.n)
        lookup = {key: i for i, key in enumerate(self.classes)}
        for i in self.I:
            for j in self.I:
                k = lookup.get((max(abs(int(i)), abs(int(j))), min(abs(int(i)), abs(int(j)))))
                if k is not None:
                    self.mapping[i+self.N, j+self.N] = k
                    self.mult[k] += 1
        self.base = np.asarray([float(Q(d['a'][f'{i},{j}'])) for i, j in self.classes]+[float(Q(d['c']))])
        self.lower = self.base-.125
        self.upper = self.base+.125
        self.lower[-1], self.upper[-1] = .9*self.base[-1], 1.1*self.base[-1]

    def objective(self, alpha, radius):
        return finite('half_mass_cost', np.r_[self.h*self.h*self.mult/2,
                                            alpha*radius**(2-alpha)/(2-alpha)])


def grid(denominator):
    points = set()
    for m in range(int(U*denominator)+1):
        u = (Q(m)+Q(1, 3))/denominator
        if u > U:
            continue
        for n in range(m+1):
            points.add((u, (Q(n)+Q(1, 3))/denominator))
    return points


def make_pools(radius):
    radius = Q(radius)
    initial_K = grid(32)
    for m in range(73):
        for n in range(m+1):
            if m or n:
                initial_K.add((Q(m, 32), Q(n, 32)))
    initial_K |= {(1-Q(j, 64), Q(j, 64)) for j in range(33)}
    validation_K = grid(64)|{(1-Q(j, 128), Q(j, 128)) for j in range(65)}
    origin = {(Q(1, 2**k)*(1-f), Q(1, 2**k)*f)
              for k in range(3, 19) for f in (Q(1, 2), Q(1, 4), Q(1, 16), Q(1, 256))}
    shift = (radius-Q(33, 20))/2
    valley = {(Q('1.541724')+shift+Q(i, 4096), Q('.132603')+shift+Q(j, 4096))
              for i in range(-2, 3) for j in range(-2, 3)}
    initial_A = grid(32)|origin|valley
    validation_A = grid(64)|origin|valley
    for destination, ks, denominator, js in ((initial_A, range(6, 15), 64, range(1, 33)),
                                            (validation_A, range(4, 16), 128, range(1, 65))):
        for sign, k, j in itertools.product((-1, 1), ks, js):
            r = radius+sign*Q(1, 2**k)
            v = r*Q(j, denominator)
            destination.add((r-v, v))
    original_sizes = [len(s) for s in (initial_K, initial_A, validation_K, validation_A)]
    initial_A = {(u, v) for u, v in initial_A if 0 < v <= u <= U and u+v != radius}
    validation_A = {(u, v) for u, v in validation_A if 0 < v <= u <= U and u+v != radius}
    groups = (initial_K, initial_A, validation_K, validation_A)
    sizes = [len(s) for s in groups]
    if any(n > cap for n, cap in zip(sizes, (6144, 4096, 12288, 12288))):
        raise ValueError('Fixed point-pool cap exceeded; do not truncate or refill.')
    if sizes != [5344, 3293, 10505, 12065]:
        raise ValueError(f'Point-pool implementation mismatch: {sizes}')
    pools = {name: np.asarray([(float(u), float(v)) for u, v in sorted(s)])
             for name, s in zip(('initial_K', 'initial_A', 'validation_K', 'validation_A'), groups)}
    return pools, {'before_A_domain_filter': original_sizes, 'actual': dict(zip(pools, sizes))}


def spline_K_rows(points, basis, bspl, guard):
    rows, cols, values = [], [], []
    for first in range(0, len(points), 256):
        guard()
        block = points[first:first+256]
        Bu = finite('Bu', bspl.Bv(block[:, 0, None]/basis.h-basis.I[None, :]))
        Bv = finite('Bv', bspl.Bv(block[:, 1, None]/basis.h-basis.I[None, :]))
        for t in range(len(block)):
            ix, iy = np.flatnonzero(Bu[t]), np.flatnonzero(Bv[t])
            c = basis.mapping[np.ix_(ix, iy)].ravel()
            a = np.outer(Bu[t, ix], Bv[t, iy]).ravel()
            mask = (c >= 0)&(a != 0)
            rows.extend([first+t]*int(mask.sum())); cols.extend(c[mask]); values.extend(a[mask])
    guard()
    out = sp.coo_matrix((values, (rows, cols)), shape=(len(points), basis.n)).tocsr()
    out.eliminate_zeros()
    finite('K_spline_rows', out.data)
    return out


def worker(index):
    start = monotonic()
    state_path = ROOT/f'pair_{index:02d}.json'
    alpha_s, radius_s = PAIRS[index]
    alpha, radius = float(Q(alpha_s)), float(Q(radius_s))
    counts = {'LP_started': 0, 'LP_completed': 0,
              'initial_K_rows_started': 0, 'initial_K_rows_completed': 0,
              'initial_A_rows_started': 0, 'initial_A_rows_completed': 0,
              'validation_K_rows_started': 0, 'validation_K_rows_completed': 0,
              'validation_A_rows_started': 0, 'validation_A_rows_completed': 0,
              'validation_K_evaluations_started': 0, 'validation_K_evaluations_completed': 0,
              'validation_A_evaluations_started': 0, 'validation_A_evaluations_completed': 0,
              'stored_initial_residual_matrix_rows_started': 0,
              'stored_initial_residual_matrix_rows_checked': 0}
    state = {'pair_index': index, 'alpha': alpha_s, 'R': radius_s, 'status': 'started',
             'counts': counts, 'solves': [], 'validation': [], 'eligible': False,
             'module_import_seconds_separate': _IMPORT_SECONDS}

    def save():
        state['worker_seconds'] = monotonic()-start
        dump(state_path, state)

    def guard():
        if monotonic()-start >= PAIR_SECONDS:
            raise StopBudget('pair_internal_time_limit')

    save()
    try:
        if not 1 < alpha < 2 or radius <= 1:
            raise PairRejected('analytic_parameter_domain')
        basis = Basis(CACHE/'baseline.json')
        bspl, singular = load_formulas()
        rindex = RADII.index(radius_s)
        z = np.load(CACHE/f'pools_{rindex}.npz')
        pools = {key: z[key] for key in z.files}
        Kbase = sp.load_npz(CACHE/'initial_K_spline.npz')
        Kvalbase = sp.load_npz(CACHE/'validation_K_spline.npz')
        guard()

        def Krows(name, base):
            points = pools[name]
            counts[name+'_rows_started'] += len(points); save(); guard()
            r = points.sum(axis=1)
            psi = finite(name+'_psi', np.maximum(r**(-alpha)-radius**(-alpha), 0.))
            matrix = sp.hstack((base, sp.csr_matrix(psi[:, None])), format='csr')
            rhs = np.where(r <= 1., 1.+2e-6, 0.)
            guard(); counts[name+'_rows_completed'] += len(points); save()
            return matrix, rhs

        def Arows(name):
            points = pools[name]
            counts[name+'_rows_started'] += len(points); save()
            if np.any(points[:, 1] <= 0) or np.any(abs(points.sum(axis=1)-radius) < 1e-12):
                raise PairRejected(name+'_invalid_A_domain')
            blocks, true_singular, effective_singular = [], [], []
            for first in range(0, len(points), 256):
                guard()
                p = points[first:first+256]; r = p.sum(axis=1)
                Bu = finite('A_Bu', bspl.Bv(p[:, 0, None]/basis.h-basis.I[None, :]))
                Bv = finite('A_Bv', bspl.Bv(p[:, 1, None]/basis.h-basis.I[None, :]))
                Lu = finite('A_Lu', bspl.LBv(p[:, 0, None]/basis.h-basis.I[None, :], alpha)*basis.h**(-alpha))
                Lv = finite('A_Lv', bspl.LBv(p[:, 1, None]/basis.h-basis.I[None, :], alpha)*basis.h**(-alpha))
                aps = finite('Apsi', singular.A_psi_vec(alpha, radius, p[:, 0], p[:, 1]))
                extra = np.zeros(len(p))
                near = r < .4
                extra[near] = .01*singular.C_alpha(alpha)*r[near]**(-2*alpha)
                effective = finite('effective_Apsi', np.minimum(aps-extra, 1e4))
                rr, cc, vv = [], [], []
                for t in range(len(p)):
                    guard()
                    ix, iy = np.flatnonzero(Bu[t]), np.flatnonzero(Bv[t])
                    c1, a1 = basis.mapping[:, iy].ravel(), np.outer(Lu[t], Bv[t, iy]).ravel()
                    c2, a2 = basis.mapping[ix, :].ravel(), np.outer(Bu[t, ix], Lv[t]).ravel()
                    for columns, values in ((c1, a1), (c2, a2)):
                        keep = (columns >= 0)&(values != 0)
                        rr.extend([t]*int(keep.sum())); cc.extend(columns[keep]); vv.extend(values[keep])
                block = sp.coo_matrix((vv, (rr, cc)), shape=(len(p), basis.n)).tocsr()
                block.eliminate_zeros(); finite('A_spline_data', block.data)
                blocks.append(block); true_singular.append(aps); effective_singular.append(effective)
                counts[name+'_rows_completed'] += len(p); save()
            guard()
            spline = sp.vstack(blocks, format='csr')
            true = np.concatenate(true_singular); effective = np.concatenate(effective_singular)
            matrix = sp.hstack((spline, sp.csr_matrix(effective[:, None])), format='csr')
            r = points.sum(axis=1)
            rhs = np.where(r <= radius, .001, np.where(r <= 2.6, .0001, 0.))
            return matrix, rhs, spline, true

        Ki, rKi = Krows('initial_K', Kbase)
        Ai, rAi, _, _ = Arows('initial_A')
        cost = basis.objective(alpha, radius)
        solver = highspy.Highs()
        for key, value in (('output_flag', False), ('threads', 1), ('solver', 'simplex'),
                           ('simplex_strategy', 1), ('primal_feasibility_tolerance', 1e-9),
                           ('dual_feasibility_tolerance', 1e-9)):
            solver.setOptionValue(key, value)
        solver.addVars(basis.nvar, basis.lower, basis.upper)
        solver.changeColsCost(basis.nvar, np.arange(basis.nvar, dtype=np.int32), cost)

        def add_rows(matrix, lower, scale_A=False):
            guard()
            if scale_A:
                scales = 1./np.maximum(1e-12, abs(matrix).max(axis=1).toarray().ravel())
                matrix = (sp.diags(scales)@matrix).tocsr(); lower = lower*scales
            solver.addRows(matrix.shape[0], lower, np.full(matrix.shape[0], highspy.kHighsInf), matrix.nnz,
                           matrix.indptr[:-1].astype(np.int32), matrix.indices.astype(np.int32), matrix.data)
            guard()

        add_rows(Ki, rKi); add_rows(Ai, rAi, True)

        def solve(stage):
            guard()
            if counts['LP_started'] >= 2:
                raise StopBudget('two_LP_limit')
            solver.setOptionValue('time_limit', min(LP_SECONDS, PAIR_SECONDS-(monotonic()-start)))
            counts['LP_started'] += 1
            record = {'stage': stage, 'status': 'started', 'same_model_basis_warm_start': stage == 2,
                      'started_worker_seconds': monotonic()-start}
            state['solves'].append(record); save()
            before = monotonic()
            solver.run()  # Exactly one call; no IPM/restart/polish fallback.
            status = solver.getModelStatus()
            record.update(status=solver.modelStatusToString(status), seconds=monotonic()-before,
                          optimal=status == highspy.HighsModelStatus.kOptimal)
            counts['LP_completed'] += 1
            solution = solver.getSolution()
            record['valid_primal'] = bool(solution.value_valid)
            infeasibility = float(solver.getInfo().max_primal_infeasibility)
            record['max_primal_infeasibility'] = infeasibility if math.isfinite(infeasibility) else None
            save(); guard()
            if not record['optimal'] or not record['valid_primal']:
                raise PairRejected('LP_nonoptimal_or_no_valid_primal')
            x = finite('LP_solution', solution.col_value)
            bound_error = max(float(np.max(basis.lower-x)), float(np.max(x-basis.upper)), 0.)
            record['bound_violation'] = bound_error
            if bound_error > 1e-8:
                raise PairRejected('LP_variable_bound_violation')
            # Stored-matrix arithmetic, not new point generation/formula rows.
            # Solver tolerances apply to scaled A rows; check original units.
            counts['stored_initial_residual_matrix_rows_started'] += Ki.shape[0]+Ai.shape[0]
            save(); guard()
            initial_K_violation = finite('stored_initial_K_residual', rKi-Ki@x)
            initial_A_violation = finite('stored_initial_A_residual', rAi-Ai@x)
            guard()
            counts['stored_initial_residual_matrix_rows_checked'] += Ki.shape[0]+Ai.shape[0]
            record['stored_initial_K_max_violation'] = float(initial_K_violation.max())
            record['stored_initial_A_max_violation'] = float(initial_A_violation.max())
            record['stored_initial_constraints_passed'] = bool(initial_K_violation.max() <= K_TOL and
                                                               initial_A_violation.max() <= A_TOL)
            save()
            if not record['stored_initial_constraints_passed']:
                raise PairRejected('unscaled_stored_initial_constraint_violation')
            record['half_mass_objective'] = float(cost@x)
            np.save(ROOT/f'pair_{index:02d}_lp{stage}.npy', x)
            save()
            return x

        x = solve(1)
        Kv, rKv = Krows('validation_K', Kvalbase)
        Av, rAv, Aval_spline, Atrue_sing = Arows('validation_A')

        def validate(x, stage):
            guard()
            counts['validation_K_evaluations_started'] += len(pools['validation_K']); save()
            kval = finite('validation_K_values', Kv@x)
            vK = rKv-kval
            counts['validation_K_evaluations_completed'] += len(kval); save(); guard()
            counts['validation_A_evaluations_started'] += len(pools['validation_A']); save()
            aval = finite('validation_effective_A', Av@x)
            true_aval = finite('validation_true_A', Aval_spline@x[:-1]+Atrue_sing*x[-1])
            vA = rAv-aval
            counts['validation_A_evaluations_completed'] += len(aval)
            badK, badA = np.flatnonzero(vK > K_TOL), np.flatnonzero(vA > A_TOL)
            record = {'stage': stage, 'K_points': len(kval), 'A_points': len(aval),
                      'bad_K': len(badK), 'bad_A': len(badA), 'max_K_violation': float(vK.max()),
                      'max_effective_A_violation': float(vA.max()), 'minimum_true_K': float(kval.min()),
                      'minimum_true_K_on_D': float(kval[pools['validation_K'].sum(axis=1) <= 1].min()),
                      'minimum_true_A': float(true_aval.min()), 'minimum_effective_A': float(aval.min()),
                      'half_mass_objective': float(cost@x)}
            state['validation'].append(record); save(); guard()
            worstK = badK[np.argsort(-vK[badK], kind='stable')[:1024]]
            worstA = badA[np.argsort(-vA[badA], kind='stable')[:1024]]
            return worstK, worstA

        worstK, worstA = validate(x, 1)
        if len(worstK) or len(worstA):
            state['separation'] = {'added_K': len(worstK), 'added_A': len(worstA),
                                   'K_validation_indices': worstK.tolist(), 'A_validation_indices': worstA.tolist()}
            save()
            if len(worstK): add_rows(Kv[worstK], rKv[worstK])
            if len(worstA): add_rows(Av[worstA], rAv[worstA], True)
            x = solve(2)
            validate(x, 2)
        else:
            state['separation'] = {'added_K': 0, 'added_A': 0, 'second_LP_skipped': True}
        guard()
        final = state['validation'][-1]
        state.update(status='complete', eligible=final['bad_K'] == 0 and final['bad_A'] == 0,
                     half_mass_objective=float(cost@x), final_solution=f'pair_{index:02d}_final.npy')
        np.save(ROOT/state['final_solution'], x)
        save()
    except StopBudget as exc:
        state.update(status='interrupted', reason=str(exc), eligible=False); save()
    except Exception as exc:
        state.update(status='rejected', reason=f'{type(exc).__name__}: {exc}', eligible=False); save()


def prepare(start):
    def guard():
        if monotonic()-start >= PREP_SECONDS:
            raise StopBudget('shared_prepare_45s_limit')
    CACHE.mkdir(exist_ok=False)
    source = ROOT.parent/'round01/scaled_certificate.json'
    if not source.is_file(): source = ROOT/'input_baseline.json'
    basis = Basis(source); (CACHE/'baseline.json').write_bytes(basis.bytes)
    bspl, _ = load_formulas(); guard()
    manifests = {}
    for index, radius in enumerate(RADII):
        guard(); pools, counts = make_pools(radius)
        np.savez(CACHE/f'pools_{index}.npz', **pools)
        manifests[radius] = counts
        if index == 0:
            for name in ('initial_K', 'validation_K'):
                matrix = spline_K_rows(pools[name], basis, bspl, guard)
                sp.save_npz(CACHE/(name+'_spline.npz'), matrix)
    arithmetic_path = ROOT.parent/'round01/scaling_arithmetic.json'
    if not arithmetic_path.is_file(): arithmetic_path = ROOT/'input_scaling_arithmetic.json'
    arithmetic = json.loads(arithmetic_path.read_text())
    baseline_upper = Q(arithmetic['new_cost_upper'])
    directory = formula_directory()
    manifest = {'baseline_sha256': sha256(basis.bytes).hexdigest(),
                'baseline_cost_upper': str(baseline_upper), 'baseline_cost_lower': arithmetic['new_cost_lower'],
                'baseline_arithmetic_sha256': sha256(arithmetic_path.read_bytes()).hexdigest(),
                'point_pools': manifests, 'formula_hashes': {name: sha256((directory/name).read_bytes()).hexdigest()
                    for name in ('bspl_closed.py', 'sing2.py', 'sing3.py')},
                'shared_K_spline_rows_constructed_once': 5344+10505,
                'prepare_seconds': monotonic()-start}
    guard(); dump(CACHE/'manifest.json', manifest)
    return basis, manifest


def main():
    path = ROOT/'result.json'
    if path.exists(): raise SystemExit('Existing result: this finite discovery must not be repeated.')
    start = monotonic()
    result = {'status': 'started', 'phase': 'shared_prepare', 'pairs': [], 'candidate_exported': False,
              'limits': {'total_seconds': TOTAL_SECONDS, 'shared_prepare_seconds': PREP_SECONDS,
                         'per_pair_seconds': PAIR_SECONDS, 'per_LP_seconds': LP_SECONDS, 'maximum_LPs': 18}}
    dump(path, result)
    basis = manifest = None
    try:
        basis, manifest = prepare(start)
        result['preparation'] = manifest; dump(path, result)
        for index, (alpha, radius) in enumerate(PAIRS):
            remaining = TOTAL_SECONDS-(monotonic()-start)
            if remaining <= 0: raise StopBudget('total_300s_limit')
            timeout = min(PAIR_SECONDS, remaining)
            result['phase'] = f'pair_{index:02d}'; dump(path, result)
            print(json.dumps({'event': 'pair_start', 'index': index, 'alpha': alpha, 'R': radius}), flush=True)
            before = monotonic()
            with (ROOT/f'pair_{index:02d}.log').open('w') as log:
                process = subprocess.Popen([sys.executable, '-I', str(Path(__file__).resolve()), '--worker', str(index)],
                                           stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = process.wait(timeout=timeout); timed_out = False
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL); process.wait(); code = process.returncode; timed_out = True
            checkpoint = ROOT/f'pair_{index:02d}.json'
            row = json.loads(checkpoint.read_text()) if checkpoint.exists() else {
                'pair_index': index, 'alpha': alpha, 'R': radius, 'status': 'no_worker_checkpoint', 'counts': {}, 'eligible': False}
            row.update(outer_seconds=monotonic()-before, process_exit_code=code, outer_timed_out=timed_out)
            if timed_out or code != 0:
                row.update(status='interrupted' if timed_out else 'process_failed', eligible=False)
            result['pairs'].append(row); dump(path, result)
            print(json.dumps({'event': 'pair_end', 'index': index, 'status': row['status'],
                              'eligible': row['eligible'], 'seconds': row['outer_seconds'],
                              'LP_started': row.get('counts', {}).get('LP_started', 0)}), flush=True)
        result['phase'] = 'unique_candidate_selection'
        if monotonic()-start >= TOTAL_SECONDS: raise StopBudget('total_300s_limit_before_export')
        eligible = [row for row in result['pairs'] if row.get('eligible')]
        if eligible:
            best = min(eligible, key=lambda r: (r['half_mass_objective'], r['pair_index']))
            result['best_numerical'] = best
            if best['half_mass_objective'] < float(Q(manifest['baseline_cost_upper']))-1e-5:
                x = np.load(ROOT/best['final_solution']); denominator = 2**50
                coefficients = {f'{i},{j}': Q(int(round(float(x[k])*denominator)), denominator)
                                for k, (i, j) in enumerate(basis.classes)}
                nonzero = {key: str(value) for key, value in coefficients.items() if value}
                c = Q(int(round(float(x[-1])*denominator)), denominator)
                if c <= 0 or not nonzero: raise PairRejected('nonpositive_c_or_empty_export_support')
                imax = max(max(map(int, key.split(','))) for key in nonzero)
                support = (imax+2)*basis.hq
                if Q(best['R']) > support: raise PairRejected('R_exceeds_actual_export_support_U')
                certificate = {'alpha': best['alpha'], 'R': best['R'], 'h': str(basis.hq), 'N': basis.N,
                               'c': str(c), 'a': nonzero}
                if monotonic()-start >= TOTAL_SECONDS: raise StopBudget('total_300s_limit_before_candidate_write')
                target = ROOT/'candidate.json'
                if target.exists(): raise RuntimeError('Refusing another candidate export.')
                dump(target, certificate)
                result.update(candidate_exported=True, candidate_sha256=sha256(target.read_bytes()).hexdigest(),
                              candidate_origin_pair=best['pair_index'], export_support_U=str(support),
                              rounding={'dyadic_denominator': str(denominator), 'zero_classes_omitted': len(coefficients)-len(nonzero),
                                        'clip_or_bump_or_feedback': False},
                              floating_solution_sha256=sha256((ROOT/best['final_solution']).read_bytes()).hexdigest())
        result['status'] = 'finite_discovery_closed'
    except StopBudget as exc:
        result.update(status='interrupted_saved_checkpoints', reason=str(exc))
    except Exception as exc:
        result.update(status='failed_saved_checkpoints', reason=f'{type(exc).__name__}: {exc}')
    counts = [row.get('counts', {}) for row in result['pairs']]
    total_LPs = sum(c.get('LP_started', 0) for c in counts)
    initial_started = sum(c.get('initial_K_rows_started', 0)+c.get('initial_A_rows_started', 0) for c in counts)
    validation_started = sum(c.get('validation_K_evaluations_started', 0)+c.get('validation_A_evaluations_started', 0) for c in counts)
    assert total_LPs <= 18 and initial_started <= 92160 and validation_started <= 442368
    result.update(elapsed_seconds=monotonic()-start, module_import_seconds_separate=_IMPORT_SECONDS,
                  actual_counts={'pairs_started': len(result['pairs']), 'LP_started': total_LPs,
                                 'initial_point_rows_started': initial_started,
                                 'validation_point_candidate_evaluations_started': validation_started,
                                 'initial_plus_validation_started': initial_started+validation_started,
                                 'stored_initial_residual_matrix_rows_checked': sum(
                                     c.get('stored_initial_residual_matrix_rows_checked', 0) for c in counts)},
                  search_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  environment={'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__,
                               'highspy': highspy.Highs().version()},
                  scope='Finite floating-point discovery only; no automatic checker/cost run, repairs, new pair, or restart.')
    dump(path, result)
    print(json.dumps({'status': result['status'], 'counts': result['actual_counts'],
                      'candidate_exported': result['candidate_exported'], 'elapsed_seconds': result['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--worker', type=int)
    args = parser.parse_args()
    if args.worker is None: main()
    elif 0 <= args.worker < len(PAIRS): worker(args.worker)
    else: raise SystemExit('Invalid fixed pair index.')
