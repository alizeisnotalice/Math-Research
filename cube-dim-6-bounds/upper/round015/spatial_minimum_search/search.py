#!/usr/bin/env python3
"""B015 bounded local candidate driver. NO native evaluation without release.json.

The optimizer is heuristic. Only validated native dyadics imply point statements.
All exceptions stop and preserve the current checkpoint; there is no auto-retry.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import json
import math
import os
import subprocess
import time
from fractions import Fraction as Q
from pathlib import Path
import prepare as prep

HERE = Path(__file__).resolve().parent
SCALE = 100000.0
LOCAL_RADIUS = Q(1, 20)
M64, M625 = Q(3, 256000), Q(3, 250000)

class StopSearch(Exception):
    pass

class SeedBudget(Exception):
    pass

def dyadic_interval(obj):
    expected = {'lower_mantissa', 'upper_mantissa', 'common_binary_exponent'}
    if set(obj) != expected:
        raise ValueError('endpoint schema mismatch')
    l, u, e = (int(obj[k]) for k in ('lower_mantissa', 'upper_mantissa', 'common_binary_exponent'))
    if l > u or abs(e) > 10000:
        raise ValueError('nonfinite/unreasonable/unordered dyadic endpoints')
    scale = Q(2)**e
    return l*scale, u*scale

def midpoint(interval):
    x = float(sum(interval)/2)
    if not math.isfinite(x):
        raise ValueError('midpoint not finite')
    return x

def validate_native(result, nums, den, t=Q(6, 125)):
    if not (result['status'] == 'generic_nonzero_supported' and
            result['precision_bits'] == 512 and result['full_group_terms'] == 967680 and
            result['computed_sign_pairs'] == 483840 and result['zero_node_terms'] == 0 and
            0 <= result['confluent_terms'] <= 967680 and result['confluent_terms'] % 2 == 0 and
            result['unique_node_values'] > 0 and len(result['gradient']) == 6):
        raise ValueError('native status/precision/full-group/count mismatch')
    k, p = dyadic_interval(result['K']), dyadic_interval(result['p_T'])
    grad = [dyadic_interval(v) for v in result['gradient']]
    euler = dyadic_interval(result['Euler_x_grad_plus5K_plusTp'])
    if not (k[0] > 0 and p[0] > 0 and euler[0] <= 0 <= euler[1]):
        raise ValueError('positivity or native Euler enclosure failure')
    for pair in [k, p, *grad]:
        if pair[1]-pair[0] > max(Q(1), abs(sum(pair)/2))*Q(2)**-256:
            raise ValueError('insufficient interval accuracy; no automatic precision retry')
    # Independent endpoint interval check of x.grad K + 5 K + T p = 0.
    el, eu = 5*k[0]+t*p[0], 5*k[1]+t*p[1]
    for n, pair in zip(nums, grad):
        a, b = Q(n, den)*pair[0], Q(n, den)*pair[1]
        el += min(a, b)
        eu += max(a, b)
    if not el <= 0 <= eu:
        raise ValueError('independent exact-dyadic Euler failure')
    return {'K_interval': [str(v) for v in k], 'p_interval': [str(v) for v in p],
            'K_midpoint_display_only': midpoint(k),
            'gradient_midpoints_optimizer_only': [midpoint(v) for v in grad],
            'rejects64': k[1] <= M64, 'rejects625': k[1] < M625,
            'at_or_below_threshold625': k[1] <= M625,
            'K_upper_minus_threshold64': str(k[1]-M64),
            'K_upper_minus_threshold625': str(k[1]-M625),
            'independent_Euler_interval': [str(el), str(eu)], 'PASS': True}

def canonical_key(frame_sha, nums, den):
    # Reduce a common numerator/denominator factor, but never merge near points.
    g = math.gcd(den, *nums)
    return json.dumps([frame_sha, '6/125', [n//g for n in nums], den//g], separators=(',', ':'))

def release_gate():
    p = HERE/'release.json'
    if not p.exists():
        raise StopSearch('NOT_RELEASED: no release.json; zero native calls')
    release = json.loads(p.read_text())
    if not (release.get('scope') == 'B015-spatial-minimum-search' and
            release.get('explicit_parent_release') is True and
            release.get('C014_final_reviewed_by_parent') is True and
            release.get('authorization_message')):
        raise StopSearch('NOT_RELEASED: incomplete explicit main authorization')
    for name in ('preparation.json', 'search.py', 'prepare.py', 'registration.md'):
        if prep.sha(HERE/name) != release['prepared_file_sha256'][name]:
            raise StopSearch('release/prepared file hash mismatch: '+name)
    report_path = Path(release['C014_final_report_path'])
    if prep.sha(report_path) != release['C014_final_report_sha256']:
        raise StopSearch('C014 frozen report mismatch')
    reg = json.loads((HERE/'preparation.json').read_text())
    if report_path != Path(reg['native_report_expected_path']):
        raise StopSearch('unregistered native report path')
    for rec in reg['inputs']:
        if prep.sha(rec['path']) != rec['sha256']:
            raise StopSearch('registered input hash mismatch: '+rec['path'])
    return reg, release

def search(reg, release):
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, minimize
    checkpoint = HERE/'checkpoint.json'
    if checkpoint.exists():
        raise StopSearch('existing checkpoint: no automatic rerun/resume of this batch')
    _, rays = prep.fan_data()
    forms = prep.all_zero_forms(rays)
    forms_sha = hashlib.sha256(json.dumps(forms, separators=(',', ':')).encode()).hexdigest()
    if forms_sha != reg['zero_forms_sha256_canonical_json']:
        raise StopSearch('zero-form preflight hash mismatch')
    state = {'status': 'RUNNING_AFTER_EXPLICIT_RELEASE', 'pid': os.getpid(),
             'start_epoch': time.time(), 'native_attempts': 0, 'callback_requests': 0,
             'native_results': [], 'proposal_records': [], 'seed_runs': [],
             'global_minimum_certified': False, 'release_sha256': prep.sha(HERE/'release.json')}
    tic = time.monotonic()
    cache = {}
    outdir = HERE/'native_outputs'
    outdir.mkdir(exist_ok=False)

    def save():
        state['elapsed_seconds'] = time.monotonic()-tic
        prep.dump(checkpoint, state)

    def remaining():
        return reg['max_seconds']-(time.monotonic()-tic)

    def hard_budget():
        if remaining() <= 0:
            raise StopSearch('TIME_BUDGET_EXHAUSTED_UNRESOLVED')
        if state['native_attempts'] >= reg['max_native_attempts']:
            raise StopSearch('NATIVE_BUDGET_EXHAUSTED_UNRESOLVED')

    def native(nums, den, proposal_record):
        key = canonical_key(reg['frame_sha256'], nums, den)
        if key in cache:
            proposal_record['cache_hit'] = True
            save()
            return cache[key]
        hard_budget()
        state['native_attempts'] += 1  # Attempt counts even if execution fails.
        index = state['native_attempts']
        command = [reg['native_binary_path'], 'full', reg['fan_path'], '6', '125', str(den), *map(str, nums)]
        record = {'index': index, 'key': key, 'numerators': nums, 'denominator': den,
                  'T': '6/125', 'command': command, 'status': 'REGISTERED_BEFORE_NATIVE',
                  'binary_sha256': prep.sha(reg['native_binary_path']),
                  'source_sha256': prep.sha(reg['native_source_path']),
                  'fan_sha256': prep.sha(reg['fan_path']), 'start_epoch': time.time()}
        registered = {x['path']: x['sha256'] for x in reg['inputs']}
        for path, field in ((reg['native_binary_path'], 'binary_sha256'),
                            (reg['native_source_path'], 'source_sha256'), (reg['fan_path'], 'fan_sha256')):
            if record[field] != registered[path]:
                raise StopSearch('source/binary/fan changed before native launch')
        state['native_results'].append(record)
        save()
        stdout_path = outdir/f'{index:02d}.stdout.json'
        stderr_path = outdir/f'{index:02d}.stderr.txt'
        started = time.monotonic()
        try:
            with stdout_path.open('w') as stdout, stderr_path.open('w') as stderr:
                process = subprocess.Popen(command, stdout=stdout, stderr=stderr, cwd=HERE)
                record['native_pid'] = process.pid
                record['status'] = 'NATIVE_RUNNING'
                save()
                try:
                    code = process.wait(timeout=max(0.001, min(120.0, remaining())))
                except subprocess.TimeoutExpired:
                    process.kill()
                    code = process.wait()
                    record['timeout'] = True
                except BaseException:
                    if process.poll() is None:
                        process.kill()
                    record['returncode'] = process.wait()
                    record['status'] = 'INTERRUPTED_CHILD_KILLED_AND_WAITED'
                    raise
            record.update({'returncode': code, 'elapsed_seconds': time.monotonic()-started,
                           'stdout_path': str(stdout_path), 'stdout_sha256': prep.sha(stdout_path),
                           'stderr_path': str(stderr_path), 'stderr_sha256': prep.sha(stderr_path)})
            if code != 0 or record.get('timeout'):
                raise StopSearch('NATIVE_ERROR_OR_TIMEOUT: no retry, whole batch stopped')
            result = json.loads(stdout_path.read_text())
            check = validate_native(result, nums, den)
            record.update({'status': 'POINT_INTERVALS_VALIDATED', 'check': check})
            cache[key] = check
            save()
            if check['rejects64']:
                state['counterexample_native_index'] = index
                raise StopSearch('POINT_OBSTRUCTION_TO_STRICT_64_STOP_IMMEDIATELY')
            if check['rejects625']:
                state.setdefault('counterexamples625', []).append(index)
                save()  # This does NOT stop the still-live 64 diagnosis.
            return check
        except BaseException:
            save()
            raise

    try:
        save()
        for seed_index, seed in enumerate(reg['seeds']):
            start_count = state['native_attempts']
            initial = np.array([float(Q(n, seed['denominator'])) for n in seed['numerators'][1:]])
            center = [Q(n, seed['denominator']) for n in seed['numerators'][1:]]
            lower = [max(Q(1, 10**8), x-LOCAL_RADIUS) for x in center]
            upper = [min(1-Q(1, 10**8), x+LOCAL_RADIUS) for x in center]
            lower_i = [math.ceil(x*prep.DEN) for x in lower]
            upper_i = [math.floor(x*prep.DEN) for x in upper]
            run = {'seed_index': seed_index, 'seed_label': seed['label'], 'status': 'RUNNING',
                   'local_radius': str(LOCAL_RADIUS), 'lower': list(map(str, lower)),
                   'upper': list(map(str, upper)), 'native_start_count': start_count}
            state['seed_runs'].append(run)
            first_request = True

            def in_box(nums, den):
                return all(a <= Q(n, den) <= b for a, n, b in zip(lower, nums[1:], upper))

            def objective(z):
                nonlocal first_request
                state['callback_requests'] += 1
                if state['callback_requests'] > reg['max_optimizer_callback_requests']:
                    raise StopSearch('CALLBACK_BUDGET_EXHAUSTED_UNRESOLVED')
                if remaining() <= 0:
                    raise StopSearch('TIME_BUDGET_EXHAUSTED_UNRESOLVED')
                if not all(math.isfinite(float(x)) for x in z):
                    raise StopSearch('NONFINITE_OPTIMIZER_PROPOSAL')
                proposed = [Q.from_float(float(x)) for x in z]
                if first_request and np.array_equal(z, initial):
                    nums, den = list(seed['numerators']), seed['denominator']
                    quantization = 'registered_exact_seed'
                else:
                    den = prep.DEN
                    free = [prep.nearest(q*den) for q in proposed]
                    free = [max(a, min(b, n)) for a, n, b in zip(lower_i, free, upper_i)]
                    previous = den
                    for i in range(5):
                        free[i] = min(free[i], previous-prep.GAP)
                        previous = free[i]
                    for i in range(4, -1, -1):
                        free[i] = max(free[i], lower_i[i], (free[i+1]+prep.GAP) if i < 4 else prep.GAP)
                    nums = [den, *free]
                    if not prep.feasible(nums, den) or not in_box(nums, den):
                        raise StopSearch('EXACT_FEASIBILITY_PROJECTION_FAILED')
                    if max(abs(Q(n, den)-q) for n, q in zip(free, proposed)) > Q(1, 10**10):
                        raise StopSearch('OPTIMIZER_PROPOSAL_TOO_INFEASIBLE_FOR_ROUNDING')
                    quantization = 'nearest_den1e14_then_exact_chain_feasibility_rounding'
                first_request = False
                record = {'seed_index': seed_index, 'proposal_binary64_hex': [float(x).hex() for x in z],
                          'quantization': quantization, 'initial_numerators': nums[:],
                          'initial_denominator': den, 'repairs': []}
                state['proposal_records'].append(record)
                pf = prep.preflight(nums, den, forms, rays)
                record['initial_preflight'] = pf
                if not pf['feasible'] or not pf['int64_safe'] or not in_box(nums, den):
                    raise StopSearch('UNSUPPORTED_EXACT_PROPOSAL')
                if pf['zero_form_count']:
                    # Each repair is a NEW exact point, never a tolerance substitution.
                    success = False
                    for rec in prep.gap_repairs(nums, den):
                        candidate = rec['numerators']
                        pfc = prep.preflight(candidate, prep.DEN, forms, rays)
                        rec.update({'preflight': pfc, 'inside_registered_local_box': in_box(candidate, prep.DEN)})
                        rec['rejection_reasons'] = [name for name, failed in (
                            ('cube_or_order', not pfc['feasible']), ('int64_range', not pfc['int64_safe']),
                            ('exact_zero_node', bool(pfc['zero_form_count'])),
                            ('outside_local_box', not rec['inside_registered_local_box'])) if failed]
                        record['repairs'].append(rec)
                        if pfc['feasible'] and pfc['int64_safe'] and not pfc['zero_form_count'] and rec['inside_registered_local_box']:
                            nums, den, success = candidate, prep.DEN, True
                            break
                    if not success:
                        save()
                        raise StopSearch('UNSUPPORTED_ZERO_NODE_AFTER_SIX_REGISTERED_REPAIRS')
                record.update({'evaluated_numerators': nums, 'evaluated_denominator': den})
                key = canonical_key(reg['frame_sha256'], nums, den)
                if key not in cache and state['native_attempts']-start_count >= reg['per_seed_new_native_cap']:
                    raise SeedBudget()
                result = native(nums, den, record)
                return result['K_midpoint_display_only']*SCALE, np.array(result['gradient_midpoints_optimizer_only'][1:])*SCALE

            # Difference constraints x2-x3,...,x5-x6 >= 1e-8.
            matrix = np.zeros((4, 5))
            for i in range(4):
                matrix[i, i] = 1
                matrix[i, i+1] = -1
            try:
                result = minimize(objective, initial, jac=True, method='SLSQP',
                                  bounds=Bounds(list(map(float, lower)), list(map(float, upper))),
                                  constraints=[LinearConstraint(matrix, 1e-8, np.inf)],
                                  options={'maxiter': 16, 'ftol': 1e-10, 'disp': False})
                run.update({'status': 'OPTIMIZER_TERMINATED', 'success_heuristic_only': bool(result.success),
                            'message': str(result.message), 'nit': int(result.nit),
                            'nfev': int(result.nfev), 'njev': int(result.njev)})
                if not result.success:
                    raise StopSearch('OPTIMIZER_FAILURE_UNRESOLVED_STOP')
            except SeedBudget:
                run['status'] = 'REGISTERED_EIGHT_POINT_SEED_CAP'
            run['new_native_attempts'] = state['native_attempts']-start_count
            save()
        state['status'] = 'FINITE_SEARCH_ENDED_GLOBAL_MINIMUM_UNRESOLVED'
    except BaseException as e:
        state['status'] = str(e) if isinstance(e, StopSearch) else 'ERROR_STOPPED_WITH_CHECKPOINT'
        state['exception_type'] = type(e).__name__
        state['exception_message'] = str(e)
        if not isinstance(e, StopSearch):
            raise
    finally:
        state['all_spawned_native_processes_waited'] = True
        save()
    print(json.dumps({'status': state['status'], 'native_attempts': state['native_attempts'],
                      'elapsed_seconds': state['elapsed_seconds'], 'global_minimum_certified': False}))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()
    if not args.run:
        print('PREPARATION ONLY: use --run only after explicit parent release; no native call.')
        return
    try:
        reg, release = release_gate()
    except StopSearch as e:
        print(str(e))
        raise SystemExit(3)
    search(reg, release)

if __name__ == '__main__':
    main()
