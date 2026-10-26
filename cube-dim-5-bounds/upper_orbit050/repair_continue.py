#!/usr/bin/env python3
"""Fixed pending sample99 arithmetic repair, then original indices100..127 only."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
from itertools import permutations, product
import ast
import hashlib
import json
import math
import statistics
import sys
import time
import numpy as np
import mpmath as mp

BASE = Path(__file__).resolve().parent
PARENT = BASE.parent/'upper_orbit049'
START = time.monotonic()
SEED, BLOCK = 20261026049, 120


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(name, data):
    path = BASE/name
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+'\n')
    tmp.replace(path)


def checkpoint(stage, count, rng_state=None, **extra):
    out = {'pid': os.getpid(), 'stage': stage, 'completed_here': count,
           'completed_total': 99+count, 'maximum_total': 128,
           'elapsed_seconds_this_run': time.monotonic()-START,
           'resume_command': 'python -B repair_continue.py', 'is_certified': False, **extra}
    if rng_state is not None:
        out['rng_state_after'] = rng_state
    write_json('checkpoint.json', out)


def mp_inf_norm(M):
    return max(mp.fsum(abs(M[i, j]) for j in range(M.cols)) for i in range(M.rows))


registration = json.loads((BASE/'registration.json').read_text())
registration_hash = sha(BASE/'registration.json')
script_hash = sha(Path(__file__))
assert sha(PARENT/'manifest.json') == registration['parent']['manifest_sha256']
parent_manifest = json.loads((PARENT/'manifest.json').read_text())
for item in parent_manifest['white_list']:
    assert sha(PARENT/item['path']) == item['sha256'], item['path']
assert sha(PARENT/'pending.json') == registration['parent']['pending_sample99_sha256']
assert sha(PARENT/'samples.jsonl') == registration['parent']['completed_99_samples_sha256']
parent_pending = json.loads((PARENT/'pending.json').read_text())
assert parent_pending['sample_id'] == 99
assert parent_pending['rng_state_after'] == registration['continuation']['rng_state_before_sample100']
old_rows = [json.loads(s) for s in (PARENT/'samples.jsonl').read_text().splitlines() if s.strip()]
assert len(old_rows) == 99 and [r['sample_id'] for r in old_rows] == list(range(99))

if (BASE/'failure.json').exists():
    print('A 050 failure already exists. No repair retry or new draw performed.')
    sys.exit(2)
continuation_path = BASE/'continuation.jsonl'
new_rows = [json.loads(s) for s in continuation_path.read_text().splitlines() if s.strip()] if continuation_path.exists() else []
assert [r['sample_id'] for r in new_rows] == list(range(99, 99+len(new_rows)))
assert len(new_rows) <= 29
assert all(r['registration_sha256'] == registration_hash and r['script_sha256'] == script_hash for r in new_rows)

# Reuse only pure numerical functions from the pinned source, not its module-level sampler.
old_source = (PARENT/'diagnostic.py').read_text()
assert sha(PARENT/'diagnostic.py') == registration['parent']['diagnostic_script_sha256']
names = {'NumericalStop', 'require', 'finite_positive', 'orbit_energy', 'evaluate'}
old_tree = ast.parse(old_source)
selected = [node for node in old_tree.body if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in names]
assert {node.name for node in selected} == names
function_hashes = {node.name: hashlib.sha256(ast.get_source_segment(old_source, node).encode()).hexdigest() for node in selected}
exec(compile(ast.Module(body=selected, type_ignores=[]), str(PARENT/'diagnostic.py')+'[pure-functions-only]', 'exec'), globals())

permutation_list = list(permutations(range(5)))
sign_list = [(1,)+e for e in product((-1, 1), repeat=4)]
group_perm = np.array([p for p in permutation_list for s in sign_list], dtype=np.int64)
group_sign = np.array([s for p in permutation_list for s in sign_list], dtype=np.float64)
V = np.array([(1,)+e for e in product((-1, 1), repeat=4)], dtype=np.float64)
runtime049 = json.loads((PARENT/'runtime.json').read_text())
C, WEIGHT_CAP = runtime049['C'], runtime049['weight_cap']
QUAD = {int(n): (np.array(row['nodes']), np.array(row['weights'])) for n, row in runtime049['quadrature'].items()}
rng = np.random.Generator(np.random.PCG64(SEED))
rng.bit_generator.state = new_rows[-1]['rng_state_after'] if new_rows else parent_pending['rng_state_after']
write_json('runtime.json', {'numpy': np.__version__, 'mpmath': mp.__version__, 'python': sys.version,
                          'reused_function_sha256': function_hashes,
                          'source_diagnostic_sha256': sha(PARENT/'diagnostic.py'),
                          'registration_sha256': registration_hash, 'script_sha256': script_hash,
                          'saved_049_quadrature_reused': True, 'new_quadrature_computed': False,
                          'source_manifest_verified': True,
                          'label': 'noncertified_arithmetic_repair_and_original_finite_continuation'})
checkpoint('ready', len(new_rows), rng.bit_generator.state)
print(json.dumps({'pid': os.getpid(), 'existing_continuation_rows': len(new_rows),
                  'registration_sha256': registration_hash, 'script_sha256': script_hash}), flush=True)


def repair_fixed_99(pending):
    raw = pending['gamma_draws']
    assert raw == parent_pending['gamma_draws']
    fractions = [float(value).as_integer_ratio() for value in raw]
    with mp.workdps(80):
        yy = [mp.mpf(n)/mp.mpf(d) for n, d in fractions]
        require(all(mp.isfinite(y) and y > 0 for y in yy), 'mp_invalid_original_gamma', {})
        rr = mp.fsum(yy)
        zz = [y/rr for y in yy]
        hh = mp.fsum(1/z for z in zz)
        pp = [(1/z)/hh for z in zz]
        vv = [(1,)+e for e in product((-1, 1), repeat=4)]
        MM = mp.matrix(5, 5)
        for i in range(5):
            for j in range(5):
                MM[i, j] = mp.fsum(p*v[i]*v[j] for p, v in zip(pp, vv))
        det_mp = mp.det(MM)
        require(mp.isfinite(det_mp) and det_mp > 0, 'mp_nonpositive_determinant', {'det': mp.nstr(det_mp, 82)})
        AA = MM**-1
        ev = mp.eigsy(MM, eigvals_only=True)
        require(all(x > 0 and mp.isfinite(x) for x in ev), 'mp_nonpositive_eigenvalue', {})
        cond_mp = ev[4]/ev[0]
        residual_mp = mp_inf_norm(MM*AA-mp.eye(5))
        asym_mp = mp_inf_norm(AA-AA.T)/mp_inf_norm(AA)
        p_sum_error = abs(mp.fsum(pp)-1)
        diagonal_error = max(abs(MM[i, i]-1) for i in range(5))
        mp_data = {
            'status': '80_digit_reconstruction_computed_not_yet_float_path_passed',
            'sample_id': 99, 'decimal_digits': 80,
            'input_exact_binary_fractions': [f'{n}/{d}' for n, d in fractions],
            'Z_80digits': [mp.nstr(z, 82) for z in zz],
            'P_80digits': [mp.nstr(p, 82) for p in pp],
            'M_80digits': [[mp.nstr(MM[i, j], 82) for j in range(5)] for i in range(5)],
            'inverse_80digits': [[mp.nstr(AA[i, j], 82) for j in range(5)] for i in range(5)],
            'determinant_80digits': mp.nstr(det_mp, 82),
            'eigenvalues_80digits': [mp.nstr(ev[i], 82) for i in range(5)],
            'condition_80digits': mp.nstr(cond_mp, 82),
            'inverse_residual_inf_80digits': mp.nstr(residual_mp, 82),
            'inverse_relative_asymmetry_80digits': mp.nstr(asym_mp, 82),
            'P_sum_error_80digits': mp.nstr(p_sum_error, 82),
            'M_diagonal_error_80digits': mp.nstr(diagonal_error, 82),
            'symmetrize_or_diagonal_reset_performed': False,
            'new_random_draw_for_repair': False,
            'source_pending_sha256': sha(PARENT/'pending.json')
        }
        write_json('repair99.json', mp_data)
        require(cond_mp < mp.mpf('1e12') and residual_mp <= mp.mpf('1e-60') and asym_mp <= mp.mpf('1e-60')
                and p_sum_error <= mp.mpf('1e-70') and diagonal_error <= mp.mpf('1e-70'),
                'mp_precision_precheck_failure', {'repair99_path': str(BASE/'repair99.json')})
        z = np.array([float(x) for x in zz])
        p = np.array([float(x) for x in pp])
        H, radial_draw = float(hh), float(rr)
        M = np.array([[float(MM[i, j]) for j in range(5)] for i in range(5)])
        A = np.array([[float(AA[i, j]) for j in range(5)] for i in range(5)])
        logd, det_value = float(mp.log(det_mp)), float(det_mp)
    require(finite_positive(z) and finite_positive(p) and math.isfinite(H) and H >= 256*(1-1e-12),
            'converted_shape_invalid', {'H': H})
    weight = C/math.sqrt(H)
    require(math.isfinite(weight) and 0 < weight <= WEIGHT_CAP*(1+1e-12), 'invalid_tilt_weight', {'weight': weight})
    eigenvalues = np.linalg.eigvalsh(M)
    sign, float_logd = np.linalg.slogdet(M)
    require(finite_positive(eigenvalues) and sign > 0 and math.isfinite(float(float_logd)) and math.isfinite(logd),
            'converted_matrix_invalid', {'eigenvalues': eigenvalues.tolist(), 'sign': float(sign)})
    condition = float(eigenvalues[-1]/eigenvalues[0])
    residual = float(np.linalg.norm(M@A-np.eye(5), ord=np.inf))
    asymmetry = float(np.linalg.norm(A-A.T, ord=np.inf)/np.linalg.norm(A, ord=np.inf))
    require(condition < 1e12, 'shape_condition_limit', {'condition': condition})
    require(np.all(np.isfinite(A)) and residual <= 1e-5 and asymmetry <= 1e-10,
            'converted_inverse_accuracy_failure', {'residual': residual, 'asymmetry': asymmetry})
    require(finite_positive(np.linalg.eigvalsh((A+A.T)/2)), 'nonpositive_inverse_symmetric_part', {})
    r16, diag16 = orbit_energy(A, logd, 16)
    factor = math.exp(-.1*logd)
    base, orbit = weight*factor, weight*factor*r16**.2
    require(finite_positive([factor, base, orbit]), 'invalid_repair_contribution', {})
    old_failure = json.loads((PARENT/'failure.json').read_text())
    old_M = np.array(old_failure['detail']['M'])
    mp_data.update({'status': 'same_sample_passed_unchanged_float_checks_and_complete_16_node_group',
                    'converted_inverse_residual_inf': residual,
                    'converted_inverse_relative_asymmetry': asymmetry,
                    'converted_condition_2': condition,
                    'converted_M_minus_old_float_M_max_abs': float(np.max(np.abs(M-old_M))),
                    'float_M_logdet_minus_mp_logdet': float(float_logd)-logd,
                    'identity_kernel_error_after_conversion': diag16['identity_kernel_maximum_error'],
                    'orbit_r16': r16, 'is_certified': False})
    write_json('repair99.json', mp_data)
    return {'sample_id': 99, 'gamma_draws': raw, 'rng_state_after': pending['rng_state_after'],
            'Z': z.tolist(), 'P': p.tolist(), 'M': M.tolist(),
            'gamma_normalizer_unused_in_weight': radial_draw,
            'H': H, 'tilt_weight': weight, 'P_sum': float(np.sum(p)),
            'determinant': det_value, 'log_determinant': logd,
            'eigenvalues': eigenvalues.tolist(), 'condition_2': condition,
            'inverse_residual_inf': residual, 'inverse_relative_asymmetry': asymmetry,
            'det_negative_tenth': factor, 'orbit_r_16': r16,
            'weighted_unimproved_contribution': base, 'weighted_orbit_contribution': orbit,
            'orbit_diagnostics_16': diag16, 'crosscheck_8': None,
            'registration_sha256': registration_hash, 'script_sha256': script_hash,
            'arithmetic': '80-digit intended shape/determinant/inverse followed by unchanged float group path',
            'repair99_sha256': sha(BASE/'repair99.json'), 'is_certified': False}


for sample_id in range(99+len(new_rows), 128):
    pending_path = BASE/'pending.json'
    if pending_path.exists():
        pending = json.loads(pending_path.read_text())
        require(pending['sample_id'] == sample_id, 'pending_resume_mismatch', {'pending': pending['sample_id'], 'expected': sample_id})
        rng.bit_generator.state = pending['rng_state_after']
    elif sample_id == 99:
        pending = dict(parent_pending)
        pending['repair_parent_pending_sha256'] = sha(PARENT/'pending.json')
        write_json('pending.json', pending)
    else:
        assert new_rows and new_rows[0]['sample_id'] == 99
        assert new_rows[0]['arithmetic'].startswith('80-digit')
        y = rng.gamma(shape=.5, scale=1.0, size=16)
        pending = {'sample_id': sample_id, 'gamma_draws': y.tolist(), 'rng_state_after': rng.bit_generator.state,
                   'seed': SEED, 'registration_sha256': registration_hash, 'script_sha256': script_hash}
        write_json('pending.json', pending)
    try:
        with np.errstate(over='raise', divide='raise', invalid='raise', under='raise'):
            row = repair_fixed_99(pending) if sample_id == 99 else evaluate(pending)
        if sample_id != 99:
            row['arithmetic'] = 'unchanged_original_049_float64_path'
    except Exception as exc:
        detail = {'stage': exc.stage, 'detail': exc.detail} if isinstance(exc, NumericalStop) else {'stage': type(exc).__name__, 'detail': str(exc)}
        failure = {'status': 'repair_or_continuation_stopped_no_retry', 'pid': os.getpid(),
                   'sample_id': sample_id, 'completed_here': len(new_rows), 'completed_total': 99+len(new_rows),
                   'pending_input': pending, 'elapsed_seconds': time.monotonic()-START,
                   'registration_sha256': registration_hash, 'script_sha256': script_hash, **detail}
        write_json('failure.json', failure)
        checkpoint('numerical_failure', len(new_rows), rng.bit_generator.state,
                   failed_sample=sample_id, failure_stage=detail['stage'])
        print(json.dumps({'status': 'failure', 'sample_id': sample_id, **detail}), flush=True)
        sys.exit(2)
    with continuation_path.open('a') as handle:
        handle.write(json.dumps(row, sort_keys=True, allow_nan=False)+'\n')
        handle.flush()
        os.fsync(handle.fileno())
    new_rows.append(row)
    pending_path.unlink()
    checkpoint('continuing_original_finite_plan', len(new_rows), rng.bit_generator.state)
    if sample_id == 99 or (sample_id+1) % 8 == 0:
        print(json.dumps({'sample_id_completed': sample_id, 'total_completed': sample_id+1,
                          'condition': row['condition_2'], 'inverse_asymmetry': row['inverse_relative_asymmetry'],
                          'elapsed_seconds': time.monotonic()-START}), flush=True)
    if time.monotonic()-START >= 110 and sample_id < 127:
        checkpoint('recoverable_time_yield', len(new_rows), rng.bit_generator.state)
        print('Saved checkpoint for the same original finite plan.', flush=True)
        sys.exit(75)


def stats(values):
    return {'mean': statistics.fmean(values), 'sample_sd': statistics.stdev(values),
            'descriptive_mc_se': statistics.stdev(values)/math.sqrt(len(values)),
            'minimum': min(values), 'maximum': max(values), 'is_confidence_certificate': False}


all_rows = old_rows+new_rows
assert len(all_rows) == 128 and [r['sample_id'] for r in all_rows] == list(range(128))
assert len(new_rows) == 29
weights = [r['tilt_weight'] for r in all_rows]
base_stats = stats([r['weighted_unimproved_contribution'] for r in all_rows])
orbit_stats = stats([r['weighted_orbit_contribution'] for r in all_rows])
c1, face, q = 1/(32*math.pi**3), math.pi**2/60, 1.25
I_base, I_orbit = c1**q*face*base_stats['mean']**q, c1**q*face*orbit_stats['mean']**q
target1, target2 = 5.36e-5, 32*(1e-4/5)**q
combined_sha = hashlib.sha256((PARENT/'samples.jsonl').read_bytes()+continuation_path.read_bytes()).hexdigest()
results = {
    'status': 'original_128_shape_plan_complete_after_fixed_sample99_arithmetic_repair',
    'seed': SEED, 'total_shapes': 128, 'unchanged_049_completed_shapes': 99,
    'same_input_repaired_shapes': 1, 'new_original_remaining_draws': 28,
    'new_seed_or_replacement_shape': False,
    'weight_statistics': stats(weights), 'weight_exact_expectation': 1,
    'weight_sample_ess': sum(weights)**2/sum(x*x for x in weights),
    'unimproved_weighted_J': base_stats, 'full_orbit_weighted_J': orbit_stats,
    'r16_statistics': stats([r['orbit_r_16'] for r in all_rows]),
    'condition_statistics': stats([r['condition_2'] for r in all_rows]),
    'max_condition_sample_id': max(all_rows, key=lambda r:r['condition_2'])['sample_id'],
    'max_identity_kernel_error': max(r['orbit_diagnostics_16']['identity_kernel_maximum_error'] for r in all_rows),
    'max_inverse_residual': max(r['inverse_residual_inf'] for r in all_rows),
    'max_inverse_relative_asymmetry': max(r['inverse_relative_asymmetry'] for r in all_rows),
    'max_node_kernel_unclipped': max(r['orbit_diagnostics_16']['maximum_node_kernel'] for r in all_rows),
    'repair_sample99': {'tilt_weight': new_rows[0]['tilt_weight'], 'condition': new_rows[0]['condition_2'],
                        'orbit_r16': new_rows[0]['orbit_r_16'],
                        'unimproved_contribution': new_rows[0]['weighted_unimproved_contribution'],
                        'orbit_contribution': new_rows[0]['weighted_orbit_contribution'],
                        'repair_record_sha256': sha(BASE/'repair99.json')},
    'first_four_8_16_crosschecks_unchanged': [{'sample_id': r['sample_id'], **r['crosscheck_8']} for r in old_rows[:4]],
    'noncertified_plugin_I': {'unimproved': I_base, 'full_orbit': I_orbit,
                             'orbit_to_unimproved_ratio': I_orbit/I_base,
                             'target_5_36e_minus5': target1, 'orbit_minus_target': I_orbit-target1,
                             'target_for_g_1e_minus4_cost32': target2, 'orbit_minus_g_target': I_orbit-target2,
                             'is_bound': False,
                             'meaning': 'Plug-in of the orbit-Minkowski expectation functional, not true face I or a certified upper bound.'},
    'determinant_calls_here': sum(r['orbit_diagnostics_16']['determinant_calls'] for r in new_rows),
    'determinant_calls_complete_main_plan': sum(r['orbit_diagnostics_16']['determinant_calls'] for r in all_rows),
    'old_crosscheck_determinant_calls': sum(r['crosscheck_8']['diagnostics_8']['determinant_calls'] for r in old_rows[:4]),
    'thresholds_changed': False, 'symmetrize_clipping_or_rejection_used': False,
    'new_certified_I_bound': False, 'new_C5_bound': False,
    'parent_manifest_sha256': sha(PARENT/'manifest.json'),
    'unchanged_prefix_sha256': sha(PARENT/'samples.jsonl'),
    'continuation_sha256': sha(continuation_path),
    'virtual_concatenation_128_jsonl_sha256': combined_sha,
    'registration_sha256': registration_hash, 'script_sha256': script_hash
}
assert results['determinant_calls_here'] == 890880
assert results['determinant_calls_complete_main_plan'] == 3932160
assert results['old_crosscheck_determinant_calls'] == 61440
# Confirm all 049 bytes still match their frozen parent manifest after continuation.
for item in parent_manifest['white_list']:
    assert sha(PARENT/item['path']) == item['sha256'], item['path']
write_json('results.json', results)
checkpoint('complete', len(new_rows), rng.bit_generator.state)
receipt = {'pid': os.getpid(), 'elapsed_seconds_this_run': time.monotonic()-START, 'exit_status': 0,
           'completed_here': len(new_rows), 'completed_total': 128, 'mpmath_digits_repair': 80,
           'results_sha256': sha(BASE/'results.json'), 'continuation_sha256': sha(continuation_path),
           'repair99_sha256': sha(BASE/'repair99.json'),
           'registration_sha256': registration_hash, 'script_sha256': script_hash,
           'all_parent_manifest_files_rechecked_unchanged': True}
write_json('run_receipt.json', receipt)
print(json.dumps(receipt, sort_keys=True), flush=True)
print(json.dumps({'weight_mean': results['weight_statistics']['mean'],
                  'J_base': base_stats['mean'], 'J_orbit': orbit_stats['mean'],
                  'I_plugin_base': I_base, 'I_plugin_orbit': I_orbit,
                  'max_condition': results['condition_statistics']['maximum'],
                  'is_certified': False}, sort_keys=True), flush=True)
