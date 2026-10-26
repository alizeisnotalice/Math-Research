#!/usr/bin/env python3
"""Only the registered 128 Dirichlet shapes; all output is noncertified diagnostic."""
import os
for thread_env in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[thread_env] = '1'
from pathlib import Path
from itertools import permutations, product
import argparse
import hashlib
import json
import math
import sys
import time
import numpy as np
import scipy
from scipy.special import roots_jacobi

BASE = Path(__file__).resolve().parent
START = time.monotonic()
parser = argparse.ArgumentParser()
parser.add_argument('--output-dir', type=Path, default=BASE)
args = parser.parse_args()
OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)
SEED, COUNT, BLOCK = 20261026049, 128, 120
Q = 1.25


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(name, data):
    path = OUT/name
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+'\n')
    tmp.replace(path)


def checkpoint(stage, completed, rng_state=None, **extra):
    data = {'pid': os.getpid(), 'stage': stage, 'completed_shapes': completed,
            'total_registered_shapes': COUNT, 'elapsed_seconds_this_run': time.monotonic()-START,
            'resume_command': 'python -B diagnostic.py', 'seed': SEED,
            'is_certified_bound': False, **extra}
    if rng_state is not None:
        data['rng_state_after_completed'] = rng_state
    write_json('checkpoint.json', data)


class NumericalStop(Exception):
    def __init__(self, stage, detail):
        self.stage, self.detail = stage, detail
        super().__init__(stage)


def require(condition, stage, detail):
    if not condition:
        raise NumericalStop(stage, detail)


def finite_positive(x):
    return bool(np.all(np.isfinite(x)) and np.all(np.asarray(x) > 0))


def load_rows():
    path = OUT/'samples.jsonl'
    if not path.exists():
        return []
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    assert [r['sample_id'] for r in rows] == list(range(len(rows)))
    assert len(rows) <= COUNT
    return rows


if (OUT/'failure.json').exists():
    print('A recorded numerical failure exists. No new draw or retry was performed.')
    sys.exit(2)
rows = load_rows()
rng = np.random.Generator(np.random.PCG64(SEED))
if rows:
    rng.bit_generator.state = rows[-1]['rng_state_after']
registration_hash = sha(BASE/'registration.json')
script_hash = sha(Path(__file__))
if rows:
    assert all(r['registration_sha256'] == registration_hash for r in rows)
    assert all(r['script_sha256'] == script_hash for r in rows)

# The whole 1920-element quotient is fixed, including the identity at index 15.
permutation_list = list(permutations(range(5)))
sign_list = [(1,)+e for e in product((-1, 1), repeat=4)]
group_perm = np.array([p for p in permutation_list for s in sign_list], dtype=np.int64)
group_sign = np.array([s for p in permutation_list for s in sign_list], dtype=np.float64)
assert len(group_perm) == 1920
assert tuple(group_perm[15]) == (0, 1, 2, 3, 4)
assert tuple(group_sign[15]) == (1, 1, 1, 1, 1)
V = np.array([(1,)+e for e in product((-1, 1), repeat=4)], dtype=np.float64)
C = 6435*math.pi/256
WEIGHT_CAP = C/16
QUAD = {}
for order in (16, 8):
    nodes, weights = roots_jacobi(order, .25, .25)
    nodes = (nodes+1)/2
    weights = weights/np.sum(weights)
    assert finite_positive(nodes) and np.all(nodes < 1) and finite_positive(weights)
    QUAD[order] = nodes, weights
write_json('runtime.json', {'seed': SEED, 'numpy': np.__version__, 'scipy': scipy.__version__,
                          'python': sys.version, 'group_representatives': 1920,
                          'group_block': BLOCK, 'C': C, 'weight_cap': WEIGHT_CAP,
                          'quadrature': {str(n): {'nodes': x.tolist(), 'weights': w.tolist(),
                                                 'weight_sum': float(np.sum(w))}
                                         for n, (x, w) in QUAD.items()},
                          'registration_sha256': registration_hash, 'script_sha256': script_hash,
                          'label': 'noncertified_float64_diagnostic'})
checkpoint('ready', len(rows), rng.bit_generator.state)
print(json.dumps({'pid': os.getpid(), 'registered_shapes': COUNT, 'resuming_at': len(rows),
                  'registration_sha256': registration_hash}), flush=True)


def orbit_energy(A, logd, order):
    nodes, weights = QUAD[order]
    energy_sum = 0.0
    minimum_rho, maximum_rho = math.inf, -math.inf
    maximum_node_kernel = -math.inf
    minimum_log_ratio = math.inf
    identity_error = None
    min_rho_index = max_rho_index = None
    det_calls = 0
    for begin in range(0, 1920, BLOCK):
        pp, ss = group_perm[begin:begin+BLOCK], group_sign[begin:begin+BLOCK]
        transformed = A[pp[:, :, None], pp[:, None, :]]*ss[:, :, None]*ss[:, None, :]
        mixed = (nodes[None, :, None, None]*A[None, None, :, :]
                 +(1-nodes)[None, :, None, None]*transformed[:, None, :, :])
        sign, ld = np.linalg.slogdet(mixed)
        bad = np.argwhere((sign <= 0) | ~np.isfinite(ld))
        if len(bad):
            k, j = map(int, bad[0])
            raise NumericalStop('mixture_nonpositive_or_nonfinite_determinant',
                                {'order': order, 'group_index': begin+k, 'node_index': j,
                                 'determinant_sign': float(sign[k, j]),
                                 'matrix': mixed[k, j].tolist()})
        log_ratio = logd+ld
        require(bool(np.all(np.isfinite(log_ratio))), 'nonfinite_log_ratio', {'order': order, 'block': begin})
        kernel = np.exp(-.5*log_ratio)
        require(finite_positive(kernel), 'nonfinite_or_zero_kernel', {'order': order, 'block': begin})
        largest = float(np.max(kernel))
        require(largest <= 1+1e-5, 'node_kernel_exceeds_declared_tolerance',
                {'order': order, 'block': begin, 'maximum': largest})
        maximum_node_kernel = max(maximum_node_kernel, largest)
        minimum_log_ratio = min(minimum_log_ratio, float(np.min(log_ratio)))
        if begin <= 15 < begin+len(pp):
            identity_error = float(np.max(np.abs(kernel[15-begin]-1)))
            require(identity_error <= 1e-5, 'identity_kernel_accuracy_failure',
                    {'order': order, 'error': identity_error})
        rho = kernel@weights
        require(finite_positive(rho), 'nonfinite_or_nonpositive_rho', {'order': order, 'block': begin})
        lo, hi = int(np.argmin(rho)), int(np.argmax(rho))
        if float(rho[lo]) < minimum_rho:
            minimum_rho, min_rho_index = float(rho[lo]), begin+lo
        if float(rho[hi]) > maximum_rho:
            maximum_rho, max_rho_index = float(rho[hi]), begin+hi
        energy_sum += float(np.sum(rho))
        det_calls += len(pp)*order
    r = energy_sum/1920
    require(bool(math.isfinite(r) and 0 < r <= 1+1e-5), 'invalid_orbit_energy', {'order': order, 'r': r})
    return r, {'quadrature_order': order, 'determinant_calls': det_calls,
               'minimum_rho': minimum_rho, 'minimum_rho_group_index': min_rho_index,
               'maximum_rho': maximum_rho, 'maximum_rho_group_index': max_rho_index,
               'maximum_node_kernel': maximum_node_kernel,
               'minimum_normalized_log_determinant': minimum_log_ratio,
               'identity_kernel_maximum_error': identity_error}


def evaluate(pending):
    i = pending['sample_id']
    y = np.array(pending['gamma_draws'], dtype=np.float64)
    require(finite_positive(y), 'invalid_gamma_draw', {'gamma_draws': pending['gamma_draws']})
    radial_draw = float(np.sum(y))
    require(finite_positive(radial_draw), 'invalid_gamma_normalizer', {'value': radial_draw})
    z = y/radial_draw
    require(finite_positive(z), 'invalid_dirichlet_coordinate', {'Z': z.tolist()})
    reciprocal = 1/z
    require(finite_positive(reciprocal), 'invalid_dirichlet_reciprocal', {'Z': z.tolist()})
    H = float(np.sum(reciprocal))
    require(math.isfinite(H) and H >= 256*(1-1e-12), 'invalid_H', {'H': H})
    p = reciprocal/H
    require(finite_positive(p), 'invalid_shape_weight', {'P': p.tolist()})
    weight = C/math.sqrt(H)
    require(math.isfinite(weight) and 0 < weight <= WEIGHT_CAP*(1+1e-12), 'invalid_tilt_weight', {'weight': weight})
    M = V.T@(p[:, None]*V)
    eigenvalues = np.linalg.eigvalsh(M)
    sign, logd = np.linalg.slogdet(M)
    require(finite_positive(eigenvalues) and sign > 0 and math.isfinite(float(logd)),
            'invalid_shape_matrix', {'M': M.tolist(), 'eigenvalues': eigenvalues.tolist(),
                                     'determinant_sign': float(sign)})
    condition = float(eigenvalues[-1]/eigenvalues[0])
    require(condition < 1e12, 'shape_condition_limit', {'condition_2': condition, 'M': M.tolist()})
    A = np.linalg.inv(M)
    residual = float(np.linalg.norm(M@A-np.eye(5), ord=np.inf))
    asymmetry = float(np.linalg.norm(A-A.T, ord=np.inf)/np.linalg.norm(A, ord=np.inf))
    require(np.all(np.isfinite(A)) and residual <= 1e-5 and asymmetry <= 1e-10,
            'inverse_accuracy_failure', {'inverse_residual_inf': residual,
                                         'inverse_relative_asymmetry': asymmetry, 'M': M.tolist()})
    require(finite_positive(np.linalg.eigvalsh((A+A.T)/2)), 'nonpositive_inverse_symmetric_part', {'M': M.tolist()})
    r16, diag16 = orbit_energy(A, float(logd), 16)
    crosscheck = None
    if i < 4:
        r8, diag8 = orbit_energy(A, float(logd), 8)
        crosscheck = {'r8': r8, 'r16': r16, 'signed_difference_r16_minus_r8': r16-r8,
                      'relative_difference_r16_minus_r8_over_r16': (r16-r8)/r16,
                      'diagnostics_8': diag8, 'is_error_certificate': False}
    det_factor = math.exp(-.1*float(logd))
    base = weight*det_factor
    orbit = base*r16**.2
    require(finite_positive([base, orbit, det_factor]), 'invalid_final_contribution',
            {'base': base, 'orbit': orbit, 'det_factor': det_factor})
    return {'sample_id': i, 'gamma_draws': pending['gamma_draws'],
            'rng_state_after': pending['rng_state_after'],
            'Z': z.tolist(), 'P': p.tolist(), 'M': M.tolist(),
            'gamma_normalizer_unused_in_weight': radial_draw,
            'H': H, 'tilt_weight': weight, 'P_sum': float(np.sum(p)),
            'determinant': math.exp(float(logd)), 'log_determinant': float(logd),
            'eigenvalues': eigenvalues.tolist(), 'condition_2': condition,
            'inverse_residual_inf': residual, 'inverse_relative_asymmetry': asymmetry,
            'det_negative_tenth': det_factor, 'orbit_r_16': r16,
            'weighted_unimproved_contribution': base, 'weighted_orbit_contribution': orbit,
            'orbit_diagnostics_16': diag16, 'crosscheck_8': crosscheck,
            'registration_sha256': registration_hash, 'script_sha256': script_hash,
            'is_certified': False}


for sample_id in range(len(rows), COUNT):
    pending_path = OUT/'pending.json'
    if pending_path.exists():
        pending = json.loads(pending_path.read_text())
        require(pending['sample_id'] == sample_id, 'resume_pending_mismatch', {'pending_id': pending['sample_id'], 'expected_id': sample_id})
        rng.bit_generator.state = pending['rng_state_after']
    else:
        y = rng.gamma(shape=.5, scale=1.0, size=16)
        pending = {'sample_id': sample_id, 'gamma_draws': y.tolist(),
                   'rng_state_after': rng.bit_generator.state, 'seed': SEED,
                   'registration_sha256': registration_hash, 'script_sha256': script_hash}
        write_json('pending.json', pending)
    try:
        with np.errstate(over='raise', divide='raise', invalid='raise', under='raise'):
            row = evaluate(pending)
    except (NumericalStop, FloatingPointError, np.linalg.LinAlgError, OverflowError) as exc:
        details = {'stage': exc.stage, 'detail': exc.detail} if isinstance(exc, NumericalStop) else {'stage': type(exc).__name__, 'detail': str(exc)}
        failure = {'status': 'stopped_without_clipping_or_replacement', 'pid': os.getpid(),
                   'sample_id': sample_id, 'completed_shapes': len(rows),
                   'pending_input': pending, **details,
                   'elapsed_seconds': time.monotonic()-START,
                   'registration_sha256': registration_hash, 'script_sha256': script_hash}
        write_json('failure.json', failure)
        checkpoint('numerical_failure', len(rows), rows[-1]['rng_state_after'] if rows else None,
                   failing_sample=sample_id, failure_stage=details['stage'])
        print(json.dumps({'status': 'numerical_failure', 'sample_id': sample_id, **details}), flush=True)
        sys.exit(2)
    with (OUT/'samples.jsonl').open('a') as handle:
        handle.write(json.dumps(row, sort_keys=True, allow_nan=False)+'\n')
        handle.flush()
        os.fsync(handle.fileno())
    rows.append(row)
    pending_path.unlink()
    checkpoint('sampling', len(rows), rng.bit_generator.state)
    if len(rows) % 16 == 0:
        print(json.dumps({'completed': len(rows), 'elapsed_seconds': time.monotonic()-START,
                          'last_condition': row['condition_2']}), flush=True)
    if time.monotonic()-START >= 110 and len(rows) < COUNT:
        checkpoint('recoverable_time_yield', len(rows), rng.bit_generator.state)
        print('Checkpoint saved. Resume exactly the same registered batch.', flush=True)
        sys.exit(75)


def describe(values):
    x = np.array(values, dtype=np.float64)
    return {'mean': float(np.mean(x)), 'sample_sd': float(np.std(x, ddof=1)),
            'descriptive_mc_se': float(np.std(x, ddof=1)/math.sqrt(len(x))),
            'minimum': float(np.min(x)), 'maximum': float(np.max(x)),
            'quantiles_0_25_50_75_100': np.quantile(x, [0, .25, .5, .75, 1]).tolist(),
            'is_confidence_certificate': False}


weights = [r['tilt_weight'] for r in rows]
base_stats = describe([r['weighted_unimproved_contribution'] for r in rows])
orbit_stats = describe([r['weighted_orbit_contribution'] for r in rows])
c1, face = 1/(32*math.pi**3), math.pi**2/60
I_base = c1**Q*face*base_stats['mean']**Q
I_orbit = c1**Q*face*orbit_stats['mean']**Q
target = 5.36e-5
target_g1e4 = 32*(1e-4/5)**Q
summary = {
    'status': 'completed_128_noncertified_shapes', 'seed': SEED, 'sample_count': COUNT,
    'weight_statistics': describe(weights), 'weight_exact_expectation': 1,
    'weight_sample_ess': float(sum(weights)**2/sum(w*w for w in weights)),
    'unimproved_weighted_det_negative_tenth': base_stats,
    'orbit_weighted_det_negative_tenth_r_fifth': orbit_stats,
    'r16_statistics': describe([r['orbit_r_16'] for r in rows]),
    'condition_statistics': describe([r['condition_2'] for r in rows]),
    'maximum_inverse_residual_inf': max(r['inverse_residual_inf'] for r in rows),
    'maximum_inverse_relative_asymmetry': max(r['inverse_relative_asymmetry'] for r in rows),
    'maximum_identity_kernel_error': max(r['orbit_diagnostics_16']['identity_kernel_maximum_error'] for r in rows),
    'maximum_node_kernel_unclipped': max(r['orbit_diagnostics_16']['maximum_node_kernel'] for r in rows),
    'first_four_8_16_crosschecks': [{'sample_id': r['sample_id'], 'condition_2': r['condition_2'], **r['crosscheck_8']}
                                 for r in rows[:4]],
    'noncertified_plugin_I': {'unimproved': I_base, 'full_orbit': I_orbit,
                             'orbit_to_unimproved_ratio': I_orbit/I_base,
                             'target_5_36e_minus5': target, 'orbit_minus_target': I_orbit-target,
                             'target_for_g_1e_minus4_and_cost_32': target_g1e4,
                             'orbit_minus_g_1e_minus4_target': I_orbit-target_g1e4,
                             'is_bound': False,
                             'meaning': 'Plug-in for the orbit-Minkowski upper-functional, not the true face integral and not a certified expectation.'},
    'determinant_calls_main': sum(r['orbit_diagnostics_16']['determinant_calls'] for r in rows),
    'determinant_calls_crosscheck': sum(r['crosscheck_8']['diagnostics_8']['determinant_calls'] for r in rows[:4]),
    'no_clipping_resampling_or_extra_seed': True,
    'new_strict_I_bound': False, 'new_C5_bound': False,
    'registration_sha256': registration_hash, 'script_sha256': script_hash,
    'samples_sha256': sha(OUT/'samples.jsonl')
}
assert summary['determinant_calls_main'] == 3932160
assert summary['determinant_calls_crosscheck'] == 61440
write_json('results.json', summary)
checkpoint('complete', len(rows), rng.bit_generator.state)
receipt = {'pid': os.getpid(), 'elapsed_seconds_this_run': time.monotonic()-START,
           'exit_status': 0, 'sample_count': len(rows), 'numpy': np.__version__, 'scipy': scipy.__version__,
           'results_sha256': sha(OUT/'results.json'), 'samples_sha256': sha(OUT/'samples.jsonl'),
           'registration_sha256': registration_hash, 'script_sha256': script_hash}
write_json('run_receipt.json', receipt)
print(json.dumps(receipt, sort_keys=True), flush=True)
print(json.dumps({'weight_mean': summary['weight_statistics']['mean'], 'J_base': base_stats['mean'],
                  'J_orbit': orbit_stats['mean'], 'I_plugin': I_orbit,
                  'I_plugin_base': I_base, 'max_condition': summary['condition_statistics']['maximum'],
                  'is_certified_bound': False}, sort_keys=True), flush=True)
