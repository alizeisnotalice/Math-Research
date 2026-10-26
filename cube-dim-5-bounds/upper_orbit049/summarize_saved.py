#!/usr/bin/env python3
"""Read-only summary of the saved stopped prefix. Draws no shapes; computes no orbit."""
from pathlib import Path
import hashlib
import json
import math
import os
import statistics
import time

BASE = Path(__file__).resolve().parent
START = time.monotonic()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stats(values):
    n = len(values)
    return {'mean': statistics.fmean(values), 'sample_sd': statistics.stdev(values),
            'descriptive_se': statistics.stdev(values)/math.sqrt(n),
            'minimum': min(values), 'maximum': max(values),
            'is_confidence_certificate': False,
            'warning': 'This is a data-dependent stopped prefix, not the completed registered 128-shape batch.'}


rows = [json.loads(x) for x in (BASE/'samples.jsonl').read_text().splitlines() if x.strip()]
failure = json.loads((BASE/'failure.json').read_text())
runtime = json.loads((BASE/'runtime.json').read_text())
assert len(rows) == 99
assert [r['sample_id'] for r in rows] == list(range(99))
assert failure['sample_id'] == 99 and failure['completed_shapes'] == 99
weights = [r['tilt_weight'] for r in rows]
base = stats([r['weighted_unimproved_contribution'] for r in rows])
orbit = stats([r['weighted_orbit_contribution'] for r in rows])
q = 1.25
c = 1/(32*math.pi**3)
face = math.pi**2/60
I_base = c**q*face*base['mean']**q
I_orbit = c**q*face*orbit['mean']**q
target = 5.36e-5
target_g = 32*(1e-4/5)**q
out = {
    'status': 'stopped_prefix_noncertified_diagnostic',
    'registered_seed': 20261026049, 'planned_shapes': 128,
    'completed_shapes': 99, 'failed_shape_id': 99,
    'drawn_shapes_including_failure': 100, 'not_drawn_shapes': 28,
    'failure_stage': failure['stage'], 'failure_detail': failure['detail'],
    'actual_batch_pid': failure['pid'], 'actual_batch_elapsed_seconds': failure['elapsed_seconds'],
    'C': runtime['C'], 'weight_cap': runtime['weight_cap'],
    'weight_statistics': stats(weights), 'exact_weight_expectation': 1,
    'weight_sample_ess': sum(weights)**2/sum(w*w for w in weights),
    'weighted_unimproved_J': base, 'weighted_orbit_J': orbit,
    'orbit_r16_statistics': stats([r['orbit_r_16'] for r in rows]),
    'condition_statistics': stats([r['condition_2'] for r in rows]),
    'maximum_completed_condition_sample_id': max(rows, key=lambda r:r['condition_2'])['sample_id'],
    'maximum_completed_inverse_residual': max(r['inverse_residual_inf'] for r in rows),
    'maximum_completed_inverse_asymmetry': max(r['inverse_relative_asymmetry'] for r in rows),
    'maximum_completed_identity_kernel_error': max(r['orbit_diagnostics_16']['identity_kernel_maximum_error'] for r in rows),
    'maximum_completed_node_kernel_unclipped': max(r['orbit_diagnostics_16']['maximum_node_kernel'] for r in rows),
    'first_four_crosschecks': [{'sample_id': r['sample_id'], 'condition': r['condition_2'], **r['crosscheck_8']}
                             for r in rows[:4]],
    'noncertified_prefix_plugin_I': {
        'unimproved': I_base, 'full_orbit': I_orbit, 'orbit_to_unimproved_ratio': I_orbit/I_base,
        'target_5_36e_minus5': target, 'orbit_minus_target': I_orbit-target,
        'target_g_1e_minus4_cost32': target_g, 'orbit_minus_g_target': I_orbit-target_g,
        'is_bound': False,
        'interpretation': 'Stopped-prefix plug-in of an upper-functional, not true I, not a complete expectation estimate, not a C5 certificate.'
    },
    'completed_determinant_calls_main': sum(r['orbit_diagnostics_16']['determinant_calls'] for r in rows),
    'completed_determinant_calls_crosscheck': sum(r['crosscheck_8']['diagnostics_8']['determinant_calls'] for r in rows[:4]),
    'no_clipping_resampling_or_threshold_change': True,
    'completed_registered_128': False, 'new_strict_I_bound': False, 'new_C5_bound': False,
    'registration_sha256': sha(BASE/'registration.json'),
    'diagnostic_script_sha256': sha(BASE/'diagnostic.py'),
    'samples_sha256': sha(BASE/'samples.jsonl'), 'failure_sha256': sha(BASE/'failure.json'),
    'pending_sha256': sha(BASE/'pending.json')
}
path = BASE/'partial_results.json'
path.write_text(json.dumps(out, indent=2, sort_keys=True, allow_nan=False)+'\n')
receipt = {'summary_pid': os.getpid(), 'summary_elapsed_seconds': time.monotonic()-START,
           'actual_batch_pid': failure['pid'], 'batch_exit_status': 2, 'summary_exit_status': 0,
           'partial_results_sha256': sha(path), 'summary_script_sha256': sha(Path(__file__)),
           'sampling_or_orbit_work_performed_by_summary': False}
(BASE/'summary_receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
print(json.dumps(receipt, sort_keys=True))
print(json.dumps({'weight_mean': out['weight_statistics']['mean'],
                  'J_base': base['mean'], 'J_orbit': orbit['mean'],
                  'I_plugin_base': I_base, 'I_plugin_orbit': I_orbit,
                  'target_g_1e_minus4': target_g,
                  'maximum_completed_condition': out['condition_statistics']['maximum'],
                  'max_crosscheck_relative_difference': max(abs(r['crosscheck_8']['relative_difference_r16_minus_r8_over_r16']) for r in rows[:4]),
                  'completed_shapes': 99, 'is_certified': False}, sort_keys=True))
