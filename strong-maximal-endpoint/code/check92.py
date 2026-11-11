"""Replay r092 scalar diagnostics and compare independent Wolfram results.

Numerical rows are diagnostics, not interval certificates or theorem claims.
Exact source replay is used for the Python rows; cross-software comparisons use
relative error without a unit-scale floor, including tiny tail probabilities.
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import re
import sys
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
from pathlib import Path

TOL = D('1e-45')
K_VALUES = (16, 32, 64, 128, 256, 512, 1024, 2048, 4096)
B_VALUES = ('1/4', '1', '4', '16')
PRECISIONS = (80, 160)
D_VALUES = (0, 1, 4)
COMMON_FIELDS = ('L', 'mu', 'sigma', 'z', 'tau', 'rank_shift',
                 'rank_shift_over_B')
FIELD_MAP = {
    'zB': 'z_B',
    'target_tail_tau_plus_BL_over_N': 'target_CDF',
    **{name: name for name in COMMON_FIELDS},
}
LINEAR_FIELD_MAP = {
    'stoploss_correction': 'stoploss_correction',
    'normalized_excess_error_over_L': 'normalized_excess_error_over_L',
    'second_moment_over_L_scale': 'curvature_scale',
}
REPLAY_FILES = (
    'python/experiments92.json', 'python/resources92.json',
    'wolfram/verification92.json', 'wolfram/resources92.json',
)
OWN_SOURCES = ('experiment92.py', 'verify92.wl', 'check92.py')
DEPENDENCIES = ('math_guard.py', 'verify85.py', 'verify84.py', 'sliver79.py')
PRIOR_FILES = (
    'rounds/r011/report.md', 'rounds/r021/report.md',
    'rounds/r022/report.md', 'rounds/r029/report.md',
    'rounds/r044/report.md', 'rounds/r089/core89.md',
    'rounds/r090/core90.md', 'rounds/r091/core91.md',
)


def dec(value):
    text = re.sub(r'`[^*]*', '', str(value)).replace('*^', 'E')
    if '/' in text and not any(character in text for character in '.Ee'):
        fraction = Q(text)
        with localcontext() as context:
            context.prec = 180
            return D(fraction.numerator)/D(fraction.denominator)
    result = D(text)
    assert result.is_finite(), value
    return result


def relative_error(a, b):
    """True relative difference, with no max(1, ...) near-zero masking."""
    x, y = dec(a), dec(b)
    scale = max(abs(x), abs(y))
    if scale == 0:
        return D(0)
    return abs(x-y)/scale


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resource(row):
    assert row['exit_code'] == 0 and row['stop_reason'] is None
    reserve = D(str(row['reserve_GiB']))
    budget = D(str(row['owned_RSS_budget_GiB']))
    peak = D(str(row['peak_owned_RSS_MiB']))
    available = D(str(row['min_system_available_MiB']))
    assert reserve >= 1
    assert available >= reserve*1024
    assert peak <= budget*1024
    assert row['samples'] > 0 and 'estimate' in row['scope'].lower()
    return row


def without_pid(rows):
    return [{key: value for key, value in row.items() if key != 'pid'}
            for row in rows]


def expected_target(tau, B, L, N):
    with localcontext() as ctx:
        ctx.prec = 180
        return dec(tau) + dec(B)*dec(L)/D(N)


def validate_scalar_rows(py_rows, wf_rows):
    expected_keys = {(k, b) for k in K_VALUES for b in B_VALUES}
    py_map = {(row['k'], row['B']): row for row in py_rows}
    wf_map = {(row['k'], row['B']): row for row in wf_rows}
    assert set(py_map) == expected_keys and set(wf_map) == expected_keys
    assert len(py_rows) == len(wf_rows) == 36
    max_cross = D(0)
    max_within = D(0)
    max_by_field = {}
    max_residual = D(0)

    # Residuals measure solver stopping error and can differ by many orders
    # between independent solvers; each is checked against its target instead
    # of comparing the two residual magnitudes as if they were observables.
    for key in sorted(expected_keys, key=lambda item: (item[0], B_VALUES.index(item[1]))):
        py_case, wf_case = py_map[key], wf_map[key]
        k, b_string = key
        assert py_case['N'] == wf_case.get('N', str(1 << k))
        assert int(wf_case['kernel_id']) > 0
        assert [r['working_precision'] for r in py_case['precision_rows']] == list(PRECISIONS)
        assert [r['working_precision'] for r in wf_case['precision_rows']] == list(PRECISIONS)
        assert len(py_case['precision_rows']) == len(wf_case['precision_rows']) == 2

        for field in FIELD_MAP:
            values = [row[field] for row in py_case['precision_rows']]
            error = relative_error(values[0], values[1])
            assert error < TOL, (key, field, error)
            max_within = max(max_within, error)

        for field, wf_field in FIELD_MAP.items():
            values_py = [row[field] for row in py_case['precision_rows']]
            values_wf = [row[wf_field] for row in wf_case['precision_rows']]
            for values in (values_py, values_wf):
                error = relative_error(values[0], values[1])
                assert error < TOL, (key, field, 'within precision', error)
                max_within = max(max_within, error)
            for py_row, wf_row in zip(py_case['precision_rows'], wf_case['precision_rows']):
                error = relative_error(py_row[field], wf_row[wf_field])
                assert error < TOL, (key, field, 'cross software', error)
                max_cross = max(max_cross, error)
                max_by_field[field] = max(max_by_field.get(field, D(0)), error)

        for py_row, wf_row in zip(py_case['precision_rows'], wf_case['precision_rows']):
            wp = py_row['working_precision']
            target = D(10) ** (-wp)
            assert dec(py_row['z_log_residual']) < target
            assert dec(py_row['zB_log_residual']) < target
            assert dec(py_row['root_log_residual_target']) == target
            assert dec(wf_row['log_tail_residual']) < target
            max_residual = max(max_residual, dec(py_row['z_log_residual']),
                               dec(py_row['zB_log_residual']),
                               dec(wf_row['log_tail_residual']))

            B = dec(b_string)
            N = 1 << k
            for row, prefix in ((py_row, 'python'), (wf_row, 'wolfram')):
                z = dec(row['z'])
                zb = dec(row['zB'] if prefix == 'python' else row['z_B'])
                tau = dec(row['tau'])
                target_prob = dec(row['target_tail_tau_plus_BL_over_N']
                                  if prefix == 'python' else row['target_CDF'])
                rank = dec(row['rank_shift'])
                rank_over_b = dec(row['rank_shift_over_B'])
                L = dec(row['L'])
                assert N == 1 << k and L > 0 and z > 0 and 0 < zb <= z
                assert 0 < tau <= target_prob < D('0.5')
                assert rank > 0 and rank_over_b > 0
                assert relative_error(target_prob, expected_target(tau, b_string, L, N)) < TOL
                assert relative_error(rank_over_b, rank/B) < TOL

            assert wf_row['pass'] is True
            assert len(py_row['linear_positive_part_rows']) == 3
            assert len(wf_row['D_rows']) == 3
            py_linear = {row['D']: row for row in py_row['linear_positive_part_rows']}
            wf_linear = {row['D']: row for row in wf_row['D_rows']}
            assert set(py_linear) == set(wf_linear) == set(D_VALUES)
            for d_value in D_VALUES:
                pl, wl = py_linear[d_value], wf_linear[d_value]
                for field, wf_field in LINEAR_FIELD_MAP.items():
                    py_values = [row['linear_positive_part_rows'][i][field]
                                 for row in py_case['precision_rows']
                                 for i in range(3)
                                 if row['linear_positive_part_rows'][i]['D'] == d_value]
                    wf_values = [row['D_rows'][i][wf_field]
                                 for row in wf_case['precision_rows']
                                 for i in range(3)
                                 if row['D_rows'][i]['D'] == d_value]
                    for values, software in ((py_values, 'python'), (wf_values, 'wolfram')):
                        within_error = relative_error(values[0], values[1])
                        assert within_error < TOL, (key, d_value, field, software,
                                                    'within precision', within_error)
                        max_within = max(max_within, within_error)
                    cross_error = relative_error(pl[field], wl[wf_field])
                    assert cross_error < TOL, (key, wp, d_value, field, cross_error)
                    max_cross = max(max_cross, cross_error)
                    max_by_field[field] = max(max_by_field.get(field, D(0)), cross_error)
                # Check positive-part outputs and their defining identities.
                for lin in (pl, wl):
                    a = dec(lin['a']) if 'a' in lin else None
                    phi = dec(lin['phi_a']) if 'phi_a' in lin else None
                    tail = dec(lin['barPhi_a']) if 'barPhi_a' in lin else None
                    psi = dec(lin['psi_a']) if 'psi_a' in lin else None
                    correction = dec(lin['stoploss_correction'])
                    excess = dec(lin['normalized_excess_error_over_L'])
                    curvature = dec(lin.get('second_moment_over_L_scale', lin.get('curvature_scale')))
                    assert correction > 0 and excess > 0 and curvature > 0
                    if a is not None:
                        assert a > 0 and phi > 0 and 0 < tail < 1 and psi > 0
                        psi_residual = abs(phi-a*tail-psi)/max(abs(phi), abs(a*tail), abs(psi))
                        assert psi_residual < D('1e-55'), (key, wp, d_value, 'psi', psi_residual)
                    assert relative_error(excess, dec(d_value)+correction/dec(py_row['L'])) < D('1e-55')

    # For each k: larger B raises the requested right-tail probability, lowers
    # its quantile z_B, and increases the induced rank shift.
    monotone_checks = 0
    for k in K_VALUES:
        ordered = [py_map[k, b]['precision_rows'][-1] for b in B_VALUES]
        targets = [dec(row['target_tail_tau_plus_BL_over_N']) for row in ordered]
        quantiles = [dec(row['zB']) for row in ordered]
        shifts = [dec(row['rank_shift']) for row in ordered]
        # B effects in target and z_B can be below the 65-digit serialization
        # resolution at large k. Permit ties there; rank_shift is separately
        # formed at guarded precision and must retain the strict ordering.
        assert all(a <= b for a, b in zip(targets, targets[1:]))
        assert all(a >= b for a, b in zip(quantiles, quantiles[1:]))
        assert all(a < b for a, b in zip(shifts, shifts[1:]))
        monotone_checks += 3

    return {
        'max_cross_relative_difference': str(max_cross),
        'max_within_precision_relative_difference': str(max_within),
        'max_relative_difference_by_field': {key: str(value) for key, value in max_by_field.items()},
        'max_solver_residual': str(max_residual),
        'monotonicity_checks': monotone_checks,
        'root_jobs': 36,
        'precision_rows': 72,
        'positive_part_rows': 216,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('output')
    parser.add_argument('--project-root', required=True)
    args = parser.parse_args()
    data_dir = Path(args.directory)
    project_root = Path(args.project_root)
    code_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(project_root/'code'))
    sys.path.insert(0, str(code_dir))

    from experiment92 import B_VALUES as EXP_B_VALUES, K_VALUES as EXP_K_VALUES, job
    from verify85 import constants85

    py, py_resource, wf, wf_resource = [json.loads((data_dir/name).read_text())
                                        for name in REPLAY_FILES]
    assert py['pass_all'] is True and wf['pass'] is True
    assert tuple(EXP_K_VALUES) == K_VALUES and tuple(EXP_B_VALUES) == B_VALUES
    assert py['working_precision'] == list(PRECISIONS)
    assert py['configured_workers'] == 12 and 1 <= py['actual_workers'] <= 12
    assert 1 <= len(set(py['worker_pids'])) <= 12
    assert py['actual_workers'] == len(set(py['worker_pids']))
    assert py['counts'] == {
        'k_values': 9, 'B_values': 4, 'root_jobs': 36,
        'precision_rows': 72, 'linear_positive_part_cases': 216,
        'high_dimensional_arrays_allocated': 0,
    }
    assert len(wf['symbolic']) == 12 and all(wf['symbolic'].values())
    assert wf['kernels_closed'] is True
    assert 1 <= wf['configured_workers'] <= 12
    assert 0 < wf['per_case_memory_limit_MiB'] <= 1024
    py_resource = resource(py_resource)
    wf_resource = resource(wf_resource)
    mu_interval = constants85()['mu_interval']
    assert py['constants']['mu_interval'] == mu_interval

    # Re-execute the actual Python worker function with spawn semantics, then
    # compare all content except process IDs exactly.
    replay_jobs = [(k, b, mu_interval) for k in K_VALUES for b in B_VALUES]
    with mp.get_context('spawn').Pool(12) as pool:
        replay = pool.map(job, replay_jobs)
    assert without_pid(py['rows']) == without_pid(replay)

    wf_rows = wf['rows']
    assert len(wf_rows) == 36
    ids = {row['kernel_id'] for row in wf_rows}
    assert all(isinstance(value, int) and value > 0 for value in ids)
    assert len(ids) <= wf['configured_workers']
    stats = validate_scalar_rows(py['rows'], wf_rows)

    own_paths = {name: code_dir/name for name in OWN_SOURCES}
    dep_paths = {name: project_root/'code'/name for name in DEPENDENCIES}
    for name, path in dep_paths.items():
        if not path.is_file() and (code_dir/name).is_file():
            dep_paths[name] = code_dir/name
    source_hashes = {name: sha(path) for name, path in {**own_paths, **dep_paths}.items()}
    input_hashes = {name: sha(data_dir/name) for name in REPLAY_FILES}
    docs = ('core92.md', 'audit92.md')
    document_hashes = {name: sha(data_dir/name) for name in docs}
    prior_hashes = {name: sha(project_root/name) for name in PRIOR_FILES}

    result = {
        'pass_all': True,
        'scope': ('Replayed deterministic scalar rank-shift and positive-part diagnostics; '
                  'does not authenticate the unknown r021 constant, infer a physical N0, '
                  'or establish an unknown C_H.'),
        'counts': py['counts'],
        'configured_python_workers': 12,
        'actual_stored_python_workers': py['actual_workers'],
        'replay_python_workers': len({row['pid'] for row in replay}),
        'configured_wolfram_workers': wf['configured_workers'],
        'wolfram_kernel_ids': sorted(ids),
        'wolfram_kernels_closed': wf['kernels_closed'],
        'symbolic_gates': wf['symbolic'],
        'working_precision': list(PRECISIONS),
        'relative_tolerance': str(TOL),
        **stats,
        'Python_resource_monitor': py_resource,
        'Wolfram_resource_monitor': wf_resource,
        'resources_are_sampled_estimates_not_OS_reservations': True,
        'source_sha256': source_hashes,
        'input_sha256': input_hashes,
        'document_sha256': document_hashes,
        'prior_dependency_sha256': prior_hashes,
        'unknown_r021_constant_authenticated': False,
        'physical_N0_inferred': False,
        'unknown_C_H_established': False,
        'new_general_C_N_coefficient': False,
        'goal_status': 'active',
        'pause_request_revoked': True,
    }
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items()
                      if 'sha256' not in key and key not in ('symbolic_gates',
                                                              'Python_resource_monitor',
                                                              'Wolfram_resource_monitor')},
                     indent=2))


if __name__ == '__main__':
    with localcontext() as context:
        context.prec = 180
        main()
