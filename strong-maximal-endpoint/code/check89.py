"""Independent acceptance checker for r089 scalar numerical experiments."""
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


def dec(x):
    if x is None:
        raise AssertionError('unexpected null numeric value')
    return D(re.sub(r'`[^*]*', '', str(x)).replace('*^', 'E'))


def finite(x):
    z = dec(x)
    assert z.is_finite(), x
    return z


def relative_error(a, b):
    a, b = finite(a), finite(b)
    if a == 0 and b == 0:
        return D(0)
    den = max(abs(a), abs(b))
    assert den > 0
    return abs(a-b)/den


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay_exact(stored, fresh):
    assert len(stored) == len(fresh)
    for old, new in zip(stored, fresh):
        assert {k: v for k, v in old.items() if k != 'pid'} == {
            k: v for k, v in new.items() if k != 'pid'
        }


def assert_precision_rows(rows):
    assert [r['working_precision'] for r in rows] == [80, 160]


def compare_fields(py_rows, wf_rows, fields, rels):
    assert_precision_rows(py_rows)
    assert [r['working_precision'] for r in wf_rows] == [80, 160]
    by_field = {f: [] for f in fields}
    for swrows in (py_rows, wf_rows):
        for field in fields:
            vals = [finite(r[field]) for r in swrows]
            err = relative_error(*vals)
            assert err < TOL, (field, 'precision', err)
            by_field[field].append(err)
    for i in range(2):
        for field in fields:
            err = relative_error(py_rows[i][field], wf_rows[i][field])
            assert err < TOL, (field, 'cross-software', i, err)
            rels.append(err)
    return {f: max(v) for f, v in by_field.items()}


def aggregate_max(target, values):
    for field, value in values.items():
        target[field] = max(target.get(field, D(0)), value)


def row_key(row, kind):
    if kind == 'root':
        return row['k'], Q(row['K_over_sigma'])
    return row['k'], row['type' if kind == 'quantile' else 'mode']


def validate_resource(z):
    assert z['exit_code'] == 0 and z['stop_reason'] is None
    reserve = D(str(z['reserve_GiB'])) * D(1024)
    budget = D(str(z['owned_RSS_budget_GiB'])) * D(1024)
    assert D(str(z['min_system_available_MiB'])) >= reserve
    assert D(str(z['peak_owned_RSS_MiB'])) <= budget
    assert 'estimate' in z['scope'].lower()
    return {
        'exit_code': z['exit_code'], 'stop_reason': z['stop_reason'],
        'reserve_GiB': z['reserve_GiB'], 'owned_RSS_budget_GiB': z['owned_RSS_budget_GiB'],
        'peak_owned_RSS_MiB': z['peak_owned_RSS_MiB'],
        'min_system_available_MiB': z['min_system_available_MiB'],
        'samples': z['samples'], 'seconds': z['seconds'],
        'scope': z['scope'], 'available_memory_is_estimate': True,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory')
    ap.add_argument('output')
    ap.add_argument('--project-root', required=True)
    a = ap.parse_args()
    d, project_root = Path(a.directory), Path(a.project_root)
    code = Path(__file__).resolve().parent
    input_names = [
        'python/experiments89.json', 'python/resources89.json',
        'wolfram/verification89.json', 'wolfram/resources89.json',
    ]
    for name in input_names:
        assert (d/name).is_file(), d/name
    py = json.loads((d/'python/experiments89.json').read_text())
    pyres = json.loads((d/'python/resources89.json').read_text())
    wf = json.loads((d/'wolfram/verification89.json').read_text())
    wfres = json.loads((d/'wolfram/resources89.json').read_text())

    sys.path.insert(0, str(code))
    sys.path.insert(0, str(project_root/'code'))
    from experiment89 import clipping_job, quantile_job, root_job

    assert py['pass_all'] is True and wf['pass'] is True
    assert len(wf['symbolic']) == 12 and all(v is True for v in wf['symbolic'].values())
    assert wf['kernels_closed'] is True
    assert 0 < wf['per_case_memory_limit_MiB'] <= 1024
    assert py['configured_workers'] == 12 and py['working_precision'] == [80, 160]
    assert 1 <= wf['configured_workers'] <= 12
    assert py['counts'] == {'root_jobs': 35, 'quantile_jobs': 28,
                            'clipping_jobs': 49, 'total_jobs': 112}
    py_resource_summary = validate_resource(pyres)
    wf_resource_summary = validate_resource(wfres)

    root_py, quant_py, clip_py = py['root_rows'], py['quantile_rows'], py['clipping_rows']
    root_wf, quant_wf, clip_wf = wf['root_rows'], wf['quantile_rows'], wf['clipping_rows']
    assert (len(root_py), len(quant_py), len(clip_py)) == (35, 28, 49)
    assert (len(root_wf), len(quant_wf), len(clip_wf)) == (35, 28, 49)
    expected_k = [16, 32, 64, 128, 256, 512, 1024]
    for rows, label in ((root_py, 'root'), (quant_py, 'quantile'), (clip_py, 'clipping')):
        assert all(r['k'] in expected_k and r['N'] == str(1 << r['k']) for r in rows), label
    for rows, label in ((root_wf, 'root'), (quant_wf, 'quantile'), (clip_wf, 'clipping')):
        assert all(r['k'] in expected_k for r in rows), label

    def all_positive_ids(rows):
        ids = {r['kernel_id'] for r in rows}
        assert ids and all(isinstance(x, int) and x > 0 for x in ids)
        assert len(ids) <= wf['configured_workers']
        return ids

    kernel_ids = all_positive_ids(root_wf+quant_wf+clip_wf)

    # Reconstruct all three Python parameter families with a fresh spawn pool.
    root_jobs = [(k, ratio) for k in expected_k for ratio in ('1', '1.25', '1.5', '1.75', '1.9')]
    quant_jobs = [(k, kind) for k in expected_k
                  for kind in ('1/N', '1/sqrt(NL)', 'K=.5sigma', 'K=.75sigma')]
    clip_jobs = [(k, mode) for k in expected_k for mode in
                 ('P=.25', 'P=1', 'P=4', 'critical', 'quadratic_c=.1',
                  'quadratic_c=.2', 'quadratic_c=.3')]
    with mp.get_context('spawn').Pool(12) as pool:
        root_replay = pool.map(root_job, root_jobs)
        quant_replay = pool.map(quantile_job, quant_jobs)
        clip_replay = pool.map(clipping_job, clip_jobs)
    replay_exact(root_py, root_replay)
    replay_exact(quant_py, quant_replay)
    replay_exact(clip_py, clip_replay)
    replay_pids = {r['pid'] for r in root_replay+quant_replay+clip_replay}

    root_fields = ('z', 'w', 'log_D', 'normalized_root')
    quant_fields = ('epsilon', 'log_epsilon', 'y', 'z', 'standardized_shift',
                    'gaussian_log_loss', 'gaussian_log_loss_over_W')
    clip_fields = ('p', 'b_parameter', 'cap_excess', 'log_s', 'log_a', 'ln_Phi_a',
                   'C0', 'v', 'b_p', 'log_high_fee_minus_q0', 'normalized_gap')
    relative_by_family = {'root': [], 'quantile': [], 'clipping': []}
    root_residual_max = quant_residual_max = clip_residual_max = D(0)
    tiny = {'python_omitted': 0, 'wolfram_omitted': 0,
            'both_omitted': 0, 'threshold_log1p_approximated': 0}
    root_rows_checked = quantile_rows_checked = clipping_rows_checked = 0
    summary_max = {'root': {}, 'quantile': {}, 'clipping': {}}

    wf_maps = {
        'root': {row_key(x, 'root'): x for x in root_wf},
        'quantile': {row_key(x, 'quantile'): x for x in quant_wf},
        'clipping': {row_key(x, 'clipping'): x for x in clip_wf},
    }
    assert len(wf_maps['root']) == 35 and len(wf_maps['quantile']) == 28 and len(wf_maps['clipping']) == 49

    for pycase in root_py:
        key = row_key(pycase, 'root')
        wfcase = wf_maps['root'][key]
        assert Q(str(pycase['K_over_sigma'])) == key[1]
        assert_precision_rows(pycase['precision_rows'])
        assert [r['working_precision'] for r in wfcase['rows']] == [80, 160]
        for pr, wr in zip(pycase['precision_rows'], wfcase['rows']):
            assert wr['pass'] is True
            residual = finite(pr['residual_over_derivative'])
            target = finite(pr['root_error_target'])
            assert residual <= target
            wfresidual = finite(wr['residual_over_derivative'])
            assert wfresidual < D(10) ** (-wr['working_precision'])
            root_residual_max = max(root_residual_max, residual, wfresidual)
        aggregate_max(summary_max['root'], compare_fields(pycase['precision_rows'], wfcase['rows'],
                                                           root_fields, relative_by_family['root']))
        root_rows_checked += 1

    for pycase in quant_py:
        key = row_key(pycase, 'quantile')
        wfcase = wf_maps['quantile'][key]
        assert_precision_rows(pycase['precision_rows'])
        assert [r['working_precision'] for r in wfcase['rows']] == [80, 160]
        for wr in wfcase['rows']:
            assert wr['pass'] is True
            res = finite(wr['log_tail_residual'])
            assert res < D(10) ** (-wr['working_precision'])
            quant_residual_max = max(quant_residual_max, res)
        aggregate_max(summary_max['quantile'], compare_fields(pycase['precision_rows'], wfcase['rows'],
                                                               quant_fields, relative_by_family['quantile']))
        quantile_rows_checked += 1

    for pycase in clip_py:
        key = row_key(pycase, 'clipping')
        wfcase = wf_maps['clipping'][key]
        pyrows, wfrows = pycase['precision_rows'], wfcase['rows']
        assert_precision_rows(pyrows)
        assert [r['working_precision'] for r in wfrows] == [80, 160]
        for pr, wr in zip(pyrows, wfrows):
            assert wr['pass'] is True
            pval = finite(pr['p']); loga = finite(pr['log_a'])
            assert pval >= 2 and D(0) < loga < D(2).ln()
            wp = wr['working_precision']
            vres = finite(pr['v_equation_residual'])
            vtarget = finite(pr['v_root_error_target'])
            assert vres <= D(100)*vtarget
            wvres = finite(wr['v_equation_residual'])
            assert wvres < D(10)**(-wp)
            clip_residual_max = max(clip_residual_max, vres, wvres)

            log_s = finite(pr['log_s'])
            py_omit = log_s < D(-10000)
            wf_log_s = finite(wr['log_s'])
            wf_omit = wf_log_s < D(-10000)
            assert py_omit == wf_omit
            assert pr['threshold_value_omitted'] is py_omit
            assert wr['threshold_value_omitted'] is wf_omit
            if py_omit:
                assert pr['s'] is None and wr['s'] is None
                tiny['python_omitted'] += 1
                tiny['wolfram_omitted'] += 1
                tiny['both_omitted'] += 1
            else:
                assert pr['s'] is not None and wr['s'] is not None
            omitted_log1p = wr['threshold_log1p_approximated']
            if omitted_log1p:
                bound = finite(wr['log_bound_on_omitted_N_times_log1p'])
                assert bound < -(D(wp+40)*D(10).ln())
                tiny['threshold_log1p_approximated'] += 1
            else:
                assert wr['log_bound_on_omitted_N_times_log1p'] is None
        aggregate_max(summary_max['clipping'], compare_fields(pyrows, wfrows, clip_fields,
                                                               relative_by_family['clipping']))
        # Compare s only when it is represented; compare its value relatively.
        for i, (pr, wr) in enumerate(zip(pyrows, wfrows)):
            if pr['s'] is None:
                assert wr['s'] is None
            else:
                err = relative_error(pr['s'], wr['s'])
                assert err < TOL, ('s', key, i, err)
                relative_by_family['clipping'].append(err)
        clipping_rows_checked += 1

    all_diffs = [v for vals in relative_by_family.values() for v in vals]
    max_rel = max(all_diffs, default=D(0))

    from verify85 import constants85
    from verify85 import sqrt_interval
    mu_iv = tuple(map(Q, constants85()['mu_interval']))
    sig_iv = sqrt_interval(mu_iv)
    critical_iv = [str(2*x) for x in sig_iv]
    assert Q(critical_iv[0]) > Q(113, 50) and Q(critical_iv[1]) < Q(227, 100)

    # Record source, inputs, and prior dependencies without asserting mathematical conclusions.
    local_sources = ['experiment89.py', 'verify89.wl', 'check89.py']
    public_sources = ['math_guard.py', 'verify85.py', 'verify86.py', 'verify84.py', 'sliver79.py']
    source_hashes = {f: sha((project_root/'code'/f) if f in public_sources else (code/f))
                     for f in local_sources+public_sources}
    input_hashes = {f: sha(d/f) for f in input_names}
    prior = [f'rounds/r{x}/report.md' for x in ('001', '011', '013', '014', '021', '022', '028', '073')]
    prior.append('rounds/r073/amplitude73.md')
    prior_hashes = {f: sha(project_root/f) for f in prior}
    out = {
        'pass_all': True,
        'scope': 'Numerical Gaussian/high-fee validation only; analytic general cap theorem not inferred from finite computations; finite-N physical onset not certified.',
        'case_counts': {'root': root_rows_checked, 'quantile': quantile_rows_checked,
                        'clipping': clipping_rows_checked, 'total': root_rows_checked+quantile_rows_checked+clipping_rows_checked},
        'dimensions_k': expected_k,
        'Python_configured_workers': py['configured_workers'],
        'Python_stored_actual_workers': py['actual_workers'],
        'Python_replay_actual_workers': len(replay_pids),
        'Python_replay_worker_pids': sorted(replay_pids),
        'Wolfram_configured_workers': wf['configured_workers'],
        'Wolfram_actual_kernel_ids': sorted(kernel_ids),
        'Wolfram_kernels_closed': wf['kernels_closed'],
        'Wolfram_symbolic_gates': wf['symbolic'],
        'working_precision': [80, 160],
        'max_relative_difference': str(max_rel),
        'max_relative_difference_by_family': {k: str(max(v, default=D(0))) for k, v in relative_by_family.items()},
        'max_root_residual_over_derivative': str(root_residual_max),
        'max_quantile_log_tail_residual': str(quant_residual_max),
        'max_v_equation_residual': str(clip_residual_max),
        'tiny_s_serialization': tiny,
        'max_precision_differences_by_field': {
            fam: {field: str(value) for field, value in fields.items()}
            for fam, fields in summary_max.items()
        },
        'Python_resource_monitor': py_resource_summary,
        'Wolfram_resource_monitor': wf_resource_summary,
        'resources_are_estimates_not_hard_OS_reservations': True,
        'new_general_C_N_coefficient': False,
        'arithmetic_critical_coefficient_interval': critical_iv,
        'constant_interval_is_not_finite_norm_certificate': True,
        'pause_request_revoked': True,
        'goal_status': 'active',
        'source_sha256': source_hashes,
        'input_sha256': input_hashes,
        'prior_dependency_sha256': prior_hashes,
    }
    Path(a.output).write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: v for k, v in out.items() if not k.endswith('sha256')}, indent=2))


if __name__ == '__main__':
    with localcontext() as ctx:
        ctx.prec = 180
        main()
