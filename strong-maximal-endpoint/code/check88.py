"""Independent final acceptance for the r088 affine-patch certificates."""
import argparse
import hashlib
import json
import multiprocessing as mp
import re
import sys
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
from pathlib import Path


def num(x):
    # Mathematica InputForm precision marks and scientific-exponent syntax.
    return D(re.sub(r'`[^*]*', '', str(x)).replace('*^', 'E'))


def qnum(x):
    z = Q(x)
    return D(z.numerator) / D(z.denominator)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite(x):
    z = num(x)
    assert z.is_finite(), x
    return z


def interval(pair):
    assert len(pair) == 2
    lo, hi = map(qnum, pair)
    assert lo.is_finite() and hi.is_finite() and lo <= hi
    return lo, hi


def replay_equal(stored, replay, ignored):
    assert len(stored) == len(replay)
    for x, y in zip(stored, replay):
        assert {k: v for k, v in x.items() if k not in ignored} == {
            k: v for k, v in y.items() if k not in ignored
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory')
    ap.add_argument('output')
    ap.add_argument('--project-root', required=True)
    a = ap.parse_args()
    d, pr = Path(a.directory), Path(a.project_root)
    code = Path(__file__).resolve().parent
    for path in (d / 'python/certificates88.json', d / 'wolfram/verification88.json',
                 d / 'wolfram/resources88.json'):
        assert path.is_file(), path
    p = json.loads((d / 'python/certificates88.json').read_text())
    w = json.loads((d / 'wolfram/verification88.json').read_text())
    resource = json.loads((d / 'wolfram/resources88.json').read_text())

    sys.path.insert(0, str(pr / 'code'))
    sys.path.insert(0, str(code))
    from experiment88 import center_job, constants as constants88, dyadic_job

    assert p['pass_all'] is True
    assert p['constants']['prior_r084'] == constants88()
    assert p['constants']['small_exp_cutoff'] == -300
    assert Q(p['constants']['small_exp_upper']) == Q(1, 10**100)
    assert p['constants']['grid_constant'] == 130
    assert p['constants']['saturation_upper_constant'] == 184
    assert p['constants']['full_upper_constant'] == 132
    assert w['pass'] is True and w['kernels_closed'] is True
    assert len(w['symbolic']) == 12 and all(v is True for v in w['symbolic'].values())
    assert resource['exit_code'] == 0 and resource['stop_reason'] is None
    assert 0 < w['per_case_memory_limit_MiB'] <= 1024
    assert w['numeric_log_cutoff'] == -400 and w['input_guard_digits'] == 80
    reserve_mib = D(str(resource['reserve_GiB'])) * D(1024)
    assert D(str(resource['min_system_available_MiB'])) >= reserve_mib
    budget_mib = D(str(resource['owned_RSS_budget_GiB'])) * D(1024)
    assert D(str(resource['peak_owned_RSS_MiB'])) <= budget_mib

    py_dyadic = p['dyadic_rows']
    py_center = p['center_rows']
    wf_dyadic = w['dyadic_rows']
    wf_center = w['center_rows']
    assert len(py_dyadic) == len(wf_dyadic) == 12
    assert len(py_center) == len(wf_center) == 24
    assert sum(len(x['rows']) for x in wf_center) == 72
    assert p['configured_workers'] == 12
    assert 1 <= w['configured_workers'] <= 12
    kernel_ids = {r['kernel_id'] for r in wf_dyadic + wf_center}
    assert kernel_ids and len(kernel_ids) <= w['configured_workers']
    assert all(isinstance(k, int) and k > 0 for k in kernel_ids)
    assert p['scope'] == ('Exact affine-patch reflected-center budget bounds and exact outward-certificate volume; '
                          'no global input allocation theorem or new C_N bound.')

    with mp.get_context('spawn').Pool(12) as pool:
        replay_dyadic = pool.map(dyadic_job, [2 ** k for k in range(10, 22)])
        centers = [(2 ** k, str(delta)) for k in range(1, 9)
                   for delta in (Q(1, (2 ** k) ** 2), Q(1, 10), Q(1, 2))]
        replay_center = pool.map(center_job, centers)
    replay_equal(py_dyadic, replay_dyadic, {'pid', 'seconds'})
    replay_equal(py_center, replay_center, {'pid', 'seconds'})
    replay_dyadic_workers = len({r['pid'] for r in replay_dyadic})
    replay_center_workers = len({r['pid'] for r in replay_center})

    dprec, dcross, dwidth = [], [], []
    cprec, ccross, cwidth = [], [], []
    finite_grid_count = tail_count = lower_count = upper_count = 0
    dyadic_dimensions = []
    peak_j_checks = []
    for py, wl in zip(py_dyadic, wf_dyadic):
        n = py['N']
        assert n == wl['N'] and n == 2 ** (10 + len(dyadic_dimensions))
        dyadic_dimensions.append(n)
        assert len(py['precision_rows']) == len(wl['rows']) == 2
        assert [x['working_precision'] for x in py['precision_rows']] == [80, 160]
        assert [x['working_precision'] for x in wl['rows']] == [80, 160]
        peak_js = [x['peak_j'] for x in py['precision_rows'] + wl['rows']]
        assert len(set(peak_js)) == 1
        peak_j_checks.append({'N': n, 'peak_j': peak_js[0], 'all_four_agree': True})
        full_lo, full_hi = interval(py['log_constructive_full_cost_interval'])
        tail_lo, tail_hi = interval(py['log_tail_upper_interval'])
        sat_lo, sat_hi = interval(py['log_saturation_upper_interval'])
        ana_lo, ana_hi = interval(py['log_full_upper_interval'])
        dwidth.extend([full_hi-full_lo, tail_hi-tail_lo, sat_hi-sat_lo, ana_hi-ana_lo])
        assert full_hi < ana_lo
        assert tail_hi < full_lo - D(300)
        for field in ('lower_saturated_log', 'lower_soft_log',
                      'upper_saturated_log', 'upper_soft_log'):
            bounds = [finite(r[field]) for r in wl['rows']]
            dprec.append(abs(bounds[0]-bounds[1]))
        for sw in (py['precision_rows'], wl['rows']):
            logs = [finite(r['finite_grid_log_sum']) for r in sw]
            tails = [finite(r['log_tail_upper']) for r in sw]
            dprec.extend([abs(logs[0]-logs[1]), abs(tails[0]-tails[1])])
            assert all(full_lo <= z <= full_hi for z in logs)
            assert all(tail_lo <= z <= tail_hi for z in tails)
            finite_grid_count += len(logs)
            tail_count += len(tails)
        dlog = [finite(x['finite_grid_log_sum']) for x in py['precision_rows']]
        wlog = [finite(x['finite_grid_log_sum']) for x in wl['rows']]
        dtail = [finite(x['log_tail_upper']) for x in py['precision_rows']]
        wtail = [finite(x['log_tail_upper']) for x in wl['rows']]
        dcross.extend([abs(dlog[i]-wlog[i]) for i in range(2)])
        dcross.extend([abs(dtail[i]-wtail[i]) for i in range(2)])
        assert all(full_lo <= z <= full_hi for z in dlog + wlog)
        assert all(tail_lo <= z <= tail_hi for z in dtail + wtail)

        py_lower = {int(x['depth_power']): interval(x['log_budget_lower_interval'])
                    for x in py['lower_rows']}
        dwidth.extend(hi-lo for lo, hi in py_lower.values())
        assert set(py_lower) == {2, 3}
        for wr in wl['rows']:
            for field, b in (('lower_saturated_log', 2), ('lower_soft_log', 3)):
                x = finite(wr[field]); lo, hi = py_lower[b]
                assert lo <= x <= hi
                lower_count += 1
            for field, iv in (('upper_saturated_log', (sat_lo, sat_hi)),
                              ('upper_soft_log', (ana_lo, ana_hi))):
                x = finite(wr[field])
                assert iv[0] <= x <= iv[1]
                upper_count += 1

    center_cases = center_volumes = 0
    for py, wl in zip(py_center, wf_center):
        n, delta = py['N'], Q(py['delta'])
        assert (n, delta) == (wl['N'], Q(wl['delta']))
        assert py['positive_terms'] == 81 and py['tail_relative_bound_certified'] is True
        assert qnum(py['relative_series_tail_upper']) <= D('1e-100')
        assert len(py['rows']) == len(wl['rows']) == 3
        for prw, wlw in zip(py['rows'], wl['rows']):
            r = Q(prw['r'])
            assert r == Q(wlw['r'])
            lo, hi = interval(prw['log_unit_rectangle_center_volume_interval'])
            cwidth.append(hi-lo)
            assert [x['working_precision'] for x in prw['precision_rows']] == [80, 160]
            assert [x['working_precision'] for x in wlw['precision_rows']] == [80, 160]
            vals_by_sw = []
            for sw in (prw['precision_rows'], wlw['precision_rows']):
                vals = [finite(x['log_center_volume']) for x in sw]
                cprec.append(abs(vals[0]-vals[1]))
                assert all(lo <= x <= hi for x in vals)
                vals_by_sw.append(vals)
            ccross.extend(abs(vals_by_sw[0][i]-vals_by_sw[1][i]) for i in range(2))
            center_volumes += 1
        center_cases += 1

    assert center_cases == 24 and center_volumes == 72
    assert max(dprec) < D('1e-45') and max(dcross) < D('1e-45')
    assert max(cprec) < D('1e-40') and max(ccross) < D('1e-40')
    max_width = max(dwidth + cwidth)
    assert max_width < D('1e-20')

    local_sources = ['experiment88.py', 'verify88.wl', 'math_guard.py', 'check88.py']
    public_sources = ['verify84.py', 'sliver79.py']
    source_files = local_sources + public_sources
    source_hashes = {f: sha((pr / 'code' / f) if f in public_sources else (code / f))
                     for f in source_files}
    inputs = ['python/certificates88.json', 'wolfram/verification88.json', 'wolfram/resources88.json']
    input_hashes = {f: sha(d / f) for f in inputs}
    prior = ['rounds/r028/report.md', 'rounds/r075/report.md', 'rounds/r081/smoothing81.md',
             'rounds/r084/obstruction84.md', 'rounds/r087/corner87.md']
    prior_hashes = {f: sha(pr / f) for f in prior}
    out = {
        'pass_all': True,
        'scope': 'Only affine patch and outward center volume; no global input allocation theorem or new C_N bound.',
        'dimensions': {'dyadic_N': dyadic_dimensions, 'center_N': sorted({x['N'] for x in py_center}),
                       'center_delta_cases': center_cases, 'physical_center_volumes': center_volumes},
        'counts': {'dyadic_dimensions': len(dyadic_dimensions), 'center_parameter_cases': center_cases,
                   'physical_center_volumes': center_volumes, 'all_symbolic_gates': len(w['symbolic']),
                   'dyadic_finite_grid_values_checked': finite_grid_count,
                   'dyadic_tail_values_checked': tail_count, 'directed_lower_values_checked': lower_count,
                   'directed_upper_values_checked': upper_count},
        'Python_configured_workers': p['configured_workers'],
        'Python_actual_dyadic_processes': p['actual_dyadic_processes'],
        'Python_actual_center_processes': p['actual_center_processes'],
        'Python_replay_dyadic_processes': replay_dyadic_workers,
        'Python_replay_center_processes': replay_center_workers,
        'Wolfram_configured_workers': w['configured_workers'],
        'Wolfram_actual_kernel_ids': sorted(kernel_ids),
        'Wolfram_kernels_closed': w['kernels_closed'],
        'Wolfram_symbolic_gates': w['symbolic'],
        'working_precision': [80, 160],
        'dyadic_peak_j_consistency': peak_j_checks,
        'max_dyadic_precision_difference': str(max(dprec)),
        'max_dyadic_cross_software_difference': str(max(dcross)),
        'max_center_precision_difference': str(max(cprec)),
        'max_center_cross_software_difference': str(max(ccross)),
        'max_directed_log_interval_width': str(max_width),
        'center_relative_series_tail_upper_max': str(max(qnum(x['relative_series_tail_upper']) for x in py_center)),
        'center_tail_certified_cases': sum(x['tail_relative_bound_certified'] is True for x in py_center),
        'tail_log_gap_minimum': str(min(interval(x['log_constructive_full_cost_interval'])[0] -
                                       interval(x['log_tail_upper_interval'])[1] for x in py_dyadic)),
        'cutoff': {'python_small_exp_cutoff': -300, 'python_small_exp_upper': '1e-100',
                   'wolfram_numeric_log_cutoff': w['numeric_log_cutoff'],
                   'wolfram_numeric_relative_omission_upper': w['numeric_relative_omission_upper']},
        'resource_monitor': resource,
        'resource_available_memory_is_estimate': True,
        'source_sha256': source_hashes,
        'input_sha256': input_hashes,
        'prior_dependency_sha256': prior_hashes,
        'goal_status': 'active',
        'pause_request_revoked': True,
    }
    Path(a.output).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: v for k, v in out.items() if not k.endswith('sha256')}, indent=2))


if __name__ == '__main__':
    with localcontext() as ctx:
        ctx.prec = 180
        main()
