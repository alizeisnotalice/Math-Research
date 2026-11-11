"""Replay and independently validate r093 exact rational certificates."""
import argparse
import hashlib
import json
import math
import multiprocessing as mp
import re
import sys
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
from pathlib import Path

import mpmath as m

TOL = D('1e-45')
N_VALUES = (2, 3, 4, 8, 16, 32, 64)
P_VALUES = ('9/10', '19/20', '49/50')
A_VALUES = ('0', '1/4', '1', '4')
PRECISIONS = (80, 160)
SIG_DIGITS = 40
REPLAY_FILES = (
    'python/certificates93.json', 'python/resources93.json',
    'wolfram/verification93.json', 'wolfram/resources93.json',
)
OWN_SOURCES = ('experiment93.py', 'verify93.wl', 'check93.py')
PRIOR_FILES = (
    'rounds/r011/report.md', 'rounds/r021/report.md',
    'rounds/r022/report.md', 'rounds/r028/report.md',
    'rounds/r079/report.md', 'rounds/r079/rational79.md',
    'rounds/r079/step_halo79.md', 'rounds/r092/core92.md',
    'rounds/r092/audit92.md',
)


def decimal_string(value):
    return re.sub(r'`[^*]*', '', str(value)).replace('*^', 'E')


def dec(value):
    text = decimal_string(value)
    if '/' in text and not any(c in text for c in '.Ee'):
        q = Q(text)
        with localcontext() as ctx:
            ctx.prec = 180
            return D(q.numerator)/D(q.denominator)
    result = D(text)
    assert result.is_finite(), value
    return result


def as_mpf(value):
    return m.mpf(decimal_string(value))


def relative_error(x, y):
    a, b = dec(x), dec(y)
    scale = max(abs(a), abs(b))
    return D(0) if scale == 0 else abs(a-b)/scale


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resource(row):
    assert row['exit_code'] == 0 and row['stop_reason'] is None
    reserve = D(str(row['reserve_GiB']))
    budget = D(str(row['owned_RSS_budget_GiB']))
    peak = D(str(row['peak_owned_RSS_MiB']))
    available = D(str(row['min_system_available_MiB']))
    assert reserve >= 1 and available >= reserve*1024
    assert peak <= budget*1024
    assert row['samples'] > 0 and 'estimate' in row['scope'].lower()
    return row


def without_pid(row):
    if isinstance(row, list):
        return [without_pid(item) for item in row]
    if isinstance(row, dict):
        return {key: without_pid(value) for key, value in row.items()
                if key not in ('pid', 'seconds')}
    return row


def _pow10(e):
    return Q(10**e) if e >= 0 else Q(1, 10**(-e))


def _compact_enclosure(value, digits=SIG_DIGITS):
    assert value >= 0
    if value == 0:
        return {'lower_significand_integer': '0', 'upper_significand_integer': '0',
                'decimal_exponent': 0, 'significand_scale_power': digits-1,
                'exact_zero': True}
    e = ((value.numerator.bit_length()-value.denominator.bit_length())*30103)//100000
    while value < _pow10(e):
        e -= 1
    while value >= _pow10(e+1):
        e += 1
    scale_e = e-digits+1
    scaled = value/_pow10(scale_e)
    lo = scaled.numerator//scaled.denominator
    hi = lo if scaled.denominator == 1 else lo+1
    return {'lower_significand_integer': str(lo),
            'upper_significand_integer': str(hi), 'decimal_exponent': e,
            'significand_scale_power': digits-1, 'exact_zero': False}


def _fraction_digest(value):
    source = f'{value.numerator}/{value.denominator}'.encode('ascii')
    return {'sha256': hashlib.sha256(source).hexdigest(),
            'numerator_bits': abs(value.numerator).bit_length(),
            'denominator_bits': value.denominator.bit_length()}


def _check_digest(stored, exact):
    assert stored == _fraction_digest(exact)


def _check_enclosure(stored, exact):
    assert stored == _compact_enclosure(exact)


def _check_lower_enclosure(stored, exact):
    expected = _compact_enclosure(exact)
    for field in ('lower_significand_integer', 'decimal_exponent',
                  'significand_scale_power', 'exact_zero'):
        assert stored[field] == expected[field]
    assert stored['exact_zero'] is False


def _check_gamma_python(case, logi):
    n, p = case['N'], Q(case['p'])
    arg = Q(2)/p-1
    blo, bhi = logi(arg)
    assert case['b_argument_2_over_p_minus_1'] == str(arg)
    assert case['b_interval'] == [str(blo), str(bhi)]
    cdf_upper = (n*bhi)**n/Q(math.factorial(n))
    margin = 1-p-cdf_upper
    assert cdf_upper < 1-p and margin > 0
    assert case['gamma_cdf_upper_formula'] == '(N*b_hi)^N/N!'
    _check_digest(case['gamma_cdf_upper_exact_fraction_sha256'], cdf_upper)
    _check_enclosure(case['gamma_cdf_upper_tight_decimal_enclosure'], cdf_upper)
    assert case['gamma_cdf_upper_log_interval'] == [str(x) for x in logi(cdf_upper)]
    assert case['threshold_1_minus_p'] == str(1-p)
    _check_digest(case['strict_margin_exact_fraction_sha256'], margin)
    _check_lower_enclosure(case['strict_margin_lower_decimal_bound'], margin)
    assert case['strict_margin_log_interval'] == [str(x) for x in logi(margin)]
    assert case['strict_margin_positive_exact'] is True
    assert case['integral_bound_uses_exp_minus_t_le_1'] is True
    return {'b_argument': arg, 'b_interval': (blo, bhi),
            'cdf_upper': cdf_upper, 'margin': margin}


def _actual_gamma_cdf_log(n, p):
    b = m.log(2/(m.mpf(p.numerator)/p.denominator)-1)
    x = n*b
    cdf = m.gammainc(n, 0, x)/m.gamma(n)
    return b, m.log(cdf)


def _gamma_checks(py_rows, wf_rows, logi):
    pmap = {(row['N'], Q(row['p'])): row for row in py_rows}
    wmap = {(row['N'], Q(row['p'])): row for row in wf_rows}
    expected_keys = {(n, Q(p)) for n in N_VALUES for p in P_VALUES}
    assert set(pmap) == expected_keys and set(wmap) == expected_keys
    assert len(py_rows) == len(wf_rows) == 21

    maximum = {}
    max_error = D(0)
    max_fraction_gap = Q(0)
    for key in sorted(expected_keys, key=lambda z: (z[0], P_VALUES.index(str(z[1])))):
        py_case, wf_case = pmap[key], wmap[key]
        n, p = key
        exact = _check_gamma_python(py_case, logi)
        assert wf_case['pass'] is True
        assert int(wf_case['kernel_id']) > 0
        linear_b = Q(wf_case['independent_log_upper'])
        wf_bound = Q(wf_case['independent_Gamma_CDF_upper'])
        wf_margin = Q(wf_case['independent_positive_margin'])
        expected_linear_b = 2*(1-p)/p
        expected_wf_bound = (n*expected_linear_b)**n/Q(math.factorial(n))
        assert linear_b == expected_linear_b and linear_b >= exact['b_interval'][1]
        assert wf_bound == expected_wf_bound
        assert wf_margin == 1-p-wf_bound and wf_margin > 0
        assert wf_bound < 1-p
        assert wf_bound >= exact['cdf_upper']
        max_fraction_gap = max(max_fraction_gap, wf_bound-exact['cdf_upper'])

        prows, wrows = py_case['precision_rows'], wf_case['precision_rows']
        assert [r['working_precision'] for r in prows] == list(PRECISIONS)
        assert [r['working_precision'] for r in wrows] == list(PRECISIONS)
        assert all(r['pass'] is True for r in wrows)
        assert len(prows) == len(wrows) == 2
        assert all(r['diagnostic_only'] is True for r in prows)
        py_actual_cdf_logs = []
        for wp, prow, wrow in zip(PRECISIONS, prows, wrows):
            with m.workdps(wp+40):
                b_num, actual_cdf_log = _actual_gamma_cdf_log(n, p)
                b_expected = m.nstr(b_num, 65)
                cdf_log_expected = m.nstr(actual_cdf_log, 65)
            py_actual_cdf_logs.append(cdf_log_expected)
            pairs = (
                ('b', prow['b'], wrow['b']),
                ('power_bound_log', prow['log_gamma_cdf_upper_at_exact_b'],
                 wrow['log_power_bound_at_exact_b']),
                ('actual_gamma_cdf_log', cdf_log_expected, wrow['log_gamma_lower_CDF']),
            )
            for name, left, right in pairs:
                err = relative_error(left, right)
                assert err < TOL, (key, wp, name, err)
                max_error = max(max_error, err)
                maximum[name] = max(maximum.get(name, D(0)), err)
            # Confirm the Python displayed CDF power bound agrees with the
            # independent Wolfram power-bound logarithm; it is not the true
            # Gamma CDF and is not substituted for the latter.
            with m.workdps(wp+40):
                wf_power = m.exp(as_mpf(wrow['log_power_bound_at_exact_b']))
                python_power = as_mpf(prow['gamma_cdf_upper_at_exact_b'])
                err = relative_error(m.nstr(wf_power, 65), m.nstr(python_power, 65))
            assert err < TOL, (key, wp, 'power bound numeric', err)
            max_error = max(max_error, err)
            maximum['power_bound_numeric'] = max(maximum.get('power_bound_numeric', D(0)), err)

        # Both precision rows independently contain high precision output, but
        # stored display is 65 significant digits. The 1e-45 checks below are
        # diagnostics and do not strengthen the exact Fraction certificate.
        for name, values in (
            ('Python b', [r['b'] for r in prows]),
            ('Python power-bound log', [r['log_gamma_cdf_upper_at_exact_b'] for r in prows]),
            ('Python power bound', [r['gamma_cdf_upper_at_exact_b'] for r in prows]),
            ('Wolfram b', [r['b'] for r in wrows]),
            ('Wolfram CDF log', [r['log_gamma_lower_CDF'] for r in wrows]),
            ('Wolfram power-bound log', [r['log_power_bound_at_exact_b'] for r in wrows]),
            ('MP actual CDF log', py_actual_cdf_logs),
        ):
            err = relative_error(values[0], values[1])
            assert err < TOL, (key, name, 'within precision', err)
            max_error = max(max_error, err)
            maximum[name] = max(maximum.get(name, D(0)), err)

    return {'gamma_cases': 21, 'max_relative_error': str(max_error),
            'max_relative_error_by_quantity': {k: str(v) for k, v in maximum.items()},
            'max_wolfram_minus_python_exact_bound': str(max_fraction_gap)}


def _polynomial_map(rows, wolfram=False):
    result = {}
    for row in rows:
        if wolfram:
            n, a, t = row['N'], Q(row['a']), Q(row['T'])
            expected_ts = (Q(0), a/2, a, a+1, a+4)
            assert 1 <= row['T_index'] <= 5 and t == expected_ts[row['T_index']-1]
        else:
            n, a, t = row['N'], Q(row['a']), Q(row['T'])
        result.setdefault((n, a, t), []).append(row)
    return result


def _polynomial_checks(py_groups, wf_rows, polynomial, fraction_digest,
                       scientific_enclosure):
    wf_groups = _polynomial_map(wf_rows, wolfram=True)
    assert set(py_groups) == set(wf_groups)
    assert len(wf_rows) == 140
    max_numeric_error = D(0)
    numeric_errors_by_branch = {}
    strict_t_a0 = 0
    strict_other = 0

    for key, py_group in py_groups.items():
        wf_group = wf_groups[key]
        assert len(py_group) == len(wf_group)
        n, a, t = key
        pt = polynomial(n, t)
        if t >= a:
            exact_ray = polynomial(n, t-a)/pt
            upper = exact_ray
            branch = 'T_ge_a_polynomial_ratio'
            assert exact_ray <= 1
            assert (exact_ray == 1) == (a == 0)
        else:
            exact_ray = None
            upper = 1/pt
            branch = 'T_lt_a_exponential_times_reciprocal'
            assert upper <= 1
            assert t-a < 0

        for py_case, wf_case in zip(py_group, wf_group):
            assert py_case['N'] == wf_case['N'] == n
            assert py_case['case'] == branch
            assert wf_case['pass'] is True and int(wf_case['kernel_id']) > 0
            assert Q(wf_case['normalized_upper']) == upper
            assert wf_case['kind'] == ('exact_polynomial_ratio' if t >= a
                                      else 'exponential_le_one_bound')
            assert py_case['exact_ray_le_one'] is True
            assert py_case['equality_at_one'] == (a == 0)
            assert py_case['exact_ray_strictly_less_than_one'] == (a != 0)
            if t >= a:
                assert py_case['exact_ray_formula'] == 'P_{N-1}(T-a)/P_{N-1}(T)'
                assert py_case['exact_fraction_le_one'] is True
                assert py_case['outside_halo_upper_fraction_sha256'] is None
                assert py_case['analytic_exp_factor_upper'] is None
                assert py_case['exact_ray_upper_fraction_sha256'] == fraction_digest(exact_ray)
                assert py_case['exact_ray_upper_tight_decimal_enclosure'] == scientific_enclosure(exact_ray)
                strict_t_a0 += int(a == 0)
                strict_other += int(a != 0)
            else:
                assert py_case['exact_ray_formula'] == 'exp(T-a)/P_{N-1}(T)'
                assert py_case['exact_ray_upper_fraction_sha256'] is None
                assert py_case['exact_ray_upper_tight_decimal_enclosure'] is None
                assert py_case['outside_halo_upper_formula'] == '1/P_{N-1}(T)'
                assert py_case['analytic_exp_factor_upper'] == '1'
                assert py_case['analytic_exp_factor_is_equality'] is False
                assert py_case['outside_halo_upper_fraction_sha256'] == fraction_digest(upper)
                assert py_case['outside_halo_upper_tight_decimal_enclosure'] == scientific_enclosure(upper)
                assert upper == 1/pt and upper <= 1
                strict_other += 1

            # Wolfram's normalized_numeric is the actual expression in both
            # branches. For T<a it is exp(T-a)/P(T), whereas Python stores only
            # the rational upper envelope 1/P(T); compare to a direct MP value,
            # never to that Python upper envelope.
            with m.workdps(160):
                tm = m.mpf(t.numerator)/t.denominator
                am = m.mpf(a.numerator)/a.denominator
                poly_m = m.fsum([tm**j/m.factorial(j) for j in range(n)])
                if t >= a:
                    numer_m = m.fsum([(tm-am)**j/m.factorial(j) for j in range(n)])
                    actual = numer_m/poly_m
                else:
                    actual = m.exp(tm-am)/poly_m
                shown = as_mpf(wf_case['normalized_numeric'])
                err = relative_error(m.nstr(actual, 65), m.nstr(shown, 65))
            assert err < TOL, (key, 'Wolfram normalized_numeric', err)
            max_numeric_error = max(max_numeric_error, err)
            numeric_errors_by_branch[branch] = max(numeric_errors_by_branch.get(branch, D(0)), err)
            assert shown > 0
            if t < a:
                with m.workdps(160):
                    upper_m = m.mpf(upper.numerator)/upper.denominator
                assert shown < upper_m

    return {'polynomial_cases': 140, 'T_ge_a_cases': sum(
                len(v) for (n, a, t), v in py_groups.items() if t >= a),
            'T_lt_a_cases': sum(len(v) for (n, a, t), v in py_groups.items() if t < a),
            'a_zero_equality_rows': strict_t_a0,
            'strict_inequality_rows': strict_other,
            'max_wolfram_numeric_relative_error': str(max_numeric_error),
            'max_numeric_relative_error_by_branch': {
                key: str(value) for key, value in numeric_errors_by_branch.items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory')
    ap.add_argument('output')
    ap.add_argument('--project-root', required=True)
    args = ap.parse_args()
    data_dir = Path(args.directory)
    project_root = Path(args.project_root)
    code_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(project_root/'code'))
    sys.path.insert(0, str(code_dir))

    from experiment93 import (A_VALUES, N_VALUES as EXP_N, P_VALUES,
                              PRECISIONS, batch1_job, batch2_job,
                              fraction_digest, polynomial, scientific_enclosure)
    from sliver79 import logi

    py, py_resource, wf, wf_resource = [json.loads((data_dir/name).read_text())
                                        for name in REPLAY_FILES]
    py_resource, wf_resource = resource(py_resource), resource(wf_resource)
    assert py['pass_all'] is True and wf['pass'] is True
    assert tuple(EXP_N) == N_VALUES
    assert tuple(str(x) for x in P_VALUES) == P_VALUES_EXPECTED
    assert tuple(str(x) for x in A_VALUES) == A_VALUES_EXPECTED
    assert py['configured_workers'] == 12 and py['actual_workers_observed_in_polynomial_batch'] == 12
    assert len(py['workers_seen']) == 12 and len(set(py['workers_seen'])) == 12
    assert py['N_values'] == list(N_VALUES)
    assert py['p_values'] == list(P_VALUES_EXPECTED)
    assert py['a_values'] == list(A_VALUES_EXPECTED)
    assert py['precision_rows'] == list(PRECISIONS)
    assert py['batch1_gamma_cdf_cases'] == 21 and py['batch2_polynomial_cases'] == 140
    assert wf['configured_workers'] == 12 and wf['kernels_closed'] is True
    assert 0 < wf['per_case_memory_limit_MiB'] <= 1024
    assert wf['pass'] is True

    # Re-run the actual Python work functions under spawn. Strip only PIDs and
    # elapsed-time fields, which are execution metadata rather than results.
    jobs1 = [(n, p) for n in N_VALUES for p in P_VALUES_EXPECTED]
    with mp.get_context('spawn').Pool(12) as pool:
        replay1 = pool.map(batch1_job, jobs1)
        replay2 = pool.map(batch2_job, N_VALUES)
    assert without_pid(py['batch1']) == without_pid(replay1)
    assert without_pid(py['batch2']) == without_pid(replay2)
    replay_workers = {row['pid'] for row in replay1+replay2}
    assert len(replay_workers) == 12

    assert len(py['batch1']) == 21 and len(wf['Gamma_rows']) == 21
    assert len(py['batch2']) == 7 and len(wf['polynomial_rows']) == 140
    wf_ids = {row['kernel_id'] for row in wf['Gamma_rows']+wf['polynomial_rows']}
    assert all(isinstance(kid, int) and kid > 0 for kid in wf_ids)
    assert len(wf_ids) == 12

    gamma_stats = _gamma_checks(py['batch1'], wf['Gamma_rows'], logi)
    py_poly_rows = [case for group in py['batch2'] for case in group['cases']]
    wf_poly_rows = wf['polynomial_rows']
    assert len(py_poly_rows) == 140
    py_poly_groups = _polynomial_map(py_poly_rows)
    poly_stats = _polynomial_checks(py_poly_groups, wf_poly_rows,
                                    polynomial, fraction_digest, scientific_enclosure)

    source_paths = {
        **{name: code_dir/name for name in OWN_SOURCES},
        'math_guard.py': project_root/'code/math_guard.py',
        'experiment88.py': project_root/'code/experiment88.py',
        'verify85.py': project_root/'code/verify85.py',
        'sliver79.py': project_root/'code/sliver79.py',
        'verify84.py': project_root/'code/verify84.py',
    }
    source_hashes = {name: sha(path) for name, path in source_paths.items()}
    input_hashes = {name: sha(data_dir/name) for name in REPLAY_FILES}
    docs = ('core93.md', 'audit93.md')
    document_hashes = {name: sha(data_dir/name) for name in docs}
    prior_hashes = {name: sha(project_root/name) for name in PRIOR_FILES}

    output = {
        'pass_all': True,
        'scope': ('Exact finite rational Gamma premises and the specified polynomial family '
                  'only; no general H_N sharp bound, no numerical certification of s_p, '
                  'and no physical N0 inference.'),
        'Python_configured_workers': py['configured_workers'],
        'Python_actual_workers': py['actual_workers_observed_in_polynomial_batch'],
        'Python_replay_workers': len(replay_workers),
        'Wolfram_configured_workers': wf['configured_workers'],
        'Wolfram_actual_kernel_ids': sorted(wf_ids),
        'Wolfram_kernels_closed': wf['kernels_closed'],
        'Python_resource_monitor': py_resource,
        'Wolfram_resource_monitor': wf_resource,
        'resources_are_sampled_estimates_not_OS_reservations': True,
        'precision_rows': list(PRECISIONS),
        'relative_tolerance': str(TOL),
        'gamma_checks': gamma_stats,
        'polynomial_checks': poly_stats,
        'source_sha256': source_hashes,
        'input_sha256': input_hashes,
        'document_sha256': document_hashes,
        'prior_dependency_sha256': prior_hashes,
        'general_H_N_sharp_bound_claimed': False,
        's_p_numerically_certified': False,
        'physical_N0_inferred': False,
        'goal_status': 'active',
        'pause_request_revoked': True,
    }
    Path(args.output).write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({key: value for key, value in output.items()
                      if 'sha256' not in key and key not in ('Python_resource_monitor',
                                                              'Wolfram_resource_monitor')},
                     indent=2))


P_VALUES_EXPECTED = ('9/10', '19/20', '49/50')
A_VALUES_EXPECTED = ('0', '1/4', '1', '4')


if __name__ == '__main__':
    with localcontext() as context:
        context.prec = 180
        main()
