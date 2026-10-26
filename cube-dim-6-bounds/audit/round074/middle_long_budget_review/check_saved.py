"""Six fixed independent checks; never execute the author's evaluator."""
import sys, time, signal, resource
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
START = time.monotonic()
CPU = time.process_time()
resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
def halt(sig, frame):
    raise TimeoutError('A074 single-attempt limit reached')
signal.signal(signal.SIGALRM, halt)
signal.setitimer(signal.ITIMER_REAL, 169)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import os, json, hashlib, traceback
import sympy as sp

H = Path(__file__).resolve().parent
read = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name, obj):
    (H / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')
def ep(v):
    lo = Q(v['lower_mantissa']) * Q(2) ** v['lower_binary_exponent']
    hi = Q(v['upper_mantissa']) * Q(2) ** v['upper_binary_exponent']
    assert 0 < lo <= hi
    return lo, hi
def pack(v):
    return {'lower': str(v[0]), 'upper': str(v[1])}
def overlap(x, y):
    return max(x[0], y[0]) <= min(x[1], y[1])
def add(x, y):
    return x[0] + y[0], x[1] + y[1]
def scale(x, q):
    assert q > 0
    return x[0] * q, x[1] * q
def mul(x, y):
    assert x[0] > 0 and y[0] > 0
    return x[0] * y[0], x[1] * y[1]

if (H / 'run_receipt.json').exists():
    raise RuntimeError('Single invocation only: no overwrite or retry')
reg = read(H / 'registration.json')
A = Path(reg['author_directory'])
run = {'status': 'RUNNING', 'pid': os.getpid(), 'attempted_configurations': 0,
       'completed_configurations': 0, 'scientific_program_invocations': 1,
       'cpu_limit_seconds': 30, 'wall_limit_including_cleanup_seconds': 180,
       'automatic_retries': 0, 'author_evaluator_calls': 0,
       'finite_panel_function_reevaluations': 0, 'old_M_integrals': 0,
       'old_short_evaluations': 0, 'body32_kernel_calls': 0}
out = {'status': 'INCOMPLETE', 'configurations': [], 'root_scale_decimal_digits': 60}
def checkpoint():
    run['wall_seconds'] = time.monotonic() - START
    run['cpu_seconds'] = time.process_time() - CPU
    save('results.json', out)
    save('run_receipt.json', run)
def begin(i):
    assert i == run['attempted_configurations'] + 1 <= 6
    run['attempted_configurations'] = i
    checkpoint()
def done(i, **data):
    out['configurations'].append({'id': i, 'pass': True, **data})
    run['completed_configurations'] = i
    checkpoint()
def pins():
    for p, expected in reg['source_sha256'].items():
        assert sha(p) == expected, p
    for p, expected in reg['author_source_pins'].items():
        assert sha(p) == expected, p
    manifest = read(A / 'manifest.json')
    for f, record in manifest['files'].items():
        assert sha(A / f) == record['sha256'], f
        assert (A / f).stat().st_size == record['bytes'], f
    binding = read(H / 'execution_binding.json')
    assert sha(H / 'registration.json') == binding['registration_sha256']
    for f, digest in binding['code_sha256'].items():
        assert sha(H / f) == digest, f
    return len(manifest['files'])

try:
    begin(1)
    manifest_files = pins()
    author_run = read(A / 'run_receipt.json')
    process = read(A / 'process_receipt.json')
    tool = read(A / 'tool_execution_receipt.json')
    R = read(A / 'results.json')
    params = read(A / 'parameters.json')
    author_reg = read(A / 'registration.json')
    author_binding = read(A / 'execution_binding.json')
    assert author_run['results_sha256'] == sha(A / 'results.json')
    assert author_run['completed_configurations'] == author_run['attempted_configurations'] == 9
    assert author_run['status'] == 'COMPLETED_STOPPED'
    assert process['status'] == 'EXITED_AND_REAPED'
    assert author_run['exit_code'] == process['exit_code'] == tool['exit_code'] == 0
    assert author_run['pid'] == process['owned_pid'] == tool['scientific_pid']
    assert 0 < author_run['elapsed_seconds'] <= process['total_wall_seconds'] < 180
    assert 0 < author_run['cpu_seconds'] < 30
    assert tool['child_exited_and_reaped'] and tool['automatic_retries'] == 0
    assert author_binding['registration_sha256'] == sha(A / 'registration.json')
    assert author_binding['theorem_sha256'] == sha(A / 'theorem.md')
    assert author_binding['parameters_sha256'] == sha(A / 'parameters.json')
    for f, digest in author_binding['code_sha256'].items():
        assert sha(A / f) == digest
    assert Q(64) / Q(1, 16) == 1024
    assert [32 * n for n in (31, 30, 29)] == [992, 960, 928]
    assert Q(10, 9) ** 3 == Q(1000, 729)
    L, ell, t, u, d = sp.symbols('L ell t u d', positive=True)
    gamma = ell * t / L + u
    extra = L ** 2 * t ** 2 + 2 * L ** 2 * t * gamma
    # Clear the positive denominator ell*gamma*d*(d+extra).
    numerator = (ell + L) * gamma * (d + extra) - ell * (gamma + t) * d
    positive_form = L * u * d + (ell + L) * gamma * extra
    assert sp.expand(numerator - positive_form) == 0
    done(1, name='source_receipt_binding_and_actual_time_denominators',
         author_manifest_files=manifest_files, author_source_pins=len(reg['author_source_pins']),
         P_divisor=1024, time_denominators=[992, 960, 928],
         omitted_parameter_loss='1000/729', average_inequality_residual='0',
         positive_difference_numerator='L*u*d+(ell+L)*gamma*(L^2*t^2+2*L^2*t*gamma)',
         author_wall_including_guard_seconds=process['total_wall_seconds'],
         author_cpu_seconds=author_run['cpu_seconds'])

    begin(2)
    z, a = sp.symbols('z a', real=True)
    dp, dm = (z + 1) ** 2 + a ** 2, (z - 1) ** 2 + a ** 2
    # Independently decompose the rational integrand, not repeat author's
    # differentiation of the complete log/atan expression.
    partial_fraction_residual = sp.expand((z + 2) * dm + (2 - z) * dp - 4 * (1 + a ** 2))
    assert partial_fraction_residual == 0
    assert sp.expand((z + 2) - ((z + 1) + 1)) == 0
    assert sp.expand((2 - z) - (-(z - 1) + 1)) == 0
    q = Q(params['q'])
    assert q == Q(6, 5) and q - 1 == Q(1, 5)
    assert 1 / q == Q(5, 6) and (q - 1) / q == Q(1, 6)
    assert Q(2048, 16) == 128
    done(2, name='independent_partial_fractions_and_interpolation_powers',
         rational_identity_residual='0',
         decomposition='[(z+2)/Dplus+(2-z)/Dminus]/[4*(1+a^2)]',
         primitive_integration='split z+2=(z+1)+1; 2-z=-(z-1)+1; log plus real atan',
         q='6/5', q_minus_one='1/5', norm_exponent='5/6', whole_pointwise_exponent='1/6')

    begin(3)
    weights = [comb(5, k - 1) for k in range(1, 6)]
    assert weights == params['weights'] == [1, 5, 10, 10, 5]
    assert sum(weights) == 31
    assert Q(params['rho']) == Q(1, 64)
    assert Q(params['epsilon']) == Q(1, 1000)
    assert Q(params['T']) == Q(1, 4)
    assert Q(params['ratio']) == Q(16, 15)
    assert params['no_refinement'] is True
    assert list(map(Q, params['whole_union_coordinate_range'])) == [Q(1, 4) - Q(1, 64), Q(3, 4) + Q(1, 64)]
    source51 = next(Path(p) for p in reg['author_source_pins'] if '/round051/' in p and p.endswith('/results.json'))
    old = read(source51)
    singleton = old['cases']['5']
    assert ep(singleton['z_lower'])[0] <= Q(1, 4) <= ep(singleton['z_lower'])[1]
    assert ep(singleton['z_upper'])[0] <= Q(3, 4) <= ep(singleton['z_upper'])[1]
    for k in range(1, 6):
        assert R['cases'][str(k)]['n'] == weights[k - 1]
        assert R['cases'][str(k)]['m'] == 6 - k
        assert R['cases'][str(k)]['old_M_recomputed'] is False
    done(3, name='whole_union_cover_weights_singletons_and_slab_structure',
         weights=weights, tube_count=31, physical_coordinate_union_count=5,
         whole_mass_candidates=2, singletons_use_exact_mean_range=True,
         marginal_time_growth_power=2, full_face_time_growth_power=1,
         overlap_assumption='positive union only; no disjointness',
         whole_mass_application='one bound for complete p_t; no restricted-branch reflection')

    begin(4)
    slabs, panels = R['time_slabs'], params['panels']
    assert len(slabs) == len(panels) == 86
    previous = Q(params['epsilon'])
    T, ratio = Q(params['T']), Q(params['ratio'])
    norm_lower_sum = norm_upper_sum = pay_lower_sum = pay_upper_sum = Q(0)
    counts = {'tube_volume_union': 0, 'tube_joint_mass_union': 0, 'one_whole_union_mass': 0}
    for row, panel in zip(slabs, panels):
        left, right = map(Q, panel)
        assert left == previous < right <= T
        assert right == min(left * ratio, T)
        assert [Q(row['left']), Q(row['right'])] == [left, right]
        width = right - left
        norm, pay = ep(row['norm_upper_rectangle']), ep(row['integrated_rectangle'])
        independent_pay = scale(norm, width)
        assert overlap(pay, independent_pay)
        norm_lower_sum += independent_pay[0]
        norm_upper_sum += independent_pay[1]
        pay_lower_sum += pay[0]
        pay_upper_sum += pay[1]
        counts[row['chosen']] += 1
        previous = right
    assert previous == T
    assert pay_upper_sum == Q(R['slab_exact_upper_sum'])
    assert overlap(ep(R['slab_total']), (pay_lower_sum, pay_upper_sum))
    assert list(counts.values()) == [46, 0, 40]
    rectangle_sum = (norm_lower_sum, norm_upper_sum)
    out['independent_width_times_saved_norm_sum'] = pack(rectangle_sum)
    out['saved_integrated_rectangle_exact_sum'] = pack((pay_lower_sum, pay_upper_sum))
    done(4, name='saved86_rows_exact_continuous_time_coverage_and_aggregation',
         saved_rows_checked=86, gaps=0, overlaps_in_time=0,
         start='1/1000', finite_end='1/4', selected_counts=counts,
         finite_panel_function_reevaluations=0,
         method='exact rational widths times saved outward norm endpoints; not sampling or quadrature')

    begin(5)
    root_records = []
    S = 10 ** reg['rational_root_scale_decimal_digits']
    def point_power(x, power, root, label):
        assert x > 0
        y = x ** power
        v = (y.numerator * S ** root) // y.denominator
        j, _ = sp.integer_nthroot(v, root)
        j = int(j)
        lower, upper = Q(j, S), Q(j + 1, S)
        # Each endpoint is justified directly by its own exact power,
        # never by containment of two independently widened intervals.
        assert lower ** root <= y < upper ** root
        root_records.append({'label': label, 'power': power, 'root': root,
                             'lower': str(lower), 'upper': str(upper),
                             'exact_power_inequalities_passed': True})
        return lower, upper
    def interval_power(iv, power, root, label):
        lower = point_power(iv[0], power, root, label + ':lower_input')[0]
        upper = point_power(iv[1], power, root, label + ':upper_input')[1]
        return lower, upper
    shared = old['shared_constants']
    K = scale(mul(ep(shared['Crem']), ep(shared['cauchy_average_factor'])), Q(1000, 729))
    aggregate = (Q(0), Q(0))
    tail_inputs = []
    for k in range(1, 6):
        case = R['cases'][str(k)]
        V, C = ep(case['V_middle']), ep(case['Ctail'])
        expected_C = add(scale(ep(old['cases'][str(k)]['D_big']), Q(1, 1024)),
                         scale(K, Q(2, 992 * 960 * 928)))
        # Rational interval propagation of these saved positive inputs is
        # exact and has no root widening. Arb must enclose that full box.
        assert C[0] <= expected_C[0] <= expected_C[1] <= C[1]
        Cp = interval_power(C, 6, 5, 'Ctail_k' + str(k))
        term = scale(mul(V, Cp), Q(case['n']))
        aggregate = add(aggregate, term)
        tail_inputs.append({'k': k, 'coefficient_interval_contains_exact_source_box': True,
                            'positive_V': True, 'positive_Ctail': True})
    coefficient = interval_power(aggregate, 5, 6, 'tail_coefficient')
    tail = scale(coefficient, 1 / (5 * T ** 5))
    assert overlap(tail, ep(R['infinite_tail']))
    out['independent_infinite_tail'] = pack(tail)
    out['root_power_certificates'] = root_records
    done(5, name='independent_fraction_complete_infinite_tail',
         tail_case_inputs=5, actual_endpoint_root_certificates=len(root_records),
         coefficient_consistency=tail_inputs, integration='integral_T^infinity t^-6 dt=1/(5*T^5)',
         saved_tail_overlap=True, old_M_integrals=0, old_short_evaluations=0)

    begin(6)
    simple = Q(R['strict_long_norm_rational_upper'])
    assert simple == Q(49091761943697, 10 ** 17)
    independent_final = add(rectangle_sum, tail)
    assert independent_final[1] < simple
    recorded_total = ep(R['long_budget_expression'])
    assert overlap(recorded_total, independent_final)
    source_ledger = next(Path(p) for p in reg['author_source_pins'] if '/primary_joint_account/' in p)
    ledger = read(source_ledger)
    allowance = Q(ledger['sufficient_remaining_D_F_B_allowance'])
    short = Q(old['union']['simple_rational_union_upper'])
    comparison = R['comparison']
    assert Q(comparison['saved_D_F_allowance']) == allowance
    assert Q(comparison['old_short_upper_read_only']) == short
    assert Q(comparison['new_same_domain_full_middle_upper']) == simple + short
    assert independent_final[0] > allowance and recorded_total[0] > allowance
    assert comparison['status'] == 'FIXED_UPPER_BUDGET_EXCEEDS_SAVED_D_PLUS_F_ALLOWANCE'
    assert comparison['formula_lower_exceeds_allowance'] is True
    for flag in ('actual_norm_lower_bound_claimed', 'model_excluded', 'new_C6_bound'):
        assert comparison[flag] is False
    for flag in ('global_I', 'global_g', 'new_C6_bound'):
        assert R[flag] is False
    assert R['body32_kernel_calls'] == 0
    pins()
    out['independent_final_enclosure_from_saved_norms_and_tail'] = pack(independent_final)
    out['strict_margin_to_claim'] = str(simple - independent_final[1])
    out['strict_long_upper'] = str(simple)
    out['same_domain_short_plus_long_upper'] = str(simple + short)
    done(6, name='strict_final_claim_saved_allowance_and_scope',
         strict_claim_verified=True, fixed_upper_expression_exceeds_allowance=True,
         actual_norm_lower_bound=False, model_excluded=False,
         short_added_before_q_power=True, full_I=False, global_g=False, C6=False,
         author_stopped_no_retry=True)
    out.update(status='PASS_VALID_LOCAL_LONG_BUDGET_NONCOMPETITIVE_ONLY',
               author_manifest_sha256=sha(A / 'manifest.json'),
               author_report_sha256=sha(A / 'report.md'),
               finite_panel_function_reevaluations=0, source_geometry_reproved=False,
               independence_scope='new algebra plus saved-certificate audit, exact sums and independent root tail; no full author replay')
    run.update(status='PASS_STOPPED', exit_code=0)
    checkpoint()
    run['results_sha256'] = sha(H / 'results.json')
    save('run_receipt.json', run)
    print(json.dumps({'status': out['status'], 'configurations': 6,
                      'strict_long_upper': str(simple), 'wall_seconds': run['wall_seconds'],
                      'cpu_seconds': run['cpu_seconds'], 'root_certificates': len(root_records),
                      'finite_panel_function_reevaluations': 0}))
except BaseException as exc:
    run.update(status='FAILED_STOP_NO_RETRY', exit_code=1,
               error=repr(exc), traceback=traceback.format_exc())
    checkpoint()
    raise
finally:
    signal.setitimer(signal.ITIMER_REAL, 0)
