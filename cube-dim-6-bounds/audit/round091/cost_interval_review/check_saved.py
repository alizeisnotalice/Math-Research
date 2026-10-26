"""Four independent symbolic/saved-certificate checks; no producer or root search."""
import sys, signal, resource, time
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
START, CPU = time.monotonic(), time.process_time()
resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
def halt(sig, frame):
    raise TimeoutError('A091 one-attempt work limit')
signal.signal(signal.SIGALRM, halt)
signal.setitimer(signal.ITIMER_REAL, 169)
from pathlib import Path
from fractions import Fraction as Q
import json, hashlib, os, traceback
import sympy as sp

H = Path(__file__).resolve().parent
read = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n, d):
    (H / n).write_text(json.dumps(d, indent=2, ensure_ascii=False) + '\n')
if (H / 'run_receipt.json').exists():
    raise RuntimeError('Single invocation only; do not overwrite evidence')
reg = read(H / 'registration.json')
A, P = Path(reg['author_directory']), Path(reg['author_project'])
inputs = reg['author_input_pins']
src = lambda key: P / inputs[key]['path']
R = read(A / 'results.json')
out = {'status': 'INCOMPLETE', 'configurations': [], 'new_root_searches': 0,
       'author_six_configuration_executions': 0, 'old_peak_evaluations': 0,
       'kernel_calls': 0, 'old_point_or_gradient_evaluations': 0}
run = {'status': 'RUNNING', 'pid': os.getpid(), 'attempted_configurations': 0,
       'completed_configurations': 0, 'scientific_program_invocations': 1,
       'CPU_limit_seconds': 30, 'wall_limit_including_cleanup_seconds': 180,
       'automatic_retries': 0}
def checkpoint():
    run['wall_seconds'], run['cpu_seconds'] = time.monotonic() - START, time.process_time() - CPU
    save('results.json', out)
    save('run_receipt.json', run)
def begin(i):
    assert i == run['attempted_configurations'] + 1 <= 4
    run['attempted_configurations'] = i
    checkpoint()
def done(i, **kw):
    out['configurations'].append({'id': i, 'pass': True, **kw})
    run['completed_configurations'] = i
    checkpoint()
def pins():
    for name in ('source_sha256', 'extra_proof_dependencies'):
        for path, digest in reg[name].items():
            assert sha(path) == digest, path
    for rec in inputs.values():
        assert sha(P / rec['path']) == rec['sha256'], rec['path']
    for f, rec in read(A / 'manifest.json')['files'].items():
        assert sha(A / f) == rec['sha256'] and (A / f).stat().st_size == rec['bytes']
    binding = read(H / 'execution_binding.json')
    assert sha(H / 'registration.json') == binding['registration_sha256']
    assert sha(H / 'proof_review.md') == binding['proof_review_sha256']
    for f, digest in binding['code_sha256'].items():
        assert sha(H / f) == digest, f
def dyadic(d):
    return tuple(Q(d[s + '_mantissa']) * Q(2) ** d[s + '_binary_exponent'] for s in ('lower', 'upper'))

try:
    begin(1)
    pins()
    z = sp.symbols('z', positive=True)
    t_of_z = z ** 6 / 32
    changed = t_of_z ** sp.Rational(1, 6) / (1 + 1024 * t_of_z ** 2) * sp.diff(t_of_z, z)
    expected = sp.Integer(32) ** (-sp.Rational(7, 6)) * 6 * z ** 6 / (1 + z ** 12)
    assert sp.simplify(changed - expected) == 0
    # At r^12=-1, Res(6z^6/(1+z^12),r)=r^-5/2.
    assert sp.cancel(6 * z ** 6 / (12 * z ** 11) - sp.Rational(1, 2) * z ** -5) == 0
    # zeta=exp(i*pi/12), Phi24(zeta)=0. The six upper roots are
    # zeta^(2k+1). Check the finite geometric residue identity exactly.
    polyvar = sp.symbols('w')
    residue_sum_twice = sum(polyvar ** ((-5 * (2 * k + 1)) % 24) for k in range(6))
    residual = sp.rem(residue_sum_twice * (polyvar ** 5 - polyvar ** 19) - 2,
                      polyvar ** 8 - polyvar ** 4 + 1, polyvar)
    assert residual == 0
    integral = sp.Integer(32) ** (-sp.Rational(7, 6)) * sp.pi / (2 * sp.sin(5 * sp.pi / 12))
    closed = sp.pi / (32 * sp.Integer(2) ** sp.Rational(1, 3) * (sp.sqrt(3) + 1))
    assert sp.simplify(integral - closed) == 0
    done(1, name='independent_all_time_scaling_and_cyclotomic_residue',
         positive_substitution_residual='0', six_upper_poles=6,
         residue_identity_mod_Phi24='0', closed_formula_residual='0',
         exact_integral='pi/[32*u^2*(v^3+1)], u^6=2,v^6=3,u,v>0',
         quadrature_calls=0, old_peaks=0)

    begin(2)
    t = sp.symbols('t', positive=True)
    denominator_per_t = (t ** 2 + (1 + 26 * t) ** 2) / t
    gap = sp.Rational(1765, 16) - denominator_per_t
    factor = (sp.Rational(1, 16) - t) * (677 * t - 16) / t
    assert sp.cancel(gap - factor) == 0
    assert Q(677, 32) - 16 == Q(165, 32) > 0
    # Reuse the accepted basis; check its inverse with scalar J^2=6J
    # coefficients, not a new basis/direction enumeration.
    assert -Q(1, 2) + 6 * Q(1, 8) - 2 * Q(1, 8) == 0
    assert abs(-Q(1, 2) + Q(1, 8)) + 5 * Q(1, 8) == 1
    coefficient = Q(1, 2 ** 26) * Q(1, 128) * Q(16, 1765) ** 6 * Q(1, 32)
    assert coefficient == Q(1, 2 ** 14 * 1765 ** 6)
    gL = Q(R['g_L_global_coarse'])
    assert gL * (Q(22, 7) ** 6) == coefficient
    assert gL > 0
    done(2, name='global_same_time_event_and_factored_time_majorant',
         retained_event_probability='2^-26', original_basis_jacobian='1/128',
         basis_inverse_row_absolute_sum='1',
         time_majorant_gap='(1/16-t)*(677*t-16)/t',
         nonnegative_on_window='[1/32,1/16]', window_length='1/32',
         g_L_global_coarse=str(gL), includes_S2=True,
         new_basis_searches=0, kernel_calls=0)

    begin(3)
    assert R['status'] == 'PASS_SIX_FIXED_CERTIFICATE_CONFIGURATIONS'
    assert R['code_sha256'] == sha(A / 'check_exact.py')
    assert R['proof_sha256'] == sha(A / 'theorem.md')
    assert R['input_sha256'] == inputs
    arun = read(A / 'run_receipt.json')
    assert arun['exit'] == 0 and arun['wall_seconds'] < 180
    root6 = read(src('roots'))['root_enclosures']
    ul, uh = map(Q, root6['u'])
    vl, vh = map(Q, root6['v'])
    assert 0 < ul < uh and 0 < vl < vh
    assert ul ** 6 < 2 < uh ** 6 and vl ** 6 < 3 < vh ** 6
    peak = read(src('peak'))['configurations'][0]
    Al, Au = dyadic(peak['A'])
    pl, pu = dyadic(peak['pi'])
    assert 0 < Al < Au and 3 < pl < pu < Q(22, 7)
    BU, IU, IL = map(Q, [R['B_U_complete_alltime_reenvelope'], R['I_U_complete'], R['I_L_disjoint_boxes']])
    gU, g0 = Q(R['g_U_full_original_point']), Q(R['conditional_g0_unproved'])
    assert BU == Au * pu / (32 * ul ** 2 * (vl ** 3 + 1))
    assert g0 == Q(1, 65000)
    targets = {
        'complete_I_upper_from_alltime036': BU ** 6,
        'structural_model_cost_floor': Q(6 ** 6),
        'model_cost_lower_from_boxes': Q(6 ** 6) * IL ** 5 / gU ** 6,
        'model_cost_upper_coarse_pair': Q(6 ** 6) * IU ** 5 / gL ** 6,
        'gL_div6_power_lower': (gL / 6) ** 6,
        'gU_div6_power_upper': (gU / 6) ** 6,
        'unproved_g0_div6_power': (g0 / 6) ** 6,
        'conditional_plan_U': Q(6 ** 6) * IU ** 5 / g0 ** 6}
    for name in ('low-primary', 'S3-primary'):
        lower = Q(read(src(name + '_box'))['certificate']['tangent_bound_lower'])
        targets[name + '_lower_power'] = lower ** 6
    records = R['fifth_root_certificates']
    assert len(records) == 10 and {r['label'] for r in records} == set(targets)
    root = {}
    direct_checks = []
    for record in records:
        label = record['label']
        lo, hi = Q(record['lower']), Q(record['upper'])
        target = targets[label]
        assert 0 <= lo < hi and lo ** 5 <= target < hi ** 5
        assert hi - lo == Q(1, 10 ** 100)
        root[label] = (lo, hi)
        direct_checks.append({'label': label, 'lower_fifth_power_le_target': True,
                              'upper_fifth_power_gt_target': True,
                              'target_fraction_sha256': hashlib.sha256(str(target).encode()).hexdigest()})
    assert root['complete_I_upper_from_alltime036'][1] == IU
    out['direct_saved_fifth_root_checks'] = direct_checks
    done(3, name='complete_source_endpoint_and_saved_root_certificate_audit',
         saved_fifth_roots_checked=10, original_sixth_root_brackets_checked=2,
         new_root_searches=0, original_512bit_A_pi_used=True,
         author_six_configuration_executions=0,
         author_wall_seconds=arun['wall_seconds'], input_from_display=False)

    begin(4)
    review_dir = P / 'audit/round061/cutoff_gradient_pilot_review'
    review = read(src('reference'))
    assert review['full_gradient_gap_closed_for_registered_cases'] is True
    assert review['complete_physical_gradient_fields'] == 12
    assert review['original_saved_boxes_audit_retained'] is True
    assert review['universal_all_points_windows_gradient_crosscheck_claim'] is False
    for f in ('report.md', 'final_summary.json', 'input_manifest.json', 'preparation_results.json',
              'full_reference/low-primary/comparison.json', 'full_reference/S3-primary/comparison.json'):
        assert sha(review_dir / f) == review['artifact_sha256'][f]
    oldpins = read(review_dir / 'input_manifest.json')['sha256']
    preparation = read(review_dir / 'preparation_results.json')['checks'][0]
    assert preparation['status'] == 'PASS_SAVED_ENDPOINTS_AND_128_RATIONAL_CORNERS'
    summaries = read(review_dir / 'final_summary.json')
    expected_window = {'kind': 'B', 'a': '1/250', 'T': '1/2'}
    assert summaries['model'] == 'body32_unit_cauchy_v1' and summaries['window'] == expected_window
    total_lower = Q(0)
    box_records = []
    for idx, name in enumerate(('low-primary', 'S3-primary')):
        box, jet = read(src(name + '_box')), read(src(name + '_jet'))
        for kind in ('box', 'jet'):
            assert inputs[name + '_' + kind]['sha256'] in oldpins.values()
        assert jet['model'] == box['certificate']['model'] == 'body32_unit_cauchy_v1'
        assert jet['window'] == box['certificate']['window'] == expected_window
        assert jet['synthetic_fixture'] is False
        assert box['certificate']['status'] == 'PASS_CONDITIONAL_ON_MATCHING_CENTER_JET'
        prep_case = next(c for c in preparation['cases'] if c['case_id'] == name)
        assert prep_case['dyadic_to_jet_mapping_exact'] and prep_case['same_window']
        prep_box = next(b for b in prep_case['boxes'] if b['width'] == '1/1000')
        assert prep_box['lower_exact'] == box['certificate']['tangent_bound_lower']
        full_case = next(c for c in review['cases'] if c['case_id'] == name)
        assert full_case['physical_fields'] == 6 and full_case['all_native_intervals_contain_independent_reference']
        intervals = [tuple(map(Q, iv)) for iv in box['face_box']]
        assert len(intervals) == 5 and all(0 <= lo < hi <= 1 for lo, hi in intervals)
        volume = Q(1)
        for lo, hi in intervals:
            volume *= hi - lo
        contribution = volume * root[name + '_lower_power'][0]
        assert R['boxes'][idx]['name'] == name
        assert R['boxes'][idx]['box'] == box['face_box']
        assert Q(R['boxes'][idx]['volume']) == volume
        assert Q(R['boxes'][idx]['I_contribution_lower']) == contribution
        total_lower += contribution
        box_records.append(intervals)
    assert box_records[0][3][1] < box_records[1][3][0]
    assert total_lower == IL > 0
    point = read(src('point'))['points'][0]
    assert gU == Q(point['check']['G_interval'][1])
    assert Q(point['check']['exact_dist_infinity_to_S2_lower_bound']) > 0
    assert all(0 <= Q(n, point['denominator']) <= 1 for n in point['numerators'])
    point_audit = read(src('point_review'))
    assert point_audit['full_old_G_chain_verified'] is True
    assert 0 < gL < gU
    assert list(map(Q, R['model_U_interval'])) == [max(root['structural_model_cost_floor'][0], root['model_cost_lower_from_boxes'][0]), root['model_cost_upper_coarse_pair'][1]]
    Sbounds = [64 * root['gL_div6_power_lower'][0] - IU, 64 * root['gU_div6_power_upper'][1] - IL]
    assert list(map(Q, R['actual_model_S_interval'])) == Sbounds
    assert list(map(Q, R['conditional_plan_U_interval'])) == list(root['conditional_plan_U'])
    assert list(map(Q, R['conditional_plan_S_interval'])) == [64 * v - IU for v in root['unproved_g0_div6_power']]
    assert Q(R['fixed_IU_even_gU_plan_S_upper']) == 64 * root['gU_div6_power_upper'][1] - IU < 0
    Ulo, Uhi = map(Q, R['model_U_interval'])
    assert Ulo < 64 < Uhi and Sbounds[0] < 0 < Sbounds[1]
    assert R['true_model_excluded'] is False and R['fixed_coarse_IU_plan_excluded'] is True
    assert R['g0_global_proved'] is False and R['new_C6_bound'] is False
    assert R['new_kernel_calls'] == R['old_producer_executions'] == R['old_point_recomputations'] == 0
    pins()
    done(4, name='same_model_box_point_and_monotone_cost_surplus_transfer',
         real_boxes=2, independent_gradient_fields_bound=12,
         boxes_strictly_separated_face_coordinate=4, gU_input='points/0/check/G_interval/1',
         essential_infimum_justification='offS2 continuity and inward positive-volume neighborhoods',
         structural_floor='radial sector argument, not a boundary-a.e. shortcut',
         true_model_interval_crosses64=True, actual_S_interval_crosses0=True,
         only_fixed_coarse_IU_plan_excluded=True, unproved_g0_kept_separate=True)
    out.update(status='PASS_SAME_MODEL_COMPLETE_COARSE_COST_INTERVAL_ONLY',
               author_manifest_sha256=sha(A / 'manifest.json'),
               new_global_g_L_valid=True, new_complete_I_U_valid=True,
               true_model_excluded=False, fixed_coarse_IU_plan_excluded=True,
               g0_global_proved=False, new_C6_bound=False,
               direct_saved_root_certificates=12, new_root_searches=0,
               no_old_source_scalar_or_kernel_recomputed=True)
    run.update(status='PASS_STOPPED', exit_code=0)
    checkpoint()
    run['results_sha256'] = sha(H / 'results.json')
    save('run_receipt.json', run)
    print(json.dumps({'status': out['status'], 'configurations': 4,
                      'cpu_seconds': run['cpu_seconds'], 'wall_seconds': run['wall_seconds'],
                      'saved_root_certificates': 12, 'root_searches': 0, 'kernel_calls': 0}))
except BaseException as exc:
    run.update(status='FAILED_STOP_NO_RETRY', exit_code=1,
               error=repr(exc), traceback=traceback.format_exc())
    checkpoint()
    raise
finally:
    signal.setitimer(signal.ITIMER_REAL, 0)
