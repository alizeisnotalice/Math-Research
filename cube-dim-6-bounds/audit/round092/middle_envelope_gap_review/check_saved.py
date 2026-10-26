"""One-shot independent audit of saved092 evidence; no F/L function calls."""
import signal
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 179.0)
import resource
resource.setrlimit(resource.RLIMIT_CPU, (15, 15))
import os, sys, time, json, hashlib, traceback
from pathlib import Path
from fractions import Fraction as Q
from math import comb
sys.set_int_max_str_digits(0)
START, CPU = time.monotonic(), time.process_time()
HOME = Path(__file__).resolve().parent
AUTHOR = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper')
MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper')
R92 = AUTHOR / 'round092/middle_envelope_gap'
rd = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n, value):
    (HOME / n).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
def ep(v):
    return (Q(v['lower_mantissa']) * Q(2)**v['lower_binary_exponent'],
            Q(v['upper_mantissa']) * Q(2)**v['upper_binary_exponent'])
def pair(v): return Q(v['lower']), Q(v['upper'])
def require(ok, label):
    if not ok: raise AssertionError(label)
    out['checks'].append(label)
def root_claim(label, radicand, degree):
    r = roots[label]
    require(Q(r['radicand']) == radicand and r['degree'] == degree,
            'root semantic binding: ' + label)
    return pair(r)
def atan_bounds(q, n):
    # n is even; n-term alternating partial sum is below atan(1/q).
    require(n % 2 == 0, 'alternating-series parity ' + str(q))
    s = sum((Q((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(n)), Q(0))
    return s, s + Q(1, (2*n+1)*q**(2*n+1))
reg = rd(HOME / 'registration.json')
with (HOME / 'attempt.json').open('x') as f:
    json.dump({'pid': os.getpid(), 'started_unix': time.time(), 'attempt': 1,
               'cpu_hard_seconds': 15, 'wall_hard_seconds': 179,
               'configuration_cap': 3, 'children': [], 'automatic_retries': 0}, f)
out = {'status': 'RUNNING', 'completed_configurations': [], 'checks': [],
       'old_time_grid_replays': 0, 'old_F_L_function_calls': 0,
       'body_native_kernel_calls': 0, 'new_E_endpoint_values': 0,
       'new_root_searches': 0}
exit_code = 1
try:
    require(sha(__file__) == reg['checker_sha256'], 'checker source pin')
    for path, expected in reg['source_sha256'].items():
        require(sha(path) == expected, 'source pin: ' + path)
    manifest = rd(R92 / 'manifest.json')
    for name, info in manifest['files'].items():
        p = R92 / name
        require(sha(p) == info['sha256'] and p.stat().st_size == info['bytes'],
                'author manifest entry: ' + name)
    r = rd(R92 / 'lower_results.json')
    p = rd(R92 / 'supplement_parameters.json')
    s51 = rd(MAIN / 'round051/all_unsigned_middles/results.json')
    g43 = rd(MAIN / 'round043/two_parameter_middle/results.json')['configurations'][2]
    z49 = rd(MAIN / 'round049/middle_joint_marginal/results.json')['configurations'][2]
    s74 = rd(AUTHOR / 'round074/middle_long_budget/results.json')
    ledger = rd(MAIN / 'round062/primary_joint_account/actual_run/results.json')
    epsilon, b = Q(p['epsilon']), Q(p['b'])
    require((epsilon, b) == (Q(1,1000), Q(1,64)), 'entire fixed time window identity')
    # Configuration1: saved sources and uniform coefficient residuals.
    a5, b5 = atan_bounds(5, 256)
    a239, b239 = atan_bounds(239, 96)
    pi_a, pi_b = 16*a5 - 4*b239, 16*b5 - 4*a239
    pi_lo, pi_hi = ep(s51['shared_constants']['pi'])
    require(0 < pi_lo <= pi_a < pi_b <= pi_hi, 'Machin pi enclosed by saved dyadic endpoints')
    crem_lo, crem_hi = ep(s51['shared_constants']['Crem'])
    avg_lo, avg_hi = ep(s51['shared_constants']['cauchy_average_factor'])
    require(crem_lo*2970*pi_b**6 <= 44743079 <= crem_hi*2970*pi_a**6,
            'Crem endpoint direction via independent pi bracket')
    require(1 < avg_lo <= avg_hi and (avg_lo-1)**2 <= Q(6,928) <= (avg_hi-1)**2,
            'Cauchy average positive root endpoints')
    K = Q(r['K_lower'])
    require(K == Q(1000,729)*crem_lo*avg_lo, 'K lower is positive product of lower endpoints')
    require(Q(r['pi_upper']) == pi_hi, 'pi denominator uses upper endpoint')
    require(1024*b*b == Q(1,4) and 2*(Q(9,16)+1024*b*b-1) == Q(-3,8),
            'D increasing c and decreasing u on full rectangle')
    dmax = (1+Q(1,16)+1024*b*b)**2 - Q(1,4)
    require(dmax == Q(377,256) == Q(r['denominator_max']), 'denominator maximum 377/256')
    cj, cw = Q(r['cJ_lower']), Q(r['cW_lower'])
    require(0 < cj*pi_hi*pi_hi*dmax <= 64, 'common marginal coefficient downward')
    require(0 < cw <= 5*cj and cw*b*pi_hi*(1+1024*b*b) <= 1,
            'both whole-mass min branches lower bounded')
    require(Q(15,64) <= Q(1,4) < Q(3,4) <= Q(49,64), 'whole interval containment')
    fs, vs = {}, {}
    for k in range(1,6):
        src = s51['cases'][str(k)]
        geom = g43 if k == 3 else src
        limits = z49 if k == 3 else src
        f = r['f_coefficients'][str(k)]
        require(ep(limits['z_lower'])[1] <= Q(1,4) and ep(limits['z_upper'])[0] >= Q(3,4),
                'common interval for k=' + str(k))
        h = {key: ep(geom[key])[1] for key in ('hv0','hv1','hw0','hw1','g')}
        require(h == {key: Q(v) for key,v in f['h_upper_inputs'].items()} and min(h.values()) > 0,
                'positive saved geometry upper endpoints k=' + str(k))
        fs[k], vs[k] = Q(f['uniform_F_lower']), Q(f['V_lower'])
        dv = (h['hv0']**2+992*b*b)*(h['hv1']**2+960*b*b)
        dw = (h['hw0']**2+992*b*b)*(h['hw1']**2+960*b*b)
        # Cross-multiplied saved coefficient inequality, never evaluates F(t).
        require(0 < fs[k]*dv*dw*(h['g']**2+928*b*b) <= K*(dv+dw),
                'uniform B denominator residual k=' + str(k))
        require(0 < vs[k] <= ep(s74['cases'][str(k)]['V_middle'])[0]
                and f['n'] == comb(5,k-1), 'saved volume lower and true multiplicity k=' + str(k))
    out['completed_configurations'].append('C1_sources_pi_and_saved_uniform_coefficient_residuals')
    # Configuration2: all saved root certificates, semantic graph, exact integral.
    roots = {v['label']: v for v in r['root_certificates']}
    require(len(roots) == len(r['root_certificates']) == 20, '20 unique saved root certificates')
    for name, v in roots.items():
        lo, hi = pair(v); rad = Q(v['radicand']); degree = v['degree']
        require(0 <= lo < hi and lo**degree <= rad < hi**degree, 'integer power enclosure: ' + name)
    fifth = {k: root_claim('f%d-fifth'%k, fs[k], 5)[0] for k in fs}
    qpower = {k: root_claim('f%d-six-fifths'%k, fs[k]**6, 5)[0] for k in fs}
    A = root_claim('volume-constant-lower', sum(comb(5,k-1)*vs[k]*qpower[k] for k in fs)**5, 6)[0]
    J = root_claim('joint-t-power-coefficient', (cj*sum(comb(5,k-1)*fifth[k] for k in fs))**5, 6)[0]
    W = root_claim('whole-t-power-coefficient', (cw*max(fifth.values()))**5, 6)[0]
    C = min(J,W)
    require([A,J,W,C] == [Q(r['branch_lower_coefficients'][key]) for key in
            ('volume_constant','joint_t5over3','whole_t5over3','minimum_t5over3')],
            'all three branches bound before minimum')
    xlo,xhi = root_claim('power-clipping-switch', (A/C)**3, 5)
    require((xlo,xhi) == pair(r['switch']) and epsilon < xlo < xhi < b,
            'unique switch strictly interior, no uncovered time')
    e83 = Q(1,10**8)
    require(e83**3 == epsilon**8, 'epsilon 8/3 exact power')
    lower = A*b - Q(5,8)*A*xhi - Q(3,8)*C*e83
    upper = A*b - Q(5,8)*A*xlo - Q(3,8)*C*e83
    require(0 < lower <= upper and (lower,upper) == pair(r['analytic_lower_integral']),
            'negative switch term uses upper endpoint for integral lower')
    require(str(lower) == manifest['strict_lower_integral'], 'reported exact integral matches seal')
    out['completed_configurations'].append('C2_all_branch_min_roots_and_complete_window_integral')
    # Configuration3: use only saved aggregate074 budget and saved062 allowance.
    allowance = Q(ledger['sufficient_remaining_D_F_B_allowance'])
    paid = Q(ledger['paid_C_E_H_I_upper'])
    short = Q(s74['comparison']['old_short_upper_read_only'])
    R = Q(s74['strict_long_norm_rational_upper'])
    Rslab = Q(s74['slab_exact_upper_sum'])
    tail = ep(s74['infinite_tail'])[1]
    g = Q(r['conditional_g0'])
    require(g == Q(1,65000) and 0 < allowance < lower < Rslab < Rslab+tail < R,
            'strict fixed-allocation failure and saved aggregate upper budget')
    ilo = paid + root_claim('new-envelope-plan-I-lower', (short+lower)**6, 5)[0]
    ihi = paid + root_claim('old-rectangle-plan-I-upper', (short+R)**6, 5)[1]
    tlo,thi = root_claim('conditional-I-target-unit', (g/6)**6, 5)
    tlo,thi = 64*tlo,64*thi
    flo,fhi = root_claim('conditional-U-factor', (6/g)**6, 5)
    blo = root_claim('envelope-plan-B-lower', ilo**5, 6)[0]
    bhi = root_claim('old-rectangle-plan-B-upper', ihi**5, 6)[1]
    for key, exact in {
        'conditional_plan_I': (ilo,ihi), 'conditional_plan_B': (blo,bhi),
        'conditional_I_target': (tlo,thi), 'conditional_plan_U': (flo*ilo,fhi*ihi),
        'conditional_plan_S_margin': (tlo-ihi,thi-ilo),
        'same_E_integral_enclosure': (lower,R),
        'analytic_lower_over_allowance': (lower/allowance,upper/allowance)
    }.items(): require(pair(r[key]) == exact, 'saved plan endpoints: ' + key)
    for key, value in {
        'saved_D_F_B_allowance': allowance, 'envelope_lower_minus_allowance': lower-allowance,
        'remaining_irrecoverable_plan_I_deficit': ilo-thi,
        'maximum_86_rectangle_error': Rslab-lower,
        'maximum_full_time_integration_gain': R-lower,
        'maximum_recoverable_S_gain': ihi-ilo,
        'maximum_recoverable_total_B_gain': bhi-blo, 'old_infinite_tail_upper': tail
    }.items(): require(Q(r[key]) == value, 'saved residual: ' + key)
    require(thi-ilo < 0 and flo*ilo > 64, 'strict S upper negative and U lower above64')
    require(all(r[key] is False for key in ('global_g0_certified','true_norm_lower_bound',
            'actual_model_cost_interval','actual_C6_bound')) and
            r['far_assumed_zero_only_for_optimistic_plan_diagnosis'] is True,
            'plan-only conclusion, unproved g and unpaid far preserved')
    out['completed_configurations'].append('C3_saved_062_allowance_and_U_S_plan_scope')
    for path, expected in reg['source_sha256'].items():
        require(sha(path) == expected, 'unchanged source pin: ' + path)
    out.update(status='PASS_SAVED_ONLY_FIXED_ENVELOPE_EXCLUSION',
               strict_integral_lower=str(lower), S_plan_upper=str(thi-ilo),
               U_plan_lower=str(flo*ilo), exact_checks=len(out['checks']))
    exit_code = 0
except BaseException as ex:
    out.update(status='FAILED_STOP_NO_RETRY', error=repr(ex), traceback=traceback.format_exc())
finally:
    # No child/process cleanup. Hard CPU/wall timers stay active through final writes and exit.
    save('results.json', out)
    receipt = {'pid': os.getpid(), 'exit_code': exit_code, 'status': out['status'],
               'wall_seconds_through_results_write': time.monotonic()-START,
               'cpu_seconds_through_results_write': time.process_time()-CPU,
               'configuration_count': len(out['completed_configurations']),
               'results_sha256': sha(HOME/'results.json'),
               'registration_sha256': sha(HOME/'registration.json'),
               'cpu_hard_seconds': 15, 'wall_hard_seconds': 179,
               'timers_active_through_receipt_and_exit': True,
               'children': [], 'active_handles': [], 'automatic_retries': 0}
    save('receipt.json', receipt)
    print(json.dumps(receipt))
sys.exit(exit_code)
