#!/usr/bin/env python3
"""Round120: exactly two new time-first spatial/Jensen diagnostics."""
import sys
sys.set_int_max_str_digits(0)
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
from pathlib import Path
import resource
import time
from datetime import datetime, timezone
from fractions import Fraction as F

HERE = Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def utc(): return datetime.now(timezone.utc).isoformat()
def req(x, label):
    if not x: raise ArithmeticError(label)
def atomic(name, data):
    p = HERE/name
    tmp = p.with_name(p.name+'.own-tmp')
    with tmp.open('w') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
        f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(tmp, p)
def dyad(m, e): return F(int(m)*2**max(e, 0), 2**max(-e, 0))
def pin_interval(p):
    return tuple(dyad(p[s+'_mantissa'], p[s+'_binary_exponent']) for s in ('lower', 'upper'))
def root5(x):
    req(x > 0, 'positive root input')
    den, lo, hi = 2**160, 0, 2**160
    for _ in range(256):
        if F(hi, den)**5 >= x: break
        hi *= 2
    else: raise ArithmeticError('registered doubling cap')
    for _ in range(512):
        if hi-lo <= 1: break
        mid = (hi+lo)//2
        if F(mid, den)**5 <= x: lo = mid
        else: hi = mid
    else: raise ArithmeticError('registered bisection cap')
    low, up = F(lo, den), F(hi, den)
    req(low**5 <= x <= up**5, 'exact fifth-root verification')
    return [str(low), str(up)]
def display(iv): return [format(float(F(s)), '.17g') for s in iv]


def cfg_geometry(inp, piL, piU):
    import sympy as sp
    p, aa = sp.symbols('pi a', positive=True)
    C = sp.Rational(3, 16)*p**sp.Rational(-5, 2)
    collapse = sp.simplify(C**6*(p**3/384)**5*(4*p)**3*aa)
    req(sp.simplify(collapse - 3*p**3*aa/2**53) == 0,
        'critical-degree angular normalization')
    R, rho = F(inp['R']), F(inp['rho'])
    req(R == 256 and rho == F(1, 64), 'fixed only radius and far region')
    tail_coefficient = 384/R
    tube_coefficient = F(496, 2)*rho**4*R**4*F(14, 3)**2/32**4
    req(tail_coefficient == F(3, 2), 'six true Cauchy marginals')
    req(tube_coefficient == F(1519, 1152), 'all496 true rank4 projected events')
    P0 = 1-tail_coefficient/3-tube_coefficient/9
    req(P0 == F(3665, 10368) and P0 > 0, 'exact coarse true-far probability')
    refined = 1-tail_coefficient/piL-tube_coefficient/piL**2
    req(refined > P0, 'same-R refined probability')
    return {'id': 'critical_angular_identity_and_true_far_event', 'status': 'PASS',
            'symbolic_J_far_fifth': str(collapse), 'R': str(R), 'rho': str(rho),
            'all_unordered_actual_pairplanes': 496, 'ordered_pair_terms_not_generated': 992,
            'tail_coefficient_over_pi': str(tail_coefficient),
            'tube_coefficient_over_pi_squared': str(tube_coefficient),
            'P_far_lower_coarse': str(P0), 'P_far_lower_refined': str(refined),
            'P_far_lower_refined_display': format(float(refined), '.17g'),
            'all_pair_overlaps': 'union upper bound for excluded tube events; no disjointness used',
            'original_time_first': True, 'original_full_law_projection_only': True}


def cfg_budget(inp, geo, piL, piU):
    P0, Pref = F(geo['P_far_lower_coarse']), F(geo['P_far_lower_refined'])
    RL, RU = F(inp['far_allowance']['lower']), F(inp['far_allowance']['upper'])
    req(0 < RL <= RU, 'full111 exact raw allocation')
    W5 = F(45, 2**86)*P0**5
    Wref5 = F(45, 2**86)*Pref**5
    req(W5 > RU**5, 'intrinsic exact-peak/exact-direction Jensen obstruction')
    intrinsic, refined = root5(W5), root5(Wref5)
    ratio = root5(W5/RU**5)
    upper5 = F(81, 2**80)/piL**3
    lower5 = F(81, 2**80)/piU**3
    full = [root5(lower5)[0], root5(upper5)[1]]
    req(F(full[1])**5 >= upper5, 'accepted108 complete actual far upper endpoint')
    from flint import arb, ctx
    ctx.prec = 384
    z = ((arb(81)/(arb(2)**80*arb.pi()**3)).log()/5).exp()
    req(bool(z > arb('0.000018')) and bool(z < arb('0.000019')),
        'new scalar Arb crosscheck')
    return {'id': 'intrinsic_Jensen_obstruction_and_complete_far_upper',
            'status': 'STRICT_OBSTRUCTION_TO_SPECIFIED_JENSEN_FAMILY',
            'J_far_lower_fifth_exact': str(W5),
            'J_far_lower_interval': intrinsic, 'J_far_lower_display': display(intrinsic),
            'J_far_refined_lower_fifth': str(Wref5),
            'J_far_refined_lower_interval': refined,
            'J_far_refined_lower_display': display(refined),
            'intrinsic_lower_over_RU_interval': ratio,
            'intrinsic_lower_over_RU_display': display(ratio),
            'residual_vs_RL_fifth': str(W5-RL**5),
            'residual_vs_RU_fifth': str(W5-RU**5),
            'accepted108_complete_actual_far_I_upper_interval': full,
            'accepted108_complete_actual_far_I_upper_display': display(full),
            'accepted108_complete_actual_far_upper_fifth_input': str(upper5),
            'arb384_complete_bound_crosscheck': str(z),
            'actual_I_far_lower_bound': False, 'actual_I_far_competitive': False,
            'exact_peak_upgrade_alone_can_rescue_this_family': False,
            'new_C6_bound': False, 'global_g0_proved': False,
            'stop': 'end both fixed configurations; no radius, cell, shape or peak refresh'}


def main():
    start, cpu = time.monotonic(), time.process_time()
    reg = json.loads((HERE/'registration.json').read_text())
    release = json.loads((HERE/'parent_release.json').read_text())
    req(release.get('authorized') is True and str(release.get('round')) == '120', 'main release')
    req(release['script_sha256'] == sha(__file__), 'released science entry')
    req(release['registration_sha256'] == sha(HERE/'registration.json'), 'released registration')
    for n, h in reg['prepared_pins'].items(): req(sha(HERE/n) == h, 'prepared pin '+n)
    sources = json.loads((HERE/'source_manifest.json').read_text())['sources']
    for s in sources: req(sha(s['path']) == s['sha256'], 'source '+s['role'])
    with (HERE/'started.json').open('x') as f:
        json.dump({'pid': os.getpid(), 'utc': utc(), 'release_sha256': sha(HERE/'parent_release.json')}, f)
        f.write('\n')
    status, done, failure = 'FAILED', [], None
    try:
        inp = json.loads((HERE/'upstream.json').read_text())
        piL, piU = pin_interval(inp['pi_saved036'])
        req(3 < piL < piU < F(22, 7), 'complete saved pi')
        done.append(cfg_geometry(inp, piL, piU))
        atomic('checkpoint.json', {'count': len(done), 'configurations': done})
        done.append(cfg_budget(inp, done[0], piL, piU))
        atomic('checkpoint.json', {'count': len(done), 'configurations': done})
        status = 'COMPLETE_TWO_FIXED_CONFIGURATIONS'
        atomic('results.json', {'status': status, 'configurations': done,
                               'new_C6_bound': False, 'full_far_competitive': False,
                               'release_sha256': sha(HERE/'parent_release.json')})
        print(json.dumps({'status': status, 'pid': os.getpid(),
                          'intrinsic_J_lower': done[1]['J_far_lower_display'],
                          'complete_far_upper': done[1]['accepted108_complete_actual_far_I_upper_display']}), flush=True)
    except BaseException as e:
        failure = type(e).__name__+': '+str(e)
        raise
    finally:
        atomic('run_receipt.json', {'status': status, 'completed_configurations': len(done),
            'failure': failure, 'pid': os.getpid(), 'utc': utc(),
            'CPU_sample_seconds': time.process_time()-cpu,
            'wall_sample_seconds_before_receipt': time.monotonic()-start,
            'scope': 'sample before final receipt write and interpreter exit; observer records complete child',
            'max_rss_bytes_macos': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'native_kernel_or_old_program_calls': 0, 'workers': 1, 'children': 0,
            'wall_limit': None, 'CPU_limit': None})


if __name__ == '__main__': main()
