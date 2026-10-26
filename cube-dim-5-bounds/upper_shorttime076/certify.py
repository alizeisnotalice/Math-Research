#!/usr/bin/env python3
"""Four new short-time bounds from stored constants; standard library only."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
CP = Path('/tmp/upper_shorttime076/checkpoint.json')
D = 2**384


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(out):
    CP.parent.mkdir(parents=True, exist_ok=True)
    tmp = CP.with_suffix('.tmp')
    tmp.write_text(json.dumps(out, indent=2) + '\n')
    tmp.replace(CP)


def timeout(_signum, _frame):
    raise TimeoutError('180-second complete-invocation limit; checkpoint preserved')


def dyadic_out(lo, hi):
    low = F((lo*D).numerator//(lo*D).denominator, D)
    high = F((hi*D).numerator//(hi*D).denominator + 1, D)
    assert 0 < low <= lo < hi < high
    return low, high


def fifth_root(q):
    n = q.numerator*D**5//q.denominator
    lo, hi = 0, 1 << ((n.bit_length()+4)//5)
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**5 <= n:
            lo = mid
        else:
            hi = mid
    low, high = F(lo, D), F(lo+1, D)
    assert low**5 <= q < high**5
    return low, high


def endpoints(row):
    lo, hi = F(row['lower_dyadic']), F(row['upper_dyadic'])
    assert 0 < lo < hi
    return lo, hi


def pair(lo, hi):
    assert lo < hi
    return {'lower':str(lo), 'upper':str(hi),
            'lower_display':float(lo), 'upper_display':float(hi)}


def main():
    reg = json.loads((HERE/'registration.json').read_text())
    assert reg['status'] == 'registered_before_computation'
    inputs = []
    for row in reg['sources']:
        assert sha(row['path']) == row['sha256']
        inputs.append(row)
    identity = {'registration_sha256':sha(HERE/'registration.json'),
                'source_sha256':sha(__file__), 'sources':inputs}
    if CP.exists():
        assert '--resume' in sys.argv, 'Inspect existing checkpoint; rerun is forbidden'
        out = json.loads(CP.read_text())
        assert all(out[k] == v for k,v in identity.items())
        assert out['status'] != 'complete_four_new_shorttime_bounds'
        out.setdefault('prior_invocations', []).append({'pid':out['pid'], 'status':out['status']})
        out['pid'] = os.getpid()
    else:
        out = {'id':reg['id'], 'status':'armed', 'pid':os.getpid(),
               'python':sys.executable, **identity, 'cases':[]}
    save(out)
    print(json.dumps({k:out[k] for k in ['status','pid','registration_sha256','source_sha256']}), flush=True)
    if '--arm' in sys.argv:
        input()
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(180)
    begin = time.monotonic()
    try:
        stored = json.loads(Path(inputs[0]['path']).read_text())['constants']['constants']
        p3l,p3u = endpoints(stored['P3'])
        p4l,p4u = endpoints(stored['P4'])
        sl,su = endpoints(stored['S'])
        fl,fu = endpoints(stored['F5'])
        pi_source = json.loads(Path(inputs[3]['path']).read_text())['pi_machin80_interval']
        pil,piu = dyadic_out(F(pi_source[0]), F(pi_source[1]))
        assert 3 < pil < piu < F(22,7)
        C3l,C3u = fl*p4l,fu*p4u
        C1l,C1u = fl*480*sl*p3l/piu,fu*480*su*p3u/pil
        old = json.loads(Path(inputs[4]['path']).read_text())
        assert old['eps'] == '1/2048'
        oldl,oldu = F(old['E_eps_Lq_norm']['lower']),F(old['E_eps_Lq_norm']['upper'])
        assert 0 < oldl < oldu
        out['derived_constants'] = {'C3':pair(C3l,C3u), 'C1':pair(C1l,C1u),
                                    'pi_outward_384_dyadic':pair(pil,piu)}
        out['old_E_063'] = {'epsilon':old['eps'], 'ball':old['E_eps_Lq_norm']['ball'], **pair(oldl,oldu)}
        out['status'] = 'running_fixed_four_cases'
        save(out)
        for index in range(len(out['cases']), 4):
            eps = F(reg['epsilons'][index])
            cap = F(reg['preselected_strict_E_caps'][index])
            multiplier = F(reg['preselected_old_E_multipliers'][index])
            L5 = F(5,6)**5*eps**6*(C3l+C1l*eps**2)/piu**4
            U5 = F(5,6)**5*eps**6*(C3u+C1u*eps**2)/pil**4
            lower,_ = fifth_root(L5)
            _,upper = fifth_root(U5)
            cap_pass = upper < cap
            relation = 'strictly_below_old_E' if upper < oldl else (
                'strictly_above_old_E' if lower > oldu else 'intervals_overlap')
            row = {'epsilon':str(eps), 'E_fifth_power':pair(L5,U5),
                   'E_norm':pair(lower,upper), 'fixed_strict_cap':str(cap),
                   'cap_verified':cap_pass, 'comparison_to_old_E':relation,
                   'ratio_new_to_old':pair(lower/oldu,upper/oldl),
                   'ratio_old_to_new':pair(oldl/upper,oldu/lower),
                   'fixed_old_multiplier':str(multiplier),
                   'new_upper_below_multiplier_times_old_lower':upper < multiplier*oldl}
            out['cases'].append(row)
            save(out)
            print(json.dumps({'epsilon':str(eps), 'E_upper_display':float(upper),
                              'strict_cap_pass':cap_pass, 'versus_old':relation}), flush=True)
        out['status'] = 'complete_four_new_shorttime_bounds'
        out['seconds_this_invocation'] = time.monotonic()-begin
        out['new_E_formulas'] = 4
        out['new_Arb_calls'] = 0
        out['new_Gamma_calls'] = 0
        out['new_Green_calls'] = 0
        out['old_baseline_recomputations'] = 0
        out['root_denominator'] = str(D)
        out['new_complete_I_g_or_C5_bound'] = False
        save(out)
        (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'status':out['status'], 'pid':out['pid'],
                          'seconds':out['seconds_this_invocation']}), flush=True)
    except Exception as error:
        out['status'] = 'stopped_with_checkpoint'
        out['error'] = repr(error)
        save(out)
        raise
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
