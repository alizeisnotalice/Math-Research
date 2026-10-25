#!/usr/bin/env python3
"""Three registered scalar Arb evaluations, then exact rational threshold tests."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import signal
import sys
import time

from flint import arb, ctx, fmpq

HERE = Path(__file__).resolve().parent
CHECKPOINT = Path('/tmp/astra_alpha1_Lp045/checkpoint.json')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ar(x):
    x = F(x)
    return arb(fmpq(x.numerator, x.denominator))


def endpoint(x, upper):
    return F(str((x.upper() if upper else x.lower()).fmpq()))


def enclosing_record(x):
    assert x.is_finite() and x > 0
    return {'ball': str(x), 'lower_dyadic': str(endpoint(x, False)),
            'upper_dyadic': str(endpoint(x, True))}


def volume(delta):
    assert 0 <= delta <= F(1, 3)
    return 108*delta**3-243*delta**4


def save_checkpoint(value):
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    temp = CHECKPOINT.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2)+'\n')
    temp.replace(CHECKPOINT)


def deadline(_signum, _frame):
    raise TimeoutError('Registered 120-second scalar checkpoint reached; no automatic retry')


def one(r, powers):
    start = time.monotonic()
    gamma = (r-1)/(2*r)
    k = 5*gamma-F(1, 2)
    assert F(5, 4) < r < F(5, 3) and 0 < k < F(1, 2)
    A, B, D = powers
    assert F(A, D) == F(5, 4)/r and F(B, D) == 1-F(5, 4)/r
    # The sole full transcendental evaluation for this registered r.
    pi = arb.pi()
    moment = (8*pi.sqrt()*arb(16)**ar(5*gamma-1) *
              ar(1-2*k).gamma()/ar(1-k).gamma() *
              (ar(1+2*gamma).gamma()/ar(1+gamma).gamma())**5 *
              ar(F(2661, 64))**ar(-gamma))
    K = ((1/(32*pi**3))**ar(r) * pi**2/60 *
         arb(5)**ar(2*r-F(5, 2)) * moment**ar(r))
    lower, upper = endpoint(K, False), endpoint(K, True)
    scale = 10**30
    cap = F((upper*scale).__floor__()+1, scale)
    floor = F((lower*scale).__floor__(), scale)
    assert 0 < floor <= lower <= upper < cap
    assert K < ar(cap) and ar(floor) < K
    rhs = F(1, 10**(6*D))
    trials = []

    def passes(n):
        delta = F(1, 2**n)
        lhs = cap**A*volume(delta)**B
        ok = lhs <= rhs
        trials.append({'n': n, 'passes': ok})
        return ok

    if passes(2):
        n = 2
    else:
        lo, hi = 2, 4
        while not passes(hi):
            lo, hi = hi, 2*hi
            if hi > 128:
                raise ArithmeticError('No pass within registered exponent bracket')
        while hi-lo > 1:
            mid = (lo+hi)//2
            if passes(mid):
                hi = mid
            else:
                lo = mid
        n = hi
    delta = F(1, 2**n)
    V = volume(delta)
    lhs = cap**A*V**B
    assert lhs <= rhs
    assert n > 2  # Otherwise the next larger dyadic is outside delta<=1/3.
    previous_V = volume(2*delta)
    previous_lhs = cap**A*previous_V**B
    assert previous_lhs > rhs
    # Certify that the preceding radius also fails for the explicit K expression,
    # not merely because of the upward rational rounding.
    previous_lower_lhs = floor**A*previous_V**B
    assert previous_lower_lhs > rhs
    return {'r': str(r), 'gamma': str(gamma), 'trace_k': str(k),
            'moment_cap': enclosing_record(moment), 'K': enclosing_record(K),
            'strict_rational_K_cap': str(cap), 'strict_rational_K_lower': str(floor),
            'integer_power_test_A_B_D': [A, B, D], 'target_rhs': str(rhs),
            'binary_exponent_trials': trials,
            'maximum_certified_dyadic_radius_exponent': n,
            'delta': str(delta), 'tube_volume': str(V),
            'exact_pass_ratio_lhs_over_rhs': str(lhs/rhs),
            'previous_delta': str(2*delta), 'previous_tube_volume': str(previous_V),
            'exact_previous_fail_ratio_lhs_over_rhs': str(previous_lhs/rhs),
            'previous_radius_fails_even_using_K_lower': True,
            'exact_previous_lower_fail_ratio': str(previous_lower_lhs/rhs),
            'improvement_factor_over_locked_2pow_minus23': str(delta/F(1, 2**23)),
            'seconds': time.monotonic()-start}


def main():
    regfile = HERE/'registration.json'
    reg = json.loads(regfile.read_text())
    assert reg['status'] == 'registered_before_computation'
    assert reg['new_r_values'] == ['8/5', '13/8', '33/20']
    for row in reg['sources']:
        assert sha(row['path']) == row['sha256'], row['path']
    old = json.loads(Path(reg['sources'][1]['path']).read_text())
    link = json.loads(Path(reg['sources'][3]['path']).read_text())
    assert old['strict_rational_K_upper'] == '1/3840'
    assert link['chosen_delta'] == '1/8388608'
    assert link['target_tube_Iq_error'] == '1/1000000'
    assert sum(w*s for w,s in zip(map(F, ['1/8','3/16','3/16','1/8','3/8']),
                                  map(F, ['81/8','27','36','36','64']))) == F(2661,64)
    ctx.prec = reg['precision_bits']
    result = {'status': 'running', 'id': reg['id'], 'pid': os.getpid(),
              'python': sys.executable, 'precision_bits': ctx.prec,
              'registration_sha256': sha(regfile), 'source_sha256': sha(__file__),
              'external_sources': reg['sources'], 'rows': [],
              'old_r_3over2_recomputed': False, 'old_delta': '1/8388608',
              'new_G_or_I_or_C5_bound': False}
    save_checkpoint(result)
    print(json.dumps({'status':'armed' if '--arm' in sys.argv else 'starting',
                      'pid':os.getpid(), 'registration_sha256':sha(regfile),
                      'source_sha256':sha(__file__), 'checkpoint':str(CHECKPOINT)}), flush=True)
    if '--arm' in sys.argv:
        input()
    signal.signal(signal.SIGALRM, deadline)
    for rs in reg['new_r_values']:
        result['active_r'] = rs
        save_checkpoint(result)
        signal.alarm(120)
        try:
            row = one(F(rs), reg['integer_power_tests'][rs])
        except Exception as error:
            result['status'] = 'stopped_with_diagnostic'
            result['diagnostic'] = repr(error)
            save_checkpoint(result)
            raise
        finally:
            signal.alarm(0)
        result['rows'].append(row)
        save_checkpoint(result)
        print(json.dumps({key:row[key] for key in ('r','strict_rational_K_cap','delta',
                   'improvement_factor_over_locked_2pow_minus23','seconds')}), flush=True)
    del result['active_r']
    result['status'] = 'three_registered_Lr_caps_and_maximal_dyadic_envelopes_certified'
    save_checkpoint(result)
    (HERE/'results.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
