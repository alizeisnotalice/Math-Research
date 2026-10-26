#!/usr/bin/env python3
"""One registered scalar batch; no Green evaluations or spatial sampling."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json
import os
import signal
import sys
import time

from flint import arb, ctx, fmpq

HERE = Path(__file__).resolve().parent
CP = Path('/tmp/upper_tube064/checkpoint.json')
D = 10**30


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ar(q):
    q = F(q)
    return arb(fmpq(q.numerator, q.denominator))


def save(out):
    CP.parent.mkdir(parents=True, exist_ok=True)
    temp = CP.with_suffix('.tmp')
    temp.write_text(json.dumps(out, indent=2) + '\n')
    temp.replace(CP)


def alarm(_signum, _frame):
    raise TimeoutError('Registered 120-second finite batch exhausted; no restart')


def integer_root(n, k):
    assert n >= 0 and k >= 2
    if n == 0:
        return 0
    lo, hi = 0, 1 << ((n.bit_length() + k - 1) // k)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**k <= n:
            lo = mid
        else:
            hi = mid
    assert lo**k <= n < (lo + 1)**k
    return lo


def root_interval(q, k):
    q = F(q)
    floor_scaled = (q.numerator * D**k) // q.denominator
    j = integer_root(floor_scaled, k)
    lo, hi = F(j, D), F(j + 1, D)
    assert lo**k <= q < hi**k
    return lo, hi


def ball(value):
    assert value.is_finite() and value > 0
    lo = F(str(value.lower().fmpq()))
    hi = F(str(value.upper().fmpq()))
    low = F((lo * D).numerator // (lo * D).denominator, D)
    high = F((hi * D).numerator // (hi * D).denominator + 1, D)
    assert 0 < low < lo < hi < high
    return {'ball': str(value), 'lower_dyadic': str(lo), 'upper_dyadic': str(hi),
            'strict_rational_lower': str(low), 'strict_rational_upper': str(high)}


def main():
    reg = json.loads((HERE / 'registration.json').read_text())
    assert reg['status'] == 'registered_before_computation'
    startup_recovery = None
    if CP.exists():
        assert '--resume-armed' in sys.argv, 'Existing checkpoint: inspect rather than restart'
        startup_recovery = json.loads((HERE / 'startup_recovery.json').read_text())
        prior = json.loads(CP.read_text())
        assert prior == startup_recovery['prior_armed_checkpoint']
        assert prior['status'] == 'armed' and 'constants' not in prior
        assert sha(__file__) == startup_recovery['resumed_source_sha256']
        assert prior['registration_sha256'] == sha(HERE / 'registration.json')
    sources = []
    for row in reg['input_files']:
        path = (HERE / row['path']).resolve()
        assert sha(path) == row['sha256']
        sources.append({'path': str(path), 'sha256': row['sha256']})
    out = {'id': reg['id'], 'status': 'armed', 'pid': os.getpid(),
           'python': sys.executable, 'registration_sha256': sha(HERE / 'registration.json'),
           'source_sha256': sha(__file__), 'input_files': sources}
    if startup_recovery is not None:
        out['startup_recovery_sha256'] = sha(HERE / 'startup_recovery.json')
        out['prior_armed_pid_no_numerical_calls'] = startup_recovery['prior_armed_checkpoint']['pid']
    save(out)
    print(json.dumps(out), flush=True)
    if '--arm' in sys.argv:
        input()
    signal.signal(signal.SIGALRM, alarm)
    signal.alarm(120)
    start = time.monotonic()
    try:
        ctx.prec = 256
        out['status'] = 'running_four_new_constants'
        save(out)
        pi = arb.pi()
        a, b, s, f = F(27, 10000), F(9, 25000), F(15, 2), F(23, 10)
        J0 = 8*pi/3 * pi.sqrt() * ar(F(3, 4)).gamma()/ar(F(5, 4)).gamma()
        m2 = (arb(2).log() + 2 - pi/2)**2 + ar(F(4, 9))
        K = ar(5)**ar(F(2, 5)) * ar(15)**ar(F(4, 5)) * ar(2)**ar(F(6, 5))
        K *= J0**ar(F(4, 5)) * (ar(a) + ar(b*s)*m2.sqrt())
        ell = ar(F(9, 4)).log()
        M5 = 512*pi*sum((ar(F(factorial(5), factorial(5-k)*3**(k+1)))*ell**(5-k)
                          for k in range(6)), arb(0))
        out['constants'] = {name: ball(value) for name, value in
                            [('J0', J0), ('m2_upper', m2), ('K', K), ('M5', M5)]}
        out['constants_seconds'] = time.monotonic() - start
        out['status'] = 'running_thirteen_exact_radius_comparisons'
        out['radius_comparisons'] = []
        save(out)
        Klo, Khi = [F(out['constants']['K'][key]) for key in
                    ['strict_rational_lower', 'strict_rational_upper']]
        Mlo, Mhi = [F(out['constants']['M5'][key]) for key in
                    ['strict_rational_lower', 'strict_rational_upper']]
        target4 = F(1, 10**24)
        for n in range(10, 23):
            delta = F(1, 2**n)
            V = 108*delta**3 - 243*delta**4
            assert 0 < delta < F(1, 3) and 0 < V
            dL, dU = root_interval(delta**6, 5)
            vL, vU = root_interval(V**4, 5)
            wL, _ = root_interval(Mlo*V**3, 5)
            _, wU = root_interval(Mhi*V**3, 5)
            SL = Klo*dL + 30*a*vL + 30*b*s*wL
            SU = Khi*dU + 30*a*vU + 30*b*s*wU
            low4, high4 = (f*SL)**5, (f*SU)**5
            succeeds = high4 <= target4
            fails = low4 > target4
            assert succeeds != fails, 'Rational enclosure did not decide the preselected comparison'
            lowI, _ = root_interval(low4, 4)
            _, highI = root_interval(high4, 4)
            out['radius_comparisons'].append({
                'n': n, 'delta': str(delta), 'tube_volume': str(V),
                'norm_bracket_lower': str(SL), 'norm_bracket_upper': str(SU),
                'integral_fourth_power_lower': str(low4),
                'integral_fourth_power_upper': str(high4),
                'target_fourth_power': str(target4),
                'integral_lower': str(lowI), 'integral_upper': str(highI),
                'integral_lower_display': float(lowI), 'integral_upper_display': float(highI),
                'passes_1e_minus_6': succeeds,
                'explicit_envelope_proven_above_target': fails})
            save(out)
        passing = [x for x in out['radius_comparisons'] if x['passes_1e_minus_6']]
        out['largest_passing_preregistered_radius'] = passing[0]['delta'] if passing else None
        out['largest_passing_n'] = passing[0]['n'] if passing else None
        out['status'] = 'complete_tube_majorant_only_no_I_g_or_C5_claim'
        out['new_Arb_expression_evaluations'] = 4
        out['new_Green_evaluations'] = 0
        out['seconds'] = time.monotonic() - start
        save(out)
        (HERE / 'results.json').write_text(json.dumps(out, indent=2) + '\n')
        print(json.dumps({'status': out['status'], 'pid': out['pid'],
                          'seconds': out['seconds'],
                          'largest_passing_radius': out['largest_passing_preregistered_radius'],
                          'comparisons': [{k: r[k] for k in ['n', 'integral_upper_display', 'passes_1e_minus_6']}
                                          for r in out['radius_comparisons']]}), flush=True)
    except Exception as error:
        out['status'] = 'stopped_with_diagnostic_no_restart'
        out['error'] = repr(error)
        save(out)
        raise
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
