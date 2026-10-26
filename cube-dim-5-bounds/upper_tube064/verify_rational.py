#!/usr/bin/env python3
"""Independent saved-evidence verification, standard library only."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import signal
import time

HERE = Path(__file__).resolve().parent
D = 10**24


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def timeout(_signum, _frame):
    raise TimeoutError('120-second independent rational verification limit')


def fifth_root(q):
    scaled = q.numerator * D**5 // q.denominator
    if scaled == 0:
        r = 0
    else:
        r = 1 << ((scaled.bit_length() + 4) // 5)
        while True:
            nxt = (4*r + scaled//r**4)//5
            if nxt >= r:
                break
            r = nxt
        while r**5 > scaled:
            r -= 1
        while (r+1)**5 <= scaled:
            r += 1
    lo, hi = F(r, D), F(r+1, D)
    assert lo**5 <= q < hi**5
    return lo, hi


def main():
    begin = time.monotonic()
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(120)
    vr = json.loads((HERE/'verification_registration.json').read_text())
    assert vr['status'] == 'registered_before_verification'
    data = json.loads((HERE/'results.json').read_text())
    assert data['registration_sha256'] == sha(HERE/'registration.json')
    assert data['source_sha256'] == sha(HERE/'certify.py')
    assert data['startup_recovery_sha256'] == sha(HERE/'startup_recovery.json')
    for row in data['input_files']:
        assert sha(row['path']) == row['sha256']
    caps = {}
    for name in ['K', 'M5']:
        row = data['constants'][name]
        lo, hi = F(row['strict_rational_lower']), F(row['strict_rational_upper'])
        assert 0 < lo < F(row['lower_dyadic']) < F(row['upper_dyadic']) < hi
        caps[name] = (lo, hi)
    a, b, s, f = F(27, 10000), F(9, 25000), F(15, 2), F(23, 10)
    rows = []
    assert [r['n'] for r in data['radius_comparisons']] == list(range(10, 23))
    for source in data['radius_comparisons']:
        delta = F(1, 2**source['n'])
        V = 27*delta**3*(4-9*delta)
        assert V == F(source['tube_volume'])
        terms = {'delta6': fifth_root(delta**6), 'volume4': fifth_root(V**4),
                 'Mlower_volume3': fifth_root(caps['M5'][0]*V**3),
                 'Mupper_volume3': fifth_root(caps['M5'][1]*V**3)}
        lo = caps['K'][0]*terms['delta6'][0] + 30*a*terms['volume4'][0]
        lo += 30*b*s*terms['Mlower_volume3'][0]
        hi = caps['K'][1]*terms['delta6'][1] + 30*a*terms['volume4'][1]
        hi += 30*b*s*terms['Mupper_volume3'][1]
        target = F(1, 10**24)
        passes = (f*hi)**5 <= target
        fails = (f*lo)**5 > target
        assert passes != fails
        assert passes == source['passes_1e_minus_6']
        assert fails == source['explicit_envelope_proven_above_target']
        rows.append({'n':source['n'], 'delta':str(delta),
                     'root_intervals':{k:[str(v[0]), str(v[1])] for k,v in terms.items()},
                     'norm_bracket_lower':str(lo), 'norm_bracket_upper':str(hi),
                     'passes':passes, 'envelope_definitely_fails':fails})
    assert rows[5]['envelope_definitely_fails'] and rows[6]['passes']
    out = {'id':vr['id'], 'status':'all_thirteen_comparisons_independently_verified',
           'pid':os.getpid(), 'registration_sha256':sha(HERE/'verification_registration.json'),
           'source_sha256':sha(__file__), 'input_results_sha256':sha(HERE/'results.json'),
           'root_denominator':str(D), 'rows':rows, 'largest_passing_delta':'1/65536',
           'new_Arb_calls':0, 'new_Green_calls':0, 'seconds':time.monotonic()-begin}
    signal.alarm(0)
    (HERE/'rational_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','pid','seconds','largest_passing_delta']}))


if __name__ == '__main__':
    main()
