#!/usr/bin/env python3
"""Eight exact checks only; no imported author code or kernel invocation.

Default mode is read-only and compares to the frozen results. --record is a
single-use writer of this audit's own results/receipt, never an author writer.
"""
import sys
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
import argparse
import datetime
from fractions import Fraction as F
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import time

HERE = Path(__file__).resolve().parent
START = time.monotonic()
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def hashes():
    m = json.loads((HERE / 'input_manifest.json').read_text())
    for p, h in m['sources'].items():
        assert sha(p) == h, ('source changed', p)
    return len(m['sources'])


def invert_det(a):
    n = len(a)
    b = [[F(v) for v in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    det = F(1)
    for j in range(n):
        pivot = next(i for i in range(j, n) if b[i][j])
        if pivot != j:
            b[pivot], b[j] = b[j], b[pivot]
            det = -det
        q = b[j][j]
        det *= q
        b[j] = [v / q for v in b[j]]
        for i in range(n):
            if i != j:
                q = b[i][j]
                b[i] = [u - q * v for u, v in zip(b[i], b[j])]
    return [row[n:] for row in b], det


def mv(a, x):
    return [sum((u * v for u, v in zip(row, x)), F(0)) for row in a]


def interval(d):
    e = int(d['common_binary_exponent'])
    scale = F(2) ** e
    lo, hi = F(int(d['lower_mantissa'])) * scale, F(int(d['upper_mantissa'])) * scale
    assert lo <= hi
    return lo, hi


def times(iv, q):
    a, b = iv[0] * q, iv[1] * q
    return min(a, b), max(a, b)


def add(ivs):
    return sum((v[0] for v in ivs), F(0)), sum((v[1] for v in ivs), F(0))


def coarse(iv, den=10**18):
    a, b = iv
    return [str(F((a * den).__floor__(), den)), str(F((b * den).__ceil__(), den))]


def variances(x):
    def w(y):
        return sum(t*t for t in y) - sum(y)**2 / len(y)
    return [w(x[:k]) + w(x[k:]) for k in range(1, 6)]


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, (list, tuple)):
        return [encode(y) for y in x]
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    return x


def compute():
    source_count = hashes()
    reg = json.loads((HERE/'registration.json').read_text())
    old = Path(reg['source_roots']['A038'])
    cases = []
    frame = [(1,) + s for s in itertools.product((-1, 1), repeat=5)]
    one = (1,)*6
    def flip(j):
        return tuple(-1 if i == j else 1 for i in range(6))
    cols = [one] + [flip(j) for j in (4, 5, 1, 2, 3)]
    a = [list(row) for row in zip(*cols)]
    inv, det = invert_det(a)
    assert abs(det) == 32
    # Independent affine coefficient equality for all b,c, not sample values.
    const = mv(inv, [1, 1, 1, 1, 0, 0])
    bc = [mv(inv, [F(i == j) for i in range(6)]) for j in (4, 5)]
    assert const == [0, F(1,2), F(1,2), 0, 0, 0]
    assert bc == [[F(1,2), F(-1,2), 0, 0, 0, 0],
                  [F(1,2), 0, F(-1,2), 0, 0, 0]]
    cases.append({'id':1, 'status':'PASS', 'determinant':det, 'inverse':inv,
                  'family_constant_coefficients':const, 'family_b_c_coefficients':bc})

    rest = [v for v in frame if v not in cols]
    transformed = [mv(inv, v) for v in rest]
    rowsums = [sum(abs(v[i]) for v in transformed) for i in range(6)]
    assert len(rest) == 26 and rowsums == [49, 15, 15, 15, 15, 15]
    cases.append({'id':2, 'status':'PASS', 'remaining_columns':26,
                  'absolute_row_sums':rowsums})

    a0 = F(1,100)**2 + F(149,100)**2
    a1 = F(1,100)**2 + F(13,20)**2
    assert a0 == F(11101,5000) and a1 == F(2113,5000)
    assert 1 + F(49,100) == F(149,100)
    assert F(1,2) + F(15,100) == F(13,20)
    prefactor = F(1, 32 * 2**26 * 226**3) / a0 / a1**2
    assert prefactor == F(5000**3, 2**31 * 226**3 * 11101 * 2113**2)
    # pi = 4 int_0^1 (1+x^2)^-1 dx < 4: no decimal pi used.
    strict_lower = prefactor / 4**6
    assert strict_lower > F(1,10**20)
    cases.append({'id':3, 'status':'PASS', 'a0':a0, 'a1':a1,
                  'constant_times_pi6':prefactor, 'strict_cstar_lower_via_pi_lt_4':strict_lower,
                  'strict_lower_minus_1e_minus20':strict_lower-F(1,10**20)})

    xstar = list(map(F, (1,1,1,1))) + [F(1,2), F(0)]
    qstar = variances(xstar)
    assert min(qstar) == F(1,8) > F(1,125000)**2
    cases.append({'id':4, 'status':'PASS', 'point':xstar,
                  'five_squared_distances':qstar, 'd2_squared':min(qstar)})

    delta = F(1,125000)
    floor = (1/delta - 40000)/5000000000
    native3 = json.loads((old/'native_outputs/03.stdout.json').read_text())
    g3 = interval(native3['G'])
    point3 = json.loads((old/'candidate_03.json').read_text())
    x3 = [F(n, point3['denominator']) for n in point3['final_numerators']]
    assert floor == F(17,10**6) > F(1,65000)
    assert F(1,65000) < g3[0] <= g3[1] < floor
    assert g3[1] < F(100033915510893,6250000000000000000)
    assert min(variances(x3)) > delta**2
    cases.append({'id':5, 'status':'PASS', 'delta':delta, 'tube_floor':floor,
                  'floor_minus_target':floor-F(1,65000), 'saved_point3_G':g3,
                  'saved_point3_G_coarse':coarse(g3), 'floor_minus_G3_upper':floor-g3[1],
                  'point3_min_q':min(variances(x3))})

    w = list(map(F, (1,1,1,1,1,0)))
    for i in (1,2,3):
        pt = json.loads((old/f'candidate_{i:02}.json').read_text())
        raw = json.loads((old/f'native_outputs/{i:02}.stdout.json').read_text())
        x = [F(n, pt['denominator']) for n in pt['final_numerators']]
        assert x[0] == 1 and all(x[j] >= x[j+1] for j in range(5)) and x[5] >= 0
        assert raw['status'] == 'generic_nonzero_supported'
        assert raw['precision_bits'] == 512 and raw['full_group_terms'] == 967680
        assert raw['computed_sign_pairs'] == 483840 and raw['zero_node_terms'] == 0
        g, grad = interval(raw['G']), [interval(d) for d in raw['gradient']]
        assert len(grad) == 6 and all(b-a < F(2)**-256 for a,b in [g]+grad)
        h = [b-a for a,b in zip(x,w)]
        assert h[0] == 0
        dot = add([times(v, u) for v,u in zip(grad,h)])
        euler = add([times(v,u) for v,u in zip(grad,x)] + [times(g,F(5))])
        assert dot[1] < 0 and euler[0] <= 0 <= euler[1]
        assert min(variances(x)) > delta**2
        cases.append({'id':5+i, 'status':'PASS', 'point':x, 'direction_to_w':h,
                      'all_six_derivative_intervals':grad, 'direction_interval':dot,
                      'direction_interval_coarse':coarse(dot), 'euler_interval':euler,
                      'conclusion':'This saved point is not a local minimum in the ordered face; no neighborhood is deleted.'})
        assert time.monotonic()-START < 900
    assert len(cases) == 8
    assert hashes() == source_count
    return encode({'status':'PASS_EIGHT_EXACT_CHECKS', 'case_count':8,
                   'source_count':source_count, 'cases':cases,
                   'new_G_H_p_native':0, 'new_numerical_kernel':0,
                   'author_writers_executed':0, 'new_C6_bound':False,
                   'global_floor_1_over_65000_proved':False})


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--record', action='store_true')
    args = p.parse_args()
    if args.record:
        assert not (HERE/'results.json').exists(), 'one registered run only'
    result = compute()
    if args.record:
        (HERE/'results.json').write_text(json.dumps(result, indent=2)+'\n')
        receipt = {'status':'COMPLETED', 'pid':os.getpid(), 'exit_code':0,
                   'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   'wall_seconds':time.monotonic()-START, 'configuration_count':8,
                   'source_sha256':sha(__file__), 'registration_sha256':sha(HERE/'registration.json'),
                   'input_manifest_sha256':sha(HERE/'input_manifest.json'),
                   'runtime':{'python':sys.version, 'executable':sys.executable,
                              'executable_sha256':sha(sys.executable),
                              'fractions_module':__import__('fractions').__file__,
                              'fractions_sha256':sha(__import__('fractions').__file__),
                              'platform':platform.platform(), 'third_party_dependencies':[]},
                   'new_G_H_p_native':0, 'child_processes':[], 'all_handles_closed':True}
        (HERE/'run_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    else:
        assert result == json.loads((HERE/'results.json').read_text()), 'saved result mismatch'
    print(json.dumps({'status':result['status'], 'pid':os.getpid(),
                      'seconds':time.monotonic()-START,
                      'record_mode':args.record, 'case_count':8,
                      'pointwise_direction_outer_bounds':[c['direction_interval_coarse'] for c in result['cases'][5:]]}))


if __name__ == '__main__':
    main()
