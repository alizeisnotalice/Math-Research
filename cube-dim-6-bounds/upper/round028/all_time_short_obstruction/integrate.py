#!/usr/bin/env python3
"""One released batch: exactly the two locked root-upper PF integrals, Arb640."""
import sys
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
import datetime, hashlib, json, math, os, signal, time, traceback
from fractions import Fraction as Q
from pathlib import Path
from flint import arb, ctx, fmpq

HERE = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(name, value):
    tmp = HERE / (name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    os.replace(tmp, HERE / name)

def exact_ball(value):
    value = Q(value)
    return arb(fmpq(value.numerator, value.denominator))

def endpoints(value):
    assert value.is_finite()
    lo, le = value.lower().man_exp()
    hi, he = value.upper().man_exp()
    return dict(lower_mantissa=str(lo), lower_binary_exponent=int(le),
                upper_mantissa=str(hi), upper_binary_exponent=int(he))

def bounds(row):
    lo = Q(int(row['lower_mantissa'])) * Q(2) ** row['lower_binary_exponent']
    hi = Q(int(row['upper_mantissa'])) * Q(2) ** row['upper_binary_exponent']
    assert lo <= hi
    return lo, hi

def table(path):
    data = json.loads(Path(path).read_text())
    assert data['normalization'] == 'p_T=(2*pi)^-6 R(T^2)'
    for row in data['zero_pole_channels']:
        if row['order'] % 2 == 0:
            assert Q(row['A']) == 0
    rows = []
    for row in data['real_z_principal']:
        s, coefficient = Q(row['b_squared']), Q(row['coefficient'])
        assert s > 0 and coefficient != 0
        n, d = math.isqrt(s.numerator), math.isqrt(s.denominator)
        assert n * n == s.numerator and d * d == s.denominator
        b, order = Q(n, d), row['order']
        assert b * b == s and type(order) is int and 1 <= order <= 3
        rows.append((s, b, order, coefficient))
    assert len(rows) == 417 and len({row[0] for row in rows}) == 359
    assert len({(row[0], row[2]) for row in rows}) == 417
    assert max(row[2] for row in rows) == 3
    return rows

def integral_at(t_exact, rows):
    """All 417 signed real-principal entries, no native, G, H, or p evaluation."""
    ctx.prec = 640
    t = exact_ball(t_exact)
    total, cache = arb(0), {}
    for s, b, order, coefficient in rows:
        if s not in cache:
            sb, bb = exact_ball(s), exact_ball(b)
            values = [None, (t / bb).atan() / bb]
            for j in (2, 3):
                values.append(t / (2 * (j - 1) * sb * (t * t + sb) ** (j - 1))
                              + (2 * j - 3) * values[j - 1] / (2 * (j - 1) * sb))
            cache[s] = values
        total += exact_ball(coefficient) * cache[s][order]
    assert len(cache) == 359
    return endpoints(total / (64 * arb.pi() ** 6))

def run():
    assert not (HERE / 'run_receipt.json').exists(), 'one batch already started'
    reg = json.loads((HERE / 'registration.json').read_text())
    assert reg['max_positive_time_PF_values'] == 2 and reg['precision_bits'] == 640
    assert reg['point'] == ['1', '1', '1', '9/10', '27/50', '0']
    for path, digest in reg['source_sha256'].items():
        assert sha(path) == digest, ('changed source', path)
    for name, digest in reg['local_code_sha256'].items():
        assert sha(HERE / name) == digest, ('changed prepared code', name)
    constants = json.loads(Path(reg['constants_path']).read_text())
    roots = json.loads(Path(reg['roots_path']).read_text())['roots']
    assert len(roots) == 2
    assert [r['T_upper'] for r in roots] == ['15536257/1073741824', '185807383/4294967296']
    assert [r['T_lower'] for r in roots] == ['62145027/4294967296', '92903691/2147483648']
    c, ell, upper_pi = map(Q, (constants['c'], constants['ell'], constants['upper_pi']))
    assert c == Q(1, 4100) and ell == Q(3141592653, 10**9) and upper_pi == Q(3141592654, 10**9)
    rho = (upper_pi / ell) ** 6
    rows = table(reg['partial_table_path'])
    start = time.perf_counter()
    receipt = dict(status='RUNNING', pid=os.getpid(), start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   registration_sha256=sha(HERE / 'registration.json'), positive_time_PF_values=0,
                   max_wall_seconds=1200, new_native=0, new_H=0, new_G_full=0,
                   extra_positive_times=0, child_processes=0)
    save('run_receipt.json', receipt)
    def timeout(signum, frame):
        raise TimeoutError('registered 1200 second batch limit')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(1200)
    result = dict(status='INCOMPLETE', precision_bits=640, point=reg['point'], c=str(c),
                  ell=str(ell), upper_pi=str(upper_pi), rho_upper=str(rho),
                  principal_entries=417, positive_rational_square_roots=359, max_order=3, evaluations=[])
    code = 1
    try:
        for index, root in enumerate(roots, 1):
            tic = time.perf_counter()
            a, b = Q(root['T_lower']), Q(root['T_upper'])
            kd = integral_at(b, rows)
            lo, hi = bounds(kd)
            assert 0 < lo <= hi
            payment = c * a - rho * hi
            payment_lower = Q((payment * 10**18).__floor__(), 10**18)
            row = dict(root_index=index, T_lower=str(a), T_upper=str(b), K_at_T_upper=kd,
                       payment_lower_exact=str(payment), payment_simple_lower=str(payment_lower),
                       payment_strictly_positive=(payment > 0), elapsed_seconds=time.perf_counter()-tic,
                       K_mid_display=float((lo + hi) / 2), payment_display=float(payment))
            result['evaluations'].append(row)
            receipt['positive_time_PF_values'] += 1
            save('results.json', result)
            save('run_receipt.json', receipt)
            print(json.dumps({k:v for k,v in row.items() if k not in ('K_at_T_upper','payment_lower_exact')}), flush=True)
            assert payment > 0, 'first payment not separated: stop, no retry or new input'
        result['status'] = 'PASS_TWO_LOCKED_CRITICAL_PAYMENTS'
        result['all_time_consequence_requires_report_proof'] = True
        result['new_C6_bound'] = False
        save('results.json', result)
        receipt['status'], code = 'COMPLETED', 0
    except BaseException as exc:
        receipt['status'], receipt['error'] = 'STOPPED_FAILURE', repr(exc)
        receipt['traceback'] = traceback.format_exc()
    finally:
        signal.alarm(0)
        receipt['elapsed_seconds'] = time.perf_counter() - start
        receipt['exit_code'] = code
        receipt['end_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        if (HERE / 'results.json').exists():
            receipt['results_sha256'] = sha(HERE / 'results.json')
        save('run_receipt.json', receipt)
    return code

if __name__ == '__main__':
    assert sys.argv[1:] == ['--run'], 'explicit --run required'
    raise SystemExit(run())
