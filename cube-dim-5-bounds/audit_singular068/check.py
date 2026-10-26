#!/usr/bin/env python3
"""Read-only exact audit. No producer imports, Arb, Gamma, geometry or Green calls."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import os
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
SCALE = 1 << 128


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def root(q, degree):
    """Greedily construct binary digits, testing the rational directly."""
    assert q >= 0 and degree >= 2
    numerator = q.numerator * SCALE ** degree
    denominator = q.denominator
    top_bit = max(0, (numerator.bit_length() - denominator.bit_length()) // degree + 2)
    acc = 0
    for bit in range(top_bit, -1, -1):
        candidate = acc + (1 << bit)
        if candidate ** degree * denominator <= numerator:
            acc = candidate
    low, high = Q(acc, SCALE), Q(acc + 1, SCALE)
    assert low ** degree <= q < high ** degree
    return low, high


def pair(row):
    lo, hi = Q(row['lower_dyadic']), Q(row['upper_dyadic'])
    assert 0 < lo < hi
    for x in (lo, hi):
        assert x.denominator & (x.denominator - 1) == 0
    if 'strict_rational_lower' in row:
        assert Q(row['strict_rational_lower']) < lo
        assert hi < Q(row['strict_rational_upper'])
    return lo, hi


def outward_decimal(q, places, upper):
    scale = 10 ** places
    z = q.numerator * scale // q.denominator
    if upper:
        z += 1
    return str(Q(z, scale))


def main():
    start = time.monotonic()
    reg = read(HERE / 'registration.json')
    assert reg['status'] == 'registered_before_exact_computation'
    recovery = read(HERE / 'schema_recovery.json')
    assert '--resume-schema' in sys.argv
    assert recovery['original_source_sha256'] == reg['audit_source_sha256']
    assert recovery['resumed_source_sha256'] == digest(__file__)
    assert read(HERE / 'results.json') == recovery['prior_full_failure_receipt']
    out = {'status': 'armed', 'pid': os.getpid(), 'python': sys.executable,
           'registration_sha256': digest(HERE / 'registration.json'),
           'audit_source_sha256': digest(__file__),
           'schema_recovery_sha256': digest(HERE / 'schema_recovery.json')}
    print(json.dumps(out), flush=True)
    if '--arm' in sys.argv:
        input()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('120-second audit budget')))
    signal.alarm(120)
    try:
        verified = []
        for item in reg['inputs']:
            assert digest(item['path']) == item['sha256'], item['path']
            verified.append(item)
        manifests = []
        for name in ('057', '064', '045'):
            path = Path(reg['roles'][name]) / 'manifest.json'
            manifest = read(path)
            for item in manifest['white_list' if name == '045' else 'files']:
                actual = path.parent / item['path']
                assert actual.stat().st_size == item['bytes'], str(actual)
                assert digest(actual) == item['sha256'], str(actual)
                manifests.append({'path': str(actual), 'sha256': item['sha256']})
        d57 = read(Path(reg['roles']['057']) / 'results.json')
        d64 = read(Path(reg['roles']['064']) / 'results.json')
        d45 = read(Path(reg['roles']['045']) / 'results.json')
        for dataset, key in ((d57, 'external_sources'), (d64, 'input_files')):
            for item in dataset[key]:
                assert digest(item['path']) == item['sha256'], item['path']
        c57 = {name: pair(value) for name, value in d57['constants']['constants'].items()}
        assert len(c57) == 8
        for name, cap in d57['constants']['strict_rational_caps'].items():
            assert c57[name][1] < Q(cap)
        c64 = {name: pair(value) for name, value in d64['constants'].items()}
        assert set(c64) == {'J0', 'm2_upper', 'K', 'M5'}
        # Rebuild the algebraic K expression from complete saved J0/m2 endpoints.
        a, b, s, f = Q(27, 10000), Q(9, 25000), Q(15, 2), Q(23, 10)
        algebraic_K = []
        for j in (0, 1):
            sqrt_m = root(c64['m2_upper'][j], 2)[j]
            radicand = 5**2 * 15**4 * 2**6 * c64['J0'][j]**4
            algebraic_K.append(root(radicand, 5)[j] * (a + b*s*sqrt_m))
        assert algebraic_K[0] < c64['K'][0] < c64['K'][1] < algebraic_K[1]
        rows = []
        assert [r['n'] for r in d64['radius_comparisons']] == list(range(10, 23))
        for saved in d64['radius_comparisons']:
            n = saved['n']
            delta = Q(1, 1 << n)
            # Directly add the five disjoint k-layer volumes, not either expanded polynomial.
            V = delta**4
            for k, choose in ((1, 4), (2, 6), (3, 4), (4, 1)):
                V += choose*delta**(4-k)*(k*(1-delta)*(2*delta)**(k-1)-(k-1)*(2*delta)**k)
            assert V == Q(saved['tube_volume'])
            rdelta = root(delta**6, 5)
            rvolume = root(V**4, 5)
            rM = (root(c64['M5'][0]*V**3, 5)[0], root(c64['M5'][1]*V**3, 5)[1])
            brackets = tuple(c64['K'][j]*rdelta[j] + 30*a*rvolume[j] + 30*b*s*rM[j] for j in (0, 1))
            powers = tuple((f*x)**5 for x in brackets)
            target = Q(1, 10**24)
            passed, failed = powers[1] < target, powers[0] > target
            assert passed != failed
            assert passed == saved['passes_1e_minus_6']
            assert failed == saved['explicit_envelope_proven_above_target']
            Ilow, Ihigh = root(powers[0], 4)[0], root(powers[1], 4)[1]
            # Audit the author's independent rational intervals by overlap, without importing its algorithm.
            assert max(Ilow, Q(saved['integral_lower'])) < min(Ihigh, Q(saved['integral_upper']))
            assert Q(saved['integral_fourth_power_lower']) == (f*Q(saved['norm_bracket_lower']))**5
            assert Q(saved['integral_fourth_power_upper']) == (f*Q(saved['norm_bracket_upper']))**5
            rows.append({'n': n, 'delta': str(delta), 'V': str(V),
                         'roots': {'delta6': list(map(str, rdelta)), 'V4': list(map(str, rvolume)), 'M5V3': list(map(str, rM))},
                         'norm_bracket': list(map(str, brackets)), 'integral_interval': [str(Ilow), str(Ihigh)],
                         'display_outward': [outward_decimal(Ilow, 22, False), outward_decimal(Ihigh, 22, True)],
                         'passes': passed, 'explicit_envelope_fails': failed})
        old = next(r for r in d45['rows'] if r['r'] == '8/5')
        oldlo, oldhi = pair(old['K'])
        assert Q(old['strict_rational_K_lower']) < oldlo < oldhi < Q(old['strict_rational_K_cap'])
        olddelta = Q(old['delta'])
        assert olddelta == Q(1, 1 << 22)
        def volume(d):
            return 108*d**3-243*d**4
        oldpass = oldhi**25 * volume(olddelta)**7 < Q(1, 10**192)
        oldfail = oldlo**25 * volume(2*olddelta)**7 > Q(1, 10**192)
        assert oldpass and oldfail
        assert Q(rows[6]['delta']) / olddelta == 64
        # Compare both explicit envelopes at the shared old radius by powers only.
        new_upper = Q(rows[-1]['integral_interval'][1])
        old_expression_lower_power32 = oldlo**25 * volume(olddelta)**7
        assert new_upper**32 < old_expression_lower_power32
        out.update({'status': 'passed_readonly_exact_audit', 'inputs': verified,
                    'manifest_entries_verified': manifests, 'endpoints057': {k:list(map(str,v)) for k,v in c57.items()},
                    'endpoints064': {k:list(map(str,v)) for k,v in c64.items()},
                    'K_reconstructed_interval': list(map(str, algebraic_K)), 'rows': rows,
                    'old045': {'r': '8/5', 'delta': str(olddelta), 'adjacent_larger_delta': str(2*olddelta),
                               'passes': oldpass, 'adjacent_envelope_fails': oldfail,
                               'radius_improvement_factor': 64, 'new_envelope_strictly_smaller_at_old_delta': True},
                    'root_algorithm': 'greedy binary digits against exact rational; denominator 2^128',
                    'new_special_function_calls': 0, 'new_Green_calls': 0, 'new_geometry_enumerations': 0,
                    'seconds': time.monotonic() - start})
        for item in reg['inputs']:
            assert digest(item['path']) == item['sha256'], item['path']
    except Exception as exc:
        out.update({'status': 'failed_stop_no_rerun', 'error': repr(exc), 'seconds': time.monotonic()-start})
        raise
    finally:
        signal.alarm(0)
        (HERE / 'results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ('status','pid','seconds','old045')}), flush=True)


if __name__ == '__main__':
    main()
