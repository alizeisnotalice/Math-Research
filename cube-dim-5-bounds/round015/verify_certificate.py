"""Minimal Fraction-only verification of the new complete upper certificate.

Python standard library only. Analytic lemmas and inherited g/UH proofs are
in primary_proof.md and their cited certificates; no decimal is an input.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json

P = Path(__file__).resolve().parent.parent


def read(name):
    return json.loads((P/name).read_bytes())


def main():
    a = read('round015/primary_upper_acceptance.json')
    for name, expected in a['source_sha256'].items():
        if name == 'results/bounds_ledger.json':
            previous = read('round014/manifest.json')['white_list']
            assert next(r['sha256'] for r in previous if r['path'] == name) == expected
            continue  # Bound to the preceding committed ledger snapshot.
        assert hashlib.sha256((P/name).read_bytes()).hexdigest() == expected, name
    s = read('upper_shorttime076/results.json')
    pi = Q(s['derived_constants']['pi_outward_384_dyadic']['lower'])
    c = Q(read('upper_facepeak134/condition_results.json')['Cbar_upper_256dyadic'])
    m = Q(read('upper_join132/portable_components.json')['objects']['analytic_stationary']['M_upper_portable_cap'])
    h = Q(read('round014/primary_integrals_acceptance.json')['complete_H512_integral_upper'])
    g = Q(read('round013/primary_globalg105_acceptance.json')['qualified_g'])
    e = Q(1, 512)
    I = Q(a['complete_I_upper'])
    quartic_increment = 625*e**6*(2*c+3*m*e**2)/(9072*pi**4)
    assert I > h > 0 and (I-h)**4 >= quartic_increment
    fourth = 3125*I**4/g**5
    lo, hi = map(Q, a['best_upper']['interval'])
    assert 0 < lo < hi < Q(799,25) < 32
    assert lo**4 <= fourth <= hi**4
    assert str(fourth) == a['best_upper']['fourth_power']
    assert Q(32)**4*g**5-3125*I**4 > 0
    print('PASS: complete exact upper certificate; no floating input or H rerun')


if __name__ == '__main__':
    main()
