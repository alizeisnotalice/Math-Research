#!/usr/bin/env python3
"""Independent interval review of the final algebraic comparison.

This intentionally does not import the production union routine.  It reads the
two stored coefficient lists and encloses u=2^(1/5), v=3^(1/5) with independent
integer bisection, then verifies candidate lower > seed upper.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib, json


def floor_root(n, k):
    lo, hi = 0, 1 << ((n.bit_length() + k - 1) // k)
    while lo + 1 < hi:
        m = (lo + hi) // 2
        if m ** k <= n: lo = m
        else: hi = m
    return lo


def enclosure(coefficients, denominator=10 ** 70):
    u0, v0 = floor_root(2 * denominator ** 5, 5), floor_root(3 * denominator ** 5, 5)
    ul, uh = Q(u0, denominator), Q(u0 + 1, denominator)
    vl, vh = Q(v0, denominator), Q(v0 + 1, denominator)
    low = high = Q(0)
    for key, value in coefficients.items():
        i, j = map(int, key.split(',')); c = Q(value)
        lo, hi = ul ** i * vl ** j, uh ** i * vh ** j
        if c >= 0:
            low += c * lo; high += c * hi
        else:
            low += c * hi; high += c * lo
    return low, high


def main():
    root = Path(__file__).resolve().parent
    certificate = root / 'rectangular_certificate.json'
    data = json.loads(certificate.read_text())
    seed = next(x for x in data['candidates'] if x['delta1'] == '0' and x['delta2'] == '0')
    candidate = next(x for x in data['candidates'] if x['delta1'] == '0' and x['delta2'] == '1/64')
    seed_lo, seed_hi = enclosure(seed['coefficients'])
    cand_lo, cand_hi = enclosure(candidate['coefficients'])
    assert cand_lo > seed_hi
    out = {
        'status': 'PASS independent algebraic comparison',
        'candidate': {'delta1': candidate['delta1'], 'delta2': candidate['delta2'], 'interval': [str(cand_lo), str(cand_hi)]},
        'seed': {'delta1': seed['delta1'], 'delta2': seed['delta2'], 'interval': [str(seed_lo), str(seed_hi)]},
        'strict_gap_lower': str(cand_lo - seed_hi),
        'certificate_sha256': hashlib.sha256(certificate.read_bytes()).hexdigest(),
        'reviewer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }
    (root / 'review_comparison.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
