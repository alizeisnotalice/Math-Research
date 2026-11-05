#!/usr/bin/env python3
"""Standalone exact verifier for the flat thin-hole band model.

It intentionally does not import the optimizer or its helpers.  All sliding
window extrema, constraints, and union-profile cells are regenerated here
from rational input.  Feasibility is checked with ``Fraction``.  The optional
``simple`` branch applies to any rational beta=1-gamma=p/q with common radius
and common height: it certifies J^q exactly.  A Decimal value of J outside
that branch is diagnostic, not a certificate.
"""
from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction as F
import hashlib, json
from pathlib import Path

getcontext().prec = 80
ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'band_rational_inputs.json'

def frac(text: str) -> F:
    return F(text)

def gamma_for_case(case: dict, default_gamma: F) -> F:
    gamma = frac(case.get('gamma', str(default_gamma)))
    if not (F(0) < gamma < F(1)):
        raise AssertionError('gamma must satisfy 0 < gamma < 1')
    return gamma

def checked_data(raw: dict) -> tuple[F, list[dict]]:
    """Parse and validate the complete rational input domain before use."""
    gamma = frac(raw['gamma'])
    if not (F(0) < gamma < F(1)):
        raise AssertionError('gamma must satisfy 0 < gamma < 1')
    cases = raw.get('cases')
    if not isinstance(cases, list) or not cases:
        raise AssertionError('cases must be a nonempty list')
    for case in cases:
        gamma_for_case(case, gamma)
        bands = [(frac(a), frac(b)) for a, b in case['bands']]
        radii = [frac(r) for r in case['radii']]
        heights = [frac(h) for h in case['heights']]
        eta = frac(case['eta'])
        if not (len(bands) == len(radii) == len(heights) and bands):
            raise AssertionError('invalid input sizes')
        if not (F(0) < eta < F(1)):
            raise AssertionError('eta must satisfy 0 < eta < 1')
        if any(a >= b for a, b in bands):
            raise AssertionError('each band must have positive length')
        # Endpoint contact is allowed; positive-length overlap is not.
        if any(bands[i][0] < bands[i - 1][0] or bands[i - 1][1] > bands[i][0]
               for i in range(1, len(bands))):
            raise AssertionError('bands must be sorted and have disjoint interiors')
        if any(r <= 0 for r in radii):
            raise AssertionError('all radii must be positive')
        if any(h < 0 for h in heights):
            raise AssertionError('all heights must be nonnegative')
    return gamma, cases

def overlap(left: F, right: F, a: F, b: F) -> F:
    return max(F(0), min(right, b) - max(left, a))

def nodes_for_band(index: int, bands: list[tuple[F,F]], radii: list[F]) -> list[F]:
    left, right = bands[index]; r = radii[index]
    endpoints = [z for ab in bands for z in ab]
    return sorted({left, right, *(z-r for z in endpoints if left <= z-r <= right),
                   *(z+r for z in endpoints if left <= z+r <= right)})

def decimal_fraction(x: F) -> Decimal:
    return Decimal(x.numerator) / Decimal(x.denominator)

def verify_case(case: dict, gamma: F) -> dict:
    bands = [(frac(a), frac(b)) for a,b in case['bands']]
    radii = [frac(r) for r in case['radii']]
    heights = [frac(h) for h in case['heights']]
    eta = frac(case['eta'])
    # checked_data has already established the input domain.  Keep this local
    # assertion so verify_case cannot silently be reused on malformed data.
    if not (F(0) < gamma < F(1) and F(0) < eta < F(1)
            and len(bands) == len(radii) == len(heights) and bands
            and all(a < b for a, b in bands) and all(r > 0 for r in radii)
            and all(h >= 0 for h in heights)
            and all(bands[i - 1][0] <= bands[i][0] and bands[i - 1][1] <= bands[i][0]
                    for i in range(1, len(bands)))):
        raise AssertionError('input is outside the verified domain')
    certificates = []
    for i, r in enumerate(radii):
        for y in nodes_for_band(i, bands, radii):
            integral = sum(overlap(y-r, y+r, a, b)*h for (a,b),h in zip(bands, heights))
            allowed = 2*(1-eta)*r*r
            if integral > allowed:
                raise AssertionError(f'window violation i={i}, y={y}: {integral}>{allowed}')
            certificates.append({'band': i, 'y': str(y), 'integral': str(integral),
                                 'allowed': str(allowed)})
    points = sorted({z + sign*r for (a,b),r in zip(bands,radii) for z in (a,b) for sign in (-1,1)})
    cells = []
    for left, right in zip(points[:-1], points[1:]):
        mid=(left+right)/2
        height=max((r for (a,b),r in zip(bands,radii) if a-r <= mid <= b+r), default=F(0))
        if height:
            cells.append({'left':str(left),'right':str(right),'radius':str(height)})
    beta = F(1)-gamma
    numerator = sum(decimal_fraction(b-a) * (decimal_fraction(h) ** (Decimal(beta.numerator)/Decimal(beta.denominator)))
                    for (a,b),h in zip(bands,heights))
    denominator = sum(decimal_fraction(F(c['right'])-F(c['left'])) *
                      (decimal_fraction(F(c['radius'])) ** (Decimal(beta.numerator)/Decimal(beta.denominator)))
                      for c in cells)
    if not cells or denominator <= 0:
        raise AssertionError('union-profile denominator must be positive')
    answer = {'name':case['name'], 'gamma':str(gamma), 'eta':str(eta),
              'feasibility':'PASS',
              'J_decimal_diagnostic':str(numerator/denominator),
              'node_count':sum(len(nodes_for_band(i,bands,radii)) for i in range(len(bands))),
              'profile_cell_count':len(cells)}
    # Here, and only here, the denominator is |union_i(B_i+[-r,r])| r^beta.
    # For beta=p/q in lowest terms, J^q is exactly rational.
    if len(set(radii)) == 1 and len(set(heights)) == 1:
        total_length = sum(b-a for a,b in bands)
        union_length = sum(F(c['right'])-F(c['left']) for c in cells)
        q = beta.denominator
        j_to_q = (total_length / union_length)**q * (heights[0] / radii[0])**beta.numerator
        expected = case.get('expect_simple_J_to_q')
        if expected not in (None, 'lt_one', 'gt_one'):
            raise AssertionError('expect_simple_J_to_q must be lt_one or gt_one')
        relation_holds = (expected is None or
                          (expected == 'lt_one' and j_to_q < 1) or
                          (expected == 'gt_one' and j_to_q > 1))
        answer.update({'simple_beta':str(beta), 'simple_power_q':q,
                       'simple_exact_J_to_q':str(j_to_q),
                       'simple_expected_relation':expected,
                       'simple_exact_relation_holds':relation_holds})
        if not relation_holds:
            raise AssertionError('requested exact simple J^q comparison failed')
    return answer

if __name__ == '__main__':
    raw=json.loads(INPUT.read_text(encoding='utf-8')); gamma,cases=checked_data(raw)
    signature=hashlib.sha256(Path(__file__).read_bytes()+INPUT.read_bytes()).hexdigest()
    output={'schema':'flat-band-fraction-verifier-v2','source_signature':signature,
            'default_gamma':str(gamma),
            'cases':[verify_case(c,gamma_for_case(c,gamma)) for c in cases]}
    for case in output['cases']:
        print(f"{case['name']}: J diagnostic={case['J_decimal_diagnostic']} nodes={case['node_count']}")
        if 'simple_exact_J_to_q' in case:
            print(f"  exact J^{case['simple_power_q']}={case['simple_exact_J_to_q']} "
                  f"expected {case['simple_expected_relation']}: {case['simple_exact_relation_holds']}")
    print('BAND FRACTION VERIFIER PASS')
