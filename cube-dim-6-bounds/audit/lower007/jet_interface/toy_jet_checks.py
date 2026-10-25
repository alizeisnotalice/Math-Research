#!/usr/bin/env python3
"""Eight rational analytic toys; no candidate-direction or old-volume scan."""
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import time


@dataclass(frozen=True, order=True)
class Jet:
    value: Q
    slope: Q = Q(0)

    def __add__(self, other):
        return Jet(self.value+other.value, self.slope+other.slope)

    def __sub__(self, other):
        return Jet(self.value-other.value, self.slope-other.slope)

    def __mul__(self, other):
        return Jet(self.value*other.value,
                   self.value*other.slope+self.slope*other.value)


ZERO, ONE = Jet(Q(0)), Jet(Q(1))


def polynomial(*coefficients):
    return tuple(Q(x) for x in coefficients)


def evaluate(poly, t):
    return sum((c*t**i for i, c in enumerate(poly)), Q(0))


def jet(poly):
    return Jet(poly[0], poly[1] if len(poly) > 1 else Q(0))


def union_volume(boxes, zero, one):
    """Same geometric definition for ordered affine germs or rational endpoints.

    Equal endpoint jets are merged for the *coefficient calculation only*.
    Positive first-order widths at value=0 are retained.
    """
    kept = tuple(sorted(set(b for b in boxes if all(b[k] < b[k+1]
                                                   for k in range(0, len(b), 2)))))
    @lru_cache(None)
    def rec(active):
        if not active:
            return zero
        if not active[0]:
            return one
        cuts = sorted({x for b in active for x in b[:2]})
        value = zero
        for lo, hi in zip(cuts, cuts[1:]):
            child = tuple(sorted({b[2:] for b in active if b[0] <= lo and hi <= b[1]}))
            value = value+(hi-lo)*rec(child)
        return value
    return rec(kept)


def main():
    started = time.monotonic()
    p = polynomial
    # Each box is a flattened tuple of endpoint polynomials (low,high,...).
    toys = [
        {'name': '1D zero-width birth', 'dimension': 1,
         'boxes': [(p(0), p(1)), (p(2), p(2, 1))], 'volume': p(1, 1)},
        {'name': '1D old touching intervals split a gap', 'dimension': 1,
         'boxes': [(p(0), p(1, -1)), (p(1, 1), p(2))], 'volume': p(2, -2)},
        {'name': '1D baseline-coincident boxes change containment', 'dimension': 1,
         'boxes': [(p(0), p(1, 1)), (p(0), p(1, -1))], 'volume': p(1, 1)},
        {'name': '2D one thin coordinate is first order', 'dimension': 2,
         'boxes': [(p(0), p(1), p(0), p(1)), (p(2), p(2, 3), p(0), p(2))],
         'volume': p(1, 6)},
        {'name': '2D two newborn widths are second order', 'dimension': 2,
         'boxes': [(p(0), p(1), p(0), p(1)), (p(2), p(2, 1), p(0), p(0, 2))],
         'volume': p(1, 0, 2)},
        {'name': '3D two newborn widths remain second order', 'dimension': 3,
         'boxes': [(p(0), p(1), p(0), p(1), p(0), p(1)),
                   (p(2), p(2, 1), p(0), p(0, 1), p(0), p(3))],
         'volume': p(1, 0, 3)},
        {'name': '1D equal jets hide a second-order gap', 'dimension': 1,
         'boxes': [(p(0), p(1)), (p(1, 0, 1), p(2))], 'volume': p(2, 0, -1)},
        {'name': '1D equal jets hide a second-order overlap', 'dimension': 1,
         'boxes': [(p(0), p(1)), (p(1, 0, -1), p(2))], 'volume': p(2)},
    ]
    assert len(toys) <= 8
    results = []
    for toy in toys:
        assert time.monotonic()-started < 1200
        boxes = toy['boxes']
        got = union_volume(tuple(tuple(jet(e) for e in b) for b in boxes), ZERO, ONE)
        expected = jet(toy['volume'])
        assert got == expected
        samples = []
        for t in (Q(1, 10), Q(1, 100)):
            exact = union_volume(tuple(tuple(evaluate(e, t) for e in b) for b in boxes), Q(0), Q(1))
            assert exact == evaluate(toy['volume'], t)
            samples.append({'t': str(t), 'exact_volume': str(exact),
                            'remainder_after_first_jet': str(exact-got.value-t*got.slope)})
        results.append({'name': toy['name'], 'dimension': toy['dimension'],
                        'certified_volume_jet': [str(got.value), str(got.slope)],
                        'analytic_volume_coefficients': [str(x) for x in toy['volume']],
                        'exact_samples': samples, 'status': 'pass'})
    # Same first jet, distinct exact t>0 volumes: this is the intended tie policy.
    assert results[6]['certified_volume_jet'] == results[7]['certified_volume_jet']
    assert results[6]['exact_samples'][0]['exact_volume'] != results[7]['exact_samples'][0]['exact_volume']
    out = {'status': 'pass', 'new_C6_bound': False, 'toy_instances': len(toys),
           'scope': 'rational 1D/2D/3D analytic toy validation only; no direction scan',
           'limits': {'toy_instances': 8, 'wall_seconds': 1200},
           'actual_candidates_evaluated': 0, 'results': results,
           'elapsed_seconds': time.monotonic()-started,
           'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'no_MATLAB': True}
    Path(__file__).with_name('toy_results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
