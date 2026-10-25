#!/usr/bin/env python3
"""Bounded floating discovery screen for staggered unit-mass columns.

Only raw witness boxes are unioned. This file produces no C6 certificate.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import ceil, floor, prod
from pathlib import Path
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve()
MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
TOL = 2e-12
RAW_CAP = 5000
CASE_CAP = 12
WALL_SECONDS = 1200


def spacings():
    u, v = 2 ** (1 / 6), 3 ** (1 / 6)
    return (1 + u**4 / 2, (u**3 + u**4 * v**3) / 4,
            *((1 + u) / 2 for _ in range(4)))


def starts(ell, radius_bound, a):
    # The actual conditions are (s+ell-1)a-r<a and sa+r>0.
    return range(floor(-radius_bound / a) - 1,
                 ceil(1 + radius_bound / a - (ell - 1)) + 2)


def check_time(started):
    if time.monotonic() - started > WALL_SECONDS:
        raise TimeoutError('20-minute global wall budget reached')


def raw_boxes(pair, delta, started):
    order = tuple(pair) + tuple(j for j in range(6) if j not in pair)
    a0 = spacings()
    a = tuple(a0[j] for j in order)
    raw = []
    for ell in product(range(1, 5), repeat=6):
        check_time(started)
        m0 = prod(ell)
        rmax = m0 ** (1 / 6) / 2
        # The enlarged symmetric j start bound is a proven superset for its
        # centered full-cell domain and all offsets of size <=delta*a_j.
        jbound = rmax + abs(delta)*a[1] + a[1]/2
        for si, sj in product(starts(ell[0], rmax, a[0]), starts(ell[1], jbound, a[1])):
            mass = m0
            radius = mass ** (1 / 6) / 2
            endpoints = []
            source_starts = [si, sj]
            for j in range(2, 6):
                source_starts.append(-(ell[j] // 2) if ell[j] % 2 else 1 - ell[j] // 2)
            valid = True
            offsets = [delta*a[1]*(1 if (si+k)%2 == 0 else -1) for k in range(ell[0])]
            for j in range(6):
                source_lo = source_starts[j] * a[j]
                source_hi = (source_starts[j] + ell[j] - 1) * a[j]
                if j == 1:
                    source_lo += min(offsets)
                    source_hi += max(offsets)
                domain_lo = -a[j]/2 if j == 1 else 0.
                domain_hi = a[j] if j == 0 else a[j] / 2
                lo, hi = max(domain_lo, source_hi - radius), min(domain_hi, source_lo + radius)
                if hi <= lo + TOL:
                    valid = False
                    break
                if any(abs(x - atom) > radius + 4e-12 for x in (lo, hi)
                       for atom in (source_lo, source_hi)):
                    raise AssertionError('cube radius endpoint check failed')
                endpoints += [lo, hi]
            if valid:
                raw.append({'lengths_reordered': ell, 'source_starts_reordered': source_starts,
                            'mass_float': mass, 'radius_float': radius, 'box': tuple(endpoints)})
                if len(raw) > RAW_CAP:
                    raise RuntimeError('5000 raw-box cap exceeded')
    return raw, order


def maximal_boxes(raw):
    unique = {tuple(item['box']): item for item in raw}
    kept = []
    for box, item in unique.items():
        # Strict float containment only: never discard on a tolerance extension.
        if any(other != box and all(other[2*j] <= box[2*j] and box[2*j+1] <= other[2*j+1]
                                    for j in range(6)) for other in unique):
            continue
        kept.append(item)
    return kept


def union_volume(boxes, started):
    @lru_cache(None)
    def rec(X):
        check_time(started)
        if not X:
            return 0.
        if len(X[0]) == 2:
            intervals = sorted(X)
            lo, hi = intervals[0]
            value = 0.
            for c, d in intervals[1:]:
                if c > hi:
                    value += hi - lo
                    lo, hi = c, d
                elif d > hi:
                    hi = d
            return value + hi - lo
        value = 0.
        cuts = sorted({x for b in X for x in b[:2]})
        for lo, hi in zip(cuts, cuts[1:]):
            if hi > lo:
                child = tuple(sorted({b[2:] for b in X if b[0] <= lo and hi <= b[1]}))
                value += (hi-lo) * rec(child)
        return value
    return rec(tuple(sorted(tuple(item['box']) for item in boxes)))


def evaluate(pair, delta_exact, started):
    case_start = time.monotonic()
    raw, order = raw_boxes(pair, float(delta_exact), started)
    kept = maximal_boxes(raw)
    volume = union_volume(kept, started)
    return {'pair_one_based': [j+1 for j in pair], 'delta': str(delta_exact),
            'coordinate_order_zero_based': order, 'raw_boxes': len(raw), 'maximal_boxes': len(kept),
            'half_cell_union_float': volume, 'ratio_float': 16 * volume,
            'elapsed_seconds': time.monotonic() - case_start}, raw, kept


def main():
    started = time.monotonic()
    ledger_bytes = (MAIN / 'results/bounds_ledger.json').read_bytes()
    frozen_bytes = (MAIN / 'coordination/phase_20261025/frozen_thresholds.json').read_bytes()
    current_upper = Fraction(json.loads(ledger_bytes)['lower']['exact_interval'][1])
    frozen_upper = Fraction(json.loads(frozen_bytes)['lower']['exact_interval'][1])
    out = {'status': 'running', 'pid': os.getpid(), 'scope': 'floating discovery only; no C6 bound',
           'normalization': 'period mass 2; full reflected union=32|U_H|; ratio=16|U_H|',
           'source_rule': 'x_i=k_i*a_i, x_j=k_j*a_j+delta*a_j*(-1)^k_i; unit masses',
           'input': {'ledger_sha256': hashlib.sha256(ledger_bytes).hexdigest(),
                     'frozen_thresholds_sha256': hashlib.sha256(frozen_bytes).hexdigest(),
                     'current_complete_upper': str(current_upper), 'frozen_complete_upper': str(frozen_upper),
                     'lengths_each_axis': [1, 2, 3, 4], 'deltas': ['1/32','1/16','1/8'],
                     'ordered_pair_representatives_one_based': [[1,2],[2,1],[1,3],[3,1]]},
           'resource_limits': {'raw_boxes_each_candidate': RAW_CAP, 'candidate_count': CASE_CAP,
                               'global_wall_seconds': WALL_SECONDS}, 'results': []}
    outfile = HERE.with_name('staggered_results.json')
    def save():
        out['elapsed_seconds'] = time.monotonic() - started
        out['code_sha256'] = hashlib.sha256(HERE.read_bytes()).hexdigest()
        outfile.write_text(json.dumps(out, indent=2) + '\n')
    save()
    try:
        cal, _, _ = evaluate((0,1), Fraction(0), started)
        cal['error_vs_current_upper_float'] = cal['ratio_float'] - float(current_upper)
        out['delta_zero_calibration'] = cal
        save()
        if abs(cal['error_vs_current_upper_float']) > 2e-9:
            out['status'] = 'calibration_mismatch_stop'
        else:
            pilot = ((0,1), Fraction(1,16))
            cases = [pilot] + [(pair, d) for pair in ((0,1),(1,0),(0,2),(2,0))
                     for d in (Fraction(1,32),Fraction(1,16),Fraction(1,8)) if (pair,d) != pilot]
            for pair, delta in cases:
                result, raw, kept = evaluate(pair, delta, started)
                result['margin_over_current_upper_float'] = result['ratio_float'] - float(current_upper)
                result['promotion_trigger'] = result['margin_over_current_upper_float'] > 1e-10
                out['results'].append(result)
                save()
                if result['promotion_trigger']:
                    HERE.with_name('promoted_candidate_boxes.json').write_text(json.dumps(
                        {'candidate': result, 'raw': raw, 'maximal': kept}, indent=2) + '\n')
                    out['status'] = 'promotion_stop_requires_exact_certificate'
                    break
            else:
                out['status'] = 'completed_numeric_screen_no_promotion'
    except (RuntimeError, TimeoutError) as exc:
        out['status'] = 'budget_stop'
        out['stop_reason'] = str(exc)
    out['best'] = max(out['results'], key=lambda r:r['ratio_float'], default=None)
    out['actual_process_status'] = 'finished, no child processes or MATLAB calls'
    save()
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
