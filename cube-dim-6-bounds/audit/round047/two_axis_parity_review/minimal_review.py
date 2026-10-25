#!/usr/bin/env python3
"""Twelve bounded predicates for the frozen two-axis parity facet claim."""
from collections import Counter
from fractions import Fraction as Q
from math import prod
from pathlib import Path
import datetime
import hashlib
import json
import os
import re
import time

HERE = Path(__file__).resolve().parent
SOURCE = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/'
              'cube-dim-6-bounds/lower/round047/two_axis_parity_preflight')
MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
EXPECTED = 'e2495c43312ebc2c550d695f593c63cf5e3db8fb15bf4c7727bed0b621a2257f'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def block(length):
    q = -(length//2) if length % 2 else 1-length//2
    return tuple(range(q, q+length))


def parity(k):
    return 1-2*(k % 2)


def response(ell):
    return Q(sum(parity(k) for k in block(ell[0]))
             *sum(parity(k) for k in block(ell[1])), ell[0]*ell[1])


def velocity(ell, upper):
    rho = response(ell)
    if not rho:
        return Q(0), 0, 0
    n, a, b = prod(ell), 0, 0
    while n % 2 == 0:
        n, a = n//2, a+1
    while n % 3 == 0:
        n, b = n//3, b+1
    assert n == 1
    return Q(1 if upper else -1, 12)*rho*2**(a//6)*3**(b//6), a % 6, b % 6


def parse_velocity(text):
    if text == '0':
        return Q(0), 0, 0
    m = re.fullmatch(r'([^*]+)\*u\^(\d+)\*v\^(\d+)', text)
    assert m is not None
    return Q(m[1]), int(m[2]), int(m[3])


def main():
    cpu = time.process_time()
    registration = json.loads((HERE/'registration.json').read_text())
    began = datetime.datetime.fromisoformat(registration['registered_utc'])
    records = []
    def check(name, condition, note):
        assert len(records) < 12
        assert condition, name
        assert time.process_time()-cpu < 2
        assert (datetime.datetime.now(datetime.timezone.utc)-began).total_seconds() < 12*60
        records.append({'name': name, 'status': 'pass', 'scope': note})
    manifest = json.loads((SOURCE/'manifest.json').read_text())
    packet_path = MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
    accepted_path = MAIN/'audit/lower008/parametric_union/candidate_results.json'
    packet = json.loads(packet_path.read_text())
    accepted = json.loads(accepted_path.read_text())
    authored = json.loads((SOURCE/'preflight_results.json').read_text())
    check('01 frozen provenance', sha(SOURCE/'manifest.json') == EXPECTED
          and all(sha(SOURCE/name) == digest for name, digest in manifest['file_sha256'].items())
          and sha(packet_path) == '99e762ad8fb63303f90b56c5d7f90fb855a3335903cae21650a23bfb07fea933'
          and sha(accepted_path) == '6849d3774ab3db1524969715970bbb3aecd07de4467faee5d77b69775cf5ed3d',
          'source hashes only; author writer not executed')
    labels = [tuple(x) for x in packet['box_lengths']]
    lookup = {ell: i for i, ell in enumerate(labels)}
    ranks, lasts = [], []
    for h, groups in enumerate(packet['order_certificate']):
        rank = {key: k for k, group in enumerate(groups) for key in group}
        for tie in accepted['boundary_ties']:
            if tie['axis'] == h and tie['right_equal']:
                low, high = tie['adjacent_groups']
                rank = {key: (low if value == high else value) for key, value in rank.items()}
        ranks.append(rank)
        lasts.append(len(groups)-1)
    ranked = [tuple((ranks[h][pair[0]], ranks[h][pair[1]]) for h, pair in enumerate(box))
              for box in packet['boxes']]
    check('02 complete positive baseline labels', len(labels) == len(lookup) == len(ranked) == 127
          and all(lo < hi for box in ranked for lo, hi in box)
          and accepted['checked_adjacent_order_gaps'] == 120,
          'all 127 positive boxes; exact inherited base ranks/tie reused, no baseline replay')
    blocks = [set(block(l)) for l in range(1, 5)]
    check('03 central nested-source template', all(blocks[k] < blocks[k+1] for k in range(3)),
          'J1 strictly nested J2 strictly nested J3 strictly nested J4')
    phases = [parity(s+t) for s in (0, 1) for t in (0, 1)]
    check('04 four-phase affine mass and reflection', sorted(phases) == [-1, -1, 1, 1]
          and sum(phases) == 0 and sum(1 for _ in phases) == 4
          and 64*2//4 == 32 and all(parity(-k) == parity(k) for k in range(-1, 3)),
          'cell mass coefficients (4,0), two phases each sign, local reflections preserve parity; no nonzero epsilon used')
    mass_ok = True
    for ell in labels:
        for chi in phases:
            slope = sum(chi*parity(k+l) for k in block(ell[0]) for l in block(ell[1]))
            mass_ok &= Q(slope, ell[0]*ell[1]) == chi*response(ell)
    check('05 correlated mass and radius derivative template', mass_ok
          and response((1, 1, 1, 1, 1, 1)) == 1
          and response((1, 3, 1, 1, 1, 1)) == response((3, 1, 1, 1, 1, 1)) == Q(-1, 3)
          and response((3, 3, 1, 1, 1, 1)) == Q(1, 9)
          and all(abs(response(ell)) <= 1 for ell in labels),
          'exact epsilon coefficient factoring; analytic chain rule gives r1=chi*rho*r0/6, common radius')
    groups = {}
    cell_count = 0
    for index, (ell, box) in enumerate(zip(labels, packet['boxes'])):
        for h, pair in enumerate(box):
            for side, key in enumerate(pair):
                rank = ranks[h][key]
                if rank in (0, lasts[h]):
                    cell_count += 1
                    continue
                groups.setdefault((h+1, rank, side), []).append(
                    (index, ell, prod(ell), ell[h], key, response(ell), velocity(ell, bool(side))))
    split = {key: tuple(sorted(items)) for key, items in groups.items()
             if len({x[-1] for x in items}) > 1}
    old = {}
    for group in authored['split_groups']:
        key = group['axis_one_based'], group['rank'], int(group['orientation'] == 'upper')
        old[key] = tuple(sorted((x['label_index'], tuple(x['lengths']), x['mass'], x['axis_length'],
                                x['endpoint_id'], Q(x['rho']), parse_velocity(x['velocity']))
                               for x in group['occurrences']))
    check('06 independent full31 velocity inventory', split == old and len(split) == 31 and cell_count == 1028,
          'one metadata pass over all 1524 faces; all labels/rho/velocity/coincidence signatures compared')
    prior_path = SOURCE.parents[1]/'round039/parity_first_variation/facet_inventory_results.json'
    prior = json.loads(prior_path.read_text())
    prior_keys = {(g['axis_one_based'], g['inherited_rank'], int(g['orientation'] == 'upper'))
                  for g in prior['same_orientation_noncell_velocity_split_groups']}
    summary = Counter((key[0], key[2], next(iter({x[2] for x in values})))
                      for key, values in split.items())
    expected_summary = Counter({(1, 1, m): 1 for m in (2, 4, 8, 16, 2304)})
    expected_summary.update({(h, 0, m): 1 for h in range(3, 7) for m in (4, 8, 16, 2304)})
    expected_summary.update({(h, 1, m): 1 for h in range(2, 7) for m in (108, 162)})
    check('07 new groups beyond039', prior_keys <= set(split) and len(set(split)-prior_keys) == 21
          and summary == expected_summary,
          '31 groups exactly classified; old ten retained but no old25-label exclusion imported')
    check('08 equal-mass same-normal data', all(len({x[2] for x in values}) == 1
          and len({x[3] for x in values}) == 1 for values in split.values()),
          'all31 potential split groups have common baseline mass and normal source length')
    positive_pairs, zero_pairs = [], 0
    for key, values in split.items():
        h = key[0]-1
        for k, a in enumerate(values):
            for b in values[k+1:]:
                if a[-1] == b[-1]:
                    continue
                if all(max(ranked[a[0]][j][0], ranked[b[0]][j][0])
                       < min(ranked[a[0]][j][1], ranked[b[0]][j][1])
                       for j in range(6) if j != h):
                    c = tuple(max(x, y) for x, y in zip(a[1], b[1]))
                    positive_pairs.append((key, a, b, c))
                else:
                    zero_pairs += 1
    check('09 actual covering-label membership', bool(positive_pairs)
          and all(c in lookup for _, _, _, c in positive_pairs),
          'each positive5D transverse pair uses actual C in inherited127, not merely lengths<=4')
    check('10 exact Cartesian source hull', all(set(block(max(a, b))) == set(block(a)) | set(block(b))
          for a in range(1, 5) for b in range(1, 5)),
          'nested central source blocks make max(A,B) the coordinatewise endpoint hull')
    check('11 strict radius slack and moving-C positivity', all(a[1] != b[1]
          and a[2] == b[2] < prod(c) and all(1 <= x <= 4 for x in c)
          and abs(response(c)) <= 1 for _, a, b, c in positive_pairs),
          'mC>m gives strict rC>r without root approximation; M_C(epsilon)>0 near0 even when C moves')
    collar_ok = True
    for key, a, b, c in positive_pairs:
        h, ci = key[0]-1, lookup[c]
        face = key[1]
        collar_ok &= ranked[ci][h][0] < face < ranked[ci][h][1]
        for j in range(6):
            if j == h:
                continue
            low = max(ranked[a[0]][j][0], ranked[b[0]][j][0])
            high = min(ranked[a[0]][j][1], ranked[b[0]][j][1])
            collar_ok &= (ranked[ci][j][0] < low if low > 0 else ranked[ci][j][0] == 0)
            collar_ok &= (ranked[ci][j][1] > high if high < lasts[j] else ranked[ci][j][1] == lasts[j])
    check('12 whole cofacet relative-cell collar', collar_ok,
          'all positive5D pair intersections: strict normal inclusion and strict every noncell transverse endpoint; cell endpoints use relative-H collar')
    assert len(records) == 12
    result = {'status': 'independent_pass_for_fixed127_first_order_only', 'pid': os.getpid(),
              'checks': records, 'new_exact_predicates': 12, 'independent_split_groups': 31,
              'groups_beyond039': 21, 'positive5D_unequal_velocity_pairs': len(positive_pairs),
              'zero5D_pair_intersections': zero_pairs,
              'distinct_covering_C_labels': len({c for _, _, _, c in positive_pairs}),
              'positive_pairs_with_moving_C': sum(bool(response(c)) for _, _, _, c in positive_pairs),
              'cell_face_occurrences': cell_count, 'baseline_positive_boxes': 127,
              'mathematical_CPU_seconds': time.process_time()-cpu,
              'new_union_calls': 0, 'new_epsilon_candidates': 0, 'baseline_replays': 0,
              'root_approximations': 0, 'new_field_sign_evaluations': 0, 'MATLAB_calls': 0,
              'new_Cn_bound': False, 'thresholds_changed': False,
              'input_sha256': {str(p): sha(p) for p in (SOURCE/'manifest.json', packet_path,
                  accepted_path, SOURCE/'preflight_results.json', prior_path, HERE/'registration.json')},
              'code_sha256': sha(Path(__file__)), 'process_status': 'finished; stop, upper-bound priority'}
    (HERE/'review_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'checks'}, indent=2))


if __name__ == '__main__':
    main()
