#!/usr/bin/env python3
"""Twelve bounded audit predicates; no union, roots, writer imports or epsilon run."""
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
              'cube-dim-6-bounds/lower/round039/parity_first_variation')
MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
EXPECTED = '8a58120ad4882974daa2051fc3fda978e3ae3e19a8e6eb85669cf84ea1dcbe49'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_velocity(mass, ell1, upper):
    if ell1 % 2 == 0:
        return (Q(0), 0, 0)
    q = -(ell1//2)
    sigma = 1-2*(q % 2)
    c = Q(sigma*(1 if upper else -1), 12*ell1)
    n, a, b = mass, 0, 0
    while n % 2 == 0:
        n, a = n//2, a+1
    while n % 3 == 0:
        n, b = n//3, b+1
    if n != 1:
        raise AssertionError('unexpected root factor')
    return c*2**(a//6)*3**(b//6), a % 6, b % 6


def parse_velocity(source):
    if source == '0':
        return Q(0), 0, 0
    m = re.fullmatch(r'([^*]+)\*u\^(\d+)\*v\^(\d+)', source)
    if not m:
        raise AssertionError('unrecognized inherited velocity format')
    return Q(m[1]), int(m[2]), int(m[3])


def endpoint_at_one(source):
    return tuple(Q(source['c0'].get(f'{i},{j}', '0'))
                 +Q(source['c1'].get(f'{i},{j}', '0')) for i in range(6) for j in range(6))


def vector(entries):
    result = [Q(0)]*36
    for (i, j), value in entries.items():
        result[6*i+j] = Q(value)
    return tuple(result)


def main():
    cpu = time.process_time()
    reg = json.loads((HERE/'registration.json').read_text())
    began = datetime.datetime.fromisoformat(reg['registered_utc'])
    records = []
    def check(name, predicate, note):
        assert len(records) < 12, '12-check cap'
        assert predicate, name
        assert time.process_time()-cpu < 2, '2-second mathematical CPU target'
        assert (datetime.datetime.now(datetime.timezone.utc)-began).total_seconds() < 15*60, '15-minute wall cap'
        records.append({'name': name, 'status': 'pass', 'scope': note})
    manifest = json.loads((SOURCE/'manifest.json').read_text())
    packet_path = MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
    accepted_path = MAIN/'audit/lower008/parametric_union/candidate_results.json'
    packet = json.loads(packet_path.read_text())
    accepted = json.loads(accepted_path.read_text())
    inherited = json.loads((SOURCE/'facet_inventory_results.json').read_text())
    check('01 frozen provenance', sha(SOURCE/'manifest.json') == EXPECTED
          and all(sha(SOURCE/name) == digest for name, digest in manifest['file_sha256'].items())
          and sha(packet_path) == '99e762ad8fb63303f90b56c5d7f90fb855a3335903cae21650a23bfb07fea933'
          and sha(accepted_path) == '6849d3774ab3db1524969715970bbb3aecd07de4467faee5d77b69775cf5ed3d',
          'hash-only integrity; no mathematical sign test')
    labels = [tuple(x) for x in packet['box_lengths']]
    ties = accepted['boundary_ties']
    ranks = []
    for axis, groups in enumerate(packet['order_certificate']):
        rank = {key: i for i, group in enumerate(groups) for key in group}
        for tie in ties:
            if tie['axis'] == axis and tie['right_equal']:
                lo, hi = tie['adjacent_groups']
                rank = {key: (lo if value == hi else value) for key, value in rank.items()}
        ranks.append(rank)
    check('02 no baseline thin boxes', len(labels) == len(set(labels)) == len(packet['boxes']) == 127
          and all(ranks[h][box[h][0]] < ranks[h][box[h][1]] for box in packet['boxes'] for h in range(6))
          and all(1 <= ell <= 4 for label in labels for ell in label),
          'reuse inherited exact ranks including the sole boundary tie; all 127 boxes strictly positive')
    caps = [vector({(0, 0): Q(1, 2), (4, 0): Q(1, 4)}),
            vector({(2, 2): Q(1, 6), (2, 4): Q(1, 12)})]
    caps += [vector({(0, 0): Q(1, 4), (1, 0): Q(1, 4)})]*4
    check('03 fixed cell faces', all(endpoint_at_one(packet['endpoints'][key]) == (Q(0),)*36
          for groups in packet['order_certificate'] for key in groups[0])
          and all(endpoint_at_one(packet['endpoints'][key]) == caps[h]
                  for h, groups in enumerate(packet['order_certificate']) for key in groups[-1]),
          'formal coefficient identities for zero and three distinct cell-cap templates, no root/sign evaluation')
    groups = {}
    cell_count = 0
    for index, (ell, box) in enumerate(zip(labels, packet['boxes'])):
        for h in range(6):
            for side, key in enumerate(box[h]):
                rank = ranks[h][key]
                if rank in (0, len(packet['order_certificate'][h])-1):
                    cell_count += 1
                    continue
                groups.setdefault((h+1, rank, side), []).append(
                    (index, ell, prod(ell), ell[0], key, canonical_velocity(prod(ell), ell[0], bool(side))))
    split = {key: tuple(sorted(items)) for key, items in groups.items()
             if len({item[-1] for item in items}) > 1}
    old = {}
    for group in inherited['same_orientation_noncell_velocity_split_groups']:
        key = group['axis_one_based'], group['inherited_rank'], int(group['orientation'] == 'upper')
        old[key] = tuple(sorted((x['label_index'], tuple(x['lengths']), x['mass'], x['l1'],
                                x['endpoint_id'], parse_velocity(x['velocity'])) for x in group['occurrences']))
    check('04 complete velocity inventory', split == old and cell_count == 1028,
          'one integer/coefficient pass over all 1524 inherited face occurrences; zero new spatial facet partition')
    check('05 all potential positive split classes', len(split) == 10
          and Counter(key[0] for key in split) == Counter({h: 2 for h in range(2, 7)})
          and all(key[2] == 1 and {x[3] for x in values} == {1, 2}
                  and {x[2] for x in values} in ({108}, {162}) for key, values in split.items()),
          'ten noncell upper groups only; excludes any missed ell1=3 or lower-face velocity split')
    fixed = {item[0]: item for values in split.values() for item in values if item[3] == 2}
    family108 = {x[1] for x in fixed.values() if x[2] == 108}
    family162 = {x[1] for x in fixed.values() if x[2] == 162}
    dominators = set()
    maps = {}
    for index, item in fixed.items():
        b = item[1]
        k = b.index(1)
        c = list(b)
        c[k] = 2
        maps[index] = k, tuple(c)
        dominators.add(tuple(c))
    check('06 static label templates', len(fixed) == 25 and len(family108) == 20 and len(family162) == 5
          and all(sorted(b[1:]) == [1, 2, 3, 3, 3] for b in family108)
          and all(sorted(b[1:]) == [1, 3, 3, 3, 3] for b in family162)
          and len(dominators) == 15 and dominators <= set(labels),
          'finite membership/bookkeeping for all 25 B and 15 existing C; not 25 numerical inclusion tests')
    check('07 u upper bound', Q(6, 5)**6 > 2, 'u<6/5')
    check('08 x upper bound', Q(4, 3)**3 > 2, 'x=u^2<4/3')
    check('09 y upper bound', Q(3, 2)**3 > 3, 'y=v^2<3/2')
    check('10 minimum dominator radius', Q(6, 5)**2 < Q(3, 2), 'r216=sqrt(6)/2>6/5')
    a2_upper = Q(4, 3)*Q(3, 2)/3+Q(4, 3)*Q(3, 2)**2/6
    tail_upper = (1+Q(6, 5))/2
    check('11 all changed-coordinate spacings', max(a2_upper, tail_upper) < Q(6, 5),
          'a2<7/6 and tails<11/10, both below r216; active a1 is never changed')
    check('12 strict normal cover template', all(item[1][key[0]-1] == 3
          and maps[item[0]][0] != key[0]-1 and maps[item[0]][0] > 0
          and maps[item[0]][1][0] == item[1][0] == 2
          and prod(maps[item[0]][1]) == 2*item[2]
          for key, values in split.items() for item in values if item[3] == 2),
          'all fixed mixed faces have unchanged normal length3, changed transverse length1->2, mass doubled, radius strictly larger')
    assert len(records) == 12
    result = {'status': 'pass_with_proof_scope', 'pid': os.getpid(), 'new_exact_check_predicates': 12,
              'checks': records, 'source_manifest_sha256': EXPECTED,
              'independent_velocity_groups': 10, 'stationary_B_labels': 25,
              'existing_stationary_C_labels': 15, 'baseline_positive_boxes': 127,
              'new_field_sign_or_root_evaluations': 0, 'new_union_calls': 0,
              'new_epsilon_candidates': 0, 'baseline_replays': 0,
              'mathematical_CPU_seconds': time.process_time()-cpu,
              'input_sha256': {str(p): sha(p) for p in (SOURCE/'manifest.json', packet_path,
                  accepted_path, SOURCE/'facet_inventory_results.json', HERE/'registration.json')},
              'code_sha256': sha(Path(__file__)), 'new_Cn_bound': False,
              'thresholds_changed': False, 'process_status': 'finished; stop, upper-bound priority'}
    (HERE/'review_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'checks'}, indent=2))


if __name__ == '__main__':
    main()
