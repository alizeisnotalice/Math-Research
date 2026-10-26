#!/usr/bin/env python3
"""Prepared finite 11330-child consumer. DO NOT RUN before primary review.

Only exact rational arithmetic on saved parent H128 U/Delta/gradient intervals.
No H, physical gradient, numerical kernel, log, exp, Arb, or new split calls.
The source qualification/prepared-input SHA must be supplied explicitly.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import os
import time

WORK = Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_curvature128')
G105 = F(21, 200000)
MODEL = 'n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
KERNEL = '128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0'
ENC = json.JSONEncoder(sort_keys=True, separators=(',', ':'))

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while b := f.read(8*1024*1024):
            h.update(b)
    return h.hexdigest()

def digest(obj):
    h = hashlib.sha256()
    for s in ENC.iterencode(obj):
        h.update(s.encode())
    return h.hexdigest()

def write(path, obj):
    with Path(path).open('x') as f:
        for s in ENC.iterencode(obj):
            f.write(s)
        f.write('\n')

def box(node, roots):
    lo = [F(x, 32) for x in roots[node[0]]]
    hi = [x+F(1, 32) for x in lo]
    for i, bit in enumerate(node[1]):
        j = i % 4
        mid = (lo[j]+hi[j])/2
        if bit == '0':
            hi[j] = mid
        else:
            assert bit == '1'
            lo[j] = mid
    return lo, hi

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reviewed-preparation-sha', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    began = time.monotonic()
    qualification_path = WORK/'source_qualification.json'
    assert sha(qualification_path) == args.reviewed_preparation_sha
    qualified = json.loads(qualification_path.read_bytes())
    assert qualified['status'] == 'source_and_ancestry_qualified_no_candidate_test'
    assert qualified['finite_child_candidates'] == 11330
    assert qualified['uniform_anchor_H64_lower'] == str(G105)
    candidate_path = Path(qualified['prepared_candidates_path'])
    assert sha(candidate_path) == qualified['prepared_candidates_sha256']
    state_lock = qualified['source_locks']['state123']
    assert sha(state_lock['path']) == state_lock['sha256']
    state = json.loads(Path(state_lock['path']).read_bytes())
    assert len(state['queue']) == 11330 and len(state['new_splits']) == 5665
    assert not state['faults'] and not state['inflight'] and not state['obstructions_H128']
    args.output.mkdir(exist_ok=False)
    registration = {'task': 'upper_curvature128_finite11330_zero_H', 'pid': os.getpid(),
        'registered_utc': datetime.now(timezone.utc).isoformat(), 'code_sha256': sha(__file__),
        'source_qualification_sha256': args.reviewed_preparation_sha,
        'prepared_candidates_sha256': qualified['prepared_candidates_sha256'],
        'state123_sha256': state_lock['sha256'], 'finite_candidates': 11330,
        'maximum_closed_corners': 181280, 'workers': 1, 'pool': None, 'wall_deadline': None,
        'target': '11/100000', 'epsilon': '1/128', 'anchor_H64_lower': str(G105),
        'new_H_gradient_kernel_log_exp_calls': 0, 'new_splits': 0,
        'rule': 'same saved parent H128 Delta and directed physical loggradient; replace Q1024 by1024-768*g105/U128;one anchor must pass all16corners'}
    write(args.output/'registration.json', registration)
    roots = [(a,b,c,d) for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
    queue = set(map(tuple, state['queue']))
    cache = state['new_cache_by_epsilon']['1/128']
    seen = set()
    passed = 0
    with candidate_path.open() as source, (args.output/'candidate_results.jsonl').open('x') as output:
        for ordinal, line in enumerate(source):
            item = json.loads(line)
            node = tuple(item['child_node'])
            assert item['candidate_ordinal'] == ordinal and node in queue and node not in seen
            seen.add(node)
            split = state['new_splits'][item['parent_split_index']]
            assert split['node'] == item['direct_parent_node']
            assert split['children'][item['side']] == item['child_node']
            cert = split['failed_certificate']
            rec = cache[item['point_key']]
            assert digest(cert) == item['saved_certificate_sha256'] and digest(rec) == item['record_sha256']
            assert cert['point_key'] == rec['pair_cache_key'][1] == item['point_key']
            assert rec['model'] == MODEL and rec['kernel_sha256'] == KERNEL
            assert rec['epsilon'] == cert['leaf_epsilon'] == '1/128'
            assert rec['qualification_target'] == cert['target'] == '11/100000'
            assert item['H64_lower'] == str(G105)
            lo, hi = box(node, roots)
            point = list(map(F, rec['point']))
            assert point[0] == 1
            hl, hu = map(F, rec['values'][0])
            assert 0 < hl <= hu and G105 <= hu
            gradients = [tuple(map(F, pair)) for pair in cert['loggradient_directed']]
            assert len(gradients) == 4
            for (gl, gu), pair in zip(gradients, rec['values'][1:]):
                a, b = map(F, pair)
                assert gl <= min(a/hl, a/hu, b/hl, b/hu) <= max(a/hl, a/hu, b/hl, b/hu) <= gu
            delta = F(cert['Delta_lower'])
            q = F(1024)-F(768)*G105/hu
            assert 256 <= q < 1024
            margins = []
            for bits in product((0,1), repeat=4):
                u = [(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)]
                margin = delta+sum(min(gl*v,gu*v) for v,(gl,gu) in zip(u,gradients))-q*sum(v*v for v in u)
                margins.append((margin,bits))
            margin, bits = min(margins)
            passed += margin > 0
            result = {'candidate_ordinal': ordinal, 'child_node': item['child_node'],
                'direct_parent_node': item['direct_parent_node'], 'point_key': item['point_key'],
                'record_sha256': item['record_sha256'], 'saved_certificate_sha256': item['saved_certificate_sha256'],
                'epsilon': '1/128', 'target': '11/100000', 'H64_lower': str(G105),
                'H64_lower_source_node': item['H64_lower_source_node'], 'H128_upper': str(hu),
                'penalty_coefficient': str(q), 'Delta_loggradient_recomputed': False,
                'minimum_log_margin': str(margin), 'worst_corner_bits': list(bits),
                'checked_points': 16, 'passed': margin > 0,
                'mechanism': 'fixed-t075_then_anchor-Jensen_time-mass',
                'lower_bound_object': 'G', 'component_function': 'H_1/128'}
            output.write(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
            if (ordinal+1) % 256 == 0:
                output.flush()
                os.fsync(output.fileno())
    assert seen == queue and len(seen) == 11330
    summary = {'status': 'finite11330_zero_H_candidate_test_complete', 'pid': os.getpid(),
        'tested': len(seen), 'passed': passed, 'pending': len(seen)-passed, 'checked_corners': 16*len(seen),
        'new_H_gradient_kernel_log_exp_calls': 0, 'new_splits': 0,
        'source_qualification_sha256': args.reviewed_preparation_sha,
        'registration_sha256': sha(args.output/'registration.json'),
        'results_sha256': sha(args.output/'candidate_results.jsonl'),
        'source_numeric_independent_and_primary_qualification_still_required': True,
        'global_g110_or_C5_claim': False, 'wall_seconds': time.monotonic()-began}
    write(args.output/'summary.json', summary)
    print(json.dumps(summary))

if __name__ == '__main__':
    main()
