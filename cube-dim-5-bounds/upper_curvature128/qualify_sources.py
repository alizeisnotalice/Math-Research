#!/usr/bin/env python3
"""Read-only source/ancestry qualification; does not evaluate real q or margins.

Writes a finite candidate pointer packet outside Git. Historical accepted H64
closed-leaf conclusions are inherited through the locked 096/098 acceptance,
not recomputed. No kernel, log, exp, Arb, gradient or H evaluation is imported.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import gzip
import hashlib
import json
import os
import time

MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
SCRATCH = Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch')
OUT = SCRATCH / 'upper_curvature128'
MODEL = 'n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
G105 = '21/200000'
TARGET = '11/100000'
KERNEL = '128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0'
LOCKS = {
    'packet_manifest': (MAIN/'upper_portable096/run098_v3/manifest.json', 'c53cef7f2626e68042b73d4bb5e80387841eb810a8fbcc1c511f9ce1d8b4780b'),
    'packet_checker': (MAIN/'upper_portable096/run098_v3/checker_summary.json', '5ad237758ac62c860e43e1977056b24e0ed9b366854dbd275dfb5aa3139c701d'),
    'source075': (MAIN/'audit_frame075/report.md', '3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225'),
    'source029': (MAIN/'upper_alpha1_029/report.md', '298f0de7277ea6a7a516b85344c5189562bf7026f3f8de43d4d3f304ff6fed81'),
    'source099': (MAIN/'audit_mixed099/report.md', '96e874223e5fbc78f9688f7eec574cdcdde184344b1fa809114bc87553d48e85'),
    'pending105': (SCRATCH/'upper_g110105/pending_gap_order.jsonl', '2e4639a521a88aec1b49074f1e38bf5eecb9b853d2321a4b8a2685420b8285d5'),
    'mapper105': (SCRATCH/'upper_g110105/source_leaf_mapper.jsonl', 'c1f95e6cdb69f096fcbc03c5dcceca15795f5c96a225c8b34d66498cb23b63c0'),
    'tasks123': (SCRATCH/'upper_continue123/planned_tasks.json', '0bf33d3fd8a4ab228252062baa2d3469cb80abb3ddb46c591392a4c63d82db72'),
    'state123': (SCRATCH/'upper_continue123/batch0002/state.json', '6c1d4c6e075053c427d797798b5ccfc20e62e46a46151511e01a330b38213301'),
    'points123_1': (SCRATCH/'upper_continue123/batch0001/point_manifest.json', '40342fbe1d0bf596dc3a31aac741237ad926b17539c8f96a5988f1e3e331b656'),
    'points123_2': (SCRATCH/'upper_continue123/batch0002/point_manifest.json', 'faf4de812b216875bebbddc802843ca8f28a9523b0f5e93a108a563374ff04e2'),
    'state107': (SCRATCH/'upper_g110107/batch0001/state.json', 'bc38ff5d5ce43e0683c14de7205e40d786f6076f0e4a3550a6aac2de50caf3e1'),
    'points107': (SCRATCH/'upper_g110107/batch0001/point_manifest.json', '71d3d488bf5ffd4e37843df98d00adfbd4d7069fb76d4ba85063f20284089ad3'),
}
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
    rid, path = node
    assert type(rid) is int and 0 <= rid < len(roots)
    lo = [F(x, 32) for x in roots[rid]]
    hi = [x+F(1, 32) for x in lo]
    for i, bit in enumerate(path):
        axis = i % 4
        mid = (lo[axis]+hi[axis])/2
        if bit == '0':
            hi[axis] = mid
        else:
            assert bit == '1'
            lo[axis] = mid
    return lo, hi

def contains(bounds, point):
    lo, hi = bounds
    assert F(point[0]) == 1
    return all(a <= F(p) <= b for a, p, b in zip(lo, point[1:], hi))

def main():
    began = time.monotonic()
    OUT.mkdir(exist_ok=True)
    assert not (OUT/'source_qualification.json').exists(), 'Preserve frozen output'
    locks = {}
    for name, (path, expected) in LOCKS.items():
        actual = sha(path)
        assert actual == expected, (name, actual, expected)
        locks[name] = {'path': str(path), 'sha256': actual, 'bytes': path.stat().st_size}
    acceptance_path = MAIN/'round013/primary_globalg105_acceptance.json'
    acceptance = json.loads(acceptance_path.read_bytes())
    assert acceptance['status'] == 'primary_accepts_complete_same_model_global_g105'
    assert acceptance['qualified_g'] == G105
    assert acceptance['full098_state_sha256'] == '73bc873f23f7b68f6f70132e590d28ebbe4cd619012ac9c81980f9fa500314af'
    for name in ('packet_manifest', 'packet_checker', 'source075', 'source029', 'source099'):
        relative = str(LOCKS[name][0].relative_to(MAIN))
        assert acceptance['input_locks'][relative] == LOCKS[name][1]
    locks['g105_primary_acceptance'] = {'path': str(acceptance_path), 'sha256': sha(acceptance_path)}
    reg = {'task': 'upper_curvature128_source_qualification', 'registered_utc': datetime.now(timezone.utc).isoformat(),
           'pid': os.getpid(), 'code_sha256': sha(__file__), 'source_locks': locks,
           'finite_objects': {'pending105': 4973, 'points107': 4096, 'points123': 6237, 'children123': 11330},
           'new_H_gradient_kernel_log_exp_calls': 0, 'new_real_q_or_corner_evaluations': 0,
           'scope': 'check source identity, epsilon, historical accepted H64 leaf membership, exact ancestry, and same-epsilon saved jet/certificate binding only'}
    write(OUT/'source_qualification_registration.json', reg)
    roots = [(a,b,c,d) for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
    assert len(roots) == 52360
    pending = [json.loads(s) for s in LOCKS['pending105'][0].open()]
    assert len(pending) == 4973
    old = {tuple(r['node']): r for r in pending}
    assert len(old) == 4973
    byroot = defaultdict(dict)
    for key, row in old.items():
        assert row['leaf_epsilon'] == '1/64'
        assert row['certificate']['target'] == TARGET
        assert row['certificate']['leaf_epsilon'] == '1/64'
        lo, hi = box(row['node'], roots)
        assert list(map(str, lo)) == row['box_lower'] and list(map(str, hi)) == row['box_upper']
        byroot[key[0]][key[1]] = row
    source_mapper_matches = 0
    for line in LOCKS['mapper105'][0].open():
        row = json.loads(line)
        key = tuple(row['node'])
        if key in old:
            pending_row = old[key]
            for field in ('source_leaf_index', 'source_leaf_sha256', 'source_certificate_sha256', 'source_record_sha256', 'box_lower', 'box_upper', 'leaf_epsilon', 'pair_cache_key'):
                assert row[field] == pending_row[field]
            source_mapper_matches += 1
    assert source_mapper_matches == 4973
    manifest = json.loads(LOCKS['packet_manifest'][0].read_bytes())
    checker = json.loads(LOCKS['packet_checker'][0].read_bytes())
    assert manifest['g'] == checker['g'] == G105
    assert checker['status'] == 'PASS_portable_compact_closed_root_packet'
    assert checker['packet_manifest_sha256'] == LOCKS['packet_manifest'][1]
    assert checker['source_locks']['state'] == acceptance['full098_state_sha256']
    portable_matches = 0
    for meta in manifest['chunks']:
        path = LOCKS['packet_manifest'][0].parent/meta['file']
        assert sha(path) == meta['sha256']
        packet = json.loads(gzip.decompress(path.read_bytes()))
        for leaf in packet['leaves']:
            key = tuple(leaf['node'])
            if key in old:
                assert leaf['epsilon'] == leaf['cache_epsilon'] == '1/64'
                assert leaf['source_penalty_coefficient'] == '256'
                assert leaf['curvature_mode'] == 'full_frame_16'
                assert leaf['source_mechanism'] in ('full_frame_16', 'full_frame_16_new_target_g105')
                assert leaf['source_quadratic'] == 'Q_full=I4/16'
                assert leaf['point_tail'] == old[key]['pair_cache_key'][1].split(',')[1:]
                portable_matches += 1
    assert portable_matches == 4973

    def ancestor(node):
        rid, path = node
        matches = [row for prefix, row in byroot[rid].items() if path.startswith(prefix)]
        assert len(matches) == 1, node
        return matches[0]

    tasks = json.loads(LOCKS['tasks123'][0].read_bytes())
    assert len(tasks) == 6237
    taskmap = {tuple(t['node']): t for t in tasks}
    assert len(taskmap) == 6237
    task_ancestors = {}
    for task in tasks:
        assert task['epsilon'] == '1/128' and task['target'] == TARGET
        inherited = ancestor(task['node'])
        assert contains(box(task['node'], roots), task['point'])
        assert contains(box(inherited['node'], roots), task['point'])
        task_ancestors[tuple(task['node'])] = inherited
    raw_rows = []
    for label, expected_count in (('points107', 4096), ('points123_1', 4096), ('points123_2', 2141)):
        pm = json.loads(LOCKS[label][0].read_bytes())
        assert pm['epsilon'] == '1/128' and pm['target'] == TARGET
        assert pm['kernel_sha256'] == KERNEL and pm['model'] == MODEL
        assert len(pm['records']) == expected_count
        for row in pm['records']:
            rec = row['record']
            assert rec['epsilon'] == '1/128' and rec['model'] == MODEL and rec['kernel_sha256'] == KERNEL
            assert rec['qualification_target'] == TARGET
            assert row['point'] == rec['point'] and row['pair_cache_key'] == rec['pair_cache_key']
            assert row['pair_cache_key'] == ['1/128', ','.join(map(str, map(F, row['point'])))]
            assert contains(box(row['node'], roots), row['point'])
            assert contains(box(ancestor(row['node'])['node'], roots), row['point'])
            if label != 'points107':
                task = taskmap[tuple(row['node'])]
                assert task['point'] == row['point'] and task['pair_cache_key'] == row['pair_cache_key']
                assert task['ordinal'] == row['global_ordinal']
                raw_rows.append(row)
    assert len(raw_rows) == 6237
    state = json.loads(LOCKS['state123'][0].read_bytes())
    assert state['status'] == 'finite_original_prefix_complete' and state['target'] == TARGET
    assert not state['faults'] and not state['inflight'] and not state['obstructions_H128']
    assert state['counts']['new_complete'] == state['counts']['new_anchor_starts'] == 6237
    assert len(state['new_splits']) == 5665 and len(state['queue']) == 11330
    assert len(state['new_passed']) == 572
    cache = state['new_cache_by_epsilon']['1/128']
    assert len(cache) == 6237
    for row in raw_rows:
        assert cache[row['pair_cache_key'][1]] == row['record']
    queue = set(map(tuple, state['queue']))
    assert len(queue) == 11330
    candidate_rows = []
    used_parents = set()
    outside = 0
    for split_index, split in enumerate(state['new_splits']):
        node = tuple(split['node'])
        assert node in taskmap and node not in used_parents
        used_parents.add(node)
        assert split['axis'] == len(node[1]) % 4
        assert split['children'] == [[node[0], node[1]+'0'], [node[0], node[1]+'1']]
        assert split['leaf_epsilon'] == '1/128' and split['target'] == TARGET
        cert = split['failed_certificate']
        assert cert['target'] == TARGET and cert['leaf_epsilon'] == '1/128'
        assert cert['point_key'] == taskmap[node]['pair_cache_key'][1]
        assert cert['penalty_coefficient'] == '1024' and cert['checked_points'] == 16
        assert not cert['passed']
        rec = cache[cert['point_key']]
        assert rec['point'] == taskmap[node]['point']
        assert len(rec['values']) == 5 and all(F(a) <= F(b) for a,b in rec['values'])
        assert F(rec['values'][0][0]) > 0
        inherited = task_ancestors[node]
        for side, child in enumerate(split['children']):
            assert tuple(child) in queue
            childbox = box(child, roots)
            outside += not contains(childbox, rec['point'])
            assert ancestor(child)['node'] == inherited['node']
            candidate_rows.append({'candidate_ordinal': len(candidate_rows), 'child_node': child,
                'direct_parent_node': list(node), 'parent_split_index': split_index, 'side': side,
                'point_key': cert['point_key'], 'record_sha256': digest(rec), 'saved_certificate_sha256': digest(cert),
                'H64_lower_source_node': inherited['node'], 'H64_lower': G105,
                'H64_source_leaf_sha256': inherited['source_leaf_sha256'],
                'H64_source_certificate_sha256': inherited['source_certificate_sha256']})
    assert set(map(lambda r: tuple(r['child_node']), candidate_rows)) == queue
    candidate_path = OUT/'prepared_candidates.jsonl'
    with candidate_path.open('x') as f:
        for row in candidate_rows:
            f.write(json.dumps(row, sort_keys=True, separators=(',', ':'))+'\n')
    result = {'status': 'source_and_ancestry_qualified_no_candidate_test', 'pid': os.getpid(),
        'source_mapper_matches': source_mapper_matches, 'accepted_g105_H64_portable_leaf_matches': portable_matches,
        'points107_in_original_H64_g105_leaves': 4096, 'points123_in_original_H64_g105_leaves': 6237,
        'tasks123_original877': sum(t['origin']=='untested_original877' for t in tasks),
        'direct_H128_parent_jets': len(used_parents), 'finite_child_candidates': len(candidate_rows),
        'child_candidates_with_parent_anchor_outside_child': outside,
        'uniform_anchor_H64_lower': G105, 'candidate_epsilon': '1/128', 'target': TARGET,
        'future_q_formula': '1024-768*(21/200000)/U128(c)',
        'saved_Delta_and_loggradient': 'same parent H128 certificate, no recomputation',
        'prepared_candidates_path': str(candidate_path), 'prepared_candidates_sha256': sha(candidate_path),
        'registration_sha256': sha(OUT/'source_qualification_registration.json'),
        'source_locks': locks, 'new_H_gradient_kernel_log_exp_calls': 0,
        'new_real_q_or_corner_evaluations': 0, 'actual_candidate_pass_count': None,
        'shared_dependencies': 'accepted075/029/099; accepted096/098 g105 H64 leaves;123 H128 raw/Delta/gradient independent and primary numerical qualification',
        'wall_seconds': time.monotonic()-began}
    write(OUT/'source_qualification.json', result)
    print(json.dumps({k:v for k,v in result.items() if k != 'source_locks'}, ensure_ascii=False))

if __name__ == '__main__':
    main()
