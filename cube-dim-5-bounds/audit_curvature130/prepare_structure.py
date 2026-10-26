"""Prepared independent structure gate. No real q or corner arithmetic.

Requires the reviewed registration SHA. Every output is exclusive creation.
Does not import or execute any author checker, kernel, or old baseline.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import time

from integer_core import (JsonView, box, contains, digest, frontier, node_key,
                          parse, roots, sha, split_partition, text)

HERE = Path(__file__).resolve().parent
TARGET = '11/100000'
G105 = '21/200000'
MODEL = 'n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
KERNEL = '128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0'


def exclusive(path,obj):
    with Path(path).open('x') as f:
        json.dump(obj,f,ensure_ascii=False,sort_keys=True,separators=(',',':'))
        f.write('\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--reviewed-registration-sha',required=True)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    start = time.monotonic()
    regpath = HERE/'preparation_registration.json'
    assert sha(regpath) == args.reviewed_registration_sha
    reg = json.loads(regpath.read_bytes())
    assert reg['new_q_and_corners_authorized'] is False
    for name,expected in reg['code_locks'].items(): assert sha(HERE/name) == expected
    paths = {}
    for name,item in reg['source_locks'].items():
        path = Path(item['path'])
        assert sha(path) == item['sha256'],name
        paths[name] = path
    assert not args.output.exists()
    args.output.mkdir()
    exclusive(args.output/'run_registration.json',{'task':'audit_curvature130_structure',
        'pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),
        'reviewed_preparation_sha256':args.reviewed_registration_sha,
        'new_H_log_gradient_kernel_calls':0,'new_q_evaluations':0,'new_corner_evaluations':0})
    R = roots()
    assert len(R) == 52360
    acceptance = json.loads(paths['g105_acceptance'].read_bytes())
    manifest = json.loads(paths['packet_manifest'].read_bytes())
    checker = json.loads(paths['packet_checker'].read_bytes())
    assert acceptance['status'] == 'primary_accepts_complete_same_model_global_g105'
    assert acceptance['qualified_g'] == manifest['g'] == checker['g'] == G105
    assert checker['status'] == 'PASS_portable_compact_closed_root_packet'
    assert checker['packet_manifest_sha256'] == sha(paths['packet_manifest'])
    assert checker['source_locks']['state'] == acceptance['full098_state_sha256']
    for label in ('packet_manifest','packet_checker'):
        relative = str(paths[label].relative_to(paths['g105_acceptance'].parent.parent))
        assert acceptance['input_locks'][relative] == sha(paths[label])
    a116 = json.loads(paths['audit116_manifest'].read_bytes())
    n116 = json.loads(paths['audit116_numeric_index'].read_bytes())
    a122 = json.loads(paths['audit122_manifest'].read_bytes())
    r122 = json.loads(paths['audit122_result'].read_bytes())
    assert a116['status'] == 'locked_first4096_H128_numeric_and_partial_cover_independently_qualified'
    assert a116['all_new_H128_lower_above_target'] is True
    assert n116['status'] == 'all4096_H128_plus4grad_strictly_verified'
    assert n116['verified_points'] == 4096
    assert a122['status'] == 'locked_all12984_single_anchor_candidates_qualified_partial'
    assert r122['status'] == 'all6492children_12984_single_anchor_candidates_independently_qualified_partial'
    assert r122['passed'] == 98007 and r122['queue'] == 6237
    assert a122['frozen_external_inputs'][str(paths['state118'])] == sha(paths['state118'])

    pending = [json.loads(line) for line in paths['pending105'].open()]
    assert len(pending) == 4973
    old = {node_key(r['node']):r for r in pending}
    assert len(old) == 4973
    ancestor_index = defaultdict(dict)
    for key,row in old.items():
        assert row['leaf_epsilon'] == '1/64'
        assert row['certificate']['target'] == TARGET
        assert row['certificate']['leaf_epsilon'] == '1/64'
        bounds = box(key,R)
        assert [text(lo) for lo,hi in bounds] == row['box_lower']
        assert [text(hi) for lo,hi in bounds] == row['box_upper']
        ancestor_index[key[0]][key[1]] = row
    mapper_matches = 0
    for line in paths['mapper105'].open():
        row = json.loads(line)
        key = node_key(row['node'])
        if key in old:
            for field in ('source_leaf_index','source_leaf_sha256','source_certificate_sha256',
                          'source_record_sha256','box_lower','box_upper','leaf_epsilon','pair_cache_key'):
                assert row[field] == old[key][field],(key,field)
            mapper_matches += 1
    assert mapper_matches == 4973

    def ancestor(node):
        rid,path = node_key(node)
        matches = [ancestor_index[rid][path[:i]] for i in range(len(path)+1)
                   if path[:i] in ancestor_index[rid]]
        assert len(matches) == 1,node
        return matches[0]

    packet_matches = set()
    packet_splits = {}
    packet_leaves = []
    for chunk in manifest['chunks']:
        path = paths['packet_manifest'].parent/chunk['file']
        assert sha(path) == chunk['sha256']
        packet = json.loads(gzip.decompress(path.read_bytes()))
        for row in packet['leaves']:
            key = node_key(row['node'])
            packet_leaves.append(key)
            if key in old:
                assert key not in packet_matches
                assert row['epsilon'] == row['cache_epsilon'] == '1/64'
                assert row['source_penalty_coefficient'] == '256'
                assert row['curvature_mode'] == 'full_frame_16'
                assert row['source_quadratic'] == 'Q_full=I4/16'
                assert ['1',*row['point_tail']] == old[key]['pair_cache_key'][1].split(',')
                packet_matches.add(key)
        for row in packet['splits']:
            key = node_key(row['node'])
            assert key not in packet_splits
            split_partition(key,row['children'],row['axis'],R)
            packet_splits[key] = (row['axis'],row['children'])
    assert packet_matches == set(old)
    source_tree = frontier(packet_leaves,list(packet_splits))
    assert source_tree['leaves'] == 100998 and source_tree['splits'] == 48638

    # Bind all retained leaf records to the three previously frozen source lists.
    expected_passed = []
    with JsonView(paths['state105']) as view:
        for row in view.array('passed'): expected_passed.append((node_key(row['node']),digest(row)))
    assert len(expected_passed) == 96025
    s107 = json.loads(paths['state107'].read_bytes())
    s118 = json.loads(paths['state118'].read_bytes())
    assert s118['source105_sha256'] == sha(paths['state105'])
    assert s118['source107_sha256'] == sha(paths['state107'])
    assert len(s107['new_passed']) == 850 and len(s118['new_passed']) == 1132
    expected_passed.extend((node_key(row['node']),digest(row)) for row in s107['new_passed'])
    expected_passed.extend((node_key(row['node']),digest(row)) for row in s118['new_passed'])
    expected_splits = dict(packet_splits)
    for row in s107['new_splits']:
        key = node_key(row['node'])
        assert key not in expected_splits
        split_partition(key,row['children'],row['axis'],R)
        expected_splits[key] = (row['axis'],row['children'])
    assert len(expected_splits) == 51884
    seed_passed = []
    seed_splits = []
    with JsonView(paths['seed123']) as view:
        assert view.get('roots') == list(map(list,R))
        assert view.get('model') == MODEL and view.get('target') == TARGET
        assert not view.get('inflight') and not view.get('resume_jobs')
        assert view.get('queue') == s118['queue']
        for row in view.array('passed'): seed_passed.append((node_key(row['node']),digest(row)))
        assert seed_passed == expected_passed
        for row in view.array('splits'):
            key = node_key(row['node'])
            split_partition(key,row['children'],row['axis'],R)
            assert expected_splits[key] == (row['axis'],row['children'])
            seed_splits.append(key)
    assert len(seed_passed) == 98007 and set(seed_splits) == set(expected_splits)
    seed_tree = frontier([key for key,h in seed_passed]+list(map(node_key,s118['queue'])),seed_splits)

    tasks = json.loads(paths['tasks123'].read_bytes())
    taskmap = {node_key(t['node']):t for t in tasks}
    assert len(taskmap) == len(tasks) == 6237
    assert set(taskmap) == set(map(node_key,s118['queue']))
    point_records = {}
    for name,count in (('points107',4096),('points123_1',4096),('points123_2',2141)):
        pm = json.loads(paths[name].read_bytes())
        assert pm['epsilon'] == '1/128' and pm['target'] == TARGET
        assert pm['model'] == MODEL and pm['kernel_sha256'] == KERNEL
        assert len(pm['records']) == count
        for row in pm['records']:
            node = node_key(row['node'])
            record = row['record']
            assert record['epsilon'] == '1/128' and record['qualification_target'] == TARGET
            assert record['model'] == MODEL and record['kernel_sha256'] == KERNEL
            assert record['point'] == row['point']
            assert record['pair_cache_key'] == row['pair_cache_key']
            assert record['pair_cache_key'] == ['1/128',','.join(text(parse(x)) for x in row['point'])]
            assert contains(box(node,R),row['point'])
            assert contains(box(ancestor(node)['node'],R),row['point'])
            if name != 'points107':
                task = taskmap[node]
                assert task['point'] == row['point'] and task['pair_cache_key'] == row['pair_cache_key']
                assert task['ordinal'] == row['global_ordinal']
                assert node not in point_records
                point_records[node] = record
    assert len(point_records) == 6237

    s123 = json.loads(paths['state123'].read_bytes())
    assert s123['seed_sha256'] == sha(paths['seed123'])
    assert s123['status'] == 'finite_original_prefix_complete' and s123['target'] == TARGET
    assert not s123['inflight'] and not s123['faults'] and not s123['obstructions_H128']
    assert s123['counts']['new_complete'] == s123['counts']['new_anchor_starts'] == 6237
    assert len(s123['new_passed']) == 572 and len(s123['new_splits']) == 5665
    queue = set(map(node_key,s123['queue']))
    assert len(queue) == len(s123['queue']) == 11330
    assert len(s123['new_cache_by_epsilon']['1/128']) == 6237
    assert set(point_records) == {node_key(r['node']) for r in s123['new_passed']+s123['new_splits']}
    for node,record in point_records.items():
        assert s123['new_cache_by_epsilon']['1/128'][record['pair_cache_key'][1]] == record
    newpass = {node_key(row['node']) for row in s123['new_passed']}
    newparents = set()
    outside = 0
    candidates = []
    for index,row in enumerate(s123['new_splits']):
        parent = node_key(row['node'])
        assert parent not in newparents and parent not in newpass
        newparents.add(parent)
        split_partition(parent,row['children'],row['axis'],R)
        cert = row['failed_certificate']
        record = point_records[parent]
        assert cert['point_key'] == record['pair_cache_key'][1]
        assert cert['target'] == row['target'] == TARGET
        assert cert['leaf_epsilon'] == row['leaf_epsilon'] == '1/128'
        assert cert['penalty_coefficient'] == '1024' and cert['checked_points'] == 16
        assert cert['passed'] is False
        inherited = ancestor(parent)
        assert contains(box(inherited['node'],R),record['point'])
        for side,child in enumerate(row['children']):
            key = node_key(child)
            assert key in queue and ancestor(child)['node'] == inherited['node']
            outside += not contains(box(key,R),record['point'])
            candidates.append({'candidate_ordinal':len(candidates),'child_node':child,
                'direct_parent_node':row['node'],'parent_split_index':index,'side':side,
                'point_key':cert['point_key'],'record_sha256':digest(record),
                'saved_certificate_sha256':digest(cert),'H64_lower_source_node':inherited['node'],
                'H64_lower':G105,'H64_source_leaf_sha256':inherited['source_leaf_sha256'],
                'H64_source_certificate_sha256':inherited['source_certificate_sha256']})
    assert newparents.union(newpass) == set(taskmap)
    assert {node_key(c['child_node']) for c in candidates} == queue and len(candidates) == 11330
    with paths['author_prepared_candidates'].open() as f:
        author_rows = [json.loads(line) for line in f]
    assert author_rows == candidates
    final_passed = [key for key,h in seed_passed]+list(newpass)
    allsplits = seed_splits+list(newparents)
    final_tree = frontier(final_passed+list(queue),allsplits)
    assert len(final_passed) == 98579 and len(allsplits) == 57549
    with (args.output/'independent_candidates.jsonl').open('x') as f:
        for row in candidates: f.write(json.dumps(row,sort_keys=True,separators=(',',':'))+'\n')
    receipt = {'task':'audit_curvature130_structure','status':'structure_and_ancestry_checked_no_q_or_corner_test',
        'pid':os.getpid(),'preparation_registration_sha256':args.reviewed_registration_sha,
        'source_mapper_matches':mapper_matches,'original_H64_portable_leaves':len(packet_matches),
        'H128_points107_in_original_H64_leaf':4096,'H128_points123_in_original_H64_leaf':6237,
        'direct_parent_jets':5665,'children_once_each':len(candidates),
        'parent_anchor_outside_child':outside,'source_tree':source_tree,'seed_tree':seed_tree,
        'current_partial_tree':final_tree,'old_and_123_passed':98579,'queue':11330,
        'candidate_packet_sha256':sha(args.output/'independent_candidates.jsonl'),
        'source_state123_sha256':sha(paths['state123']),'new_q_evaluations':0,'new_corner_evaluations':0,
        'new_H_gradient_log_kernel_calls':0,'numeric_source123_final_receipt_required':True,
        'global_g110_qualified':False,'Cn_C5_I_updated':False,'wall_seconds':time.monotonic()-start}
    exclusive(args.output/'structure_receipt.json',receipt)
    print(json.dumps(receipt,ensure_ascii=False))


if __name__ == '__main__': main()
