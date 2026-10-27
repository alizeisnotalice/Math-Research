"""Prepared158 source and fixed-pair gates. No q/log/corner evaluation."""
import gzip
import hashlib
import json
from pathlib import Path
from rational168 import Q


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while b := handle.read(1<<20):h.update(b)
    return h.hexdigest()


def digest(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def read(path):
    return json.loads(Path(path).read_bytes())


def lines(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as handle:
        for line in handle:yield json.loads(line)


def source_gates(plan):
    P = plan['paths'];fixed = read(P['finite158'])
    old_qualification = read(P['result157'])
    assert old_qualification['status'] == 'finite155_newq_log_jointcorners_and_complete110_cover_independently_qualified'
    assert old_qualification['objects'] == 107530 and old_qualification['leaves'] == 111201
    assert old_qualification['pending115'] == 1240 and old_qualification['passed115'] == 109961
    old_binding = read(P['binding157']);logref = fixed['existing155_log']
    assert logref['sha256'] == old_binding['common_log_sha256'] == sha(logref['path'])
    log_upper = read(logref['path'])['upper']
    old_log_check = read(P['logcheck157'])
    assert old_log_check['upper'] == log_upper
    assert old_log_check['status'] == 'positive_polynomial_integral_and_positive_tail_bound_passed'
    assert logref['new_log_evaluations'] == 0
    selected = {r['canonical_id']:r for r in lines(P['needed_sources158'])}
    assert len(selected) == 4769 and len({s['new_core_sha256'] for s in selected.values()}) == 4769
    old_pending = {r['leaf_index']:r for r in read(P['pending155'])};assert len(old_pending) == 1240
    failed_canonical = {r['canonical_id'] for r in old_pending.values()};failed_signatures = {}
    needed_aliases = {}
    for cid,source in selected.items():
        assert source['canonical_id'] == min(a['canonical_id'] for a in source['aliases'])
        assert source['new_core_sha256'] == digest(source['core'])
        assert source['source_common_log_sha256'] == logref['sha256']
        assert source['new_q_or_Delta_recalculated'] is False and source['new_global_floor_not_used_for_mass']
        for alias in source['aliases']:
            key = alias['source_path'],alias['source_line_index']
            assert key not in needed_aliases;needed_aliases[key] = (cid,alias)
    matched = set();full_old = {}
    for i in range(15):
        qpath = Path(P['output155'])/f'shard{i:02d}/q_updates.jsonl'
        for ordinal,event in enumerate(lines(qpath)):
            key = str(qpath),ordinal
            body = event['body']
            if key not in needed_aliases and body['canonical_id'] not in failed_canonical:continue
            expected = {'epsilon':body['epsilon'],'point':body['point'],'target':'23/200000',
                        'Delta_lower':body['Delta115_lower'],'loggradient_directed':body['loggradient_directed'],
                        'q_new':body['q_new'],'q_rule':'153_time_shape_uniform_or_original_g105_mass'}
            if body['canonical_id'] in failed_canonical:failed_signatures[body['canonical_id']] = digest(expected)
            if key not in needed_aliases:continue
            cid,alias = needed_aliases[key];source = selected[cid]
            assert expected == source['core']
            assert body['canonical_id'] == alias['canonical_id'] and event['hash'] == alias['source_event_hash']
            assert digest(body) == alias['source_body_sha256'] and body['raw_record_sha256'] == alias['raw_record_sha256']
            assert body['source_penalty_sha256'] == source['source153_penalty_sha256']
            if body['canonical_id'] == cid:
                assert source['saved_Delta110_lower'] == body['Delta110_lower']
                assert source['Delta110_identity_source_body_sha256'] == digest(body)
            matched.add(key)
        for event in lines(Path(P['output155'])/f'shard{i:02d}/leaf_results.jsonl'):
            row = event['body'];index = row['leaf_index'];assert index not in full_old
            assert Q(row['margin110_new']) > Q(0)
            full_old[index] = {k:row[k] for k in ['leaf_index','node','canonical_id','epsilon','margin110_new','margin115','passed115']}
    assert matched == set(needed_aliases) and set(full_old) == set(range(111201))
    assert set(failed_signatures) == failed_canonical
    pooled = 0;seen_selected = set()
    for row in lines(P['pool_index158']):
        assert row['source_index'] == pooled;pooled += 1
        cid = row['canonical_id']
        if cid in selected:
            assert all(selected[cid][k] == v for k,v in row.items());seen_selected.add(cid)
    assert pooled == fixed['source_library_unique_newcores'] and seen_selected == set(selected)
    children = read(P['children158']);pairs = list(lines(P['pairs158']))
    assert len(children) == 1240 and len(pairs) == 21494
    receiver_indices = set();pair_cursor = 0;seen_pairs = set();used_sources = set()
    leaves149 = {}
    for row in lines(P['leaves149']):
        if row['leaf_index'] in old_pending:leaves149[row['leaf_index']] = row
    assert set(leaves149) == set(old_pending)
    for i,child in enumerate(children):
        assert child['child_index'] == i and child['leaf_index'] not in receiver_indices
        index = child['leaf_index'];receiver_indices.add(index)
        assert child['prior155_failed115'] == old_pending[index]
        original = full_old[index];leaf = leaves149[index]
        assert child['prior155_failed115']['canonical_id'] == original['canonical_id'] == leaf['canonical_id']
        assert child['failed155_newcore_sha256'] == failed_signatures[original['canonical_id']]
        assert child['old_full110_remains_qualified']
        assert not original['passed115'] and original['margin110_new'] == old_pending[index]['margin110_new']
        assert original['margin115'] == old_pending[index]['margin115']
        assert child['node'] == original['node'] == leaf['node']
        assert child['box_lower'] == leaf['box_lower'] and child['box_upper'] == leaf['box_upper'] and leaf['closed_box']
        assert child['pair_start'] == pair_cursor and child['pair_count'] == len(child['selected_sources'])
        assert 0 < child['pair_count'] <= 32
        block = pairs[pair_cursor:pair_cursor+child['pair_count']];pair_cursor += child['pair_count']
        assert len({r['canonical_id'] for r in block}) == len(block)
        assert [r['canonical_id'] for r in block] == [r['canonical_id'] for r in child['selected_sources']]
        for pair in block:
            cid = pair['canonical_id'];source = selected[cid]
            used_sources.add(cid)
            assert pair['pair_index'] not in seen_pairs;seen_pairs.add(pair['pair_index'])
            assert pair['child_index'] == i and pair['leaf_index'] == index and pair['node'] == child['node']
            assert pair['source_index'] == source['source_index'] and pair['new_core_sha256'] == source['new_core_sha256']
            assert pair['new_core_sha256'] != child['failed155_newcore_sha256']
    assert receiver_indices == set(old_pending) and pair_cursor == 21494 and seen_pairs == set(range(21494))
    assert used_sources == set(selected)
    seen = set();owned = set()
    for spec in fixed['shards']:
        data = read(spec['input_path']);sid = spec['shard_id']
        assert len(data['pairs']) == spec['pairs'] and len(data['sources']) == spec['sources'] and len(data['children']) == spec['children']
        for cid,source in data['sources'].items():assert source == selected[cid]
        for text,child in data['children'].items():
            i = int(text);assert i%15 == sid and i not in owned and child == children[i];owned.add(i)
        for pair in data['pairs']:
            i = pair['pair_index'];assert i not in seen and pair == pairs[i];seen.add(i)
    assert owned == set(range(1240)) and seen == set(range(21494))
    projected = {}
    for part in fixed['full_old_margin_parts']:
        assert sha(part['absolute_path']) == part['sha256']
        for row in lines(part['absolute_path']):
            i = row['leaf_index'];assert i not in projected
            assert all(row[k] == full_old[i][k] for k in ['leaf_index','node','canonical_id','margin110_new','margin115','passed115'])
            projected[i] = full_old[i]
    assert set(projected) == set(full_old)
    return fixed,children,full_old,log_upper
