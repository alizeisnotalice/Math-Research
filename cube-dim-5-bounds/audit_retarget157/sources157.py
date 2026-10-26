"""Prepared original-source/closed-support/forest gates, no q/log/corners."""
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
from integer157 import bounds,inside,roots
from json_read157 import Document
from rational157 import Q


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while chunk := handle.read(8*1024*1024):
            h.update(chunk)
    return h.hexdigest()


def digest(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def load(path):
    return json.loads(Path(path).read_bytes())


def lines(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as handle:
        for line in handle:
            yield json.loads(line)


def source_gates(plan):
    P = plan['paths']
    manifest151 = load(P['manifest151'])
    def dataset(kind):
        for part in manifest151['datasets'][kind]:
            yield from lines(Path(P['portable151'])/part['path'])
    objects = {}
    for i,row in enumerate(lines(P['objects149'])):
        c = row['core'];cid = row['canonical_id']
        assert i == row['canonical_index'] and digest(c) == cid and cid not in objects
        assert c['epsilon'] in ('1/64','1/128') and c['target110'] == '11/100000'
        assert len(c['point']) == 5 and Q(c['point'][0]) == 1 and len(c['loggradient_directed']) == 4
        assert row['typed_cache_key'] == [c['epsilon'],','.join(c['point'])]
        objects[cid] = {k:row[k] for k in ['canonical_index','H_interval','typed_cache_key','raw_record_sha256','source_alias_indices']}
    assert len(objects) == 107530
    jets_needed = defaultdict(list);coefficient_seen = set()
    for row in dataset('coefficients'):
        oid = row['canonical_id'];obj = objects[oid]
        assert oid not in coefficient_seen;coefficient_seen.add(oid)
        assert digest(row['core']) == oid and row['raw_record_sha256'] == obj['raw_record_sha256']
        jets_needed[row['jet_id']].append(oid)
    assert coefficient_seen == set(objects)
    jet_seen = set()
    for row in dataset('jets'):
        expected = row.pop('minimal_projection_sha256')
        assert digest(row) == expected and row['jet_id'] not in jet_seen
        jet_seen.add(row['jet_id'])
        assert row['model'] == 'n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
        point_key = ','.join(str(Q(x)) for x in row['point'])
        assert digest({'epsilon':row['epsilon'],'point_key':point_key,
                       'raw_record_sha256':row['raw_record_sha256']}) == row['jet_id']
        for oid in jets_needed.get(row['jet_id'],[]):
            obj = objects[oid]
            assert obj['typed_cache_key'] == [row['epsilon'],point_key]
            assert obj['H_interval'] == row['values'][0] and obj['raw_record_sha256'] == row['raw_record_sha256']
    assert set(jets_needed) <= jet_seen
    aliases = {}
    for row in lines(P['provenance149']):
        i = row['source_alias_index'];assert i not in aliases
        aliases[i] = {'canonical_id':row['canonical_id'],'node':row['node'],'hash':digest(row),
                      'proof_node':None if row['curvature_proof'] is None else row['curvature_proof']['node'],
                      'margin':row['saved_g110_margin'],'raw_record_sha256':row['raw_record_sha256']}
    assert len(aliases) == 111201
    supports = {tuple(row['node']):row for row in dataset('supports')}
    assert len(supports) == 3812
    accepted105 = load(P['accepted105'])
    assert accepted105['status'] == 'primary_accepts_complete_same_model_global_g105'
    assert accepted105['full098_state_sha256'] == plan['input_locks'][P['state098']]
    found = set()
    with Document(P['state098']) as doc:
        fields = dict(doc.object_spans())
        for row in doc.array(fields['passed']):
            node = tuple(row['node'])
            if node not in supports:
                continue
            s = supports[node];c = s['certificate'];original = row['certificate']
            assert s['source_state_sha256'] == plan['input_locks'][P['state098']]
            assert digest(row) == s['source_leaf_sha256']
            assert original == s['original_certificate'] and digest(original) == s['original_certificate_sha256']
            assert row['leaf_epsilon'] == '1/64' and row['target'] == '21/200000'
            assert original['passed'] and Q(original['minimum_log_margin']) > Q(0)
            assert c['passed'] and c['checked_points'] == 16 and c['epsilon'] == '1/64'
            assert c['target'] == '21/200000' and c['penalty_coefficient'] == '256'
            assert c['minimum_log_margin'] == original['minimum_log_margin']
            assert s['typed_cache_key'][0] == '1/64'
            found.add(node)
    assert found == set(supports)
    table = roots();all_nodes = set();used_objects = set();byroot = defaultdict(list)
    leaf_data = {}
    for i,row in enumerate(lines(P['leaves149'])):
        assert row['leaf_index'] == i and tuple(row['node']) not in all_nodes and row['closed_box']
        all_nodes.add(tuple(row['node']));byroot[row['node'][0]].append(row['node'][1])
        oid = row['canonical_id'];alias = aliases[row['source_alias_index']]
        assert row['source_alias_index'] in objects[oid]['source_alias_indices']
        assert alias['canonical_id'] == oid and alias['node'] == row['node']
        assert alias['margin'] == row['saved_g110_margin'] and Q(alias['margin']) > Q(0)
        assert alias['raw_record_sha256'] == objects[oid]['raw_record_sha256']
        leaf_data[i] = digest({k:row[k] for k in ['leaf_index','node','canonical_id','source_alias_index','epsilon','box_lower','box_upper','closed_box','saved_g110_margin']})
        for (a,b,d),lo,hi in zip(bounds(row['node'],table),row['box_lower'],row['box_upper']):
            assert Q(a,d) == Q(lo) and Q(b,d) == Q(hi)
    assert len(all_nodes) == 111201 and len(byroot) == 52360
    ancestors = set()
    for rid,path in all_nodes:
        ancestors.update((rid,path[:j]) for j in range(len(path)))
    assert len(ancestors) == 58841 and not ancestors & all_nodes
    for rid in range(52360):
        paths = sorted(byroot[rid]);assert all(not b.startswith(a) for a,b in zip(paths,paths[1:]))
        depth = max(map(len,paths));assert sum(1<<(depth-len(p)) for p in paths) == 1<<depth
    shard_plan = load(P['finite155']);seen_leaves = set();seen_sources = set();used_supports = set()
    for spec in shard_plan['shards']:
        data = load(spec['input_path']);sid = spec['shard_id']
        assert data['shard_id'] == sid and len(data['sources']) == spec['objects'] and len(data['leaves']) == spec['leaves']
        for source in data['sources']:
            cid = source['canonical_id'];obj = objects[cid];c = source['old_core'];alias = aliases[source['old_proof_alias_index']]
            assert cid not in seen_sources and obj['canonical_index']%15 == sid
            seen_sources.add(cid)
            assert digest(c) == cid and source['canonical_index'] == obj['canonical_index']
            assert source['old_H_interval'] == obj['H_interval'] and source['raw_record_sha256'] == obj['raw_record_sha256']
            assert source['typed_cache_key'] == obj['typed_cache_key']
            assert source['old_proof_alias_index'] in obj['source_alias_indices']
            assert source['old_proof_alias_sha256'] == alias['hash'] and alias['canonical_id'] == cid
            mass = source['mass_qualification']
            if c['curvature_rule_class'] == '129_anchor_fixed_Jensen':
                assert c['epsilon'] == '1/128' and mass is not None
                node = tuple(mass['support_node']);s = supports[node]
                assert alias['proof_node'] == list(node) and inside(c['point'],list(node),table)
                assert mass['lower'] == '21/200000' and mass['lower_epsilon'] == '1/64' and mass['upper_epsilon'] == '1/128'
                assert mass['upper'] == obj['H_interval'][1] and mass['point'] == c['point']
                assert mass['support_source_leaf_sha256'] == s['source_leaf_sha256'] and mass['support_source_state_sha256'] == s['source_state_sha256']
                assert mass['support_certificate_sha256'] == digest(s['certificate']) and mass['original_support_certificate_sha256'] == s['original_certificate_sha256']
                assert mass['source_alias_index'] == source['old_proof_alias_index'] and mass['source_alias_sha256'] == alias['hash']
                assert mass['point_in_support_closedbox'] and mass['no_new_raw_lower_lookup'] and mass['global_H128_floor_not_used_as_H64']
                used_supports.add(node)
            else:
                assert c['curvature_rule_class'] == '075_full_frame' and mass is None
        for leaf in data['leaves']:
            i = leaf['leaf_index'];assert i not in seen_leaves and digest(leaf) == leaf_data[i]
            assert objects[leaf['canonical_id']]['canonical_index']%15 == sid
            assert leaf['epsilon'] == objects[leaf['canonical_id']]['typed_cache_key'][0]
            seen_leaves.add(i);used_objects.add(leaf['canonical_id'])
    assert len(seen_sources) == len(used_objects) == 107530 and seen_sources == set(objects)
    assert seen_leaves == set(range(111201)) and used_supports == set(supports)
    accepted = load(P['acceptance017'])
    assert accepted['status'] == 'primary_accepts_complete_global_gstar_and_C5_upper'
    assert accepted['closed_forest'] == {'roots':52360,'passed':111201,'splits':58841,'pending':0}
    return shard_plan['shards'], {'objects':107530,'leaves':111201,'roots':52360,'splits':58841,
                                  'original_H64_supports':3812,'old_closed110_domain_bound':True}
