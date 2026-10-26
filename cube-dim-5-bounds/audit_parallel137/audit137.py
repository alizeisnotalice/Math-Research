"""PREPARED: only1153 successful certificates,18448 corners; no execution yet.

All89426 journal rows are checked as metadata/hash-chain/first-pass protocol.
Failed arithmetic is never repeated. No imports of130/131/137 author code.
"""
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import time

from exact137 import (box,canonical_roots,check_child,eq,f,full_tree,in_box,key,
                      le,new_q_check,positive_certificate)
from json_read137 import Document

HERE=Path(__file__).resolve().parent


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        while b:=handle.read(8*1024*1024):h.update(b)
    return h.hexdigest()


def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def load(p):return json.loads(Path(p).read_bytes())


def lines(p):
    with Path(p).open() as handle:
        for line in handle:yield json.loads(line)


def write_new(path,x):
    with Path(path).open('x') as handle:
        json.dump(x,handle,sort_keys=True,separators=(',',':'));handle.write('\n')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reviewed-plan-sha',required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();start=time.monotonic()
    planpath=HERE/'audit_plan.json';assert sha(planpath)==args.reviewed_plan_sha
    plan=load(planpath)
    for name,h in plan['code_locks'].items():assert sha(HERE/name)==h
    for path,h in plan['all_input_locks'].items():assert sha(path)==h,path
    P={name:Path(path) for name,path in plan['paths'].items()}
    assert not args.output.exists();args.output.mkdir()
    write_new(args.output/'run_registration.json',{'task':'audit_parallel137',
        'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),
        'reviewed_plan_sha256':args.reviewed_plan_sha,'maximum_pass_certificates':1153,
        'maximum_new_integer_corners':18448,'metadata_journal_rows':89426,
        'failed_corner_tests_repeated':0,'new_H_log_gradient_kernel_calls':0,'new_splits':0})
    children=load(P['children']);pairs=list(lines(P['pairs']));specplan=load(P['shard_plan'])
    master=load(P['master_state']);summary=load(P['master_summary'])
    assert len(children)==5427 and len(pairs)==89426 and len(specplan['shards'])==15
    assert master['status']==summary['status']=='complete_fixed89426_objects'
    assert master['faults'] is None and summary['master_fault'] is None and not master['new_splits']
    accepted105=load(P['accepted_g105'])
    assert accepted105['status']=='primary_accepts_complete_same_model_global_g105'
    assert accepted105['qualified_g']=='21/200000'
    assert accepted105['full098_state_sha256']==plan['all_input_locks'][str(P['state098'])]
    a123=load(P['audit123_manifest']);a130=load(P['audit130_manifest'])
    assert a123['status']=='both_original_finite_new_value_and_structure_batches_independently_qualified'
    assert a123['source_lock_sha256']==sha(P['audit123_source_lock'])
    assert a130['status']=='finite11330_new_curvature_component_independently_qualified_partial'
    assert a130['new_passed']==5903 and a130['pending']==5427
    R=canonical_roots();needed={p['source_id'] for p in pairs};sources={};poolids=set()
    assert len(needed)==10549
    for i,row in enumerate(lines(P['source_pool'])):
        sid=row['source_id'];assert sid==digest(row['core']) and sid not in poolids
        assert row['source_index']==i;poolids.add(sid)
        if sid in needed:sources[sid]=row
    assert len(poolids)==111331 and set(sources)==needed
    for i,child in enumerate(children):
        assert child['child_index']==i;check_child(child,R)
        block=pairs[child['pair_start']:child['pair_start']+child['pair_count']]
        assert len(block)==child['pair_count'] and 0<len(block)<=32
        assert [p['source_id'] for p in block]==[r['source_id'] for r in child['selected_sources']]
        for pair in block:
            assert pair['child_index']==i and pair['child_node']==child['node']
            assert sources[pair['source_id']]['source_index']==pair['source_index']
            assert pair['curvature_rule']==sources[pair['source_id']]['core']['curvature_rule']
    assert [p['pair_index'] for p in pairs]==list(range(89426))
    # Read all fixed journals once; do not recalculate failed margins.
    qualified={};best={};events=set();count=Counter();durable=[]
    for spec in specplan['shards']:
        sid=spec['shard_id'];data=load(spec['input_path']);folder=P['run_dir']/f'shard{sid:02d}'
        cp=load(folder/'checkpoint.json');assert cp['status']=='complete'
        fingerprint=cp['fingerprint']
        assert fingerprint['shard_id']==sid and fingerprint['input_sha256']==spec['input_sha256']
        assert fingerprint['core_sha256']==specplan['core_sha256']
        assert fingerprint['core_sha256']==sha(P['core131'])
        assert fingerprint['worker_code_sha256']==sha(P['worker137'])
        assert fingerprint['source_pool_sha256']==sha(P['source_pool'])
        assert fingerprint['child_inputs_sha256']==sha(P['children'])
        assert fingerprint['candidate_pairs_sha256']==sha(P['pairs'])
        assert fingerprint['prepared_manifest_sha256']==master['prepared_manifest_sha256']
        assert len(data['pairs'])==spec['pair_count'] and len(data['children'])==spec['child_count']
        for k,c in data['children'].items():assert c==children[int(k)] and int(k)%15==sid
        for k,s in data['sources'].items():assert s==sources[k]
        previous='0'*64;offset=ev=sk=0;first={}
        with (folder/'results.jsonl').open('rb') as handle:
            for cursor,raw in enumerate(handle):
                assert raw.endswith(b'\n');offset+=len(raw);row=json.loads(raw);h=row.pop('hash')
                assert row['seq']==cursor+1 and row['previous']==previous and digest(row)==h
                previous=h;pair=data['pairs'][cursor];gi=row['global_pair_index'];ci=row['child_index']
                assert pair==pairs[gi] and gi not in events and ci%15==sid;events.add(gi)
                assert gi==pair['pair_index'] and ci==pair['child_index'] and row['source_id']==pair['source_id']
                if row['kind']=='skipped':
                    assert str(ci) in first and row['reason']=='child already certified by one full16corner source';sk+=1
                else:
                    assert row['kind']=='evaluated' and str(ci) not in first;ev+=1
                    cert=row['certificate'];source=sources[row['source_id']]
                    assert cert['source_id']==source['source_id'] and cert['source_pool_object_sha256']==digest(source)
                    assert cert['child_node']==children[ci]['node'] and cert['checked_points']==16
                    for field in ('epsilon','point','target','Delta_lower','loggradient_directed','curvature_rule','penalty_coefficient'):
                        assert cert[field]==source['core'][field]
                    margin=f(cert['minimum_log_margin']);assert cert['passed']==(margin[0]>0)
                    if ci not in best or not le(cert['minimum_log_margin'],best[ci]['minimum_log_margin']):best[ci]=cert
                    if cert['passed']:first[str(ci)]=cert
                count[row['kind']]+=1
        assert cursor+1==spec['pair_count']==cp['next_pair']
        assert cp['last_event_hash']==previous and cp['valid_prefix_bytes']==offset
        assert cp['evaluated']==ev and cp['skipped']==sk and cp['qualified_children']==first
        assert not set(qualified)&set(first);qualified.update(first)
        worker=summary['worker_summaries'][sid]
        assert worker['checkpoint_sha256']==sha(folder/'checkpoint.json')
        assert worker['results_sha256']==sha(folder/'results.jsonl')
        assert cp['fingerprint']==load(folder/'worker_registration.json')['fingerprint']
        durable.append({'shard_id':sid,'rows':cursor+1,'evaluated':ev,'skipped':sk,'first_passes':len(first)})
    assert events==set(range(89426)) and count=={'evaluated':71812,'skipped':17614}
    assert qualified==master['qualified_children'] and len(qualified)==1153
    pending=[c for c in children if str(c['child_index']) not in qualified]
    assert [c['node'] for c in pending]==master['pending_children'] and len(pending)==4274
    # Bind each selected source to its saved ORIGINAL certificate and raw record.
    aliases=defaultdict(list)
    for source in sources.values():
        new_q_check(source)
        for alias in source['aliases']:aliases[(alias['origin'],alias['index'])].append((source,alias))
    original_records={};matched=set();raw_needed=defaultdict(list)
    def match(source,a,node,cert,record):
        assert node==a['source_node'] and digest(cert)==a['source_certificate_sha256']
        assert digest(record)==a['source_record_sha256']==source['raw_record_sha256']
        c=source['core'];assert record['point']==c['point']
        assert record['values'][0]==source['H_interval']
        assert cert['point_key']==','.join(c['point']) and cert['leaf_epsilon']==c['epsilon']
        assert cert['target']==c['target']=='11/100000'
        assert eq(cert['penalty_coefficient'],a['source_original_penalty'])
        assert eq(cert['Delta_lower'],c['Delta_lower'])
        assert all(eq(x,y) for rp,cp in zip(cert['loggradient_directed'],c['loggradient_directed']) for x,y in zip(rp,cp))
        if 'epsilon' in record:assert eq(record['epsilon'],c['epsilon'])
        original_records[source['source_id']]=True;matched.add((a['origin'],a['index']))
    with Document(P['state105']) as doc:
        fields=dict(doc.object_spans())
        for i,row in enumerate(doc.array(fields['passed'])):
            for source,a in aliases.get(('105_passed',i),[]):
                assert row['node']==a['source_node'] and row['leaf_epsilon']==source['core']['epsilon']
                raw_needed[(row['leaf_epsilon'],row['certificate']['point_key'])].append((source,a,row['node'],row['certificate']))
    for i,row in enumerate(lines(P['pending105'])):
        for source,a in aliases.get(('105_failed',i),[]):
            raw_needed[('1/64',row['certificate']['point_key'])].append((source,a,row['node'],row['certificate']))
    states={}
    for project in ('107','123'):
        state=load(P['state'+project]);states[project]=state
        proofs={r['certificate']['point_key']:(r['node'],r['certificate']) for r in state['new_passed']}
        for r in state['new_splits']:proofs[r['failed_certificate']['point_key']]=(r['node'],r['failed_certificate'])
        cache=state['new_cache_by_epsilon']['1/128']
        for i,(point,record) in enumerate(sorted(cache.items())):
            for source,a in aliases.get((project+'_H128',i),[]):
                node,cert=proofs[point];match(source,a,node,cert,record)
    ancestry=defaultdict(list)
    for source in sources.values():
        if source['H64_g105_ancestry']:ancestry[key(source['H64_g105_ancestry']['node'])].append(source)
    h64nodes=set();ancestor_matches=set();base_splits=[]
    with Document(P['state098']) as doc:
        fields=dict(doc.object_spans());assert doc.load_span(fields['roots'])==R
        for row in doc.array(fields['passed']):
            node=key(row['node'])
            if row['leaf_epsilon']=='1/64':h64nodes.add(node)
            for source in ancestry.get(node,[]):
                a=source['H64_g105_ancestry'];assert row['leaf_epsilon']=='1/64'
                assert row['target']==row['certificate']['target']=='21/200000'
                assert row['certificate']['passed'] and f(row['certificate']['minimum_log_margin'])[0]>0
                assert digest(row)==a['source_leaf_sha256']
                assert a['source_state_sha256']==plan['all_input_locks'][str(P['state098'])]
                assert a['point_inside_ancestor'] is True and in_box(source['core']['point'],row['node'],R)
                for alias in source['aliases']:
                    rid,path=key(alias['source_node']);assert rid==node[0] and path.startswith(node[1])
                ancestor_matches.add(source['source_id'])
        for row in doc.array(fields['splits']):base_splits.append(row)
        eps_spans=dict(doc.object_spans(fields['cache_by_epsilon'][0]))
        for epsilon,span in eps_spans.items():
            for point,value_span in doc.object_spans(span[0]):
                use=raw_needed.get((epsilon,point))
                if use:
                    record=doc.load_span(value_span)
                    for source,a,node,cert in use:match(source,a,node,cert,record)
    assert set(original_records)==needed and matched==set(aliases)
    for source in sources.values():
        c=source['core'];a=source['H64_g105_ancestry']
        if a:assert source['source_id'] in ancestor_matches
        elif c['epsilon']=='1/128':
            assert all(x['origin']=='105_passed' for x in source['aliases'])
            for x in source['aliases']:
                rid,path=key(x['source_node']);assert not any((rid,path[:j]) in h64nodes for j in range(len(path)+1))
    # Only the1153 successful new certificates get actual new corner arithmetic.
    with (args.output/'winning_integer_checks.jsonl').open('x') as out:
        for ci,cert in sorted(qualified.items(),key=lambda z:int(z[0])):
            source=sources[cert['source_id']];answer=positive_certificate(source,children[int(ci)],cert,R)
            out.write(json.dumps({'child_index':int(ci),'node':cert['child_node'],'source_id':cert['source_id'],**answer},sort_keys=True,separators=(',',':'))+'\n')
    # Exact residual and gap order; failed margin values are read, never recomputed.
    ordered=list(lines(P['residual']));assert len(ordered)==4274
    residual_ids=set()
    def compare(a,b):
        x,y=a['best_crossanchor_failed_certificate']['minimum_log_margin'],b['best_crossanchor_failed_certificate']['minimum_log_margin']
        if not eq(x,y):return -1 if le(x,y) else 1
        da,db=len(a['node'][1]),len(b['node'][1])
        if da!=db:return -1 if da<db else 1
        return (key(a['node'])>key(b['node']))-(key(a['node'])<key(b['node']))
    for i,row in enumerate(ordered):
        ci=row['child_index'];assert ci not in residual_ids;residual_ids.add(ci)
        assert str(ci) not in qualified and row['gap_rank']==i
        assert row['best_crossanchor_failed_certificate']==best[ci] and not best[ci]['passed']
        assert row['node']==children[ci]['node'] and row['box_lower']==children[ci]['box_lower'] and row['box_upper']==children[ci]['box_upper']
        assert row['prior_direct_parent_newq_failure']==children[ci]['prior_direct_parent_newq_failure']
        assert eq(row['physical_volume'],f'1/{32**4*(1<<len(row["node"][1]))}')
        if i:assert compare(ordered[i-1],row)<=0
    assert residual_ids=={c['child_index'] for c in pending}
    # Rebuild full forest without reading the giant123 raw-cache seed.
    mapper=load(P['accepted_mapper']);s123lock=load(P['audit123_source_lock'])
    assert s123lock['input_locks'][str(P['accepted_mapper'])]==sha(P['accepted_mapper'])
    assert len(mapper)==98007 and [r['seed_passed_index'] for r in mapper]==list(range(98007))
    oldpass=[r['node'] for r in mapper]+[r['node'] for r in states['123']['new_passed']]
    result128=list(lines(P['results128']));assert len(result128)==11330
    assert len({key(r['child_node']) for r in result128})==11330
    assert {key(r['child_node']) for r in result128}==set(map(key,states['123']['queue']))
    assert {key(r['child_node']) for r in result128 if not r['passed']}=={key(c['node']) for c in children}
    passed=oldpass+[r['child_node'] for r in result128 if r['passed']]+[r['child_node'] for r in qualified.values()]
    splits=base_splits+states['107']['new_splits']+states['123']['new_splits']
    for split in splits:
        rid,path=key(split['node']);assert split['axis']==len(path)%4
        assert split['children']==[[rid,path+'0'],[rid,path+'1']]
    assert len(passed)==105635 and len(splits)==57549
    tree=full_tree(passed+[r['node'] for r in pending],[r['node'] for r in splits])
    for path,h in plan['all_input_locks'].items():assert sha(path)==h,path
    result={'task':'audit_parallel137','status':'all1153_pass_certificates_and_full_protocol_structure_qualified_partial',
        'pid':os.getpid(),'reviewed_plan_sha256':args.reviewed_plan_sha,'winning_certificates':1153,
        'new_integer_corners':18448,'failed_corner_tests_repeated':0,'journal_rows':89426,
        'evaluated_receipts':71812,'skip_receipts':17614,'all_first_pass_skip_hash_chains_checked':True,
        'selected_sources_raw_certificate_record_point_bound':len(sources),
        'winning_integer_checks_sha256':sha(args.output/'winning_integer_checks.jsonl'),
        'pending':4274,'passed':105635,'tree':tree,'durable_shards':durable,
        'all_source_hashes_unchanged':True,'globalg110_qualified':False,'S_credit':0,
        'new_H_log_gradient_kernel_calls':0,'new_splits':0,'Cn_C5_I_updated':False,
        'wall_seconds':time.monotonic()-start}
    write_new(args.output/'audit_result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
