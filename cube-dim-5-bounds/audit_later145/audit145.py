"""PREPARED finite named-stage checker; no arithmetic has been executed.

Counts, paths, raw providers, qualification receipts and previous frontiers are
explicit frozen configuration. Never discover or run unknown future stages.
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


def rows_from_layer(layer,states,base_splits,plan):
    if layer['kind']=='raw098_splits':return base_splits
    path=Path(layer['path']);assert str(path) in plan['all_input_locks']
    if layer['kind']=='jsonl':rows=list(lines(path))
    else:
        data=load(path)
        for name in layer.get('fields',[]):data=data[name]
        rows=list(data.values()) if layer['kind']=='mapping_values' else data
    if 'passed_filter' in layer:rows=[r for r in rows if r['passed']==layer['passed_filter']]
    return rows


def nodes_from_layer(layer,states,base_splits,plan):
    rows=rows_from_layer(layer,states,base_splits,plan)
    return [r[layer['node_field']] for r in rows] if layer.get('node_field') else rows


def split_rows_from_layer(layer,states,base_splits,plan):
    return rows_from_layer(layer,states,base_splits,plan)


def audit_stage(cfg,plan,output,frontier):
    start=time.monotonic();C=cfg['counts'];P={k:Path(v) for k,v in cfg['paths'].items()}
    assert C['corners']==16*C['passes'] and C['passes']+C['pending']==C['children']
    assert C['pairs']==C['evaluated']+C['skipped']
    assert not output.exists();output.mkdir()
    write_new(output/'run_registration.json',{'stage_id':cfg['id'],'pid':os.getpid(),
        'maximum_new_corners':C['corners'],'maximum_pass_certificates':C['passes'],
        'metadata_rows':C['pairs'],'new_H_log_gradient_kernel_calls':0,'failed_corner_tests_repeated':0})
    children=load(P['children']);pairs=list(lines(P['pairs']));specplan=load(P['shard_plan'])
    master=load(P['master_state']);summary=load(P['master_summary'])
    assert len(children)==C["children"] and len(pairs)==C["pairs"] and len(specplan['shards'])==cfg['workers']
    assert master['status']==summary['status']==cfg['master_status']
    assert master['faults'] is None and summary['master_fault'] is None and not master['new_splits']
    accepted105=load(P['accepted_g105'])
    assert accepted105['status']=='primary_accepts_complete_same_model_global_g105'
    assert accepted105['qualified_g']=='21/200000'
    assert accepted105['full098_state_sha256']==plan['all_input_locks'][str(P['state098'])]
    R=canonical_roots();needed={p['source_id'] for p in pairs};sources={};poolids=set()
    assert len(needed)==C["selected_sources"]
    pool_refs=load(P['source_pool_refs'])
    assert pool_refs['whole_library']==cfg['library_count']
    pool_count=0
    for pool in cfg['pools']:
        path=Path(pool['path']);ref=pool_refs[pool['ref_key']]
        assert ref['path']==str(path) and ref['sha256']==plan['all_input_locks'][str(path)]
        assert ref['objects']==pool['objects'] and pool['index_start']==pool_count
        block_count=0
        for row in lines(path):
            sid=row['source_id'];assert sid==digest(row['core']) and sid not in poolids
            assert row['source_index']==pool_count;pool_count+=1;block_count+=1;poolids.add(sid)
            if sid in needed:sources[sid]=row
        assert block_count==pool['objects']
    assert pool_count==len(poolids)==cfg['library_count'] and set(sources)==needed
    pair_cursor=0
    for i,child in enumerate(children):
        assert child['child_index']==i;check_child(child,R)
        block=pairs[pair_cursor:pair_cursor+child['pair_count']];pair_cursor+=child['pair_count']
        assert len(block)==child['pair_count'] and 0<len(block)<=32
        assert [p['source_id'] for p in block]==child['selected_sources']
        assert child['selected_sources'][0]==child['parent_source_id']
        for pair in block:
            assert pair['child_index']==i and pair['child_node']==child['node']
            assert sources[pair['source_id']]['source_index']==pair['source_index']
            assert pair['curvature_rule']==sources[pair['source_id']]['core']['curvature_rule']
            pool=next(p for p in cfg['pools'] if p['index_start']<=pair['source_index']<p['index_start']+p['objects'])
            assert pair['source_pool']==pool['pair_label']
    assert pair_cursor==C["pairs"] and [p['pair_index'] for p in pairs]==list(range(C["pairs"]))
    # Read all fixed journals once; do not recalculate failed margins.
    qualified={};best={};events=set();count=Counter();durable=[]
    for spec in specplan['shards']:
        assert plan['all_input_locks'][spec['input_path']]==spec['input_sha256']
        sid=spec['shard_id'];data=load(spec['input_path']);folder=P['run_dir']/f'shard{sid:02d}'
        cp=load(folder/'checkpoint.json');assert cp['status']=='complete'
        fingerprint=cp['fingerprint']
        assert fingerprint['shard_id']==sid and fingerprint['input_sha256']==spec['input_sha256']
        assert fingerprint['core_sha256']==specplan['core_sha256']
        assert fingerprint['core_sha256']==plan['all_input_locks'][str(P['core131'])]
        assert fingerprint['worker_code_sha256']==plan['all_input_locks'][str(P['worker137'])]
        assert fingerprint['source_pool_sha256']==plan['all_input_locks'][str(P['source_pool_refs'])]
        assert fingerprint['child_inputs_sha256']==plan['all_input_locks'][str(P['children'])]
        assert fingerprint['candidate_pairs_sha256']==plan['all_input_locks'][str(P['pairs'])]
        assert fingerprint['prepared_manifest_sha256']==master['prepared_manifest_sha256']
        assert len(data['pairs'])==spec['pair_count'] and len(data['children'])==spec['child_count']
        for k,c in data['children'].items():assert c==children[int(k)] and int(k)%cfg['workers']==sid
        for k,s in data['sources'].items():assert s==sources[k]
        previous='0'*64;offset=ev=sk=0;first={}
        with (folder/'results.jsonl').open('rb') as handle:
            for cursor,raw in enumerate(handle):
                assert raw.endswith(b'\n');offset+=len(raw);row=json.loads(raw);h=row.pop('hash')
                assert row['seq']==cursor+1 and row['previous']==previous and digest(row)==h
                previous=h;pair=data['pairs'][cursor];gi=row['global_pair_index'];ci=row['child_index']
                assert pair==pairs[gi] and gi not in events and ci%cfg['workers']==sid;events.add(gi)
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
    assert events==set(range(C["pairs"])) and count=={'evaluated':C["evaluated"],'skipped':C["skipped"]}
    assert qualified==master['qualified_children'] and len(qualified)==C["passes"]
    pending=[c for c in children if str(c['child_index']) not in qualified]
    assert [c['node'] for c in pending]==master['pending_children'] and len(pending)==C["pending"]
    # Bind each selected source to its saved ORIGINAL certificate and raw record.
    aliases=defaultdict(list);jet_uses=defaultdict(list);expected_aliases=set()
    providers={p['origin']:p for p in cfg['jet_providers']};provider_for_source={}
    for source in sources.values():
        new_q_check(source)
        for alias in source['aliases']:
            expected_aliases.add(digest(alias))
            if alias['origin'] in providers:
                jet_uses[alias['origin']].append((source,alias))
                provider_for_source[source['source_id']]=providers[alias['origin']]
            else:aliases[(alias['origin'],alias['index'])].append((source,alias))
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
        original_records[source['source_id']]=True;matched.add(digest(a))
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
    jet_states={}
    for provider in cfg['jet_providers']:
        path=Path(provider['state_path']);state=load(path);jet_states[provider['id']]=state
        assert state['status']=='finite_original_prefix_complete' and not state['faults'] and not state['inflight'] and not state['obstructions_H128']
        E=provider['state_counts']
        assert len(state['new_passed'])==E['passed'] and len(state['new_splits'])==E['splits'] and len(state['queue'])==E['pending']
        cache=state['new_cache_by_epsilon']['1/128'];assert len(cache)==E['records']
        proofs={r['certificate']['point_key']:(r['node'],r['certificate']) for r in state['new_passed']}
        for row in state['new_splits']:
            k=row['failed_certificate']['point_key'];assert k not in proofs
            proofs[k]=(row['node'],row['failed_certificate'])
        assert set(proofs)==set(cache)
        for source,alias in jet_uses[provider['origin']]:
            point=','.join(source['core']['point']);record=cache[point];node,cert=proofs[point]
            assert record['epsilon']=='1/128' and record['qualification_target']=='11/100000'
            assert record['pair_cache_key']==['1/128',point]
            match(source,alias,node,cert,record)
            ancestry=source['H64_g105_ancestry']
            assert ancestry['proof']==cert[provider['certificate_proof_field']]
            assert ancestry[provider['ancestry_state_sha_field']]==plan['all_input_locks'][str(path)]
            assert ancestry['node']==cert['H64_lower_source_node']
            assert in_box(source['core']['point'],node,R)
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
                assert row['target']=='21/200000'
                if 'target' in row['certificate']:
                    assert row['certificate']['target']==row['target']
                assert row['certificate']['passed'] and f(row['certificate']['minimum_log_margin'])[0]>0
                if 'source_leaf_sha256' in a:
                    assert digest(row)==a['source_leaf_sha256']
                    assert a['source_state_sha256']==plan['all_input_locks'][str(P['state098'])]
                else:
                    provider=provider_for_source[source['source_id']]
                    assert all(x['origin']==provider['origin'] for x in source['aliases'])
                    proof=a['proof']
                    assert proof['node']==row['node'] and proof['lower']=='21/200000'
                    assert proof['child_is_descendant'] and proof['newpoint_inside_child']
                    assert a[provider['ancestry_state_sha_field']]==plan['all_input_locks'][provider['state_path']]
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
    assert set(original_records)==needed and matched==expected_aliases
    for source in sources.values():
        c=source['core'];a=source['H64_g105_ancestry']
        if a:assert source['source_id'] in ancestor_matches
        elif c['epsilon']=='1/128':
            assert all(x['origin']=='105_passed' for x in source['aliases'])
            for x in source['aliases']:
                rid,path=key(x['source_node']);assert not any((rid,path[:j]) in h64nodes for j in range(len(path)+1))
    # Exact residual and gap order; failed margin values are read, never recomputed.
    ordered=list(lines(P['residual']));assert len(ordered)==C["pending"]
    residual_ids=set()
    def compare(a,b):
        x,y=a[cfg['residual_certificate_field']]['minimum_log_margin'],b[cfg['residual_certificate_field']]['minimum_log_margin']
        if not eq(x,y):return -1 if le(x,y) else 1
        da,db=len(a['node'][1]),len(b['node'][1])
        if da!=db:return -1 if da<db else 1
        return (key(a['node'])>key(b['node']))-(key(a['node'])<key(b['node']))
    for i,row in enumerate(ordered):
        ci=row['child_index'];assert ci not in residual_ids;residual_ids.add(ci)
        assert str(ci) not in qualified and row['gap_rank']==i
        assert row[cfg['residual_certificate_field']]==best[ci] and not best[ci]['passed']
        assert row['node']==children[ci]['node'] and row['box_lower']==children[ci]['box_lower'] and row['box_upper']==children[ci]['box_upper']
        for field,value in children[ci].items():assert row[field]==value
        assert eq(row['physical_volume'],f'1/{32**4*(1<<len(row["node"][1]))}')
        if i:assert compare(ordered[i-1],row)<=0
    assert residual_ids=={c['child_index'] for c in pending}
    # Explicit initial layers or the previous named stage's returned frontier.
    if frontier is None:
        passed=[];splits=[]
        for layer in cfg['initial_frontier']['passed_layers']:
            passed.extend(nodes_from_layer(layer,states,base_splits,plan))
        for layer in cfg['initial_frontier']['split_layers']:
            splits.extend(split_rows_from_layer(layer,states,base_splits,plan))
        pending_before=nodes_from_layer(cfg['initial_frontier']['pending_layer'],states,base_splits,plan)
    else:passed,splits,pending_before=frontier
    B=cfg['frontier_before_jet']
    assert (len(passed),len(splits),len(pending_before))==(B['passed'],B['splits'],B['pending'])
    full_tree(passed+pending_before,[r['node'] for r in splits])
    phase=jet_states[cfg['current_jet_provider']]
    parents=[r['node'] for r in phase['new_passed']+phase['new_splits']]
    assert len(parents)==len(pending_before) and set(map(key,parents))==set(map(key,pending_before))
    assert len(set(map(key,phase['queue'])))==C['children']
    assert {key(c['node']) for c in children}==set(map(key,phase['queue']))
    split_phase={key(r['node']):r for r in phase['new_splits']}
    cache=phase['new_cache_by_epsilon']['1/128']
    current_provider=next(p for p in cfg['jet_providers'] if p['id']==cfg['current_jet_provider'])
    for child in children:
        parent=split_phase[key(child['direct_parent_node'])]
        assert child['node'] in parent['children']
        assert child['parent_failed_margin']==parent['failed_certificate']['minimum_log_margin']
        source=sources[child['parent_source_id']]
        assert all(a['origin']==current_provider['origin'] for a in source['aliases'])
        assert source['aliases'][0]['source_node']==parent['node']
        assert source['raw_record_sha256']==digest(cache[parent['failed_certificate']['point_key']])
    passed=passed+[r['node'] for r in phase['new_passed']]+[r['child_node'] for r in qualified.values()]
    splits=splits+phase['new_splits']
    for split in splits:
        rid,path=key(split['node']);assert split['axis']==len(path)%4
        assert split['children']==[[rid,path+'0'],[rid,path+'1']]
    assert len(passed)==C['passed'] and len(splits)==C['splits']
    tree=full_tree(passed+[r['node'] for r in pending],[r['node'] for r in splits])
    pass_eps=Counter(c['epsilon'] for c in qualified.values())
    outside=sum(not c['anchor_inside_child'] for c in qualified.values())
    assert pass_eps==C['pass_eps'] and outside==C['outside']
    # Config/source/protocol/tree gates precede every fresh corner calculation.
    with (output/'winning_integer_checks.jsonl').open('x') as out:
        for ci,cert in sorted(qualified.items(),key=lambda z:int(z[0])):
            source=sources[cert['source_id']];answer=positive_certificate(source,children[int(ci)],cert,R)
            out.write(json.dumps({'child_index':int(ci),'node':cert['child_node'],'source_id':cert['source_id'],**answer},sort_keys=True,separators=(',',':'))+'\n')
    for path,h in plan['all_input_locks'].items():assert sha(path)==h,path
    result={'task':'audit_later145','stage_id':cfg['id'],
        'status':'finite_stage_passes_sources_protocol_fulltree_qualified_partial',
        'pid':os.getpid(),'winning_certificates':C['passes'],'new_integer_corners':C['corners'],
        'journal_rows':C['pairs'],'evaluated_receipts':C['evaluated'],'skip_receipts':C['skipped'],
        'selected_sources_bound':len(sources),'selected_new_sources_by_origin':{k:len(v) for k,v in jet_uses.items()},
        'pass_epsilon_counts':dict(pass_eps),'outside_anchor_passes':outside,'tree':tree,
        'passed':C['passed'],'pending':C['pending'],'durable_shards':durable,
        'winning_integer_checks_sha256':sha(output/'winning_integer_checks.jsonl'),
        'all_source_hashes_unchanged':True,'globalg110_qualified':False,'S_credit':0,
        'new_H_log_gradient_kernel_calls':0,'failed_corner_tests_repeated':0,'new_splits':0,
        'wall_seconds':time.monotonic()-start}
    write_new(output/'stage_result.json',result)
    return result,(passed,splits,[r['node'] for r in pending])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True)
    ap.add_argument('--reviewed-plan-sha',required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();assert sha(args.plan)==args.reviewed_plan_sha;plan=load(args.plan)
    assert plan['execution_ready'] is True,'Missing frozen numerical qualification/configuration'
    assert args.output.is_relative_to(HERE) and not args.output.exists()
    for name,h in plan['code_locks'].items():assert sha(HERE/name)==h
    for path,h in plan['all_input_locks'].items():assert sha(path)==h,path
    stage_file=Path(plan['stage_list']['path']);assert sha(stage_file)==plan['stage_list']['sha256']
    config=load(stage_file);stages=config['stages'];ids=[s['id'] for s in stages]
    assert len(ids)==len(set(ids)) and not set(ids)&set(config['already_qualified_stage_ids'])
    assert all(i and not set(i)-set('abcdefghijklmnopqrstuvwxyz0123456789_-') for i in ids)
    assert sum(s['counts']['corners'] for s in stages)==plan['maximum_new_corners']
    assert sum(s['counts']['passes'] for s in stages)==plan['maximum_new_pass_certificates']
    for q in plan['inherited_qualifications']:
        assert plan['all_input_locks'][q['path']]==q['sha256']
        m=load(q['path']);assert m['status']==q['status']
        for k,v in q['required_fields'].items():assert m[k]==v
    # Every new raw-jet provider needs a named FINAL, hash-locked numerical receipt.
    for stage in stages:
        for provider in stage['jet_providers']:
            q=provider['numeric_qualification'];assert q['ready'] is True
            assert plan['all_input_locks'][q['manifest_path']]==q['manifest_sha256']
            m=load(q['manifest_path']);assert m['status']==q['status']
            for k,v in q['required_fields'].items():assert m[k]==v
            lock=load(q['source_lock_path']);assert m['source_lock_sha256']==sha(q['source_lock_path'])
            assert lock['input_locks'][provider['state_path']]==plan['all_input_locks'][provider['state_path']]
        for p in stage['paths'].values():
            if p!=stage['paths']['run_dir']:assert p in plan['all_input_locks']
    args.output.mkdir();write_new(args.output/'execution_registration.json',{
        'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),
        'reviewed_plan_sha256':args.reviewed_plan_sha,'stage_ids':ids,
        'maximum_new_corners':plan['maximum_new_corners'],'new_H_log_gradient_kernel_calls':0})
    frontier=None;results=[];previous=None
    for stage in stages:
        assert stage['previous_stage_id']==previous
        result,frontier=audit_stage(stage,plan,args.output/stage['id'],frontier)
        results.append(result);previous=stage['id']
    for path,h in plan['all_input_locks'].items():assert sha(path)==h,path
    write_new(args.output/'execution_result.json',{'status':'explicit_finite_stage_list_complete_partial',
        'stage_ids':ids,'new_corners':sum(r['new_integer_corners'] for r in results),
        'new_passes':sum(r['winning_certificates'] for r in results),'globalg110_qualified':False,
        'S_credit':0,'new_H_log_gradient_kernel_calls':0,'failed_corner_tests_repeated':0})


if __name__=='__main__':main()
