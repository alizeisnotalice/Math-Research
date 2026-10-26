#!/usr/bin/env python3
"""Finite migration audit; no production functions and no new H/log calls."""
import sys
sys.dont_write_bytecode=True
import json,gzip,os,time
from pathlib import Path
from fractions import Fraction as F
from math import lcm
from interfaces import REG,HERE,common as C


def write(path,obj):
    temporary=Path(str(path)+'.tmp');temporary.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n');temporary.replace(path)


def choose(node,roots,vertices):
    lo,hi=C.physical_box(roots,node);mid=[(a+b)/2 for a,b in zip(lo,hi)];eta=min(b-a for a,b in zip(lo,hi))/4096
    candidates=[mid]+[[x+eta*beta**j for j,x in enumerate(mid)] for beta in (F(1,8),F(1,16),F(1,32))];tests=[]
    for i,tail in enumerate(candidates):
        point=[F(1),*tail];D=lcm(*(x.denominator for x in point));nums=[int(x*D) for x in point]
        q=[sum(a*b for a,b in zip(nums,vertex)) for vertex in vertices]
        test={'candidate':i,'common_denominator':str(D),'distinct_count':len(set(q)),'zero_count':q.count(0),'eligible':len(set(q))==280 and 0 not in q};tests.append(test)
        if test['eligible']:return list(map(str,point)),i,tests
    raise AssertionError('no generic anchor')


def wrapped(row):
    return {**row,'leaf_epsilon':'1/64','lower_bound_object':'G','component_lower_function':'H_epsilon','target':'21/200000'}


def main():
    began=time.monotonic();rpath=HERE/'source_gate_registration.json';registration=json.loads(rpath.read_text())
    for p,h in {**REG['protected_old_locks'],**REG['preparation_input_locks'],**registration['code_locks']}.items():C.require(C.file_sha(p)==h,'frozen input changed: '+p)
    old=json.loads(Path(REG['source089']).read_text());C.require(old['status']=='counterexample' and old['faults']==[] and len(old['queue'])==190 and len(old['inflight'])==4,'089 stop state differs')
    expected_cache={};accepted=[];mapper=[];split_hashes=[];split_nodes=[];roots=None
    stream=C.Stream(REG['source086'])
    while True:
        field=stream.key()
        if field is None:break
        if field=='cache':
            for key,value in stream.items(True):C.require(key==','.join(value['point']),'source point key differs');expected_cache[key]=C.digest(value)
        elif field=='passed':
            for i,row in enumerate(stream.items()):
                cert=row['certificate'];C.require(cert['passed'] and F(cert['minimum_log_margin'])>0 and cert['epsilon']=='1/64' and cert['target']=='21/200000','source accepted qualification metadata differs')
                accepted.append((C.node(row['node']),C.digest(wrapped(row))))
                mapper.append({'certificate_sha256':C.digest(cert),'leaf_epsilon':'1/64','new_passed_index':len(accepted)-1,'node':row['node'],'source':'source086_passed','source_index':i,'source_leaf_sha256':C.digest(row)})
        elif field=='splits':
            for row in stream.items():split_hashes.append(C.digest(row));split_nodes.append(C.node(row['node']))
        elif field=='roots':roots=stream.value()
        else:stream.value()
    stream.close();C.require(len(expected_cache)==180040 and len(accepted)==96104 and len(split_nodes)==45600,'source finite counts differ')
    for key,row in old['new_cache'].items():C.require(key not in expected_cache,'old cache overlap');expected_cache[key]=C.digest(row)
    C.require(len(expected_cache)==184792,'positive source cache count differs')
    for i,row in enumerate(old['new_passed']):
        cert=row['certificate'];C.require(cert['passed'] and F(cert['minimum_log_margin'])>0,'089 accepted record not strict')
        accepted.append((C.node(row['node']),C.digest(wrapped(row))))
        mapper.append({'certificate_sha256':C.digest(cert),'leaf_epsilon':'1/64','new_passed_index':len(accepted)-1,'node':row['node'],'source':'source089_new_passed','source_index':i,'source_leaf_sha256':C.digest(row)})
    for row in old['new_splits']:split_hashes.append(C.digest(row));split_nodes.append(C.node(row['node']))
    pending=[]
    for i,node in enumerate(old['queue']):pending.append({'node':node,'source089_queue_index':i,'source_kind':'queue','leaf_epsilon_for_next_evaluation':'1/128'})
    counter_keys=set()
    for ordinal,entry in old['inflight'].items():
        got=entry['returned'];row=got['record'];key=','.join(row['point'])
        C.require(got['kind']=='counterexample' and 0<F(row['values'][0][0])<=F(row['values'][0][1])<F(21,200000),'H64 obstruction lost')
        C.require(key not in expected_cache,'counter already positive cached');expected_cache[key]=C.digest(row);counter_keys.add(key)
        pending.append({'node':entry['node'],'source089_inflight_ordinal':ordinal,'source_kind':'retained_counter_inflight','leaf_epsilon_for_next_evaluation':'1/128','H64_record_sha256':C.digest(row),'H64_point':row['point'],'H64_interval':row['values'][0]})
    pending.sort(key=lambda x:(0 if x['source_kind']=='retained_counter_inflight' else 1,F(x['H64_interval'][1]) if x['source_kind']=='retained_counter_inflight' else F(old['hardness'][str(x['node'][0])+':'+x['node'][1]][0]),*x['node']))
    for i,row in enumerate(pending):
        lo,hi=C.physical_box(roots,C.node(row['node']));row.update(new_queue_index=i,box_lower=list(map(str,lo)),box_upper=list(map(str,hi)))
    root=Path(REG['source_initial']).parent
    C.require(mapper==json.loads((root/'accepted_H64_mapper.json').read_text()),'accepted mapper not exact source identity')
    C.require(pending==json.loads((root/'pending_transfer_mapper.json').read_text()),'194 migration mapper differs')
    stream=C.Stream(REG['source_initial']);metadata={};seen_cache=set();seen_passed=[];seen_splits=[];queue=None
    while True:
        field=stream.key()
        if field is None:break
        if field=='cache_by_epsilon':
            stream.char('{');C.require(stream.value()=='1/128','wrong H128 namespace');stream.char(':');C.require(stream.value()=={},'initial H128 bucket must be empty')
            stream.char(',');C.require(stream.value()=='1/64','wrong H64 namespace');stream.char(':')
            for key,row in stream.items(True):
                C.require(key not in seen_cache and key in expected_cache and expected_cache[key]==C.digest(row),'H64 typed bucket changed source record')
                C.require(key==','.join(row['point']),'typed cache point mismatch');seen_cache.add(key)
            stream.char('}')
        elif field=='passed':
            for i,row in enumerate(stream.items()):
                C.require(i<len(accepted) and (C.node(row['node']),C.digest(row))==accepted[i],'retained H64 leaf changed or epsilon relabelled');seen_passed.append(C.node(row['node']))
        elif field=='splits':
            for i,row in enumerate(stream.items()):
                node=C.node(row['node']);C.require(i<len(split_hashes) and C.digest(row)==split_hashes[i],'split history changed')
                C.require(row['axis']==len(node[1])%4 and row['children']==[[node[0],node[1]+'0'],[node[0],node[1]+'1']],'closed child geometry differs');seen_splits.append(node)
        elif field=='roots':C.require(stream.value()==roots,'root geometry changed')
        elif field=='queue':queue=stream.value()
        else:metadata[field]=stream.value()
    stream.close()
    C.require(seen_cache==set(expected_cache) and len(seen_cache)==184796,'H64 cache omission/addition')
    C.require(len(seen_passed)==100804 and len(seen_splits)==48638 and queue==[x['node'] for x in pending] and len(queue)==194,'mixed frontier counts/queue differ')
    C.require(metadata['target']=='21/200000' and metadata['pending_evaluation_epsilon']=='1/128' and metadata['lower_bound_object']=='G' and metadata['inflight']==metadata['resume_jobs']=={} and metadata['new_H128_records']==0 and metadata['g105_closed_cover'] is False,'mixed object/epsilon/status differs')
    fingerprint=C.digest({'preparation_registration_sha256':C.file_sha(root/'preparation_registration.json'),'accepted_mapper_sha256':C.file_sha(root/'accepted_H64_mapper.json'),'pending_mapper_sha256':C.file_sha(root/'pending_transfer_mapper.json'),'source086_sha256':REG['source086_sha256'],'source089_stop_sha256':REG['source089_sha256']})
    C.require(metadata['fingerprint098']==fingerprint,'migration fingerprint differs')
    tree=C.validate_tree(roots,seen_passed+list(map(C.node,queue)),set(seen_splits),'mixed098initial')
    fan=json.loads(gzip.decompress(Path(REG['fan']).read_bytes()));requests=json.loads((root/'planned_H128_requests.json').read_text());tasks=json.loads((root/'execution_tasks_194.json').read_text())
    C.require(len(requests)==len(tasks)==194,'request count differs');hit64=0
    for i,(request,task,node) in enumerate(zip(requests,tasks,queue)):
        point,candidate,tests=choose(C.node(node),roots,fan['vertices']);key=','.join(point);same=key in expected_cache
        expected={'ordinal':i,'node':node,'epsilon':'1/128','point':point,'pair_cache_key':['1/128',key],'candidate':candidate,'geometry_tests':tests,'H128_cache_hit':False,'H64_same_point_present':same,'target':'21/200000'}
        C.require(request==expected,'planned epsilon-key/anchor differs');hit64+=same
        if same:
            C.require(key in counter_keys,'unexpected H64 hit');row=next(x['returned']['record'] for x in old['inflight'].values() if ','.join(x['point'])==key);expected['old_H64_interval']=row['values'][0]
        expected.update(point_checkpoint=str(root/'batch0001'/('point%04d.json'%i)),stop_file=str(root/'STOP'))
        C.require(task==expected,'execution task differs from typed request')
    C.require(hit64==4 and len({tuple(x['pair_cache_key']) for x in requests})==194,'epsilon cache isolation differs')
    for p,h in {**REG['protected_old_locks'],**REG['preparation_input_locks'],**registration['code_locks']}.items():C.require(C.file_sha(p)==h,'input changed during finite audit')
    result={'task':'audit_mixed100_source_gate','status':'mixed098_initial_source_migration_typed_cache_and_tree_qualified_partial','pid':os.getpid(),'wall_seconds':time.monotonic()-began,'registration_sha256':C.file_sha(rpath),'accepted_H64_leaves':100804,'retained_H64_cache_records':184796,'empty_H128_cache':True,'H64_counters_retained':4,'pending_from_queue':190,'pending_from_counter':4,'pending_H128':194,'all194_tasks_independently_reconstructed':True,'H64_same_point_observations_not_reused':4,'fingerprint098':fingerprint,'tree':tree,'old_protected_files_unchanged':len(REG['protected_old_locks']),'new_H_or_gradient_or_log_calls':0,'global_g105_qualified':False}
    path=HERE/'source_gate_result.json';write(path,result);print(json.dumps(result),flush=True);print(json.dumps({'result_sha256':C.file_sha(path)}),flush=True)


if __name__=='__main__':main()
