#!/usr/bin/env python3
"""Frozen107 first4096 source/guard/replay/jet/corner and partial-tree audit."""
import sys
sys.dont_write_bytecode=True
import json,gzip,os,time
from pathlib import Path
from fractions import Fraction as F
from math import lcm
from interfaces import REG,HERE,common as C
from certificate_check import certificate


def write(path,obj):
    p=Path(path);q=Path(str(p)+'.tmp');q.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n');q.replace(p)


def geometry(point,vertices):
    p=list(map(F,point));D=lcm(*(x.denominator for x in p));nums=[int(x*D) for x in p];values=[sum(x*y for x,y in zip(v,nums)) for v in vertices]
    return {'common_denominator':str(D),'distinct_count':len(set(values)),'eligible':len(set(values))==280 and 0 not in values,'zero_count':values.count(0)}


def source_gate():
    needed={row['pair_cache_key'][1] for row in jsonlines(Path(REG['source105']).parent/'pending_gap_order.jsonl')};cache={};old128=set();roots=None;splits=[];reader=C.Stream(REG['source098'])
    while True:
        field=reader.key()
        if field is None:break
        if field=='cache_by_epsilon':
            reader.char('{')
            for i,epsilon in enumerate(('1/128','1/64')):
                if i:reader.char(',')
                C.require(reader.value()==epsilon,'source typed epsilon differs');reader.char(':')
                for key,record in reader.items(True):
                    if epsilon=='1/128':old128.add(key)
                    elif key in needed:cache[key]={'sha256':C.digest(record),'H_interval':record['values'][0]}
            reader.char('}')
        elif field=='roots':roots=reader.value()
        elif field=='splits':
            for row in reader.items():
                n=C.node(row['node']);C.require(row['axis']==len(n[1])%4 and row['children']==[[n[0],n[1]+'0'],[n[0],n[1]+'1']],'source split geometry differs');splits.append(n)
        elif field=='passed':
            for _ in reader.items():pass
        else:reader.value()
    reader.close();C.require(len(old128)==194 and len(cache)==4973 and len(splits)==48638,'source typed records/tree counts differ')
    reader=C.Stream(REG['source105']);accepted=[];metadata={};old194=0
    while True:
        field=reader.key()
        if field is None:break
        if field=='passed':
            for row in reader.items():
                C.require(row['certificate']['passed'] is True and row['certificate']['target']=='11/100000','unqualified source newtarget leaf')
                accepted.append(C.node(row['node']));old194+=row['leaf_epsilon']=='1/128'
        else:metadata[field]=reader.value()
    reader.close();C.require(len(accepted)==96025 and old194==194 and len(metadata['queue'])==4973 and not metadata['inflight'] and not metadata['resume_jobs'],'source105 qualification/queue differs')
    tree=C.validate_tree(roots,accepted+list(map(C.node,metadata['queue'])),set(splits),'source105NewtargetPartial')
    return cache,old128,roots,splits,accepted,metadata,tree


def jsonlines(path):
    with Path(path).open() as f:
        for line in f:yield json.loads(line,object_pairs_hook=C.no_duplicates)


def requests_gate(cache,old128,roots,source):
    root=Path(REG['batch_directory']).parent;pending=list(jsonlines(Path(REG['source105']).parent/'pending_gap_order.jsonl'));planned=json.loads((root/'planned_requests.json').read_text());tasks=json.loads((root/'production_tasks_first4096.json').read_text())
    fan=json.loads(gzip.decompress(Path(REG['fan']).read_bytes()));reserved={}
    for label in ('095','103'):
        path=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-alt-cont/cube-dim-5-bounds')/('upper_adaptive'+label)/'jobs.json.gz'
        reserved['upper_adaptive'+label]={','.join(x['point']) for x in json.loads(gzip.decompress(path.read_bytes()))}
    seen=set();changed=0;guards=0
    C.require(len(planned)==len(pending)==4973 and source['queue']==[r['node'] for r in planned] and len(tasks)==4096,'planned set is not exact source pending set')
    for i,(row,request) in enumerate(zip(pending,planned)):
        n=C.node(row['node']);lo,hi=C.physical_box(roots,n);C.require(row['box_lower']==list(map(str,lo)) and row['box_upper']==list(map(str,hi)),'source pending physical box differs')
        mid=[(a+b)/2 for a,b in zip(lo,hi)];eta=min(b-a for a,b in zip(lo,hi))/4096;candidates=[[F(1),*mid]]+[[F(1),*[x+eta*beta**j for j,x in enumerate(mid)]] for beta in (F(1,8),F(1,16),F(1,32))]
        oldpoint=list(map(F,row['pair_cache_key'][1].split(',')));oldindex=candidates.index(oldpoint);tests=[];picked=None
        for j in range(oldindex,4):
            point=list(map(str,candidates[j]));key=','.join(point);generic=geometry(point,fan['vertices']);collisions=[name for name,points in reserved.items() if key in points]
            test={'candidate':j,'point':point,'inside_actual_box':all(a<=x<=b for a,x,b in zip(lo,candidates[j][1:],hi)),'generic':generic,'planned_H128_collision_sources':collisions,'already_new107':key in seen,'completed_H128_cache_hit':key in old128};tests.append(test)
            if generic['eligible'] and not collisions and key not in seen:picked=(j,point,key);break
        C.require(picked is not None,'no independently eligible fresh candidate');j,point,key=picked;seen.add(key)
        expected={'H128_cache_hit':key in old128,'H64_source_failed_g110_certificate':row['certificate'],'H64_source_record_sha256':row['source_record_sha256'],'box_lower':row['box_lower'],'box_upper':row['box_upper'],'candidate':j,'candidate_tests':tests,'epsilon':'1/128','node':row['node'],'ordinal':i,'original_H64_candidate':oldindex,'original_H64_point':list(map(str,oldpoint)),'pair_cache_key':['1/128',key],'point':point,
            'selection_reason':'avoid already registered H128 mathematical object,not evidence original anchor failed H128' if j!=oldindex else 'same physical anchor but new epsilon is a distinct mathematical object','source_gap_rank':i,'source_leaf_index':row['source_leaf_index'],'target':'11/100000'}
        C.require(request==expected and expected['H128_cache_hit'] is False,'planned source/typed candidate/reason differs')
        oldkey=','.join(map(str,oldpoint));C.require(cache[oldkey]['sha256']==row['source_record_sha256'],'source H64 guard record differs')
        if i<4096:
            same=point==list(map(str,oldpoint));guard={'epsilon':'1/64','point':point,'H_interval':cache[oldkey]['H_interval']} if same else None
            expected.update(geometry_tests=[{'candidate':x['candidate'],**x['generic']} for x in tests],same_point_H64_guard=guard,point_checkpoint=str(root/'batch0001'/('point%04d.json'%i)),stop_file=str(root/'STOP'))
            C.require(tasks[i]==expected,'execution task/same-point guard differs');changed+=not same;guards+=same
    C.require(changed==125 and guards==3971 and len(seen)==4973 and not seen.intersection(old128) and all(not seen.intersection(x) for x in reserved.values()),'changed-point guard/dedup classification differs')
    return planned,tasks,{'changed_first4096':changed,'same_point_guards':guards,'all4973_typed_requests_unique':True,'old194_H128_intersection':0,'frozen095_103_intersections':0,'reserved_jobs_counts':{name:len(points) for name,points in reserved.items()}}


def numeric_gate(pm):
    root=HERE/'numeric_batch0001';index=json.loads((root/'index.json').read_text());C.require(index['status']=='all4096_H128_plus4grad_strictly_verified' and index['verified_points']==index['global_verified_prefix']==4096,'new numeric gate incomplete')
    for p,h in index['input_locks'].items():C.require(C.file_sha(p)==h,'numeric inputs changed')
    count=0
    for plan in index['plans']:
        state=json.loads((root/('block%04d.json'%plan['block'])).read_text());C.require(state['status']=='complete' and len(state['points'])==plan['stop_exclusive']-plan['start'],'numeric prefix incomplete')
        C.require(state['input_locks']==index['input_locks'] and C.file_sha(root/('block%04d.json'%plan['block']))==index['block_progress'][plan['block']]['sha256'],'numeric durable block/code lock changed')
        for proof in state['points']:
            entry=pm['records'][proof['index']];C.require(proof['point']==entry['point'] and proof['raw_checkpoint_sha256']==entry['point_checkpoint_sha256'] and proof['H_multiplier']==16 and proof['gradient_multiplier']==32 and proof['independent_EPS64_input']==[str(2*F(x)) for x in entry['point']],'scaled numeric point association differs')
            C.require(len(proof['independent_H128_values_512_dyadic'])==5,'independent jet missing components')
            for ball,production in zip(proof['independent_H128_values_512_dyadic'],entry['record']['values']):
                a,b=map(F,ball);L,U=map(F,production);C.require(L<a<=b<U,'independent scaled ball outside production')
            count+=1
    C.require(count==4096 and all(x['exitcode']==0 and x['is_alive'] is False for x in index['pool_exits']),'numeric finite count/pool differs')
    return index


def main():
    began=time.monotonic();rp=HERE/'structure_registration.json';reg=json.loads(rp.read_text())
    for p,h in {**REG['input_locks'],**REG['protected_old_locks'],**reg['input_locks']}.items():C.require(C.file_sha(p)==h,'frozen input changed: '+p)
    cache,old128,roots,splits,accepted,source,source_tree=source_gate();planned,tasks,request_result=requests_gate(cache,old128,roots,source)
    b=Path(REG['batch_directory']);initial=json.loads((b/'initial_overlay.json').read_text());final=json.loads((b/'state.json').read_text());br=json.loads((b/'registration.json').read_text());pm=json.loads((b/'point_manifest.json').read_text());numeric=numeric_gate(pm)
    C.require(initial['queue']==source['queue'] and initial['target']=='11/100000' and initial['source098_sha256']==REG['source098_sha256'] and initial['source105_sha256']==REG['source105_sha256'],'initial immutable source/target/queue differs')
    C.require(initial['new_cache_by_epsilon']=={'1/128':{}} and initial['new_passed']==initial['new_splits']==initial['faults']==initial['obstructions_H128']==[] and not initial['inflight'],'initial overlay contains values or retired leaves')
    expected=json.loads(json.dumps(initial));queue=list(initial['queue']);dispatch={};received={};seq=0;previous='0'*64;merge_count=0;certs=[];calls=0
    for event in jsonlines(b/'journal.jsonl'):
        body=event['body'];kind=event['kind'];C.require(event['seq']==seq+1 and event['previous']==previous and C.digest({k:v for k,v in event.items() if k!='hash'})==event['hash'],'journal hash chain differs');seq+=1;previous=event['hash']
        if kind=='dispatch':
            task=body['task'];i=task['ordinal'];C.require(i==len(dispatch) and task==tasks[i] and queue.pop(0)==task['node'],'dispatch is not exact first4096 finite prefix');dispatch[i]=task
        elif kind=='received':
            i=body['ordinal'];got=body['returned'];task=dispatch[i];entry=pm['records'][i];C.require(i not in received and got['kind']=='complete' and got['epsilon']=='1/128' and got['pair_cache_key']==task['pair_cache_key'],'fault/counter/duplicate/wrong epsilon returned')
            C.require(entry['ordinal']==i and entry['node']==task['node'] and entry['point']==task['point'] and entry['target']=='11/100000' and entry['stage']=='complete' and entry['point_checkpoint']==task['point_checkpoint'],'point manifest association differs')
            C.require(C.file_sha(entry['point_checkpoint'])==entry['point_checkpoint_sha256'],'raw point SHA differs');raw=json.loads(Path(entry['point_checkpoint']).read_text());record=got['record']
            C.require(raw['stage']=='complete' and raw['partial'] is None and raw['record']==entry['record']==record,'complete partial=None/raw record differs')
            C.require(raw['epsilon']==record['epsilon']=='1/128' and raw['target']==record['qualification_target']=='11/100000' and raw['pair_cache_key']==record['pair_cache_key']==task['pair_cache_key'],'raw/record target/epsilon/pair differs')
            C.require(record['point']==task['point'] and record['candidate']==task['candidate'] and record['geometry_tests']==task['geometry_tests'] and record['model']==pm['model'] and raw['kernel_sha256']==record['kernel_sha256']==pm['kernel_sha256'],'anchor/model/kernel differs')
            C.require(raw['pid']==entry['worker_pid']==got['pid'] and raw['calls']==got['calls']==1 and record['bits']==256 and raw['precision_index']==0 and record['attempts']==got['attempts']==raw['attempts'] and len(record['attempts'])==1 and record['attempts'][0]['all_cones']==3504 and record['attempts'][0]['numeric_quality'] is True,'jet completeness/precision/PID differs')
            C.require(F(record['values'][0][0])>F(11,100000),'source H128 anchor not above target');guard=task['same_point_H64_guard']
            if guard:C.require(guard['point']==record['point'] and F(record['values'][0][1])>=F(guard['H_interval'][0]),'wrong-point/monotonicity guard fault')
            received[i]=got;calls+=got['calls']
        elif kind=='merge':
            i=body['ordinal'];C.require(i==merge_count and i in received,'merge ordering differs');merge_count+=1;task=dispatch[i];record=received[i]['record'];cert=body['certificate'];n=C.node(task['node']);proof=certificate(n,cert,record,roots);certs.append(proof)
            key=','.join(record['point']);C.require(key not in old128 and key not in expected['new_cache_by_epsilon']['1/128'],'typed cache overwrote old/new mathematical object');expected['new_cache_by_epsilon']['1/128'][key]=record
            if cert['passed']:expected['new_passed'].append({'node':list(n),'leaf_epsilon':'1/128','lower_bound_object':'G','target':'11/100000','certificate':cert})
            else:
                kids=[[n[0],n[1]+'0'],[n[0],n[1]+'1']];expected['new_splits'].append({'node':list(n),'axis':len(n[1])%4,'children':kids,'leaf_epsilon':'1/128','target':'11/100000','failed_certificate':cert});queue.extend(kids)
        elif kind=='finish':
            pool=body['pool_evidence'];C.require(body['status']=='first4096_complete' and len(dispatch)==len(received)==merge_count==4096 and len(pool)==4 and all(x['exitcode']==0 and x['alive'] is False for x in pool),'finite pool/finish differs');C.require({x['pid'] for x in received.values()}<=set(x['pid'] for x in pool),'returned PID outside pool');expected.update(pool_evidence=pool,status=body['status'],wall_seconds=body['wall_seconds'])
        else:raise AssertionError('unknown journal event')
    expected.update(queue=queue,next_merge=merge_count,journal_seq=seq,journal_hash=previous,counts={'cache_hits':0,'kernel_segment_calls':calls,'new_anchor_starts':4096,'new_complete':4096})
    C.require(expected==final and seq==12289 and len(final['new_passed'])==850 and len(final['new_splits'])==3246 and len(queue)==7369 and final['inflight']=={} and final['faults']==final['obstructions_H128']==[],'exact overlay or finite counts differ')
    C.require(queue[:877]==source['queue'][4096:],'untested original877 not retained exactly')
    tree=C.validate_tree(roots,accepted+[C.node(x['node']) for x in final['new_passed']]+list(map(C.node,queue)),set(splits)|{C.node(x['node']) for x in final['new_splits']},'107First4096PartialFinal')
    mapper=[];root=b.parent
    for i,request in enumerate(planned[4096:]):mapper.append({'H128_evaluation_performed':False,'frozen_planned_point':request['point'],'kind':'untested_original_H128_request','node':request['node'],'pending_evaluation_epsilon':'1/128','pending_index':i,'source105_leaf_index':request['source_leaf_index'],'source_H64_failed_certificate_sha256':C.digest(request['H64_source_failed_g110_certificate']),'source_gap_rank':request['source_gap_rank'],'target':'11/100000'})
    inside_count=0
    for split in final['new_splits']:
        parent=split['node'];cert=split['failed_certificate'];key=cert['point_key'];record=final['new_cache_by_epsilon']['1/128'][key]
        for child in split['children']:
            lo,hi=C.physical_box(roots,C.node(child));inside=all(a<=F(x)<=b for a,x,b in zip(lo,record['point'][1:],hi));inside_count+=inside
            mapper.append({'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi)),'child_certificate_evaluated':False,'kind':'closed_child_of_failed_H128_parent','new_child_H_evaluation_performed':False,'node':child,'parent':parent,'parent_H128_interval':record['values'][0],'parent_anchor_inside_child':inside,'parent_failed_certificate_sha256':C.digest(cert),'parent_record_sha256':C.digest(record),'pending_evaluation_epsilon':'1/128','pending_index':len(mapper),'stored_parent_jet_pair_cache_key':['1/128',key],'target':'11/100000'})
    C.require(mapper==list(jsonlines(root/'pending_mapper_batch0001.jsonl')) and [x['node'] for x in mapper]==queue and inside_count==3248,'complete pending mapper/parent-box geometry differs')
    pids=[br['master_pid'],*[x['pid'] for x in final['pool_evidence']],numeric['master_pid'],*[x['pid'] for x in numeric['pool_exits']]]
    for pid in pids:
        try:os.kill(pid,0)
        except ProcessLookupError:continue
        raise AssertionError('production/audit PID still present or reused')
    for p,h in {**REG['input_locks'],**REG['protected_old_locks'],**reg['input_locks']}.items():C.require(C.file_sha(p)==h,'input changed after finite audit')
    result={'task':'audit_g110116_structure','status':'first4096_H128_values_source_journal_and_partial_tree_independently_qualified','pid':os.getpid(),'wall_seconds':time.monotonic()-began,'registration_sha256':C.file_sha(rp),'source_gate':source_tree,'request_gate':request_result,'new_H128_points':4096,'new_independent_components':20480,'new_Delta_checks':4096,'new_integer_corners':65536,'new_passed':850,'new_geometry_split':3246,'passed':96875,'splits':51884,'queue':7369,'untested_original':877,'closed_children':6492,'inflight':0,'faults':0,'H128_obstructions':0,'old_typed_cache_records':184990,'logical_total_typed_cache':189086,'source96025_passed_and194_H128_preserved':True,'all_new_H128_lower_above110u':True,'geometry_failure_not_precision_failure':True,'journal_events':12289,'exact_overlay_reconstructed':True,'tree':tree,'certificates':certs,'child_parent_anchor_inside':3248,'child_parent_anchor_outside':3244,'actual_PIDs_absent':pids,'protected_old_files_unchanged':204,'old_H_values_recomputed':False,'new118_or877_H_dispatched':False,'global_g110_qualified':False,'Cn_C5_I_updated':False}
    p=HERE/'structure_result.json';write(p,result);print(json.dumps({k:v for k,v in result.items() if k not in ['certificates','actual_PIDs_absent'] and not isinstance(v,dict)}),flush=True);print(json.dumps({'result_sha256':C.file_sha(p)}),flush=True)


if __name__=='__main__':main()
