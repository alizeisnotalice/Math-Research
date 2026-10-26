#!/usr/bin/env python3
"""Frozen194 jet/epsilon/certificate/replay/fullstate and exact closed tree gate."""
import sys
sys.dont_write_bytecode=True
import json,os,time
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import lcm
from interfaces import REG,HERE,common as C,arith as A


def write(path,value):
    q=Path(str(path)+'.tmp');q.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n');q.replace(path)


def snapshot(path):
    stream=C.Stream(path);answer={'metadata':{},'passed_nodes':[],'passed_hashes':[],'split_nodes':[],'split_hashes':[],'cache':{}}
    while True:
        field=stream.key()
        if field is None:break
        if field=='cache_by_epsilon':
            stream.char('{')
            for j,expected_epsilon in enumerate(('1/128','1/64')):
                if j:stream.char(',')
                epsilon=stream.value();C.require(epsilon==expected_epsilon,'typed cache epsilon key differs');stream.char(':');records={}
                for key,row in stream.items(True):
                    C.require(key==','.join(row['point']),'cache point key differs')
                    if epsilon=='1/128':C.require(row['epsilon']==epsilon,'H128 record epsilon differs')
                    records[key]=C.digest(row)
                answer['cache'][epsilon]=records
            stream.char('}')
        elif field=='passed':
            for row in stream.items():
                C.require(row['leaf_epsilon'] in ('1/64','1/128') and row['lower_bound_object']=='G' and row['component_lower_function']=='H_epsilon' and row['target']=='21/200000','leaf epsilon/object qualification differs')
                answer['passed_nodes'].append(C.node(row['node']));answer['passed_hashes'].append(C.digest(row))
        elif field=='splits':
            for row in stream.items():
                node=C.node(row['node']);C.require(row['axis']==len(node[1])%4 and row['children']==[[node[0],node[1]+'0'],[node[0],node[1]+'1']],'closed physical partition differs')
                answer['split_nodes'].append(node);answer['split_hashes'].append(C.digest(row))
        else:answer['metadata'][field]=stream.value()
    stream.close();return answer


def certificate(node,cert,record,roots):
    C.require(record['epsilon']==cert['leaf_epsilon']=='1/128' and cert['target']=='21/200000' and cert['lower_bound_object']=='G' and cert['component_function']=='H_epsilon','certificate epsilon/model binding differs')
    C.require(cert['mechanism']=='075_full_frame_16' and cert['quadratic']=='I4/(16epsilon^2)' and F(cert['penalty_coefficient'])==1024 and cert['checked_points']==16,'EPS128 curvature differs')
    C.require(cert['point_key']==','.join(record['point']),'certificate anchor differs')
    values=[list(map(F,pair)) for pair in record['values']];C.require(len(values)==5 and all(a<=b for a,b in values) and values[0][0]>F(21,200000),'H128/target interval fails')
    for lo,hi in values:C.require((hi-lo)*2**40<=max(abs(lo),abs(hi),F(1,10**6)),'actual production width condition fails')
    gradients=[list(map(F,pair)) for pair in cert['loggradient_directed']];C.require(len(gradients)==4 and all(lo<=hi for lo,hi in gradients),'gradient interval order differs')
    for (gl,gu),(dl,du) in zip(gradients,values[1:]):
        for h in values[0]:
            for d in (dl,du):C.require(gl*h<=d<=gu*h,'same-epsilon physical derivative/H containment fails')
    delta=F(cert['Delta_lower']);ratio=values[0][0]/F(21,200000);lower,upper=A.log_bounds(ratio.numerator,ratio.denominator)
    C.require(0<delta<=F(lower,A.S),'fresh H128 Delta not below independent integer log')
    lo,hi=C.physical_box(roots,node);point=list(map(F,record['point']))
    C.require(point[0]==1 and all(a<=x<=b for a,x,b in zip(lo,point[1:],hi)),'anchor outside independently constructed closed box')
    offsets=[(a-x,b-x) for a,x,b in zip(lo,point[1:],hi)];vd=lcm(*(x.denominator for pair in offsets for x in pair));gd=lcm(delta.denominator,*(x.denominator for pair in gradients for x in pair))
    zz=[tuple(int(x*vd) for x in pair) for pair in offsets];gg=[tuple(int(x*gd) for x in pair) for pair in gradients];di=int(delta*gd);rows=[]
    for bits in product((0,1),repeat=4):
        v=[zz[j][bits[j]] for j in range(4)];num=di*vd*vd+vd*sum(min(a*x,b*x) for (a,b),x in zip(gg,v))-1024*gd*sum(x*x for x in v);rows.append((num,bits))
    minimum,bits=min(rows);margin=F(minimum,gd*vd*vd)
    C.require(minimum>0 and cert['passed'] is True and F(cert['minimum_log_margin'])==margin and tuple(cert['worst_corner_bits'])==bits,'all16 exact H128 corners do not strictly pass')
    return {'node':list(node),'margin':str(margin),'certificate_sha256':C.digest(cert),'record_sha256':C.digest(record),'closed_corners':16,'penalty':1024}


def main():
    began=time.monotonic();rp=HERE/'closed_gate_registration.json';registration=json.loads(rp.read_text())
    for p,h in {**registration['input_locks'],**REG['protected_old_locks']}.items():C.require(C.file_sha(p)==h,'frozen/protected input changed: '+p)
    srcgate=json.loads((HERE/'source_gate_result.json').read_text());C.require(srcgate['status']=='mixed098_initial_source_migration_typed_cache_and_tree_qualified_partial','source gate incomplete')
    numeric_dir=HERE/'numeric128_batch0001';numeric=json.loads((numeric_dir/'index.json').read_text());C.require(numeric['status']=='all194_H128_plus4grad_strictly_verified' and numeric['verified_points']==numeric['global_verified_prefix']==194,'new numeric gate incomplete')
    for p,h in numeric['input_locks'].items():C.require(C.file_sha(p)==h,'numeric proof input changed')
    b=Path(registration['batch_directory']);br=json.loads((b/'registration.json').read_text());initial=json.loads((b/'initial_overlay.json').read_text());final=json.loads((b/'state.json').read_text());tasks=json.loads((b.parent/'execution_tasks_194.json').read_text());pm=json.loads((b/'point_manifest.json').read_text())
    C.require(br['source_sha256']==initial['source_sha256']==final['source_sha256']==REG['source_initial_sha256'] and br['epsilon']==final['epsilon']=='1/128' and br['target']==final['target']=='21/200000','run source/epsilon/target differs')
    C.require(br['finite_selected_nodes']==initial['queue']==[t['node'] for t in tasks] and C.digest(initial['queue'])==br['finite_selected_sha256'] and len(tasks)==194,'registered queue/tasks differ')
    C.require(initial['new_cache_H128']=={} and initial['new_passed']==initial['new_splits']==initial['faults']==initial['obstructions_H128']==[] and initial['inflight']=={},'initial overlay contains unqualified values')
    source=snapshot(REG['source_initial']);roots=source['metadata']['roots']
    C.require(source['metadata']['queue']==initial['queue'] and len(source['passed_nodes'])==100804 and len(source['cache']['1/64'])==184796 and source['cache']['1/128']=={},'frozen source migration differs')
    proofs={}
    for plan in numeric['plans']:
        state=json.loads((numeric_dir/('block%04d.json'%plan['block'])).read_text());C.require(state['status']=='complete' and len(state['points'])==plan['stop_exclusive']-plan['start'],'numeric successful block incomplete')
        for proof in state['points']:
            index=proof['index'];entry=pm['records'][index];C.require(proof['point']==entry['point'] and proof['raw_checkpoint_sha256']==entry['point_checkpoint_sha256'] and proof['H_multiplier']==16 and proof['gradient_multiplier']==32,'numeric source association differs')
            C.require(proof['independent_EPS64_input']==[str(2*F(x)) for x in entry['point']],'independent scaling input differs')
            C.require(len(proof['independent_H128_values_512_dyadic'])==5,'numeric proof missing components')
            for ball,original in zip(proof['independent_H128_values_512_dyadic'],entry['record']['values']):
                lo,hi=map(F,ball);L,U=map(F,original);C.require(L<lo<=hi<U,'numeric enclosure not strictly inside original')
            C.require(index not in proofs,'duplicate numeric proof');proofs[index]=proof
    C.require(set(proofs)==set(range(194)),'numeric proof set differs')
    expected=json.loads(json.dumps(initial));queue=list(initial['queue']);dispatch={};received={};seq=0;previous='0'*64;merges=0;certs=[];calls=0
    for line in (b/'journal.jsonl').read_text().splitlines():
        event=json.loads(line);body=event['body'];kind=event['kind'];C.require(event['seq']==seq+1 and event['previous']==previous and C.digest({k:v for k,v in event.items() if k!='hash'})==event['hash'],'journal chain differs');seq+=1;previous=event['hash']
        if kind=='dispatch':
            task=body['task'];i=task['ordinal'];C.require(i==len(dispatch) and task==tasks[i] and task['node']==queue.pop(0) and task['epsilon']=='1/128' and task['H128_cache_hit'] is False,'dispatch not registered finite typed task')
            dispatch[i]=task;expected['inflight'][str(i)]=dict(task)
        elif kind=='received':
            i=body['ordinal'];got=body['returned'];task=dispatch[i];entry=pm['records'][i]
            C.require(i not in received and got['kind']=='complete' and got['epsilon']=='1/128' and got['pair_cache_key']==task['pair_cache_key'],'counter/fault/duplicate/wrong epsilon returned')
            C.require(entry['ordinal']==i and entry['node']==task['node'] and entry['point']==task['point'] and entry['epsilon']=='1/128' and entry['stage']=='complete','point manifest association differs')
            C.require(entry['point_checkpoint']==task['point_checkpoint'] and C.file_sha(task['point_checkpoint'])==entry['point_checkpoint_sha256'],'raw checkpoint SHA/path differs')
            raw=json.loads(Path(task['point_checkpoint']).read_text());record=got['record'];C.require(raw['stage']=='complete' and raw['partial'] is None and raw['record']==entry['record']==record,'raw/return record differs')
            C.require(raw['epsilon']==record['epsilon']=='1/128' and raw['pair_cache_key']==task['pair_cache_key'] and raw['point']==record['point']==task['point'],'raw epsilon/point differs')
            C.require(raw['kernel_sha256']==record['kernel_sha256']==pm['kernel_sha256']==br['kernel_sha256'],'kernel binding differs')
            C.require(raw['pid']==entry['worker_pid']==got['pid'] and raw['calls']==got['calls']==1 and raw['precision_index']==0 and record['bits']==256,'worker/precision/calls differs')
            C.require(record['attempts']==raw['attempts']==got['attempts'] and len(record['attempts'])==1 and record['attempts'][0]['all_cones']==3504 and record['attempts'][0]['numeric_quality'] is True,'incomplete/poor jet')
            C.require(record['candidate']==task['candidate'] and record['geometry_tests']==task['geometry_tests'],'physical anchor generic qualification differs')
            if task['old_H64_interval']:C.require(F(record['values'][0][1])>=F(task['old_H64_interval'][0]),'cross epsilon monotonicity consistency fault')
            received[i]=got;calls+=got['calls']
        elif kind=='merge':
            i=body['ordinal'];C.require(i==merges and i in received,'merge order differs');merges+=1;task=dispatch[i];record=received[i]['record'];cert=body['certificate'];node=C.node(task['node'])
            certs.append(certificate(node,cert,record,roots));key=','.join(record['point']);C.require(key not in expected['new_cache_H128'],'H128 cache key repeated');expected['new_cache_H128'][key]=record
            expected['new_passed'].append({'node':list(node),'certificate':cert,'leaf_epsilon':'1/128','component_lower_function':'H_epsilon','lower_bound_object':'G','target':'21/200000'});expected['inflight'].pop(str(i))
        elif kind=='finish':
            C.require(len(dispatch)==len(received)==merges==194 and body['status']=='finite194_complete','finite194 not completely resolved')
            pool=body['pool_evidence'];C.require(len(pool)==12 and all(p['exitcode']==0 and p['alive'] is False for p in pool),'production pool exit differs');C.require({g['pid'] for g in received.values()}<=set(x['pid'] for x in pool),'point PID not in pool')
            expected.update(pool_evidence=pool,status=body['status'],wall_seconds=body['wall_seconds'])
        else:raise AssertionError('unknown journal event')
    expected.update(queue=queue,next_merge=merges,journal_seq=seq,journal_hash=previous,counts={'H128_cache_hits':0,'kernel_segment_calls':calls,'new_anchor_starts':194,'new_complete':194})
    C.require(expected==final and seq==583 and queue==[] and final['inflight']=={} and final['faults']==final['obstructions_H128']==[] and final['new_splits']==[],'exact overlay reconstruction/closed status differs')
    full=snapshot(registration['final_state']);C.require(full['cache']['1/64']==source['cache']['1/64'] and full['cache']['1/128']=={k:C.digest(v) for k,v in final['new_cache_H128'].items()},'typed final cache changed or cross-epsilon values substituted')
    C.require(full['passed_hashes']==source['passed_hashes']+[C.digest(row) for row in final['new_passed']] and full['split_hashes']==source['split_hashes'],'final old records not retained or new records differ')
    metadata={**source['metadata'],'task':'upper_mixed098','status':'mixed_G105_closed_candidate','queue':[],'inflight':{},'resume_jobs':{},'accepted_leaf_epsilon':'per leaf,1/64 or1/128','g105_closed_cover':True,'new_H128_records':194,'production_authorized':True,
        'final098_source_sha256':REG['source_initial_sha256'],'final098_overlay_sha256':C.file_sha(b/'state.json'),'final098_registration_sha256':C.file_sha(b/'registration.json'),
        'final098_scope':'same G mixed-epsilon candidate;await main and100;no I or C5 conclusion',
        'fingerprint098_final':C.digest({'source':REG['source_initial_sha256'],'overlay':C.file_sha(b/'state.json'),'journal':C.file_sha(b/'journal.jsonl'),'point_manifest':C.file_sha(b/'point_manifest.json')})}
    C.require(full['metadata']==metadata,'final source metadata/fingerprint migration differs')
    tree=C.validate_tree(roots,full['passed_nodes'],set(full['split_nodes']),'mixedG105ClosedFinal')
    C.require(len(full['passed_nodes'])==100998 and len(full['split_nodes'])==48638,'closed frontier finite counts differ')
    pids=[br['master_pid'],*[p['pid'] for p in final['pool_evidence']],numeric['master_pid'],*[p['pid'] for p in numeric['pool_exits']]]
    for pid in pids:
        try:os.kill(pid,0)
        except ProcessLookupError:continue
        raise AssertionError('production/audit PID still present or reused')
    for p,h in {**registration['input_locks'],**REG['protected_old_locks']}.items():C.require(C.file_sha(p)==h,'frozen input changed after audit')
    result={'task':'audit_mixed100_closed_gate','status':'mixed_H64_H128_complete_closed_cover_numeric_and_structural_gates_passed_for_main',
        'pid':os.getpid(),'wall_seconds':time.monotonic()-began,'registration_sha256':C.file_sha(rp),'final_state_sha256':C.file_sha(registration['final_state']),
        'retained_H64_passed':100804,'new_H128_passed':194,'total_passed':100998,'new_H128_all5_values':970,'fresh_H128_Delta_checks':194,'fresh_H128_closed_corners':3104,'penalty':1024,
        'minimum_new_corner_margin':str(min(F(x['margin']) for x in certs)),'certificates':certs,'journal_events':583,'overlay_exact':True,'final_every_field_exact':True,
        'cache_H64_records_retained':184796,'cache_H128_records':194,'typed_cache_total':184990,'old_H64_counters_retained_in_history':4,'queue':0,'inflight':0,'faults':0,'H128_obstructions':0,
        'tree':tree,'all_production_and_audit_pids_os_absent':pids,'old_protected_files_unchanged':133,'old_H64_values_recomputed':False,
        'mixed_closed_face_G105_bound_under_locked075_099':True,'global_extension':'accepted083/029: full G symmetry and degree-4 only; main final proof acceptance retained','same_epsilon_H_homogeneity_used':False,'Cn_C5_I_updated':False}
    path=HERE/'closed_gate_result.json';write(path,result);print(json.dumps({k:result[k] for k in ['pid','status','wall_seconds','total_passed','fresh_H128_closed_corners','queue','minimum_new_corner_margin']}),flush=True);print(json.dumps({'result_sha256':C.file_sha(path)}),flush=True)


if __name__=='__main__':main()
