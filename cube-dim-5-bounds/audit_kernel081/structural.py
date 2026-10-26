#!/usr/bin/env python3
"""Frozen070 structural and directed-certificate audit, independent of producer."""
import argparse,json,gzip,os,time,resource
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from itertools import product,zip_longest
from math import lcm
from collections import Counter
from common import require,file_sha,digest,canonical,Stream,node,physical_box,validate_tree
from arith import log_bounds,S

HERE=Path(__file__).resolve().parent


def write(name,obj):
    path=HERE/name;temp=Path(str(path)+'.tmp');temp.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n');temp.replace(path)


def geometry(point,vertices):
    p=list(map(F,point));D=lcm(*(x.denominator for x in p));integer=[int(x*D) for x in p]
    projections=[sum(a*b for a,b in zip(integer,v)) for v in vertices]
    return {'eligible':0 not in projections and len(set(projections))==280,'zero_count':projections.count(0),'distinct_count':len(set(projections)),'common_denominator':str(D)}


def select_point(n,roots,vertices):
    lo,hi=physical_box(roots,n);mid=[(a+b)/2 for a,b in zip(lo,hi)];eta=min(b-a for a,b in zip(lo,hi))/4096
    tails=[mid]+[[x+eta*beta**j for j,x in enumerate(mid)] for beta in (F(1,8),F(1,16),F(1,32))]
    tests=[]
    for index,tail in enumerate(tails):
        point=['1',*map(str,tail)];test={'candidate':index,**geometry(point,vertices)};tests.append(test)
        if test['eligible']:return point,index,tests
    raise AssertionError('no eligible finite generic candidate')


def certify(n,certificate,record,roots):
    require(certificate['mechanism']=='full_frame_16' and certificate['quadratic']=='Q_full=I4/16' and certificate['checked_points']==16,'wrong certificate mechanism')
    key=','.join(record['point']);require(certificate['anchor']==key,'certificate/cache anchor mismatch')
    values=[[F(a),F(b)] for a,b in record['values']]
    require(len(values)==5 and all(a<=b for a,b in values),'unordered value intervals')
    require(record['bits']==256 and len(record['attempts'])==1,'unexpected precision path')
    attempt=record['attempts'][0]
    require(attempt['bits']==256 and attempt['all_cones']==3504 and attempt['numeric_quality'] is True,'incomplete or bad production attempt')
    for a,b in values:require((b-a)*2**40<=max(abs(a),abs(b),F(1,10**6)),'actual interval width fails numeric-quality condition')
    require(values[0][0]>F(1,10000),'new H lower endpoint not above target')
    gradient=[[F(a),F(b)] for a,b in certificate['loggradient_directed']]
    require(len(gradient)==4 and all(a<=b for a,b in gradient),'unordered physical gradient ratios')
    for (gl,gu),(dl,du) in zip(gradient,values[1:]):
        for h in values[0]:
            for derivative in (dl,du):require(gl*h<=derivative<=gu*h,'directed physical derivative/H ratio fails')
    delta=F(certificate['Delta_lower']);ratio=values[0][0]*10000
    lower,upper=log_bounds(ratio.numerator,ratio.denominator)
    require(delta>0 and delta<=F(lower,S),'saved directed Delta is not below independent fixed-point log lower bound')
    lo,hi=physical_box(roots,n);point=list(map(F,record['point']))
    require(point[0]==1 and all(a<=x<=b for a,x,b in zip(lo,point[1:],hi)),'anchor outside independently reconstructed box')
    offsets=[(a-x,b-x) for a,x,b in zip(lo,point[1:],hi)]
    vd=lcm(*(v.denominator for pair in offsets for v in pair));gd=lcm(delta.denominator,*(v.denominator for pair in gradient for v in pair))
    zz=[tuple(int(x*vd) for x in pair) for pair in offsets];gg=[tuple(int(x*gd) for x in pair) for pair in gradient];di=int(delta*gd)
    rows=[]
    for bits in product((0,1),repeat=4):
        z=[zz[j][bits[j]] for j in range(4)]
        numerator=di*vd*vd+vd*sum(min(a*x,b*x) for (a,b),x in zip(gg,z))-256*gd*sum(x*x for x in z)
        rows.append((numerator,bits))
    minimum=min(x[0] for x in rows);margin=F(minimum,gd*vd*vd)
    require(certificate['passed']==(minimum>0) and F(certificate['minimum_log_margin'])==margin,'exact16corner margin or verdict differs')
    require((minimum,tuple(certificate['worst_corner_bits'])) in rows,'saved worst corner is not a minimizer')
    return {'node':list(n),'anchor':key,'passed':minimum>0,'minimum_log_margin':str(margin),'checked_points':16,'independent_Delta_log_lower':str(F(lower,S))}


def baseline(path,counts):
    stream=Stream(path);arrays={k:[] for k in ('passed','queue','splits')};roots=None;metadata={}
    cache_hash=sha256();cache_count=0
    while True:
        key=stream.key()
        if key is None:break
        if key=='cache':
            cache_hash.update(b'{');first=True
            for k,v in stream.items(True):
                if not first:cache_hash.update(b',')
                first=False;cache_hash.update(canonical(k));cache_hash.update(b':');cache_hash.update(canonical(v));cache_count+=1
            cache_hash.update(b'}')
        elif key in arrays:
            arrays[key]=[node(v if key=='queue' else v['node']) for v in stream.items()]
        else:
            value=stream.value()
            if key=='roots':roots=value
            elif key in ('model','epsilon','target','active','inflight','resume_jobs','fingerprint'):metadata[key]=value
    stream.close()
    for k in arrays:require(len(arrays[k])==counts[k], 'source array count mismatch')
    require(cache_count==counts['cache']==178722,'source cache count mismatch')
    require(roots==[[i,j,k,l] for i in range(32) for j in range(i+1) for k in range(j+1) for l in range(k+1)],'root enumeration differs')
    require(metadata['model']=='alpha1' and metadata['epsilon']=='1/64' and metadata['target']=='1/10000','wrong model qualification')
    require(metadata['active'] is None and metadata['inflight']=={} and metadata['resume_jobs']=={},'source contains active/resume work')
    tree=validate_tree(roots,arrays['passed']+arrays['queue'],set(arrays['splits']),'source077')
    return {'arrays':arrays,'roots':roots,'metadata':metadata,'tree':tree,'cache_sha256':cache_hash.hexdigest()}


def compare_final(source_path,final_path,last,reg):
    old,new=Stream(source_path),Stream(final_path);a,b=old.key(),new.key();extra={};missing=object()
    direct={'old_cache':0,'new_cache':0,'old_passed':0,'old_splits':0,'new_passed':0,'new_splits':0,'old_top_fields':0}
    additions={'task070','status070','fingerprint070','conversion070'}
    while a is not None or b is not None:
        if b in additions and (a is None or b<a):extra[b]=new.value();b=new.key();continue
        require(a==b and a is not None,'unapproved final schema change')
        if a=='cache':
            olditer=old.items(True);entry=next(olditer,missing);seen=set()
            for key,record in new.items(True):
                if entry is not missing and key==entry[0]:
                    require(record==entry[1] and key not in last['new_cache'],'old cache entry changed or overlaps new cache');direct['old_cache']+=1;entry=next(olditer,missing)
                else:
                    require(key in last['new_cache'] and key not in seen and (entry is missing or key<entry[0]),'unexpected/missing final cache key')
                    require(record==last['new_cache'][key],'new cache differs from overlay');seen.add(key);direct['new_cache']+=1
            require(entry is missing and seen==set(last['new_cache']),'final old/new cache coverage mismatch')
        elif a in ('passed','splits'):
            target=new.items()
            for record in old.items():require(next(target,missing)==record,'source record prefix changed');direct['old_'+a]+=1
            for expected in last['new_'+a]:require(next(target,missing)==expected,'new final record differs from overlay');direct['new_'+a]+=1
            require(next(target,missing) is missing,'extra final array record')
        elif a=='queue':
            old.value();require(list(new.items())==[],'final queue is nonempty')
        else:
            require(old.value()==new.value(),'source top-level field changed: '+a);direct['old_top_fields']+=1
        a,b=old.key(),new.key()
    old.close();new.close();require(set(extra)==additions,'wrong final added keys')
    require(extra['task070']=='upper_cover070' and extra['status070']=='closed_root_cover_candidate_waiting_main_review','wrong final status/task')
    final_overlay=Path(reg['last_overlay']);runreg=Path(reg['run_registration'])
    expected={'source_state':reg['source'],'source_state_sha256':reg['source_sha256'],'last_overlay':str(final_overlay),'last_overlay_sha256':file_sha(final_overlay),
              'new_counts':last['counts'],'new_passed':len(last['new_passed']),'new_splits':len(last['new_splits']),'math':'075 Q_full=I4/16','closed_root_cover_candidate':True,'global_g_or_C5_claim':False}
    require(extra['conversion070']==expected,'final conversion metadata mismatch')
    require(extra['fingerprint070']==digest({'run_registration':file_sha(runreg),'last_overlay':file_sha(final_overlay)}),'final fingerprint migration differs')
    return {'direct_object_comparisons':direct,'new_top_fields':extra}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');args=parser.parse_args()
    reg=json.loads((HERE/'registration.json').read_text())
    if not args.run:print(json.dumps({'registered':True,'new_points':1318,'wall_deadline':None}));return
    started=time.monotonic();pid=os.getpid();print(json.dumps({'pid':pid,'stage':'immutable_locks'}),flush=True)
    require(file_sha(reg['manifest'])==reg['manifest_sha256'],'manifest mismatch');manifest=json.loads(Path(reg['manifest']).read_text())
    inputs={reg['manifest']:reg['manifest_sha256'],reg['source']:reg['source_sha256'],reg['final_state']:reg['final_state_sha256'],reg['fan']:reg['fan_sha256']}
    inputs.update(manifest['source_locks']);inputs.update(manifest['external_files'])
    production_root=Path(reg['manifest']).parent
    inputs[str(production_root/'run.py')]=manifest['production_code_sha256']
    for item in manifest['batches']:
        for name,h in item['files'].items():inputs[str(Path(item['directory'])/name)]=h
    for p,h in {**inputs,**reg['protected_old_locks']}.items():require(file_sha(p)==h,'immutable input/protected hash mismatch: '+p)
    fan=json.loads(gzip.decompress(Path(reg['fan']).read_bytes()))
    base=baseline(reg['source'],reg['source_counts']);require(base['cache_sha256']=='f3b174e2bab3fcaed2c291090b378ddac7414393a6fb7a7bceffedca0ffaccac','source cache canonical mismatch')
    result={'task':'audit_kernel081_structural','pid':pid,'status':'pending_values','code_locks':{f:file_sha(HERE/f) for f in ['registration.json','common.py','arith.py','structural.py']},'inputs':inputs,
            'source_tree':base['tree'],'source_cache_sha256':base['cache_sha256'],'batches':[],'new_certificates':[],'new_value_independence':'pending numeric_result all1318 H+4grad','wall_deadline':None,'MATLAB_calls':0,'author_modules_imported':False}
    previous=None;total_events=0;allpoints=set();allworkers=set()
    remaining_path=next(p for p in manifest['source_locks'] if p.endswith('remaining_queue.json'));remaining=json.loads(Path(remaining_path).read_text())['records']
    for item in manifest['batches']:
        root=Path(item['directory']);br=json.loads((root/'registration.json').read_text());initial=json.loads((root/'initial_overlay.json').read_text());final=json.loads((root/'state.json').read_text())
        require(br['source_state_sha256']==reg['source_sha256'] and br['source_state']==reg['source'] and br['wall_deadline'] is None and br['workers']==12 and br['max_new_anchor_starts']==4096,'batch scope mismatch')
        run_registration=Path(item['directory']).parent/'registration.json'
        require(br['code_sha256']==manifest['production_code_sha256'] and br['run_registration_sha256']==file_sha(run_registration) and br['locks']==manifest['source_locks'],'registered code/math/source links differ')
        require(br['base_cache_canonical_sha256']==base['cache_sha256'],'batch source cache reference differs')
        require(br['selected_nodes']==initial['queue'][:4096] and digest(br['selected_nodes'])==br['selected_sha256'],'registered finite selection mismatch')
        require(len(br['selected_nodes'])==len(initial['queue'])==item['selected'],'this finite batch did not consume exactly the complete registered frontier')
        if previous is None:
            require(br['parent_overlay'] is None and br['parent_overlay_sha256'] is None,'first overlay has a parent')
            require(set(map(node,initial['queue']))==set(base['arrays']['queue']) and initial['queue']==[r['node'] for r in remaining],'initial difficult ordering differs from frozen077')
            for key in ('new_cache','new_passed','new_splits'):require(not initial[key],'first overlay contains unregistered new data')
            require(initial['counts']=={'cache_hits':0,'kernel_segment_calls':0,'new_anchor_starts':0,'new_complete':0},'first numerical counters do not start at zero')
        else:
            prior_item=manifest['batches'][item['batch']-2]
            require(br['parent_overlay']==str(Path(prior_item['directory'])/'state.json') and br['parent_overlay_sha256']==prior_item['state_sha256'],'overlay parent hash/path mismatch')
            for key in ('new_cache','new_passed','new_splits','queue','hardness','counts'):require(initial[key]==previous[key],'cumulative overlay inheritance changed: '+key)
        expected=json.loads(json.dumps(initial));expected['new_cache']=dict(initial['new_cache']);expected['new_passed']=list(initial['new_passed']);expected['new_splits']=list(initial['new_splits'])
        dispatched={};received={};merges=[];queue=list(initial['queue']);hard=dict(initial['hardness']);pointmap={r['path']:r for r in json.loads((root/'point_manifest.json').read_text())['point_checkpoints']}
        prefix='0'*64;sequence=0
        for line in (root/'journal.jsonl').read_text().splitlines():
            event=json.loads(line);body=event['body'];kind=event['kind'];claimed=event['hash'];payload={k:v for k,v in event.items() if k!='hash'}
            require(event['previous']==prefix and event['seq']==sequence+1 and digest(payload)==claimed,'journal hash chain differs');prefix=claimed;sequence+=1
            if kind=='dispatch':
                ordinal=body['ordinal'];require(ordinal==len(dispatched) and body['node']==br['selected_nodes'][ordinal] and queue.pop(0)==body['node'],'dispatch ordinal/frontier differs')
                require(body['cached'] is False and ordinal not in dispatched,'unexpected old cache evaluation or repeated dispatch')
                point,candidate,tests=select_point(node(body['node']),base['roots'],fan['vertices'])
                require(body['point']==point and body['candidate']==candidate and body['geometry_tests']==tests,'independent generic anchor/guard differs')
                dispatched[ordinal]=body
            elif kind=='received':
                ordinal=body['ordinal'];got=body['result'];require(ordinal in dispatched and ordinal not in received and got['kind']=='complete' and got['calls']==1,'received fault/counter/retry or duplicate')
                task=dispatched[ordinal];entry=pointmap[task['point_checkpoint']];require(file_sha(entry['path'])==entry['sha256'],'raw point checkpoint hash differs')
                raw=json.loads(Path(entry['path']).read_text());record=got['record']
                require(entry['point']==task['point'] and entry['bits']==256 and entry['stage']=='complete','point manifest qualification fields differ')
                require(raw['task_sha256']==digest(task) and raw['stage']=='complete' and raw['partial'] is None and raw['precision_index']==0 and raw['kernel_segment_calls']==1,'point task/precision/completeness differs')
                require(raw['record']==record and raw['point']==task['point'] and raw['attempts']==record['attempts'] and raw['pid']==got['pid']==entry['pid'],'raw point/record/worker mismatch')
                require(record['point']==task['point'] and record['candidate']==task['candidate'] and record['geometry_tests']==task['geometry_tests'],'returned anchor fields differ')
                key=','.join(record['point']);require(key not in allpoints and key not in expected['new_cache'],'new point/key repeats');allpoints.add(key)
                expected['new_cache'][key]=record;received[ordinal]=got
            elif kind=='merge':
                ordinal=body['ordinal'];require(ordinal==len(merges) and ordinal in received,'merge order/received result differs')
                task=dispatched[ordinal];n=node(task['node']);record=received[ordinal]['record'];cert=body['certificate'];proof=certify(n,cert,record,base['roots']);result['new_certificates'].append(proof);merges.append(ordinal)
                hard.pop(str(n[0])+':'+n[1],None)
                if cert['passed']:expected['new_passed'].append({'node':list(n),'certificate':cert})
                else:
                    children=[[n[0],n[1]+'0'],[n[0],n[1]+'1']];expected['new_splits'].append({'node':list(n),'axis':len(n[1])%4,'children':children,'failed_certificate':cert});queue.extend(children)
                    for child in children:hard[str(child[0])+':'+child[1]]=[record['values'][0][0],cert['minimum_log_margin']]
            elif kind=='finish':
                require(len(merges)==len(dispatched)==len(received)==len(br['selected_nodes']),'batch ends with incomplete task sets')
                pool=body['pool_evidence'];require(pool['close_called'] is True and pool['join_returned'] is True and pool['master_pid']==br['master_pid'] and pool['terminated_on_actual_fault'] is False,'pool did not close cleanly')
                require(len(pool['worker_pids'])==12 and sorted(pool['worker_pids'])==sorted(e['pid'] for e in pool['exits']) and all(e['exitcode']==0 and e['is_alive'] is False for e in pool['exits']),'12worker exit receipts differ')
                require(pool['master_pid']==item['master_pid'] and pool['worker_pids']==item['worker_pids'] and {r['pid'] for r in received.values()}<=set(pool['worker_pids']),'pool/actual result PID identities differ')
                allworkers.update(pool['worker_pids']);allworkers.add(pool['master_pid']);expected['pool_evidence']=pool;expected['status']=body['status']
            else:raise AssertionError('unregistered journal event kind')
        queue.sort(key=lambda n:(*map(F,hard[str(n[0])+':'+n[1]]),*n))
        expected['queue']=queue;expected['hardness']=hard;expected['inflight']={};expected['resume_jobs']={};expected['next_merge']=len(merges);expected['journal_hash']=prefix;expected['journal_seq']=sequence
        expected['counts']={**initial['counts'],'new_anchor_starts':initial['counts']['new_anchor_starts']+len(dispatched),'new_complete':initial['counts']['new_complete']+len(received),'kernel_segment_calls':initial['counts']['kernel_segment_calls']+len(received)}
        require(expected==final,'independent journal-derived overlay differs from frozen state')
        require(final['counterexample'] is None and final['faults']==[] and not final['inflight'] and not final['resume_jobs'],'batch has counter/fault/live work')
        frontier=base['arrays']['passed']+[node(r['node']) for r in final['new_passed']]+list(map(node,queue))
        splits=set(base['arrays']['splits'])|{node(r['node']) for r in final['new_splits']}
        tree=validate_tree(base['roots'],frontier,splits,'batch'+str(item['batch']))
        result['batches'].append({'batch':item['batch'],'selected':len(dispatched),'journal_events':sequence,'overlay_exactly_reconstructed':True,'tree':tree,'cumulative_new_passed':len(final['new_passed']),'cumulative_new_splits':len(final['new_splits']),'queue':len(queue),'all12_workers_closed_joined_exit0':True})
        total_events+=sequence;previous=final;write('structural_result.json',result)
        print(json.dumps({'pid':pid,'structural_batch':item['batch'],'queue':len(queue),'verified_new_certificates':len(result['new_certificates'])}),flush=True)
    require(len(allpoints)==1318 and len(previous['new_cache'])==1318 and len(previous['new_passed'])==1001 and len(previous['new_splits'])==317 and not previous['queue'],'final finite counts/closed cover differ')
    reg_for_compare={**reg,'last_overlay':str(Path(manifest['batches'][-1]['directory'])/'state.json'),'run_registration':str(Path(manifest['batches'][0]['directory']).parent/'registration.json')}
    result['materialization']=compare_final(reg['source'],reg['final_state'],previous,reg_for_compare)
    require(total_events==3959,'journal event count differs');result['total_journal_events']=total_events
    result['new_certificates_verified']=len(result['new_certificates']);result['new_corner_points_verified']=16*len(result['new_certificates'])
    result['new_passed_minimum_exact_margin']=str(min(F(r['minimum_log_margin']) for r in result['new_certificates'] if r['passed']))
    result['new_cached_minimum_exact_H_lower']=str(min(F(r['values'][0][0]) for r in previous['new_cache'].values()))
    result['final_counts']={'passed':len(base['arrays']['passed'])+len(previous['new_passed']),'splits':len(base['arrays']['splits'])+len(previous['new_splits']),'queue':0,'cache':reg['source_counts']['cache']+len(previous['new_cache']),'roots':52360}
    absent=[]
    for pid_to_check in sorted(allworkers):
        try:os.kill(pid_to_check,0)
        except ProcessLookupError:absent.append(pid_to_check)
    result['production_master_worker_pids_os_absent']=absent;require(len(absent)==65,'a production PID is still present or reused')
    for p,h in {**inputs,**reg['protected_old_locks']}.items():require(file_sha(p)==h,'input/old protection changed during audit')
    result['all_input_old_locks_unchanged']=True;result['status']='structural_and_all16corner_passed_pending_values';result['global_g_qualified']=False
    result['wall_seconds']=time.monotonic()-started;result['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    write('structural_result.json',result)
    print(json.dumps({k:result[k] for k in ('pid','status','new_certificates_verified','new_corner_points_verified','final_counts','wall_seconds','peak_rss_bytes')}),flush=True)
    print(json.dumps({'structural_result_sha256':file_sha(HERE/'structural_result.json')}),flush=True)


if __name__=='__main__':main()
