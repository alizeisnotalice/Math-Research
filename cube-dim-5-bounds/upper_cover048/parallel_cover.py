#!/usr/bin/env python3
"""One finite eight-process cover batch; import frozen046 API, single writer."""
from __future__ import annotations
import argparse,copy,gzip,hashlib,importlib.util,json,multiprocessing as mp,os,resource,sys,time
from collections import deque
from fractions import Fraction as F
from pathlib import Path
from flint import ctx

if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,x):
    q=p.with_suffix(p.suffix+'.tmp')
    with q.open('wb') as f:f.write(canon(x)+b'\n');f.flush();os.fsync(f.fileno())
    os.replace(q,p)
def load_api(path):
    spec=importlib.util.spec_from_file_location('locked046',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
API=None;KERNEL=None
def worker_init(api_path,api_sha,table_path,table_sha,epsilon,ready):
    global API,KERNEL
    assert sha(api_path)==api_sha and sha(table_path)==table_sha
    API=load_api(api_path);table=json.loads(gzip.decompress(Path(table_path).read_bytes()));KERNEL=API.FirstJetKernel(table,F(epsilon))
    ready.put({'pid':os.getpid(),'initialized':True})
def run_job(task):
    began=time.monotonic();attempts=[];calls=0;pi=0;partial=None
    try:
        while pi<3:
            bits=(256,384,512)[pi]
            if time.monotonic()>=task['deadline']:
                return {'kind':'partial','partial':partial,'precision_index':pi,'attempts':attempts,'calls':calls,'pid':os.getpid(),'seconds':time.monotonic()-began}
            got=KERNEL.evaluate(task['point'],bits,partial,task['deadline'],128);calls+=1
            if not got['complete']:
                return {'kind':'partial','partial':got['partial'],'precision_index':pi,'attempts':attempts,'calls':calls,'pid':os.getpid(),'seconds':time.monotonic()-began}
            values=got['values'];packed=list(map(API.pack,values));quality=API.numeric_quality(values)
            attempt={k:got[k] for k in ('bits','seconds','cancellation_ratio','minimum_integer_projection_gap','common_denominator')}
            attempt['numeric_quality']=quality;attempts.append(attempt)
            counter=values[0].is_finite() and API.endpoint(values[0],True)<F(task['target'])
            separated=values[0].is_finite() and API.endpoint(values[0])>F(task['target'])
            if counter or (quality and separated):
                record={'point':task['point'],'values':packed,'bits':bits,'candidate':task['candidate'],
                        'geometry_tests':task['geometry_tests'],'attempts':attempts}
                return {'kind':'counterexample' if counter else 'complete','record':record,'calls':calls,
                    'pid':os.getpid(),'seconds':time.monotonic()-began,'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
            attempts[-1]['values']=packed;pi+=1;partial=None
        return {'kind':'numerical_diagnostic','attempts':attempts,'calls':calls,'pid':os.getpid(),'seconds':time.monotonic()-began}
    except BaseException as exc:
        return {'kind':'worker_error','exception':repr(exc),'attempts':attempts,'calls':calls,'pid':os.getpid(),'seconds':time.monotonic()-began}

def as_json(state):return {**state,'queue':list(state['queue'])}
def partition(state,api):
    proxy={**state,'queue':list(state['queue'])+[entry['node'] for entry in state['inflight'].values()]}
    api.check_tree(proxy)
def apply_event(state,event,api):
    kind=event['kind'];body=event['body'];active=state['inflight']
    if kind=='dispatch':
        assert state['queue'].popleft()==body['node'];assert str(body['ordinal']) not in active
        active[str(body['ordinal'])]=body;state['batch_new_starts']+=not body['cached']
        state['total_new_anchor_starts']+=not body['cached']
        state['geometry_candidates_rejected']+=sum(not row['eligible'] for row in body['geometry_tests'])
    elif kind=='received':
        entry=active[str(body['ordinal'])];assert 'returned' not in entry;entry['returned']=body['result']
        result=body['result'];state['total_kernel_attempts_completed']+=len(result.get('record',{}).get('attempts',result.get('attempts',[])))
    elif kind=='merge':
        key=str(body['ordinal']);entry=active[key];result=entry['returned'];assert body['ordinal']==state['next_merge']
        state['next_merge']+=1;state['batch_merges']+=1
        if result['kind']=='complete':
            record=result['record'];anchor=api.pointkey(record['point'])
            if entry['cached']:state['batch_cache_hits']+=1;state['cache_hits']+=1;assert state['cache'][anchor]==record
            else:assert anchor not in state['cache'];state['cache'][anchor]=record;state['batch_new_completed']+=1
            cert=body['certificate'];node=entry['node']
            if cert['passed']:state['passed'].append({'node':node,'certificate':cert})
            else:
                box=api.box_of(node,state['roots']);axis=max(range(4),key=lambda j:box[1][j]-box[0][j]);children=[[node[0],node[1]+'0'],[node[0],node[1]+'1']]
                state['splits'].append({'node':node,'axis':axis,'children':children,'failed_certificate':cert});state['queue'].extend(children)
            del active[key]
        else:
            entry['merged_unresolved']=True
            if result['kind']=='counterexample':
                state['counterexample']={'point':entry['point'],'H_interval':result['record']['values'][0],
                    'node':entry['node'],'target':state['target'],'scope':'fixed-epsilon H target only; not fullG or C5 failure'}
    elif kind=='requeue':
        ordinals=body['ordinals'];assert ordinals==sorted(map(int,active))
        for ordinal in reversed(ordinals):
            entry=active.pop(str(ordinal));assert entry['merged_unresolved'];state['queue'].appendleft(entry['node'])
            state['resume_jobs'][api.nodekey(entry['node'])]=entry
    elif kind=='finish':
        state['status']=body['status'];state['pool_evidence']=body['pool_evidence']
        state['segments'].append({'task':'upper_cover048','pid':body['pool_evidence']['master_pid'],
            'new_anchor_starts':state['batch_new_starts'],'new_anchors_completed':state['batch_new_completed'],
            'seconds':body['seconds'],'stop':body['status']})
    else:raise ValueError(kind)
    state['journal_seq']=event['seq'];state['journal_hash']=event['hash']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source046',type=Path,required=True);parser.add_argument('--upstream-root',type=Path,required=True)
    parser.add_argument('--parent-state',type=Path,required=True);parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--recover-only',action='store_true');parser.add_argument('--recovery-output',type=Path)
    args=parser.parse_args();began=time.monotonic();base=Path(__file__).resolve().parent
    reg=json.loads((base/'registration.json').read_text());sources={'registration.json':sha(base/'registration.json'),'parallel_cover.py':sha(Path(__file__))}
    for rel,h in reg['source046_hashes'].items():assert sha(args.source046/rel)==h
    assert sha(args.parent_state)==reg['parent_state_sha256'];oldreg=json.loads((args.source046/'registration.json').read_text())
    for rel,h in oldreg['upstream_read_only_inputs'].items():assert sha(args.upstream_root/rel)==h
    ref=args.source046.parent/'upper_harnack043'
    for rel,h in oldreg['reference043_read_only_inputs'].items():assert sha(ref/rel)==h
    api=load_api(args.source046/'arb_cover.py');parent=json.loads(args.parent_state.read_text());api.check_tree(parent)
    assert parent['fingerprint']==reg['parent_fingerprint'] and parent['active'] is None
    assert digest({'sources':parent['source_hashes'],'inputs':oldreg['upstream_read_only_inputs'],'reference':oldreg['reference043_read_only_inputs']})==parent['fingerprint']
    selected=parent['queue'][:4096];assert len(selected)==reg['selected_queue_prefix']['count'] and digest(selected)==reg['selected_queue_prefix']['sha256']
    fingerprint=digest({'sources':sources,'source046':reg['source046_hashes'],'parent':reg['parent_state_sha256'],
        'inputs':oldreg['upstream_read_only_inputs'],'reference':oldreg['reference043_read_only_inputs'],'selection':digest(selected)})
    work=args.work_dir.resolve();assert not any((p/'.git').exists() for p in (work,*work.parents));work.mkdir(parents=True,exist_ok=True)
    statepath=work/'state.json';initialpath=work/'initial_state.json';journalpath=work/'journal.jsonl'
    if args.recover_only:
        state=json.loads(initialpath.read_text());assert state['fingerprint']==fingerprint;state['queue']=deque(state['queue'])
        expected_seq=0;previous='0'*64
        for line in journalpath.read_bytes().splitlines():
            event=json.loads(line);hashed=event.pop('hash');assert event['seq']==expected_seq+1 and event['previous']==previous and digest(event)==hashed
            event['hash']=hashed;apply_event(state,event,api);expected_seq+=1;previous=hashed
        partition(state,api);target=args.recovery_output or work/'recovered_state.json';atomic(target,as_json(state))
        print(json.dumps({'replay_events':expected_seq,'output':str(target),'sha256':sha(target),'density_evaluations':0,'inflight':len(state['inflight'])},indent=2));return
    assert not statepath.exists() and not journalpath.exists(),'one finite batch; recover-only does not dispatch'
    state=copy.deepcopy(parent);state['queue']=deque(state['queue']);parentcachehash=digest(state['cache']);parentpassed=len(state['passed']);parentsplits=len(state['splits'])
    state.update({'fingerprint':fingerprint,'parent_fingerprint':parent['fingerprint'],'parent_state_sha256':reg['parent_state_sha256'],
        'status':'initialized048','source048_hashes':sources,'selected_nodes':selected,'selected_sha256':digest(selected),'inflight':{},'resume_jobs':{},
        'batch_new_starts':0,'batch_new_completed':0,'batch_merges':0,'batch_cache_hits':0,'next_merge':0,'journal_seq':0,'journal_hash':'0'*64,
        'model':'alpha1','epsilon':reg['epsilon'],'target':reg['target'],'transition':'046 frozen parent -> only active048 handle'})
    partition(state,api);atomic(initialpath,as_json(state));atomic(statepath,as_json(state));initialsha=sha(initialpath)
    journal=journalpath.open('ab');compactions=1;completed_since_compact=0;events=0
    def emit(kind,body,sync=False):
        nonlocal events
        event={'seq':state['journal_seq']+1,'previous':state['journal_hash'],'kind':kind,'body':body};event['hash']=digest(event)
        journal.write(canon(event)+b'\n');journal.flush()
        if sync:os.fsync(journal.fileno())
        apply_event(state,event,api);events+=1
    def compact():
        nonlocal compactions,completed_since_compact
        journal.flush();os.fsync(journal.fileno());partition(state,api);atomic(statepath,as_json(state));compactions+=1;completed_since_compact=0
    tablepath=args.upstream_root/'astra_vertex017/boundary_triangulation.json.gz'
    table=json.loads(gzip.decompress(tablepath.read_bytes()));kernel=api.FirstJetKernel(table,F(reg['epsilon']))
    context=mp.get_context('spawn');ready=context.Queue();pool=context.Pool(reg['workers'],worker_init,
        (str(args.source046/'arb_cover.py'),reg['source046_hashes']['arb_cover.py'],str(tablepath),oldreg['upstream_read_only_inputs']['astra_vertex017/boundary_triangulation.json.gz'],reg['epsilon'],ready))
    processes=list(pool._pool);readies=[ready.get(timeout=30) for _ in range(reg['workers'])];pids=sorted(x['pid'] for x in readies)
    print(json.dumps({'master_pid':os.getpid(),'worker_pids':pids,'workers_initialized':len(readies),'selected_nodes':4096,
        'fingerprint':fingerprint,'parent_state_sha256':reg['parent_state_sha256']},sort_keys=True),flush=True)
    deadline=began+reg['seconds'];pending={};dispatch=0;stop_reason=None;cachekeys_inflight=set();worker_results=[]
    try:
        while True:
            if time.monotonic()>=deadline and stop_reason is None:stop_reason='time_checkpoint'
            while stop_reason is None and len(pending)<reg['workers'] and dispatch<len(selected):
                if state['batch_new_starts']>=reg['max_new_anchor_starts']:stop_reason='anchor_checkpoint';break
                node=selected[dispatch];assert state['queue'][0]==node;box=api.box_of(node,state['roots']);point,candidate,tests=api.choose_anchor(box,kernel)
                task={'ordinal':dispatch,'node':node,'point':None if point is None else list(map(str,point)),
                    'candidate':candidate,'geometry_tests':tests,'deadline':deadline,'target':reg['target'],'cached':False}
                if point is not None:
                    key=api.pointkey(point);task['cached']=key in state['cache'];assert key not in cachekeys_inflight
                emit('dispatch',task)
                if point is None:
                    emit('received',{'ordinal':dispatch,'result':{'kind':'no_generic_anchor','calls':0}},True);stop_reason='no_generic_anchor'
                elif task['cached']:emit('received',{'ordinal':dispatch,'result':{'kind':'complete','record':state['cache'][key],'calls':0}},True)
                else:cachekeys_inflight.add(key);pending[dispatch]=pool.apply_async(run_job,(task,))
                dispatch+=1
            for ordinal,result in list(pending.items()):
                if not result.ready():continue
                try:got=result.get()
                except BaseException as exc:got={'kind':'worker_error','exception':repr(exc),'calls':0}
                emit('received',{'ordinal':ordinal,'result':got},True);worker_results.append(got);del pending[ordinal]
                cachekeys_inflight.discard(api.pointkey(state['inflight'][str(ordinal)]['point']))
                if got['kind'] not in ('complete','partial') and stop_reason is None:stop_reason=got['kind']
                if got['kind']=='partial' and stop_reason is None:stop_reason='time_checkpoint'
            while str(state['next_merge']) in state['inflight'] and 'returned' in state['inflight'][str(state['next_merge'])]:
                ordinal=state['next_merge'];entry=state['inflight'][str(ordinal)];got=entry['returned'];body={'ordinal':ordinal}
                if got['kind']=='complete':
                    ctx.prec=got['record']['bits'];body['certificate']=api.certificate(api.box_of(entry['node'],state['roots']),got['record'],F(reg['epsilon']),F(reg['target']))
                emit('merge',body);completed_since_compact+=1
                if completed_since_compact>=64:
                    compact();print(f"merged {state['batch_merges']}; starts {state['batch_new_starts']}; pass+{len(state['passed'])-parentpassed}; split+{len(state['splits'])-parentsplits}; inflight {len(state['inflight'])}; {time.monotonic()-began:.3f}s",flush=True)
            if dispatch==len(selected) and stop_reason is None:stop_reason='anchor_checkpoint'
            if not pending and (stop_reason is not None or dispatch==len(selected)):break
            time.sleep(0.002)
    finally:
        pool.close();pool.join();ready.close();ready.join_thread()
    exits=[{'pid':p.pid,'exitcode':p.exitcode,'is_alive':p.is_alive()} for p in processes]
    assert len(exits)==reg['workers'] and all(p['exitcode']==0 and not p['is_alive'] for p in exits)
    unresolved=sorted(map(int,state['inflight']))
    if unresolved:emit('requeue',{'ordinals':unresolved},True)
    pool_evidence={'master_pid':os.getpid(),'worker_pids':pids,'close_called':True,'join_returned':True,'exits':exits}
    emit('finish',{'status':stop_reason or 'partial_cover','pool_evidence':pool_evidence,'seconds':time.monotonic()-began},True);compact();journal.close()
    assert digest({key:state['cache'][key] for key in parent['cache']})==parentcachehash
    assert sha(args.parent_state)==reg['parent_state_sha256'];snapshot=work/'batch001_state.json';atomic(snapshot,as_json(state));assert sha(snapshot)==sha(statepath)
    newpasses=state['passed'][parentpassed:];newsplits=state['splits'][parentsplits:];margins=[F(x['certificate']['minimum_log_margin']) for x in newpasses]
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result={'status':state['status'],'model':'alpha1','epsilon':reg['epsilon'],'target':reg['target'],'global_g_or_C5_claim':False,
      'parent_state_sha256':reg['parent_state_sha256'],'parent_passed_reused':parentpassed,'selected_count':len(selected),'selected_sha256':digest(selected),
      'new_anchor_starts':state['batch_new_starts'],'new_complete_cached_anchors':state['batch_new_completed'],'cache_hits':state['batch_cache_hits'],
      'new_passed':len(newpasses),'new_split':len(newsplits),'total_passed':len(state['passed']),'total_split':len(state['splits']),
      'remaining_queue':len(state['queue']),'remaining_inflight':len(state['inflight']),'resume_jobs':len(state['resume_jobs']),
      'worst_new_passed_log_margin':str(min(margins)) if margins else None,
      'counterexample':state.get('counterexample'),'worker_outcome_histogram':{kind:sum(r['kind']==kind for r in worker_results) for kind in set(r['kind'] for r in worker_results)},
      'worker_kernel_calls':sum(r.get('calls',0) for r in worker_results),'worker_complete_kernel_seconds_sum':sum(a['seconds'] for r in worker_results for a in r.get('record',{}).get('attempts',r.get('attempts',[]))),
      'pool':pool_evidence,'runtime':{'wall_seconds':time.monotonic()-began,'main_peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,
          'max_reported_worker_peak_rss_bytes':max((r.get('peak_rss',0) for r in worker_results),default=0),'full_state_compactions':compactions,'journal_events':events},
      'state':{'path':str(statepath),'snapshot':str(snapshot),'sha256':sha(statepath),'bytes':statepath.stat().st_size,'fingerprint':fingerprint},
      'journal':{'path':str(journalpath),'sha256':sha(journalpath),'bytes':journalpath.stat().st_size,'final_chain_hash':state['journal_hash']},
      'initial_state':{'path':str(initialpath),'sha256':initialsha},
      'tree_sha256':digest({'roots':state['roots'],'splits':state['splits'],'passed':state['passed'],'queue':list(state['queue'])}),
      'remaining_queue_sha256':digest(list(state['queue'])),'source_hashes':sources,'source046_hashes':reg['source046_hashes']}
    args.output_dir.mkdir(parents=True,exist_ok=True);out=args.output_dir/'results.json';atomic(out,result)
    print(json.dumps({k:result[k] for k in ('status','new_anchor_starts','new_passed','new_split','remaining_queue','resume_jobs','pool','runtime','state')},indent=2),flush=True)


if __name__=='__main__':main()
