"""C114 fixed14 private BFS roots,16 slots each/224 total; no time ceiling."""
import sys
sys.dont_write_bytecode=True
import concurrent.futures as cf
import copy,datetime,hashlib,importlib.util,json,os,threading,time,traceback
from pathlib import Path
from fractions import Fraction as Q
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def require(x,s):
    if not x:raise ValueError(s)
def atomic(path,x):
    path=Path(path);tmp=path.with_name(path.name+'.pending')
    with tmp.open('wb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)
    fd=os.open(path.parent,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)
def new(path,x):
    with Path(path).open('xb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

rp=P/'parent_release.json';require(rp.exists(),'C114 PREPARED_NOT_RELEASED; no native/math')
release=read(rp);reg=read(P/'registration.json');inputs=read(P/'inputs.json')
require(release.get('authorized')is True and release.get('round')=='114'
    and release.get('script_sha256')==sha(Path(__file__))
    and release.get('registration_sha256')==sha(P/'registration.json'),'Exact114 release required')
require(release.get('max_roots')==14 and release.get('max_native_per_root')==16
    and release.get('max_native_total')==224 and release.get('max_parallel_native')==14,'Fixed14x16 release')
require(release.get('continuation_owner')=='C114'and release.get('112_frontier_retired_from_dispatch')is True,
    'Parent must transfer/retire112 residual dispatch roots')
require(not any((P/n).exists()for n in ['started.json','checkpoint.json','forest.json','run_receipt.json','roots']),
    'One114 batch only; observe original handles, never restart unknown PIDs')
for name,key in [('inputs.json','inputs_sha256'),('point_box_adapter.py','adapter_sha256'),
                 ('source_manifest.json','source_manifest_sha256')]:
    require(sha(P/name)==reg[key],'Prepared114 source changed: '+name)
start=time.monotonic();startutc=datetime.datetime.now(datetime.timezone.utc).isoformat()
new(P/'started.json',{'round':'114','pid':os.getpid(),'started_utc':startutc,'release_sha256':sha(rp),
    'max_roots':14,'max_native_per_root':16,'max_native_total':224,'program_time_limit':None})
lock=threading.RLock();checker_lock=threading.Lock();live={};root_states={};seen={};launched=set();gaps=[]
def current_head():
    require(sha(inputs['112_forest_path'])==inputs['112_forest_sha256'],'112 forest changed')
    files=sorted(Path(inputs['old_state_dir']).glob('snapshot_????????.json'))
    require(files and sha(files[-1])==inputs['old_head_sha256'],'Old068 head changed')
    return files[-1],sha(files[-1])
def publish(rid,event=None):
    with lock:
        state=root_states[rid];folder=P/'roots'/rid
        if event is not None:
            with(folder/'events.jsonl').open('ab')as f:f.write(canon(event));f.flush();os.fsync(f.fileno())
        atomic(folder/'checkpoint.json',state)
        slim={k:{n:v[n]for n in ['status','slots_used','native_calls','queue','passed_leaves','blocked_leaves','inflight']}
              for k,v in root_states.items()}
        atomic(P/'checkpoint.json',{'round':'114','owner_PID':os.getpid(),'roots':slim,
            'native_launches':len(launched),'max_native_total':224,'release_sha256':sha(rp),
            '112_forest_sha256':inputs['112_forest_sha256'],'global_g_certificate':False})
def observe(rid,node,**fields):
    with lock:
        s=root_states[rid]
        s['last_job_observation']={k:v for k,v in fields.items()if k in ['status','token','native_pid','verification_sha256','error']}
        pid=fields.get('native_pid')
        if pid is not None:
            launched.add(fields['token']if'token'in fields else s['inflight']['token'])
            if s['inflight']is not None:s['inflight']['native_pid']=pid
            require(len(launched)<=224,'Global native launch cap')
        publish(rid)
def gap_record(rid,result,folder):
    row={'root_id':rid,'node_id':result['old_node_id'],'token':result['token'],'status':result['status']}
    if 'center_interval'in result:
        lo,hi=map(Q,result['center_interval']);row.update(center_lower_margin=str(lo-Q(1,65000)),
            center_interval_relative_width=str((hi-lo)/lo))
        cp=folder/'box_certificate.json'
        if cp.exists():
            c=read(cp);worst=min(c['corners'],key=lambda x:Q(x['certificate_lower']))
            b=Q(worst['gradient_dot_lower']);f=Q(worst['minimizing_center_value'])
            curvature=Q(worst['quadratic'])/Q('1/32')**2;grad=b/f
            row.update(whole_box_lower_margin=str(Q(c['tangent_bound_lower'])-Q(1,65000)),
                whole_box_to_center_lower_ratio=str(Q(c['tangent_bound_lower'])/lo),
                worst_corner_bits=worst['corner_bits'],curvature_exponent_penalty=str(curvature),
                gradient_exponent_term=str(grad),robust_prefactor_to_center_lower_ratio=str(f/lo),
                curvature_ge_absolute_gradient_term=curvature>=abs(grad),
                interpretation='certificate-envelope loss only, not actual G or model obstruction')
    with lock:gaps.append(row)
    return row

def run_root(spec):
    rid=spec['id'];folder=P/'roots'/rid;adapter=None;popped_node_id=None
    try:
        adapter=load('c114_adapter_'+rid.replace('.','_'),P/'point_box_adapter.py')
        state=root_states[rid]
        adapter.P=folder;adapter.base=base;adapter.M=M;adapter.checker=checker
        adapter.checker_lock=checker_lock;adapter.lock=lock;adapter.live=live
        adapter.pins=pins;adapter.reg=reg;adapter.release_path=rp;adapter.current_head=current_head
        adapter.plan={'old_attempts_dir':inputs['old_attempts_dir'],'old_head_sha256':inputs['old_head_sha256']}
        adapter.old={'state':state};adapter.status=lambda node,**kw:observe(rid,node,**kw)
        publish(rid,{'kind':'ROOT_STARTED','root_id':rid,'source112_node':spec['source112_id']})
        while True:
            with lock:
                if not state['queue']or state['slots_used']>=16 or state['blocked_leaves']:break
                popped_node_id=state['queue'][0]
                node_id=state['queue'].pop(0);node=state['nodes'][node_id]
                require(node['status']=='queued'and not node['attempts'],'No repeated node attempt')
                expected=base.expected_input(state,node);key=tuple(expected['center'])
                if key in seen:
                    node['status']='blocked';state['blocked_leaves'].append(node_id)
                    state['status']='BLOCKED_DUPLICATE_CENTER_NO_RELAUNCH'
                    publish(rid,{'kind':'DUPLICATE_CENTER_BLOCKED','node_id':node_id,'prior':seen[key]});break
                serial=state['slots_used']+1;token=f'{node_id}.a0.q{serial:03d}'
                seen[key]={'root_id':rid,'node_id':node_id,'token':token,'reserved_not_reissuable':True}
                state['slots_used']+=1;node['status']='inflight'
                state['inflight']={'node_id':node_id,'token':token,'native_pid':None}
                enc=base.encode_center(expected['center'],**pins['fan_bounds'])
                command=[pins['paths']['binary'],'full',pins['paths']['fan'],'1','32','1','2',
                    str(enc['denominator']),*map(str,enc['numerators'])]
                job={'old_node_id':node_id,'continuation_node_id':node_id,'depth':node['depth'],
                    'token':token,'input':expected,'command':command}
                node['attempts'].append({'token':token,'input_sha256':digest(expected),'command':command})
                publish(rid,{'kind':'RESERVE_ONCE','node_id':node_id,'token':token,'slot':serial,
                    'input_sha256':digest(expected),'queue_after_pop':state['queue']})
            result=adapter.run_job(job);jobfolder=folder/'jobs'/token
            posterior=gap_record(rid,result,jobfolder)
            with lock:
                state['native_calls']+=result.get('native_calls',0)
                require(state['native_calls']<=state['slots_used']<=16,'Per-root cap')
                state['inflight']=None;node['attempts'][-1]['result_status']=result['status']
                if result['status']=='WHOLE_BOX_PASSED':
                    node['status']='passed';node['proof']={'token':token,
                        'verification_sha256':sha(jobfolder/'verification.json'),
                        'box_certificate_sha256':sha(jobfolder/'box_certificate.json')}
                    state['passed_leaves'].append(node_id)
                elif result['status']=='WHOLE_BOX_UNRESOLVED_SPLIT':
                    node['status']='split';node['split']=result['split'];node['children']=[n['id']for n in result['children']]
                    for child in result['children']:
                        require(child['id']not in state['nodes']and child['parent']==node_id
                            and child['depth']==node['depth']+1,'New exact child identity')
                        state['nodes'][child['id']]=child
                        if child['status']=='empty':state['empty_leaves'].append(child['id'])
                        else:state['queue'].append(child['id'])
                    require(len(state['nodes'])<=33,'Finite sixteen-call root geometry cap')
                    base.queue_sort(state)
                else:
                    node['status']='blocked';state['blocked_leaves'].append(node_id)
                    state['status']='BLOCKED_CENTER_OR_ERROR_NO_RETRY'
                publish(rid,{'kind':'IMPORT_ORIGINAL_RESULT_ONCE','node_id':node_id,'token':token,
                    'status':result['status'],'native_calls':state['native_calls'],'slots_used':state['slots_used'],
                    'children':node['children'],'posterior_gap':posterior})
        with lock:
            if not state['blocked_leaves']:
                state['status']='ROOT_FULLY_PASSED'if not state['queue']else'ROOT_SLOT_CAP_REACHED'
            state['root_complete']=not state['queue']and not state['blocked_leaves']and state['inflight']is None
            state['residual_frontier']=list(state['queue'])+list(state['blocked_leaves'])
            require(all(n['status']in['split','empty','passed']for n in state['nodes'].values())
                if state['root_complete']else True,'Root completion requires every nonempty leaf passed')
            publish(rid,{'kind':'ROOT_STOP','status':state['status'],'root_complete':state['root_complete'],
                'native_calls':state['native_calls'],'slots_used':state['slots_used'],'residual_frontier':state['residual_frontier']})
        return state
    except BaseException:
        error=traceback.format_exc()
        with lock:
            active=root_states[rid]['inflight']
            proc=None if active is None else live.get(active['node_id'])
        if adapter is not None:adapter.stop_owned(proc)
        with lock:
            state=root_states[rid];state['status']='ROOT_ERROR_STOP_NO_RETRY';state['error']=error
            if state['inflight']is not None:
                node=state['nodes'][state['inflight']['node_id']];node['status']='blocked'
                if node['id']not in state['blocked_leaves']:state['blocked_leaves'].append(node['id'])
            state['inflight']=None
            if popped_node_id in state['nodes']:
                popped=state['nodes'][popped_node_id]
                represented=popped_node_id in state['queue']or popped_node_id in state['blocked_leaves']
                if popped['status']not in['passed','split','empty']and not represented:
                    popped['status']='blocked';state['blocked_leaves'].append(popped_node_id)
            state['root_complete']=False
            state['residual_frontier']=list(state['queue'])+list(state['blocked_leaves']);publish(rid)
        return state

outcome='PREFLIGHT_FAILURE_STOP';failure=None;results=[];pool=None
try:
    for item in read(P/'source_manifest.json')['files']:
        require(sha(item['path'])==item['sha256'],'Source pin changed: '+item['path'])
    current_head();f112=read(inputs['112_forest_path']);handoff=read(inputs['112_handoff_path'])
    require(f112['new_native_calls']==7 and f112['residual_frontier']==handoff['new_residual_frontier'],
        '112 authoritative frontier identity')
    require([s['source112_id']for s in inputs['roots']]==f112['residual_frontier'],'Exactly14 saved child roots')
    require([s['source_node']for s in inputs['roots']]==f112['new_child_nodes'],'Saved112 child geometry unchanged')
    olddir=Path(inputs['old_batch_dir']);pins=read(olddir/'input_manifest.json')
    base=load('c114_frozen_geometry',pins['paths']['C063_core'])
    M=load('c114_frozen_native_contract',olddir/'math_contract.py')
    checker=load('c114_frozen_A057',pins['paths']['A057_checker'])
    import flint
    require(flint.__FLINT_VERSION__=='3.6.0','Pinned FLINT runtime')
    (P/'roots').mkdir()
    for item in inputs['past_centers']:seen[tuple(item['center'])]=item
    for spec in inputs['roots']:
        rid=spec['id'];node=copy.deepcopy(spec['source_node']);node.update(id=rid,parent=None,children=[],attempts=[],status='queued')
        node['source112_id']=spec['source112_id'];(P/'roots'/rid).mkdir()
        root_states[rid]={'root_id':rid,'source112_id':spec['source112_id'],'status':'READY','slots_used':0,
            'native_calls':0,'nodes':{rid:node},'queue':[rid],'passed_leaves':[],'empty_leaves':[],
            'blocked_leaves':[],'inflight':None,'input_manifest_sha256':inputs['original_input_manifest_sha256']}
    pool=cf.ThreadPoolExecutor(max_workers=14)
    futures=[pool.submit(run_root,s)for s in inputs['roots']]
    for fut in cf.as_completed(futures):results.append(fut.result())
    pool.shutdown(wait=True);pool=None;current_head()
    with lock:
        calls=sum(s['native_calls']for s in root_states.values());slots=sum(s['slots_used']for s in root_states.values())
        require(calls<=slots<=224 and len(launched)<=224 and len(results)==14,'Global224/14 invariant')
        require(calls==len(launched),'Observed unique launch ledger disagrees with completed results')
        forest={'round':'114','112_forest_sha256':inputs['112_forest_sha256'],'old068_head_sha256':inputs['old_head_sha256'],
            'inherited_passed_leaf':'r000','root_mapping':{s['source112_id']:s['id']for s in inputs['roots']},
            'roots':root_states,'new_native_calls':calls,'slots_used':slots,'max_native_per_root':16,'max_native_total':224,
            'whole_root_complete':all(s['root_complete']for s in root_states.values()),
            'residual_frontier':[n for s in root_states.values()for n in s['residual_frontier']],
            'closed_coverage_union_preserved':True,'global_g_certificate':False,'new_C6_bound':False}
        new(P/'forest.json',forest);new(P/'posterior_gaps.json',{'round':'114','jobs':gaps,
            'scope':'exact posterior envelope decomposition, no true-model exclusion'})
    outcome='COMPLETED_FOURTEEN_FINITE_ROOTS'
except BaseException:
    failure=traceback.format_exc()
finally:
    with lock:owned=list(live.values())
    for proc in owned:
        if proc.poll()is None:
            proc.terminate()
            try:proc.wait(timeout=5)
            except Exception:proc.kill();proc.wait()
    if pool is not None:pool.shutdown(wait=True)
    receipt={'round':'114','status':outcome,'pid':os.getpid(),'started_utc':startutc,
        'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-start,
        'program_time_limit':None,'observed_unique_native_launches':len(launched),'roots_returned':len(results),
        'native_calls_by_root':{k:s['native_calls']for k,s in root_states.items()},
        'slots_by_root':{k:s['slots_used']for k,s in root_states.items()},'automatic_retries':0,
        '112_or_old_native_calls_replayed':0,'release_sha256':sha(rp),'failure':failure,'global_g_certificate':False}
    atomic(P/'run_receipt.json',receipt)
    print(json.dumps({'round':'114','status':outcome,'pid':os.getpid(),'native_launches':len(launched),
        'roots_returned':len(results),'elapsed_seconds':receipt['elapsed_seconds'],'global_g_certificate':False}))
if outcome!='COMPLETED_FOURTEEN_FINITE_ROOTS':raise SystemExit(1)
