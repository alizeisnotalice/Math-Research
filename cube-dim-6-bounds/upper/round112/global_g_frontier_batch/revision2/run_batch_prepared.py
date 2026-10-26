"""C112: one reviewed seven-job continuation forest, no old state mutations."""
import sys
sys.dont_write_bytecode=True
import concurrent.futures as cf
import copy,datetime,hashlib,importlib.util,json,math,os,signal,subprocess,threading,time,traceback
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def canon(x):return (json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode()
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
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

release_path=P/'parent_release.json'
require(release_path.exists(),'C112 PREPARED_NOT_RELEASED; refusing native/math')
release=read(release_path);reg=read(P/'registration.json');plan=read(P/'inputs.json')
require(release.get('authorized')is True and release.get('round')=='112'
    and release.get('script_sha256')==sha(Path(__file__))
    and release.get('registration_sha256')==sha(P/'registration.json'), 'C112 release mismatch')
require(release.get('continuation_owner')=='C112'and release.get('old_queue_retired_from_dispatch')is True,
    'Parent must own/retire seven old queue dispatch IDs before new launch')
require(release.get('max_new_nodes')==7 and release.get('max_native_calls')==7,'Fixed7/7 release required')
require(not any((P/n).exists()for n in ['started.json','checkpoint.json','forest.json','run_receipt.json']),
    'C112 one invocation only; observe original job handles, never relaunch')
require(sha(P/'inputs.json')==reg['inputs_sha256']and sha(P/'source_manifest.json')==reg['source_manifest_sha256'],
    'Prepared input bytes changed')
require(sha(P/'old_tree_geometry.json')==reg['old_tree_geometry_sha256'],'Prepared old coverage geometry changed')
start=time.monotonic();startutc=datetime.datetime.now(datetime.timezone.utc).isoformat()
new(P/'started.json',{'round':'112','pid':os.getpid(),'started_utc':startutc,
    'release_sha256':sha(release_path),'max_jobs':7,'program_time_limit':None})
lock=threading.RLock();checker_lock=threading.Lock();live={};states={j['old_node_id']:{'status':'NOT_STARTED','token':j['token']}for j in plan['jobs']}
def checkpoint():
    with lock:atomic(P/'checkpoint.json',{'round':'112','owner_PID':os.getpid(),'jobs':copy.deepcopy(states),
        'old_head_sha256':plan['old_head_sha256'],'release_sha256':sha(release_path),
        'native_launches':sum('native_pid'in x for x in states.values()),'global_g_certificate':False})
def status(node,**kw):
    with lock:states[node].update(kw);checkpoint()
def current_head():
    files=sorted(Path(plan['old_state_dir']).glob('snapshot_????????.json'))
    require(bool(files),'Old state disappeared')
    return files[-1],sha(files[-1])
def stop_owned(proc):
    if proc is not None and proc.poll()is None:
        proc.terminate()
        try:proc.wait(timeout=5)
        except subprocess.TimeoutExpired:proc.kill();proc.wait()

def validate_native(m,expected):
    # Same accepted pointjet physics as C068; no obsolete timing assertion.
    require(m.get('status')=='all_real_nodes_supported'and m.get('precision_bits')==512,'Native status/precision')
    require(m.get('full_group_terms')==967680 and m.get('computed_sign_pairs')==483840,'Incomplete body32 group')
    for k in ['confluent_terms','zero_node_terms']:
        v=m.get(k);require(type(v)is int and 0<=v<=967680 and v%2==0,'Native count '+k)
    require(type(m.get('unique_node_values'))is int and 1<=m['unique_node_values']<=2903040,'Distinct-node count')
    require(isinstance(m.get('kernel_seconds'),(int,float))and not isinstance(m['kernel_seconds'],bool)
        and math.isfinite(m['kernel_seconds'])and m['kernel_seconds']>=0,'Native elapsed metadata')
    require(isinstance(m.get('gradient'),list)and len(m['gradient'])==6,'Six physical derivatives')
    F=M.dyadic(m['band']);g=[M.dyadic(d)for d in m['gradient']]
    pa=M.dyadic(m['p_a']);pt=M.dyadic(m['p_T']);ne=M.dyadic(m['Euler_xgrad_plus5B_plusTpT_minusapa'])
    require(F[0]>0 and pa[0]>0 and pt[0]>0,'Positive band/endpoints')
    require(ne[0]<=0<=ne[1],'Native Euler')
    def scale(i,q):return tuple(sorted(v*q for v in i))
    terms=[scale(F,Q(5)),scale(pt,Q('1/2')),scale(pa,-Q('1/32'))]
    terms += [scale(i,Q(x))for i,x in zip(g,expected['center'])]
    e=tuple(sum(t[k]for t in terms)for k in [0,1]);require(e[0]<=0<=e[1],'Exact-rational Euler')
    return F,g,pa,pt,ne,e

def run_job(job):
    node=job['old_node_id'];folder=P/'jobs'/job['token'];proc=None;jstart=time.monotonic()
    try:
        # Detect an orphan legacy launch even when not yet imported into a snapshot.
        legacy=list(Path(plan['old_attempts_dir']).glob(node+'.a*'))
        if legacy:
            out={'status':'SKIPPED_LEGACY_ATTEMPT_EXISTS','legacy_attempts':list(map(str,legacy)),
                'node_id':node,'token':job['token'],'native_calls':0,'global_g_certificate':False}
            status(node,**out);return out
        require(current_head()[1]==plan['old_head_sha256'],'Old head changed before unique launch')
        folder.mkdir(parents=True,exist_ok=False)
        expected=base.expected_input(old['state'],old['state']['nodes'][node])
        require(expected==job['input'],'Prepared midpoint/box/window differ from old exact geometry')
        encoded=base.encode_center(expected['center'],**pins['fan_bounds'])
        command=[pins['paths']['binary'],'full',pins['paths']['fan'],'1','32','1','2',
            str(encoded['denominator']),*map(str,encoded['numerators'])]
        require(command==job['command'],'Prepared native argv mismatch')
        intent={'round':'112','token':job['token'],'old_node_id':node,'continuation_node_id':job['continuation_node_id'],
            'input':expected,'input_sha256':digest(expected),'old_head_sha256':plan['old_head_sha256'],
            'release_sha256':sha(release_path),'command':command,'wrapper_pid':os.getpid(),
            'source_binary_fan_runtime_hashes':reg['native_pins'],'program_time_limit':None}
        new(folder/'intent.json',intent);status(node,status='INTENT_SAVED',intent_sha256=sha(folder/'intent.json'))
        with (folder/'native.json').open('xb')as out,(folder/'stderr.txt').open('xb')as err:
            proc=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=out,stderr=err,cwd=folder)
            with lock:live[node]=proc
            new(folder/'original_process.json',{'token':job['token'],'native_pid':proc.pid,'wrapper_pid':os.getpid(),
                'command':command,'input_sha256':digest(expected),'intent_sha256':sha(folder/'intent.json')})
            status(node,status='NATIVE_RUNNING',native_pid=proc.pid,command=command)
            code=proc.wait()
        with lock:live.pop(node,None)
        receipt={'schema':'C112-original-native-receipt-v1','token':job['token'],'old_node_id':node,
            'input_sha256':digest(expected),'intent_sha256':sha(folder/'intent.json'),
            'release_sha256':sha(release_path),'owned_native_pid':proc.pid,'owned_wrapper_pid':os.getpid(),
            'original_process_sha256':sha(folder/'original_process.json'),'command':command,
            'native_exit_code':code,'native_sha256':sha(folder/'native.json'),'stderr_sha256':sha(folder/'stderr.txt'),
            'stderr_empty':(folder/'stderr.txt').stat().st_size==0,'source_binary_fan_runtime_hashes':reg['native_pins'],
            'elapsed_to_native_exit':time.monotonic()-jstart,'timeout_seconds':None,'retries':0}
        new(folder/'native_receipt.json',receipt)
        require(code==0 and receipt['stderr_empty'],'Native error/stderr; stop this job')
        m=read(folder/'native.json');F,g,pa,pt,ne,e=validate_native(m,expected)
        jet={**{k:expected[k]for k in ['model','model_definition','window','center']},
            'F_interval':list(map(str,F)),'face_gradient_intervals':[list(map(str,a))for a in g[1:]],
            'source_proof':f'C112_ORIGINAL_NATIVE:{folder}/native.json#sha256={receipt["native_sha256"]};receipt={sha(folder/"native_receipt.json")}',
            'provenance':{'kind':'C016_NATIVE_SAME_WINDOW_POSITIVE_BAND','precision_bits':512,'full_group_terms':967680,
                'computed_sign_pairs':483840,'synthetic_fixture':False,'native_sha256':receipt['native_sha256'],
                'receipt_sha256':sha(folder/'native_receipt.json'),'gradient_gate_sha256':pins['F061_manifest_sha256'],
                'independent_full_gradient_crosschecked_for_this_new_attempt':False}}
        base.validate_jet(jet,expected,test_only=False);new(folder/'jet.json',jet)
        result={'round':'112','old_node_id':node,'continuation_node_id':job['continuation_node_id'],
            'token':job['token'],'input_sha256':digest(expected),'center_interval':list(map(str,F)),
            'six_physical_gradient_intervals':[list(map(str,v))for v in g],
            'endpoint_intervals':{'p_a':list(map(str,pa)),'p_T':list(map(str,pt))},
            'native_Euler':list(map(str,ne)),'independent_exact_Euler':list(map(str,e)),
            'native_counts':{k:m[k]for k in ['full_group_terms','computed_sign_pairs','confluent_terms','zero_node_terms','unique_node_values']},
            'native_pid':proc.pid,'native_calls':1,'children':[],'global_g_certificate':False,
            'new_C6_bound':False,'independent_gradient_crosscheck_for_this_new_jet':False}
        if F[0]>base.TARGET:
            # Native processes are parallel; Python/Arb corner checks serialize.
            with checker_lock:
                cert=checker.certify(base.certificate_input(jet,expected,test_only=False),precision_bits=384)
                passed=base.check_saved_box(jet,expected,cert)
            new(folder/'box_certificate.json',cert)
            result.update(center_classification='CENTER_STRICTLY_ABOVE',corners_checked=32,
                box_certificate_sha256=sha(folder/'box_certificate.json'),whole_box_lower=cert['tangent_bound_lower'],
                whole_box_upper=cert['tangent_bound_upper'])
            if passed:result['status']='WHOLE_BOX_PASSED'
            elif any(lo!=hi for lo,hi in base.box(expected['face_box'])):
                split=base.split_box(expected['face_box']);children=[]
                for suffix,raw in [('0',split['lower_raw_box']),('1',split['upper_raw_box'])]:
                    child=base.make_node(job['continuation_node_id']+suffix,job['continuation_node_id'],4,raw,0)
                    child['old_frontier_root']=node;children.append(child)
                result.update(status='WHOLE_BOX_UNRESOLVED_SPLIT',split=split,children=children)
            else:result['status']='WHOLE_BOX_UNRESOLVED_SINGLETON_BLOCKED'
        else:
            result.update(status='FIXED_WINDOW_CENTER_BLOCKED',corners_checked=0,
                center_classification='CENTER_STRICTLY_BELOW'if F[1]<base.TARGET else'CENTER_OVERLAP_OR_EQUALITY')
        new(folder/'verification.json',result)
        status(node,**result,verification_sha256=sha(folder/'verification.json'))
        return result
    except Exception:
        stop_owned(proc)
        with lock:live.pop(node,None)
        failure={'status':'OWN_JOB_ERROR_STOP_NO_RETRY','old_node_id':node,'token':job['token'],
            'error':traceback.format_exc(),'native_pid':None if proc is None else proc.pid,
            'native_exit_code':None if proc is None else proc.returncode,'global_g_certificate':False,
            'native_calls':int(proc is not None),'elapsed_seconds':time.monotonic()-jstart}
        if folder.exists():atomic(folder/'failure.json',failure)
        status(node,**failure);return failure

outcome='PREFLIGHT_FAILURE_STOP';failure=None;results=[]
try:
    # No controller/old writer is called. Only accepted read-only typed recovery.
    for item in read(P/'source_manifest.json')['files']:
        require(sha(item['path'])==item['sha256'],'Source pin changed: '+item['path'])
    olddir=Path(plan['old_batch_dir']);sys.path.insert(0,str(olddir))
    V=load('c112_accepted_c068_core',olddir/'core_v2.py')
    C=load('c112_accepted_c068_context',olddir/'batch_context.py')
    M=load('c112_accepted_native_math',olddir/'math_contract.py')
    context=C.Context();old=V.recover(olddir/'state',context);base=context.base;checker=context.checker;pins=context.pins
    require(old['sequence']==25 and old['sha256']==plan['old_head_sha256'],'Old typed head changed')
    require(old['state']['inflight']is None and old['state']['halted']is None,'Old batch inflight/halt exists')
    require(old['state']['batch']['native_reservations']==8 and old['state']['production_kernel_calls']==8,'Old8-call ledger changed')
    require(old['state']['queue']==plan['old_frontier']and old['state']['passed_leaves']==['r000'],'Old leaf partition changed')
    require(all(not old['state']['nodes'][n]['attempts']for n in plan['old_frontier']),'Previously started frontier; stop without refreshing slots')
    geom=read(P/'old_tree_geometry.json')
    keys=['id','parent','children','depth','status','raw_box','box','tightening','cutoff_index','split']
    require({n['id']:n for n in geom['nodes']}=={node_id:{k:v for k,v in n.items()if k in keys}for node_id,n in old['state']['nodes'].items()},
        'Prepared inherited coverage geometry differs from authenticated old tree')
    require(geom['passed_leaf_proof']==old['state']['nodes']['r000']['proof'],'Inherited r000 proof changed')
    centers=[tuple(j['input']['center'])for j in plan['jobs']]
    old_centers={tuple(a['input']['center'])for n in old['state']['nodes'].values()for a in n['attempts']}
    require(len(set(centers))==7 and not set(centers)&old_centers,'New centers repeat each other or old8 calls')
    new(P/'old_recovery.json',{'sequence':old['sequence'],'sha256':old['sha256'],'inflight':None,
        'queue':old['state']['queue'],'passed_leaves':old['state']['passed_leaves'],
        'old_native_calls':8,'old_native_replayed':0,'typed_recovery_saved_math_only':True})
    (P/'jobs').mkdir(exist_ok=False);checkpoint()
    with cf.ThreadPoolExecutor(max_workers=7)as pool:
        futures=[pool.submit(run_job,j)for j in plan['jobs']]
        for fut in cf.as_completed(futures):results.append(fut.result())
    require(current_head()[1]==plan['old_head_sha256'],'Old head changed during continuation; do not merge')
    results.sort(key=lambda r:r.get('old_node_id',r.get('node_id','')))
    passed=['r000'];residual=[];new_children=[]
    for r in results:
        node=r.get('old_node_id',r.get('node_id'))
        if r['status']=='WHOLE_BOX_PASSED':passed.append('C112.'+node)
        elif r['status']=='WHOLE_BOX_UNRESOLVED_SPLIT':
            new_children+=r['children'];residual += [c['id']for c in r['children']if c['status']!='empty']
        else:residual.append('C112.'+node)
    native_calls=sum(r['native_calls']for r in results)
    require(native_calls<=7 and len(results)==7,'Fixed seven-slot ledger violation')
    forest={'round':'112','old_head_sha256':plan['old_head_sha256'],'old_sequence':25,
        'old_tree_geometry_sha256':sha(P/'old_tree_geometry.json'),'old_root_domain':'[0,1]^5 intersect descending order',
        'inherited_passed_leaf':'r000','old_frontier_to_continuation':{j['old_node_id']:j['continuation_node_id']for j in plan['jobs']},
        'jobs':results,'new_child_nodes':new_children,'passed_leaves':passed,'residual_frontier':residual,
        'old_native_calls':8,'new_native_calls':native_calls,'all_nonempty_root_leaves_passed':not residual,
        'coverage_union_preserved_by_exact_closed_splits':True,'global_g_certificate':False,'new_C6_bound':False,
        'parent_must_take_over_frontier_before_future_dispatch':True}
    new(P/'forest.json',forest);outcome='COMPLETED_SEVEN_FIXED_SLOTS'
except BaseException:
    failure=traceback.format_exc()
finally:
    with lock:owned=list(live.values())
    for proc in owned:stop_owned(proc)
    receipt={'round':'112','status':outcome,'pid':os.getpid(),'started_utc':startutc,
        'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-start,
        'program_time_limit':None,'native_launches':sum('native_pid'in x and x['native_pid']is not None for x in states.values()),
        'results_received':len(results),'states':states,'release_sha256':sha(release_path),'failure':failure,
        'automatic_retries':0,'old_native_calls_replayed':0,'old_state_modified':False,'global_g_certificate':False}
    atomic(P/'run_receipt.json',receipt)
    print(json.dumps({'round':'112','status':outcome,'pid':os.getpid(),'native_launches':receipt['native_launches'],
        'results_received':len(results),'elapsed_seconds':receipt['elapsed_seconds'],'global_g_certificate':False}))
if outcome!='COMPLETED_SEVEN_FIXED_SLOTS':raise SystemExit(1)
