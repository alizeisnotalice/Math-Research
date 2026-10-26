"""C117:96 private roots,8 new single-window kernels/root,16 nodes/root,14 workers."""
import sys
sys.dont_write_bytecode=True
import concurrent.futures as cf
import copy,datetime,hashlib,importlib.util,json,os,subprocess,threading,time,traceback
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def require(x,s):
    if not x:raise ValueError(s)
def load(n,p):
    s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def resolve(ref):
    x=read(ref['path'])
    for k in ref.get('selector',[]):x=x[k]
    return x
rp=P/'parent_release.json';require(rp.exists(),'C117 PREPARED_NOT_RELEASED')
rel=read(rp);reg=read(P/'registration.json');inp=read(P/'inputs.json');cat=read(P/'cache_catalog.json')
require(rel.get('authorized')is True and rel.get('round')=='117'and rel.get('script_sha256')==sha(Path(__file__))
    and rel.get('registration_sha256')==sha(P/'registration.json'),'Exact117 release')
require([rel.get(k)for k in['roots','max_new_native_per_root','max_nodes_per_root','max_new_native_total','parallel']]==[96,8,16,768,14], 'Fixed96/8/16/768/14')
require(rel.get('continuation_owner')=='C117'and rel.get('old_frontier_retired_from_dispatch')is True,'Parent dispatch takeover')
require(not any((P/n).exists()for n in['started.json','roots','new_native','cache_index.json','run_receipt.json']), 'Never restart original/unknown117 work')
for n,k in [('inputs.json','inputs_sha256'),('cache_catalog.json','catalog_sha256'),('source_manifest.json','source_manifest_sha256'),('math_helpers.py','helpers_sha256'),('native_validator.py','validator_sha256')]:require(sha(P/n)==reg[k],'Prepared bytes changed:'+n)
H=load('c117_helpers',P/'math_helpers.py')
tic=time.monotonic();utc=datetime.datetime.now(datetime.timezone.utc).isoformat()
H.new(P/'started.json',{'round':'117','PID':os.getpid(),'utc':utc,'release_sha256':sha(rp),'time_limit':None})
lock=threading.RLock();cert_guard=threading.Lock();check_lock=threading.Lock();cache={};certcache={};sumcache={};live={};roots={};new_reserved=0;started=set();new_cert_count=0
def nk(center,w):return hashlib.sha256(canon({'center':[str(Q(x))for x in center],'a':str(Q(w['a'])),'T':str(Q(w['T'])),'engineSHA':cat['engineSHA']})).hexdigest()
def ck(k,box):return hashlib.sha256(canon({'native_key':k,'box':[[str(Q(x))for x in b]for b in box],'target':'1/65000','quadratic':'frame32'})).hexdigest()
def sk(keys,box):return hashlib.sha256(canon({'keys':keys,'box':[[str(Q(x))for x in b]for b in box]})).hexdigest()
def command(e):
    z=base.encode_center(e['center'],**inp['fan_bounds']);a=Q(e['window']['a']);T=Q(e['window']['T'])
    return[inp['native_paths']['binary'],'full',inp['native_paths']['fan'],str(a.numerator),str(a.denominator),str(T.numerator),str(T.denominator),str(z['denominator']),*map(str,z['numerators'])]
def expected(node,w):return {'input_schema':'C117-general-positive-band-point-v1','model':base.MODEL,'model_definition':base.MODEL_DEFINITION,'center':base.midpoint(node['box']),'face_box':node['box'],'window':w,'target':inp['target'],'quadratic':'frame32'}
def flush_cache():
    H.atomic(P/'cache_index.json',{'engineSHA':cat['engineSHA'],'new_reserved':new_reserved,'new_started':len(started),
        'entries':{k:{a:b for a,b in v.items()if not a.startswith('_')}for k,v in cache.items()},'no_reissue_unknown':True})
def publish(rid,event=None):
    with lock:
        s=roots[rid];d=P/'roots'/rid
        if event is not None:
            with(d/'events.jsonl').open('ab')as f:f.write(H.canon(event));f.flush();os.fsync(f.fileno())
        H.atomic(d/'checkpoint.json',s)
        H.atomic(P/'checkpoint.json',{'round':'117','PID':os.getpid(),'new_reserved':new_reserved,'new_started':len(started),'new_certificates':new_cert_count,
            'roots':{r:{k:s[k]for k in['status','new_reserved','nodes_processed','queue','passed_leaves','blocked_leaves','inflight']}for r,s in roots.items()},'global_g_certificate':False})
def pins_of(r):
    if 'source_binary_fan_runtime_hashes'in r:return r['source_binary_fan_runtime_hashes']
    return {inp['native_paths']['source']:r['source_sha256'],inp['native_paths']['binary']:r['binary_sha256'],inp['native_paths']['fan']:r['fan_sha256'],**r['runtime_sha256']}
def authenticate(v):
    e=v['expected'];r=read(v['receipt_path']);i=read(v['intent_path']);n=read(v['native_path'])
    require(i['input']==e and r['input_sha256']==H.digest(e)and r['command']==command(e),'Cache input/command')
    require(r['native_exit_code']==0 and Path(v['stderr_path']).read_bytes()==b''and r['native_sha256']==sha(v['native_path'])and r['stderr_sha256']==sha(v['stderr_path']),'Cache original bytes/exit')
    require(pins_of(r)==reg['native_pins'],'Cache engine/fan/runtime pins')
    data=V.validate_native(n,e);jet=resolve(v['jet_ref'])
    require(jet['F_interval']==list(map(str,data[0]))and jet['face_gradient_intervals']==[list(map(str,g))for g in data[1][1:]],'Cache jet/raw mismatch')
    base.validate_jet(jet,e,test_only=False);return data,jet
def get_data(rid,node,e):
    global new_reserved
    k=nk(e['center'],e['window']);owner=False
    with lock:
        if k not in cache:
            if roots[rid]['new_reserved']>=8 or new_reserved>=768:return None,'NATIVE_BUDGET_PARTIAL_PREFIX',k
            roots[rid]['new_reserved']+=1;new_reserved+=1
            token='C117.'+k;event=threading.Event()
            cache[k]={'key':k,'origin':'117','status':'RESERVED_LAUNCH_UNKNOWN','owner_root':rid,'token':token,'expected':e,'_event':event}
            flush_cache();owner=True
        v=cache[k];event=v['_event']
        roots[rid]['inflight']={'node_id':node['id'],'cache_key':k,'kind':'NEW_NATIVE'if owner else'CACHE_REFERENCE','owner_root':v.get('owner_root'),'native_pid':v.get('native_pid')}
        publish(rid,{'kind':'POINT_WINDOW_LOOKUP','node':node['id'],'key':k,'new_owner':owner})
    if not owner:
        event.wait()
        with lock:return (v if v['status']=='AVAILABLE'else None),('AVAILABLE'if v['status']=='AVAILABLE'else'FAILED_OR_UNKNOWN_CACHE_NO_RETRY'),k
    proc=None;d=P/'new_native'/k
    try:
        d.mkdir();cmd=command(e);H.new(d/'intent.json',{'round':'117','token':v['token'],'input':e,'input_sha256':H.digest(e),'command':cmd,'release_sha256':sha(rp),'owner_root':rid,'source_binary_fan_runtime_hashes':reg['native_pins']})
        with(d/'native.json').open('xb')as out,(d/'stderr.txt').open('xb')as err:
            proc=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=out,stderr=err,cwd=d)
            with lock:
                live[k]=proc;started.add(k);v.update(status='NATIVE_RUNNING',native_pid=proc.pid)
                H.new(d/'original_process.json',{'token':v['token'],'native_pid':proc.pid,'wrapper_pid':os.getpid(),'command':cmd,'intent_sha256':sha(d/'intent.json')});flush_cache();publish(rid)
            code=proc.wait()
        with lock:live.pop(k,None)
        r={'schema':'C117-one-window-native-v1','token':v['token'],'input_sha256':H.digest(e),'command':cmd,'native_exit_code':code,'owned_native_pid':proc.pid,'owned_wrapper_pid':os.getpid(),'native_sha256':sha(d/'native.json'),'stderr_sha256':sha(d/'stderr.txt'),'stderr_empty':(d/'stderr.txt').stat().st_size==0,'intent_sha256':sha(d/'intent.json'),'source_binary_fan_runtime_hashes':reg['native_pins'],'release_sha256':sha(rp),'one_window_kernel_evaluations':1,'time_limit':None}
        H.new(d/'native_receipt.json',r);require(code==0 and r['stderr_empty'],'Native error, no retry')
        data=V.validate_native(read(d/'native.json'),e)
        jet=H.jet_from(e,data[0],data[1],f'C117_ORIGINAL_NATIVE:{d}/native.json#sha256={r["native_sha256"]};receipt={sha(d/"native_receipt.json")}',{'kind':'C016_NATIVE_SAME_WINDOW_POSITIVE_BAND','precision_bits':512,'full_group_terms':967680,'computed_sign_pairs':483840,'synthetic_fixture':False,'native_sha256':r['native_sha256'],'receipt_sha256':sha(d/'native_receipt.json'),'gradient_gate_sha256':inp['F061_manifest_sha256'],'independent_full_gradient_crosschecked_for_this_new_attempt':False})
        base.validate_jet(jet,e,test_only=False);H.new(d/'jet.json',jet)
        with lock:
            v.update(status='AVAILABLE',native_path=str(d/'native.json'),receipt_path=str(d/'native_receipt.json'),stderr_path=str(d/'stderr.txt'),intent_path=str(d/'intent.json'),process_path=str(d/'original_process.json'),jet_ref={'path':str(d/'jet.json'),'selector':[]},_data=data,_jet=jet)
            flush_cache();event.set();publish(rid)
        return v,'AVAILABLE',k
    except BaseException:
        error=traceback.format_exc();H.stop_owned(proc)
        with lock:
            live.pop(k,None);v.update(status='FAILED_NO_RETRY',error=error,native_pid=None if proc is None else proc.pid);event.set()
            if d.exists():H.atomic(d/'failure.json',{'key':k,'error':error,'native_pid':v['native_pid'],'exit':None if proc is None else proc.returncode})
            flush_cache();publish(rid)
        return None,'NATIVE_ERROR_PARTIAL_PREFIX',k
def prefix(keys,e):
    ds=[cache[k]['_data']for k in keys];F=ds[0][0];g=ds[0][1];inter=[]
    for j in range(1,len(ds)):
        b=(max(ds[j][3][0],ds[j-1][2][0]),min(ds[j][3][1],ds[j-1][2][1]));require(b[0]<=b[1],'Adjacent density mismatch')
        inter.append(list(map(str,b)));F=H.add(F,ds[j][0]);g=[H.add(a,b)for a,b in zip(g,ds[j][1])]
    ee=copy.deepcopy(e);ee['window']={'kind':'B','a':cache[keys[-1]]['expected']['window']['a'],'T':'1/2'}
    E=H.euler(F,g,ds[-1][2],ds[0][3],ee)
    return {'keys':keys.copy(),'F_interval':list(map(str,F)),'six_gradients':[list(map(str,v))for v in g],'outer_window':ee['window'],'Euler_interval':list(map(str,E)),'internal_density_intersections':inter},F
def certify_band(k,node):
    global new_cert_count
    ckey=ck(k,node['box'])
    with cert_guard:
        if ckey in certcache:
            v=certcache[ckey];c=resolve(v['certificate_ref']);return c,{'cache':'REUSED','key':ckey,**v['certificate_ref']}
        v=cache[k];e=copy.deepcopy(v['expected']);e['face_box']=node['box'];c=H.box_certificate(e,v['_jet'])
        d=P/'new_certificates'/ckey[:2];d.mkdir(exist_ok=True);path=d/(ckey+'.json');H.new(path,c)
        ref={'path':str(path),'selector':[]};certcache[ckey]={'key':ckey,'native_key':k,'expected':e,'certificate_ref':ref,'origin':'117'}
        with lock:new_cert_count+=1
        return c,{'cache':'NEW','key':ckey,**ref}
def whole_box(keys,node):
    if sk(keys,node['box'])in sumcache:
        s=sumcache[sk(keys,node['box'])];return tuple(map(Q,s['interval'])),{'cache':'REUSED115_BANDSET','source_path':s['source_path'],'keys':keys}
    lower=Q(0);upper=Q(0);refs=[]
    for k in keys:
        c,ref=certify_band(k,node);lower+=Q(c['tangent_bound_lower']);upper+=Q(c['tangent_bound_upper']);refs.append(ref)
    return(lower,upper),{'cache':'PER_BAND_SUM','certificates':refs,'keys':keys}
def run_root(spec):
    rid=spec['id'];popped=None;s=roots[rid]
    try:
        while s['queue']and s['nodes_processed']<16 and not s['blocked_leaves']:
            with lock:
                popped=s['queue'][0];node=s['nodes'][s['queue'].pop(0)];s['nodes_processed']+=1;node['status']='inflight'
                s['inflight']={'node_id':node['id'],'kind':'PREFIX'};publish(rid,{'kind':'PROCESS_NODE','node':node['id'],'ordinal':s['nodes_processed']})
            keys=[];finished=False
            for wi,w in enumerate(inp['windows']):
                e=expected(node,w);v,reason,k=get_data(rid,node,e)
                if v is None:
                    with lock:node.update(status='blocked',stop_reason=reason,next_window=w);s['blocked_leaves'].append(node['id']);s['status']=reason
                    finished=True;break
                keys.append(k);pre,F=prefix(keys,e)
                with lock:node['partial_prefix']=pre;publish(rid,{'kind':'PREFIX_AVAILABLE','node':node['id'],'keys':keys,'center_class':('ABOVE'if F[0]>Q(inp['target'])else'BELOW'if F[1]<Q(inp['target'])else'OVERLAP')})
                if F[0]>Q(inp['target']):
                    bound,refs=whole_box(keys,node)
                    with lock:
                        node['whole_box_interval']=list(map(str,bound));node['certificate_refs']=refs
                        if bound[0]>=Q(inp['target']):
                            node['status']='passed';s['passed_leaves'].append(node['id'])
                        else:
                            split=base.split_box(node['box']);kids=[base.make_node(node['id']+suffix,node['id'],node['depth']+1,raw,0)for suffix,raw in [('0',split['lower_raw_box']),('1',split['upper_raw_box'])]]
                            node.update(status='split',split=split,children=[c['id']for c in kids])
                            for c in kids:
                                require(c['id']not in s['nodes'],'Duplicate geometric child');s['nodes'][c['id']]=c
                                if c['status']=='empty':s['empty_leaves'].append(c['id'])
                                else:s['queue'].append(c['id'])
                            base.queue_sort(s)
                        finished=True
                    break
                if F[1]<Q(inp['target'])and wi+1<len(inp['windows']):continue
                with lock:
                    reason='MINIMUM_WINDOW_CENTER_BELOW'if F[1]<Q(inp['target'])else'CENTER_OVERLAP_OR_EQUALITY'
                    node.update(status='blocked',stop_reason=reason);s['blocked_leaves'].append(node['id']);s['status']=reason
                finished=True;break
            require(finished,'Unrepresented node prefix')
            with lock:s['inflight']=None;publish(rid,{'kind':'NODE_STOP','node':node['id'],'status':node['status'],'prefix_keys':keys})
        with lock:
            if not s['blocked_leaves']:s['status']='ROOT_FULLY_PASSED'if not s['queue']else'NODE_PROCESS_CAP_REACHED'
            s['root_complete']=not s['queue']and not s['blocked_leaves'];s['residual_frontier']=s['queue']+s['blocked_leaves'];s['inflight']=None
            publish(rid,{'kind':'ROOT_STOP','status':s['status'],'root_complete':s['root_complete']})
        return s
    except BaseException:
        with lock:
            s['status']='ROOT_ERROR_STOP';s['error']=traceback.format_exc()
            for entry in cache.values():
                if entry.get('owner_root')==rid and entry['status']!='AVAILABLE'and not entry['_event'].is_set():
                    entry.update(status='FAILED_NO_RETRY',error=s['error'],unknown_original_token_retained=True)
                    entry['_event'].set()
            if popped in s['nodes']:
                n=s['nodes'][popped]
                if n['status']not in['passed','split','empty']:
                    n['status']='blocked'
                    if popped in s['queue']:s['queue'].remove(popped)
                    if popped not in s['blocked_leaves']:s['blocked_leaves'].append(popped)
            s.update(inflight=None,root_complete=False,residual_frontier=s['queue']+s['blocked_leaves']);flush_cache();publish(rid)
        return s
status='PREFLIGHT_FAILURE_STOP';failure=None;results=[];pool=None
try:
    for f in read(P/'source_manifest.json')['files']:require(sha(f['path'])==f['sha256'],'Source changed:'+f['path'])
    require(len(inp['roots'])==96 and len(inp['inherited_passed_leaves'])==42,'Root/pass inventory')
    require(sha(inp['source114_forest_path'])==inp['source114_forest_sha256'],'114 forest changed')
    require(sha(sorted(Path(inp['old068_state_dir']).glob('snapshot_????????.json'))[-1])==inp['old068_head_sha256'],'Old068 head changed')
    base=H.load('c117_base',inp['base_path']);M=H.load('c117_dyadics',inp['math_contract_path']);checker=H.load('c117_A057',inp['checker_path']);V=H.load('c117_validator',P/'native_validator.py');V.M=M
    H.base=base;H.checker=checker;H.checker_lock=check_lock
    import flint
    require(flint.__FLINT_VERSION__=='3.6.0','FLINT runtime')
    for row in cat['native_entries']:
        require(row['key']==nk(row['expected']['center'],row['expected']['window'])and row['key']not in cache,'Cache key duplicate/mismatch')
        data,jet=authenticate(row);ev=threading.Event();ev.set();cache[row['key']]={**row,'status':'AVAILABLE','_event':ev,'_data':data,'_jet':jet}
    require(len(cache)==183,'Exactly183 old native windows')
    for row in cat['band_certificates']:
        require(row['key']==ck(row['native_key'],row['expected']['face_box']),'Certificate key');c=resolve(row['certificate_ref'])
        base.check_saved_box(cache[row['native_key']]['_jet'],row['expected'],c);certcache[row['key']]=row
    for row in cat['bandset_certificates115']:sumcache[sk(row['keys'],row['box'])]=row
    (P/'roots').mkdir();(P/'new_native').mkdir();(P/'new_certificates').mkdir();flush_cache()
    for spec in inp['roots']:
        rid=spec['id'];d=P/'roots'/rid;d.mkdir();n=copy.deepcopy(spec['source_node']);n.update(id=rid,parent=None,status='queued',attempts=[],children=[])
        roots[rid]={'root_id':rid,'source114_id':spec['source114_id'],'status':'READY','new_reserved':0,'nodes_processed':0,'nodes':{rid:n},'queue':[rid],'passed_leaves':[],'empty_leaves':[],'blocked_leaves':[],'inflight':None}
    pool=cf.ThreadPoolExecutor(max_workers=14);fs=[pool.submit(run_root,s)for s in inp['roots']]
    for f in cf.as_completed(fs):results.append(f.result())
    pool.shutdown(wait=True);pool=None
    require(new_reserved<=768 and len(started)<=new_reserved and len(results)==96,'Finite kernel/root cap')
    H.new(P/'certificate_cache_index.json',{'entries':list(certcache.values()),'new_certificates':new_cert_count,'old_exp_recomputed':0})
    H.new(P/'forest.json',{'round':'117','inherited_passed_leaves':inp['inherited_passed_leaves'],'source114_forest_sha256':inp['source114_forest_sha256'],'roots':roots,'new_kernel_reservations':new_reserved,'new_kernel_evaluations_started':len(started),'old_native_replayed':0,'new_certificates':new_cert_count,'residual_frontier':[n for r in roots.values()for n in r['residual_frontier']],'whole_root_complete':all(r['root_complete']for r in roots.values()),'global_g_certificate':False,'new_C6_bound':False})
    status='COMPLETED_96_FINITE_ROOTS'
except BaseException:failure=traceback.format_exc()
finally:
    with lock:owned=list(live.values())
    for p in owned:H.stop_owned(p)
    if pool is not None:pool.shutdown(wait=True)
    H.atomic(P/'run_receipt.json',{'round':'117','status':status,'PID':os.getpid(),'started_utc':utc,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-tic,'new_kernel_reserved':new_reserved,'new_kernel_started':len(started),'roots_returned':len(results),'time_limit':None,'automatic_retries':0,'old_native_replayed':0,'failure':failure,'release_sha256':sha(rp),'global_g_certificate':False})
    print(json.dumps({'round':'117','status':status,'PID':os.getpid(),'new_kernel_started':len(started),'roots_returned':len(results),'elapsed_seconds':time.monotonic()-tic,'global_g_certificate':False}))
if status!='COMPLETED_96_FINITE_ROOTS':raise SystemExit(1)
