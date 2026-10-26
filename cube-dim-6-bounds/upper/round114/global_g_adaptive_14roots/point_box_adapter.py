"""C114 adapter extracted from pinned112v2 functions; no top-level launch."""
import copy,json,hashlib,math,os,subprocess,time,traceback
from pathlib import Path
from fractions import Fraction as Q

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
        intent={'round':'114','token':job['token'],'old_node_id':node,'continuation_node_id':job['continuation_node_id'],
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
        receipt={'schema':'C114-original-native-receipt-v1','token':job['token'],'old_node_id':node,
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
            'source_proof':f'C114_ORIGINAL_NATIVE:{folder}/native.json#sha256={receipt["native_sha256"]};receipt={sha(folder/"native_receipt.json")}',
            'provenance':{'kind':'C016_NATIVE_SAME_WINDOW_POSITIVE_BAND','precision_bits':512,'full_group_terms':967680,
                'computed_sign_pairs':483840,'synthetic_fixture':False,'native_sha256':receipt['native_sha256'],
                'receipt_sha256':sha(folder/'native_receipt.json'),'gradient_gate_sha256':pins['F061_manifest_sha256'],
                'independent_full_gradient_crosschecked_for_this_new_attempt':False}}
        base.validate_jet(jet,expected,test_only=False);new(folder/'jet.json',jet)
        result={'round':'114','old_node_id':node,'continuation_node_id':job['continuation_node_id'],
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
                    child=base.make_node(job['continuation_node_id']+suffix,job['continuation_node_id'],job['depth']+1,raw,0)
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
