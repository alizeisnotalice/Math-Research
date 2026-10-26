"""C115: exactly2 new low-band natives; saved high jets, six local certificates."""
import sys
sys.dont_write_bytecode=True
import concurrent.futures as cf
import copy,datetime,hashlib,importlib.util,json,os,subprocess,threading,time,traceback
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def require(x,s):
    if not x:raise ValueError(s)
def new(p,x):
    with Path(p).open('xb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
def atomic(p,x):
    p=Path(p);t=p.with_name(p.name+'.pending')
    with t.open('wb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
    os.replace(t,p)
def load(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def stop_owned(proc):
    if proc is not None and proc.poll()is None:
        proc.terminate()
        try:proc.wait(timeout=5)
        except subprocess.TimeoutExpired:proc.kill();proc.wait()
rp=P/'parent_release.json';require(rp.exists(),'C115 PREPARED_NOT_RELEASED; no mathematics/native')
rel=read(rp);reg=read(P/'registration.json');inputs=read(P/'inputs.json')
require(rel.get('authorized')is True and rel.get('round')=='115'
    and rel.get('script_sha256')==sha(Path(__file__))
    and rel.get('registration_sha256')==sha(P/'registration.json')
    and rel.get('max_native_calls')==2 and rel.get('max_parallel_native')==2,'Exact finite115 release')
for n,k in [('inputs.json','inputs_sha256'),('native_validator.py','validator_sha256'),('source_manifest.json','source_manifest_sha256')]:
    require(sha(P/n)==reg[k],'Prepared source changed: '+n)
require(not any((P/n).exists()for n in ['started.json','jobs','results.json','run_receipt.json','checkpoint.json']),
    'One115 pilot only; observe original handles, never reissue')
start=time.monotonic();startutc=datetime.datetime.now(datetime.timezone.utc).isoformat()
new(P/'started.json',{'round':'115','pid':os.getpid(),'started_utc':startutc,'release_sha256':sha(rp),
    'max_native':2,'program_time_limit':None})
lock=threading.RLock();checker_lock=threading.Lock();live={};states={x['id']:{'status':'NOT_STARTED'}for x in inputs['points']}
def update(pid,**data):
    with lock:
        states[pid].update(data);atomic(P/'checkpoint.json',{'round':'115','owner_PID':os.getpid(),'points':states,
            'release_sha256':sha(rp),'global_g_certificate':False})
def command(expected):
    e=base.encode_center(expected['center'],**inputs['fan_bounds']);a=Q(expected['window']['a']);T=Q(expected['window']['T'])
    return[inputs['native_paths']['binary'],'full',inputs['native_paths']['fan'],str(a.numerator),str(a.denominator),
        str(T.numerator),str(T.denominator),str(e['denominator']),*map(str,e['numerators'])]
def add(i,j):return(i[0]+j[0],i[1]+j[1])
def scale(i,c):return tuple(sorted(c*x for x in i))
def euler(F,g,pa,pt,expected):
    a=Q(expected['window']['a']);T=Q(expected['window']['T'])
    terms=[scale(F,Q(5)),scale(pt,T),scale(pa,-a)]+[scale(v,Q(x))for v,x in zip(g,expected['center'])]
    out=tuple(sum(t[k]for t in terms)for k in[0,1]);require(out[0]<=0<=out[1],'Merged exact Euler');return out
def jet_from(expected,F,g,source,provenance):
    return{**{k:expected[k]for k in['model','model_definition','window','center']},'F_interval':list(map(str,F)),
        'face_gradient_intervals':[list(map(str,v))for v in g[1:]],'source_proof':source,'provenance':provenance}
def box_certificate(expected,jet):
    # Explicit general-window schema; sum jet is not forged native metadata.
    data={**{k:expected[k]for k in['model','model_definition','window','center','face_box','target','quadratic']},'jet':jet}
    with checker_lock:
        c=checker.certify(data,precision_bits=384);base.check_saved_box(jet,expected,c)
    return c

def run_point(spec):
    pid=spec['id'];folder=P/'jobs'/pid;proc=None;tic=time.monotonic()
    try:
        high=Path(spec['high_source_folder']);eh=spec['high_input'];el=spec['low_input'];em=spec['merged_input']
        # Authenticate the actual saved high native/jet, never execute its command.
        hr=read(high/'native_receipt.json');hi=read(high/'intent.json');hn=read(high/'native.json');hj=read(high/'jet.json')
        require(hi['input']==eh and hr['input_sha256']==digest(eh)and hr['command']==command(eh),'Saved high command/input')
        require(hr['native_exit_code']==0 and (high/'stderr.txt').read_bytes()==b''
            and hr['native_sha256']==sha(high/'native.json')and hr['source_binary_fan_runtime_hashes']==reg['native_pins'],
            'Saved high original bytes/native pins')
        FH,gH,paH,ptH,neH,eH=validator.validate_native(hn,eh)
        require(list(map(str,FH))==hj['F_interval']and[list(map(str,v))for v in gH[1:]]==hj['face_gradient_intervals'],
            'Saved high jet differs from complete native intervals')
        base.validate_jet(hj,eh,test_only=False)
        folder.mkdir(parents=True,exist_ok=False);cmd=command(el);require(cmd==spec['low_command'],'Prepared low argv')
        intent={'round':'115','token':spec['token'],'input':el,'input_sha256':digest(el),'command':cmd,
            'release_sha256':sha(rp),'source_binary_fan_runtime_hashes':reg['native_pins'],'owned_wrapper_pid':os.getpid()}
        new(folder/'intent.json',intent);update(pid,status='LOW_INTENT_SAVED',token=spec['token'])
        with(folder/'native.json').open('xb')as out,(folder/'stderr.txt').open('xb')as err:
            proc=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=out,stderr=err,cwd=folder)
            with lock:live[pid]=proc
            new(folder/'original_process.json',{'token':spec['token'],'owned_native_pid':proc.pid,'owned_wrapper_pid':os.getpid(),
                'input_sha256':digest(el),'command':cmd,'intent_sha256':sha(folder/'intent.json')})
            update(pid,status='LOW_NATIVE_RUNNING',native_pid=proc.pid);code=proc.wait()
        with lock:live.pop(pid,None)
        receipt={'schema':'C115-general-positive-band-native-v1','token':spec['token'],'input_sha256':digest(el),
            'window':el['window'],'center':el['center'],'command':cmd,'native_exit_code':code,'owned_native_pid':proc.pid,
            'owned_wrapper_pid':os.getpid(),'native_sha256':sha(folder/'native.json'),'stderr_sha256':sha(folder/'stderr.txt'),
            'stderr_empty':(folder/'stderr.txt').stat().st_size==0,'intent_sha256':sha(folder/'intent.json'),
            'process_record_sha256':sha(folder/'original_process.json'),'source_binary_fan_runtime_hashes':reg['native_pins'],
            'release_sha256':sha(rp),'program_time_limit':None,'elapsed_to_native_exit':time.monotonic()-tic,'retries':0}
        new(folder/'native_receipt.json',receipt);require(code==0 and receipt['stderr_empty'],'Low native error/stderr')
        ln=read(folder/'native.json');FL,gL,paL,ptL,neL,eL=validator.validate_native(ln,el)
        common=(max(ptL[0],paH[0]),min(ptL[1],paH[1]));require(common[0]<=common[1],'Internal1/32 density mismatch')
        FM=add(FL,FH);gM=[add(l,h)for l,h in zip(gL,gH)];eM=euler(FM,gM,paL,ptH,em)
        jl=jet_from(el,FL,gL,f'C115_ORIGINAL_LOW_NATIVE:{folder}/native.json#sha256={receipt["native_sha256"]};receipt={sha(folder/"native_receipt.json")}',
            {'kind':'C016_NATIVE_SAME_WINDOW_POSITIVE_BAND','precision_bits':512,'full_group_terms':967680,
             'computed_sign_pairs':483840,'synthetic_fixture':False,'native_sha256':receipt['native_sha256'],
             'receipt_sha256':sha(folder/'native_receipt.json'),'gradient_gate_sha256':inputs['F061_manifest_sha256'],
             'independent_full_gradient_crosschecked_for_this_new_attempt':False})
        base.validate_jet(jl,el,test_only=False)
        jm=jet_from(em,FM,gM,f'C115_EXACT_ADJACENT_BAND_SUM:low={receipt["native_sha256"]};high={sha(high/"native.json")}',
            {'kind':'EXACT_SUM_OF_TWO_ADJACENT_POSITIVE_BANDS','low_window':el['window'],'high_window':eh['window'],
             'low_native_sha256':receipt['native_sha256'],'high_native_sha256':sha(high/'native.json'),
             'single_merged_native_claim':False})
        new(folder/'low_jet.json',jl);new(folder/'merged_jet.json',jm)
        cl=box_certificate(el,jl);ch=box_certificate(eh,hj);cm=box_certificate(em,jm)
        for name,c in [('low',cl),('high_from_saved_jet',ch),('merged',cm)]:new(folder/(name+'_box_certificate.json'),c)
        sep=(Q(cl['tangent_bound_lower'])+Q(ch['tangent_bound_lower']),Q(cl['tangent_bound_upper'])+Q(ch['tangent_bound_upper']))
        merged=(Q(cm['tangent_bound_lower']),Q(cm['tangent_bound_upper']))
        ordering=('SEPARATE_STRICTLY_LARGER_EXPRESSION'if sep[0]>merged[1]else
                  'MERGED_STRICTLY_LARGER_EXPRESSION'if merged[0]>sep[1]else'EXPRESSION_INTERVALS_OVERLAP')
        best=max(sep[0],merged[0]);target=Q(el['target'])
        result={'round':'115','id':pid,'source114_node':spec['source114_node'],'token':spec['token'],'native_pid':proc.pid,
            'low_native_calls':1,'high_or_merged_native_calls':0,'center':el['center'],'face_box':el['face_box'],
            'low_F_interval':list(map(str,FL)),'high_F_interval':list(map(str,FH)),'merged_F_interval':list(map(str,FM)),
            'low_six_gradients':[list(map(str,v))for v in gL],'high_six_gradients':[list(map(str,v))for v in gH],
            'merged_six_gradients':[list(map(str,v))for v in gM],'internal_density_intersection':list(map(str,common)),
            'internal_density_dyadics_bit_identical':ln['p_T']==hn['p_a'],'internal_endpoint_cancelled_analytically':True,
            'low_Euler':list(map(str,eL)),'high_saved_Euler':list(map(str,eH)),'merged_outer_endpoint_Euler':list(map(str,eM)),
            'separate_band_wholebox_expression_interval':list(map(str,sep)),
            'merged_jet_wholebox_expression_interval':list(map(str,merged)),'expression_ordering':ordering,
            'separate_over_merged_expression_interval':[str(sep[0]/merged[1]),str(sep[1]/merged[0])],
            'best_local_wholebox_lower':str(best),'local_wholebox_target_pass':best>=target,
            'merged_center_lower_above_target':FM[0]>target,'merged_center_upper_below_target':FM[1]<target,
            'corners':96,'global_g_certificate':False,'new_C6_bound':False,'status':'COMPLETED_POINT_PILOT'}
        new(folder/'verification.json',result);update(pid,status='COMPLETED_POINT_PILOT',native_pid=proc.pid,
            local_box_pass=result['local_wholebox_target_pass'],verification_sha256=sha(folder/'verification.json'))
        return result
    except BaseException:
        stop_owned(proc)
        with lock:live.pop(pid,None)
        fail={'round':'115','id':pid,'status':'OWN_POINT_ERROR_STOP_NO_RETRY','native_pid':None if proc is None else proc.pid,
            'low_native_calls':int(proc is not None),'error':traceback.format_exc(),'global_g_certificate':False}
        if folder.exists():atomic(folder/'failure.json',fail)
        update(pid,**fail);return fail

status='FIRST_FAILURE_STOP';failure=None;results=[];pool=None
try:
    for item in read(P/'source_manifest.json')['files']:require(sha(item['path'])==item['sha256'],'Source pin changed: '+item['path'])
    require(len(inputs['points'])==2,'Exactly two fixed point pilots')
    base=load('c115_frozen_geometry',inputs['base_path']);M=load('c115_frozen_dyadics',inputs['math_contract_path'])
    checker=load('c115_frozen_A057',inputs['checker_path']);validator=load('c115_general_native_validator',P/'native_validator.py');validator.M=M
    import flint
    require(flint.__FLINT_VERSION__=='3.6.0','Pinned FLINT runtime')
    (P/'jobs').mkdir();pool=cf.ThreadPoolExecutor(max_workers=2)
    futures=[pool.submit(run_point,s)for s in inputs['points']]
    for future in cf.as_completed(futures):results.append(future.result())
    pool.shutdown(wait=True);pool=None
    require(len(results)==2 and all(r['status']=='COMPLETED_POINT_PILOT'for r in results),'Point failure retained; no replacement')
    require(sum(r['low_native_calls']for r in results)==2,'Fixed2 low-native budget')
    new(P/'results.json',{'round':'115','status':'PASS_TWO_LOCAL_BAND_COMPARISON_PILOTS_PROOF_REVIEW_SEPARATE',
        'points':sorted(results,key=lambda r:r['id']),'low_native_calls':2,'high_native_calls':0,'union_native_calls':0,
        'labelled_corners':192,'global_g_certificate':False,'new_C6_bound':False});status='COMPLETED_TWO_POINT_PILOT'
except BaseException:failure=traceback.format_exc()
finally:
    with lock:owned=list(live.values())
    for proc in owned:stop_owned(proc)
    if pool is not None:pool.shutdown(wait=True)
    receipt={'round':'115','status':status,'pid':os.getpid(),'started_utc':startutc,
        'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-start,
        'program_time_limit':None,'states':states,'low_native_calls':sum(r['low_native_calls']for r in results),
        'saved_high_or_union_native_replays':0,'automatic_retries':0,'release_sha256':sha(rp),'failure':failure,
        'global_g_certificate':False}
    atomic(P/'run_receipt.json',receipt)
    print(json.dumps({'round':'115','status':status,'pid':os.getpid(),'low_native_calls':receipt['low_native_calls'],
        'points_returned':len(results),'elapsed_seconds':receipt['elapsed_seconds'],'global_g_certificate':False}))
if status!='COMPLETED_TWO_POINT_PILOT':raise SystemExit(1)
