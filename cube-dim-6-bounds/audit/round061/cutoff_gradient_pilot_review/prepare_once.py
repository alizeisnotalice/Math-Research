"""Nine fixed F061 configurations; one build, two 128-tuple partials, no full run."""
import signal
import time
START = time.monotonic()
def timeout(signum, frame):
    raise TimeoutError('180-second preparation program wall limit')
signal.signal(signal.SIGALRM, timeout)
signal.setitimer(signal.ITIMER_REAL, 180)

from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
from math import factorial
import hashlib
import json
import os
import subprocess
import sys
import tarfile

P = Path(__file__).resolve().parent
REG = json.loads((P/'registration.json').read_text())
PIN = json.loads((P/'input_manifest.json').read_text())
SRC = Path(REG['author_directory'])
RECEIPT = P/'preparation_receipt.json'
if RECEIPT.exists() or (P/'preparation_results.json').exists():
    raise RuntimeError('Preparation is single-attempt; no replay')
receipt = {'status':'STARTED','pid':os.getpid(),'compile_attempts':0,'native_scalar_processes':0,
           'benchmark_processes':0,'partial_tuples_completed':0,'partial_tuple_field_contributions':0,
           'full_reference_cases':0,'author_native_calls':0,'old_baseline_reruns':0,
           'automatic_retries':0,'wall_timeout_seconds':180,'processes':[]}
child = None

def write_receipt():
    receipt['wall_seconds']=time.monotonic()-START
    temp=P/'preparation_receipt.json.tmp'
    temp.write_text(json.dumps(receipt,indent=2)+'\n')
    temp.replace(RECEIPT)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save_new(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')
def checkpoint(checks):
    temp=P/'preparation_checkpoint.json.tmp'
    temp.write_text(json.dumps({'completed_configurations':checks,'elapsed_seconds':time.monotonic()-START},indent=2)+'\n')
    temp.replace(P/'preparation_checkpoint.json')
def dyadic(raw):
    assert set(raw)=={'lower_mantissa','upper_mantissa','common_binary_exponent'}
    assert all(isinstance(x,str) for x in raw.values())
    scale=F(2)**int(raw['common_binary_exponent'])
    a,b=int(raw['lower_mantissa'])*scale,int(raw['upper_mantissa'])*scale
    assert a<=b
    return a,b

def process(command,label):
    global child
    out=P/(label+'_stdout.txt');err=P/(label+'_stderr.txt')
    entry={'label':label,'command':command,'started_elapsed_seconds':time.monotonic()-START}
    receipt['processes'].append(entry)
    with out.open('xb') as stdout,err.open('xb') as stderr:
        child=subprocess.Popen(command,cwd=P,stdout=stdout,stderr=stderr,start_new_session=True)
        entry['owned_pid']=child.pid;write_receipt()
        code=child.wait(timeout=max(0.001,180-(time.monotonic()-START)))
    entry.update(exit_code=code,elapsed_seconds=time.monotonic()-START-entry['started_elapsed_seconds'],
                 stdout_sha256=sha(out),stderr_sha256=sha(err))
    child=None;write_receipt()
    if code!=0:raise RuntimeError(f'{label} failed once: {err.read_text()}')
    return out

try:
    write_receipt()
    for path,digest in PIN['sha256'].items():assert sha(path)==digest,path
    manifest=read(SRC/'pilot_manifest.json')
    assert sha(SRC/'pilot_manifest.json')==REG['author_pilot_manifest_sha256']
    for name,digest in manifest['artifact_sha256'].items():assert sha(SRC/name)==digest,name
    for path,digest in read(SRC/'input_manifest.json')['source_sha256'].items():assert sha(path)==digest,path
    checks=[]

    # 1: exact saved-artifact comparison and all 128 saved corners; no kernel replay.
    saved=[]
    for case in ('low-primary','S3-primary'):
        folder=SRC/'production'/case
        native,jet,rec,result=(read(folder/name) for name in ('native.json','jet.json','receipt.json','result.json'))
        assert rec['native_sha256']==sha(folder/'native.json') and rec['jet_sha256']==sha(folder/'jet.json')
        assert rec['status']=='COMPLETED_NO_RETRY' and rec['native_exit_code']==0 and rec['stderr_empty'] is True
        assert native['status']=='all_real_nodes_supported' and native['precision_bits']==512
        assert native['full_group_terms']==967680 and native['computed_sign_pairs']==483840
        assert rec['window']==jet['window']==REG['window'] and not jet['synthetic_fixture']
        assert rec['model']==jet['model']==REG['model'] and rec['model_definition']==jet['model_definition']
        center=list(map(F,jet['center']));assert center==list(map(F,rec['center']))
        command=rec['command'];assert command[3:7]==['1','250','1','2']
        assert [F(int(v),int(command[7])) for v in command[8:]]==center
        flo,fhi=dyadic(native['band']);assert [flo,fhi]==list(map(F,jet['F_interval']))
        gradient=[dyadic(x) for x in native['gradient']];assert len(gradient)==6
        assert [list(g) for g in gradient[1:]]==[list(map(F,g)) for g in jet['face_gradient_intervals']]
        assert flo>F(1,65000) and dyadic(native['p_a'])[0]>0 and dyadic(native['p_T'])[0]>0
        assert dyadic(native['Euler_xgrad_plus5B_plusTpT_minusapa'])[0]<=0<=dyadic(native['Euler_xgrad_plus5B_plusTpT_minusapa'])[1]
        boxrows=[]
        for width_index in (1,2):
            item=read(folder/f'box_{width_index}.json');cert=item['certificate'];width=F(1,250*2**width_index)
            box=[[max(F(0),c-width),min(F(1),c+width)] for c in center[1:]]
            assert item['halfwidth_index']==width_index and F(item['halfwidth_exact'])==width
            assert box==[list(map(F,pair)) for pair in item['face_box']]
            assert cert['model']==REG['model'] and cert['window']==REG['window'] and cert['quadratic']=='frame32'
            assert cert['conditional_on_jet_source']==jet['source_proof'] and cert['precision_bits']==384
            assert cert['target']=='1/65000' and len(cert['corners'])==32
            seen=set();low=[];high=[]
            for row in cert['corners']:
                bits=tuple(row['corner_bits']);assert bits not in seen and bits in tuple(product((0,1),repeat=5));seen.add(bits)
                d=[box[j][bits[j]]-center[j+1] for j in range(5)]
                b=sum(min(g[0]*v,g[1]*v) for g,v in zip(gradient[1:],d));f=min(fhi,max(flo,b));q=sum(v*v for v in d)/32;e=b/f-q/F(1,250)**2
                assert d==list(map(F,row['displacement'])) and b==F(row['gradient_dot_lower'])
                assert f==F(row['minimizing_center_value']) and q==F(row['quadratic']) and e==F(row['exponent'])
                lo,hi=F(row['certificate_lower']),F(row['certificate_upper']);assert lo<=hi
                assert row['passes_target']==(lo>=F(1,65000));low.append(lo);high.append(hi)
            assert min(low)==F(cert['tangent_bound_lower']) and min(high)==F(cert['tangent_bound_upper'])
            expected='PASS_CONDITIONAL_ON_MATCHING_CENTER_JET' if min(low)>=F(1,65000) else 'UNRESOLVED'
            assert cert['status']==expected and (width_index==2)==(expected.startswith('PASS'))
            boxrows.append({'width':str(width),'corners':32,'status':expected,'lower_exact':str(min(low))})
        saved.append({'case_id':case,'dyadic_to_jet_mapping_exact':True,'same_window':True,'six_physical_fields':6,'boxes':boxrows})
    checks.append({'configuration':1,'status':'PASS_SAVED_ENDPOINTS_AND_128_RATIONAL_CORNERS','cases':saved,'new_exponential_evaluations':0})
    checkpoint(checks)

    import sympy as s
    import flint
    from flint import arb,ctx,fmpq
    ctx.prec=1024
    def ball(q):q=F(q);return arb(fmpq(q.numerator,q.denominator))
    def ends(v):return F(str(v.lower().fmpq())),F(str(v.upper().fmpq()))
    z=s.symbols('z',real=True);t=s.symbols('t',positive=True);beta=s.symbols('beta',real=True)
    # 2: independently verify every time primitive by symbolic differentiation.
    primitive_rows=[]
    for r in range(1,8):
        integrand=s.expand_complex(s.I**(r-6)/(t-s.I*beta)**r).as_real_imag()[0]
        primitive=-s.atan(beta/t) if r==1 else s.expand_complex(-s.I**(r-6)*(t-s.I*beta)**(1-r)/(r-1)).as_real_imag()[0]
        assert s.cancel(s.diff(primitive,t)-integrand)==0
        primitive_rows.append({'pole_order':r,'derivative_identity':True,'zero_pole_integrand':str(s.simplify(integrand.subs(beta,0)))})
    checks.append({'configuration':2,'status':'PASS_TIME_PRIMITIVES','orders':primitive_rows})
    checkpoint(checks)

    # Only one compile attempt, pinned headers/runtime, no author executable.
    build=read(PIN['build_record_path']);archive=Path(build['headers_archive'])
    assert sha(archive)==build['headers_sha256']
    scratch=P/'scratch';scratch.mkdir()
    with tarfile.open(archive) as stream:stream.extractall(scratch,filter='data')
    lib=Path(next(path for path in build['runtime_dylib_sha256'] if 'libflint.' in path))
    binary=scratch/'physical_pf_reference'
    command=['/usr/bin/clang++','-O3','-std=c++20','-I',str(scratch/'include'),str(P/'physical_pf_reference.cpp'),str(lib),'-Wl,-rpath,'+str(lib.parent),'-o',str(binary)]
    receipt['compile_attempts']=1;write_receipt();process(command,'compile')
    process(['/usr/bin/install_name_tool','-change','/DLC/flint/.dylibs/libflint.24.0.dylib','@rpath/libflint.24.0.dylib',str(binary)],'binary_patch')
    receipt['binary_sha256']=sha(binary)
    receipt['native_scalar_processes']=1;write_receipt()
    native_scalar=read(process([str(binary),'scalars'],'scalar_native'))
    assert native_scalar['precision_bits']==768 and len(native_scalar['rows'])==5

    # 3--7: separate exact SymPy apart oracle for all six physical directions.
    sets=[list(map(F,['-3/8','-1/9','0','2/7','3/5','5/6'])),list(map(F,['0','0','1/7','1/7','-2/5','3/4'])),[F(0)]*6,[F(2,11)]*6,list(map(F,['-2/13','-2/13','-2/13','0','0','5/17']))]
    a,T=F(1,250),F(1,2);primitive_cache={}
    def phi(b,r):
        key=(b,r)
        if key not in primitive_cache:
            if r==1:primitive_cache[key]=(ball(b)/ball(a)).atan()-(ball(b)/ball(T)).atan()
            else:
                aa,tt,bb=map(lambda q:s.Rational(q.numerator,q.denominator),(a,T,b))
                exact=s.simplify(s.expand_complex(s.I**(r-6)*((aa-s.I*bb)**(1-r)-(tt-s.I*bb)**(1-r))/(r-1)).as_real_imag()[0])
                primitive_cache[key]=ball(F(str(exact)))
        return primitive_cache[key]
    for fixture,base in enumerate(sets):
        count=Counter(base);symbase=[s.Rational(b.numerator,b.denominator) for b in base]
        denominator=s.prod(z-b for b in symbase);comparisons=[]
        for k in range(6):
            H=s.cancel(sum(s.Rational((-1)**(j+k),j+k+1)/(z-b) for j,b in enumerate(symbase))/denominator)
            apart=s.apart(H,z);assert s.cancel(apart-H)==0
            value=arb(0);coefficients=[]
            for b,m in sorted(count.items()):
                sb=s.Rational(b.numerator,b.denominator)
                analytic=s.cancel(apart*(z-sb)**(m+1))
                for r in range(1,m+2):
                    n=m+1-r
                    coefficient=s.cancel(s.diff(analytic,z,n).subs(z,sb)/factorial(n))
                    q=F(str(coefficient));value+=ball(q)*phi(b,r)
                    if q:coefficients.append([str(b),r,str(q)])
            lo,hi=ends(value);nl,nu=dyadic(native_scalar['rows'][fixture]['gradient'][k]);assert nl<=lo<=hi<=nu,(fixture,k)
            comparisons.append({'physical_field':k+1,'reference_lower_exact':str(lo),'reference_upper_exact':str(hi),'native768_contains_reference1024':True,'exact_PF_coefficients':coefficients})
        checks.append({'configuration':fixture+3,'status':'PASS_NEW_SCALAR_SIX_FIELDS','nodes':list(map(str,base)),'synthetic_only':True,'fields':comparisons})
        checkpoint(checks)

    # 8--9: exactly 128 real tuples per case. Save reusable normalized partial sums.
    fan_path=read(SRC/'registration.json')['fan_path'];benchmarks=[]
    for configuration,case in ((8,'low-primary'),(9,'S3-primary')):
        receipt['benchmark_processes']+=1;write_receipt()
        result=read(process([str(binary),'benchmark',case,fan_path],'benchmark_'+case))
        assert result['computed_tuples']==128 and result['physical_fields']==6 and result['tuple_field_contributions']==768
        assert not result['full_reference'] and result['precision_bits']==768
        for raw in result['gradient_partial_normalized']:dyadic(raw)
        save_new(P/('benchmark_'+case+'.json'),result)
        receipt['partial_tuples_completed']+=128;receipt['partial_tuple_field_contributions']+=768;write_receipt()
        estimate=result['kernel_seconds']*967552/128
        benchmarks.append({'case_id':case,'computed_tuples':128,'tuple_field_contributions':768,'seconds':result['kernel_seconds'],'linear_remaining_seconds_estimate':estimate,'not_a_runtime_bound':True,'zero_base_tuples':result['zero_base_tuples'],'repeated_base_tuples':result['repeated_base_tuples'],'reuse_saved_partial':True})
        checks.append({'configuration':configuration,'status':'PASS_REGISTERED_PARTIAL_BENCHMARK','summary':benchmarks[-1]})
        checkpoint(checks)
    assert len(checks)==9 and receipt['partial_tuple_field_contributions']==1536
    for path,digest in PIN['sha256'].items():assert sha(path)==digest,path
    out={'status':'PASS_NINE_PREPARATION_CONFIGURATIONS_NO_FULL_REFERENCE_RELEASE','checks':checks,'compile_attempts':1,'native_scalar_processes':1,'new_scalar_configurations':5,'synthetic_scalar_fields':30,'real_partial_tuples':256,'real_partial_tuple_field_contributions':1536,'complete_reference_cases':0,'complete_reference_components':0,'partial_real_kernel_work':True,'author_native_calls':0,'old_reference_reruns':0,'precision_bits_native':768,'precision_bits_scalar_oracle':1024,'benchmarks':benchmarks,'binary_sha256':sha(binary),'source_sha256':sha(P/'physical_pf_reference.cpp'),'python':sys.version,'sympy':s.__version__,'python_flint':flint.__version__,'wall_seconds':time.monotonic()-START,'source_hashes_unchanged':True,'full_gradient_gap_closed':False,'global_g_certificate':False,'new_C6_bound':False}
    save_new(P/'preparation_results.json',out)
    receipt.update(status='PASS_PREPARATION_ONLY_WAITING_FULL_REFERENCE_RELEASE',results_sha256=sha(P/'preparation_results.json'))
    print(json.dumps({k:out[k] for k in ['status','real_partial_tuples','real_partial_tuple_field_contributions','complete_reference_cases','wall_seconds','benchmarks']},indent=2))
except BaseException as error:
    signal.setitimer(signal.ITIMER_REAL,0)
    if child is not None and child.poll() is None:
        try:os.killpg(child.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        child.wait();receipt['killed_and_waited_only_owned_process_group']=True
    receipt.update(status='FAILED_OR_TIMEOUT_STOP_NO_RETRY',error=repr(error))
    raise
finally:
    signal.setitimer(signal.ITIMER_REAL,0)
    write_receipt()
