"""E108: one independent saved-certificate check; no scientific replay."""
import datetime, hashlib, json, os, sys, time, traceback
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
P = Path(__file__).resolve().parent
start, cpu = time.monotonic(), time.process_time()
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def save(n,x): (P/n).write_text(json.dumps(x,indent=2)+'\n')
if any((P/n).exists() for n in ['started.json','review_results.json','run_receipt.json']):
    raise SystemExit('one attempt only')
reg=read(P/'registration.json')
assert sha(Path(__file__))==reg['script_sha256']
assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
save('started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'registration_sha256':sha(P/'registration.json'),'script_sha256':sha(Path(__file__))})
checks=0; completed=[]; status='FIRST_FAILURE_STOP'; failure=None
def check(b,label):
    global checks
    assert b,label
    checks+=1
def ends(d):
    z=F(2)**int(d['common_binary_exponent'])
    out=(F(d['lower_mantissa'])*z,F(d['upper_mantissa'])*z)
    check(out[0]<=out[1],'ordered saved dyadics')
    return out
def const(a): return (F(a),F(a))
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def sub(a,b): return (a[0]-b[1],a[1]-b[0])
def mul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)
def div(a,b):
    check(b[0]>0,'positive divisor')
    return mul(a,(1/b[1],1/b[0]))
def power(a,n):
    check(a[0]>=0,'nonnegative power')
    return a[0]**n,a[1]**n
def root_exact(a,n):
    check(a>=0,'nonnegative root')
    S=1<<96; N=a.numerator*S**n//a.denominator
    lo,hi=0,1<<((N.bit_length()+n-1)//n)
    while lo<hi:
        m=(lo+hi+1)//2
        if m**n<=N: lo=m
        else: hi=m-1
    out=F(lo,S),F(lo+1,S)
    check(out[0]**n<=a<=out[1]**n,'exact outward root residual')
    return out
def root(a,n): return root_exact(a[0],n)[0],root_exact(a[1],n)[1]
def coarse(a):
    S=1<<96
    return F((a[0]*S).__floor__(),S),F((a[1]*S).__ceil__(),S)
def includes(a,b): return a[0]<=b[0]<=b[1]<=a[1]
def json_intervals(d): return {k:list(map(str,v)) for k,v in d.items()}
def recursive_equal(o,target):
    if isinstance(o,dict): return o==target or any(recursive_equal(v,target) for v in o.values())
    if isinstance(o,list): return any(recursive_equal(v,target) for v in o)
    return False
try:
    src=read(P/'source_manifest.json')
    pre={}
    for row in src['files']:
        path=Path(row['resolved_path']); h=sha(path)
        check(h==row['sha256'],f'source pin {path}')
        check(path.stat().st_size==row['bytes'],f'source size {path}')
        pre[str(path)]=h
    author=Path(reg['author_directory'])
    manifest=read(author/'manifest.json')
    check(sha(author/'manifest.json')==reg['author_manifest_sha256'],'author final manifest')
    for row in manifest['file_whitelist']:
        check(sha(author/row['path'])==row['sha256'],'entire author whitelist')
        check((author/row['path']).stat().st_size==row['bytes'],'whitelist size')
    amap=read(author/'administrative_source_snapshots.json')
    upstream=read(author/'source_manifest.json')
    for row in amap['files']:
        old=next(x for x in upstream['files'] if x['path']==row['original_path'])
        check(old['sha256']==row['asof_sha256']==sha(author/row['snapshot']),'mutable asof map')
    result=read(author/'results.json'); a,b,c=result['configurations']
    receipt=read(author/'run_receipt.json'); tool=read(author/'actual_tool_receipt.json')
    release=read(author/'parent_release.json')
    check(result['status']=='PASS_THREE_CONFIGURATIONS_PROOF_AUDIT_SEPARATE','producer status')
    check(receipt['status']=='COMPLETED' and receipt['failure'] is None,'producer completion')
    check(receipt['completed_configurations']==3 and len(result['configurations'])==3,'producer count')
    check(receipt['pid']==tool['observed_PID']==50992 and tool['exit_code']==0,'producer PID exit')
    check(receipt['children_created']==receipt['retries']==0,'producer no child/retry')
    check(receipt['release_sha256']==sha(author/'parent_release.json'),'producer release binding')
    check(result['source_code_sha256']==sha(author/'run_prepared.py')==release['script_sha256'],'producer code')
    check(result['registration_sha256']==sha(author/'registration.json')==release['registration_sha256'],'producer registration')
    cert=read(author/'certificate.json')
    for k,v in cert['exact_outward_saved_objects'].items(): check(v==b[k],'certificate result binding')
    data=read(author/'saved104_inputs.json')
    old104=read(Path(reg['source104_results']))
    oldconstants=read(Path(reg['source104_constants']))
    check(sha(Path(reg['source104_results']))==data['source_results_sha256'],'104 results origin')
    check(sha(Path(reg['source104_constants']))==data['source_upstream_constants_sha256'],'104 constants origin')
    for name in ['core_full_norm','core_ten_integral','complete_middle_integral','old051short074long_middle_expression']:
        check(recursive_equal(old104,data[name]),'verbatim104 saved dyadic '+name)
    check(recursive_equal(oldconstants,data['paid_CEH_I_upper']),'verbatim paidCEH rational')
    completed.append('all author/upstream bytes and historical mutable map; producer evidence and complete104 inputs')

    # Verify the saved polynomial against the independent parity derivation.
    # No32-sign enumeration and no author configuration/program replay.
    expected={}
    for i in range(6):
        for j in range(i,6):
            expected[((i,i),(j,j))]=F(1 if i==j else 2)
    for i in range(6):
        for j in range(i+1,6): expected[((i,j),(i,j))]=F(4)
    stored={tuple(map(tuple,row['entries'])):F(row['coefficient']) for row in a['nonzero_quadratic_moment_coefficients']}
    check(stored==expected and len(stored)==len(a['nonzero_quadratic_moment_coefficients']),'saved quadratic parity polynomial')
    check(a['body32_frame_gram']==[[32*int(i==j) for j in range(6)] for i in range(6)],'saved frame identity')
    hist=a['zero_projection_witness_histogram']
    check(hist=={'0':4,'2':12,'4':12,'6':4},'binomial zero witness')
    check(sum(int(k)*v for k,v in hist.items())==96 and sum(int(k)**2*v for k,v in hist.items())==384,'saved witness moments')
    for r,row in enumerate(c['dimensions'],1):
        M=F(2)+F(2*r,3)
        check(row['dimension']==r and F(row['M_r'])==M,'rank coefficient parameter')
        check(F(row['exact_square_of_pi_free_peak_coefficient'])==M**r/F(32)**(2*r),'rank coefficient exact square')
        cr=cert['r1to6_exact_coefficients'][r-1]
        check(F(cr['square_of_pi_free_peak_coefficient'])==M**r/F(32)**(2*r),'certificate rank square')
    c6=F(216,2**30)
    check(F(c['full6_pi_free_coefficient'])==c6,'rank6 coefficient')
    check(c6/F(1,160)==F(c['full6_ratio_to_saved001_det160_peak'])==F(135,4194304),'rank6 oldpeak ratio')
    completed.append('saved parity coefficients, zero witness and all-r rational coefficients')

    # Independent coarse exact intervals, deliberately not a384-digit audit.
    pi=tuple(F(s) for s in read(Path(reg['pi_certificate']))['pi_interval'])
    check(3<pi[0]<pi[1]<F(22,7),'inherited exact pi bounds')
    lam=div(const(8),mul(const(3),mul(root(const(3),2),pi)))
    sharp=root(const(F(8,27)),2)
    fac=root(power(lam,6),5)
    qnew=div(const(1),mul(const(4096),power(pi,3)))
    independent={'Q3_new':qnew,'ratio_to098':lam,
        'conditional_sharpKhinchine_ratio_to098':sharp,
        'ratio_to_conditional_sharpKhinchine':div(lam,sharp),
        'integral_scaling_lambda_power6over5':fac}
    check(0<lam[0]<=lam[1]<F(1,2)<sharp[0]<1,'new vs old and sharp coefficient')
    base=root(const(F(1,390000)**6),5)
    target=mul(const(64),base); paid=const(F(data['paid_CEH_I_upper']))
    independent['one_core_full_G_norm']=mul(ends(data['core_full_norm']),lam)
    core=mul(ends(data['core_ten_integral']),fac)
    middle=mul(ends(data['complete_middle_integral']),fac)
    independent['ten_core_I_allocation']=core
    independent['complete31_middle_I_allocation']=middle
    independent['paidCEH_plus_ten_core_plan']=add(paid,core)
    independent['conditional_core_remaining_allowance']=sub(target,add(paid,core))
    independent['paidCEH_plus_complete_middle_plan']=add(paid,middle)
    independent['conditional_complete_middle_S_far_zero']=sub(target,add(paid,middle))
    independent['conditional_complete_middle_U_far_zero']=div(add(paid,middle),base)
    independent['new_over_old051short074long_middle_expression']=div(middle,ends(data['old051short074long_middle_expression']))
    independent={k:coarse(v) for k,v in independent.items()}
    for k,v in independent.items(): check(includes(v,ends(b[k])),'author saved result inside independently coarse enclosure '+k)
    for r,row in enumerate(c['dimensions'],1):
        coefficient=div(root(const((F(2)+F(2*r,3))**r),2),mul(const(32**r),power(pi,r)))
        check(includes(coarse(coefficient),ends(row['peak_coefficient_including_pi'])),'rank peak coarse interval')
    check(includes(coarse(div(const(c6),power(pi,6))),ends(c['full6_peak_including_pi'])),'full6 saved interval')
    check(independent['conditional_core_remaining_allowance'][0]>F(5,1000000),'core allowance strict positive')
    check(independent['conditional_complete_middle_S_far_zero'][1]<-F(2,10000),'middle allocation strict negative')
    check(independent['conditional_complete_middle_U_far_zero'][0]>1100,'conditional middle Uplan exceeds64')
    check(independent['new_over_old051short074long_middle_expression'][0]>2,'new middle still larger than old expression')
    check(b['g0_status']=='1/65000 unproved globally' and result['Cn_improved'] is False,'scope qualifications')
    completed.append('independent rational root bounds, all-r peaks and complete saved104 rescaling/signs')
    for path,h in pre.items(): check(sha(Path(path))==h,'source unchanged after review')
    output={'status':'PASS_INDEPENDENT_ANALYTIC_REVIEW_AND_ONE_SAVED_EXACT_AUDIT',
        'assertions':checks,'completed_checks':completed,'independent_root_bits':96,
        'fine_384_bit_digits_independently_certified':False,
        'independent_exact_intervals':json_intervals(independent),
        'source_pins_verified_pre_and_post':len(pre),'source_hashes':pre,
        'scientific_replays_kernel_native_grids_MATLAB_children':0,
        'theorem_proof_review':'proof_review.md','global_g_fullI_far_S_unpaid':True,
        'new_Cn_or_lower_project_work':False,'CPU_seconds_before_result_save':time.process_time()-cpu,
        'wall_seconds_before_result_save':time.monotonic()-start}
    save('review_results.json',output); status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    save('run_receipt.json',{'status':status,'pid':os.getpid(),'assertions':checks,
        'completed_checks':completed,'failure':failure,'attempts':1,'automatic_retries':0,
        'children':0,'current_user_runtime_ceiling':None,
        'wall_seconds_before_receipt_write':time.monotonic()-start,
        'CPU_seconds_since_start_before_receipt_write':time.process_time()-cpu,
        'CPU_process_total_before_receipt_write':time.process_time(),
        'timing_scope':'one Python process, before final receipt write/stdout/normal exit; outer tool stopwatch recorded separately',
        'source_script_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json')})
    print(json.dumps({'status':status,'assertions':checks,'completed_checks':len(completed),'pid':os.getpid()}))
if status!='COMPLETED': raise SystemExit(1)
