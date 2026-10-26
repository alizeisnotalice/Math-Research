"""One E123 analytic Gamma scalar configuration; no MC or old program."""
import sys,os
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import datetime,hashlib,json,math,time,traceback
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
start,cpu=time.monotonic(),time.process_time()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
if any((P/n).exists() for n in ['started.json','review_results.json','run_receipt.json']):raise SystemExit('one attempt only')
reg=read(P/'registration.json')
assert sha(Path(__file__))==reg['script_sha256']
assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
save('started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'registration_sha256':sha(P/'registration.json')})
count=0;status='FIRST_FAILURE_STOP';stage='source';failure=None
def ck(ok,label):
    global count
    assert ok,label
    count+=1
def root_bounds(x,n,bits=256):
    ck(x>0 and n>0,'positive exact root')
    D=1<<bits;N=x.numerator*D**n//x.denominator
    lo,hi=0,1<<((N.bit_length()+n-1)//n)
    while lo<hi:
        m=(lo+hi+1)//2
        if m**n<=N:lo=m
        else:hi=m-1
    a,b=Q(lo,D),Q(lo+1,D)
    ck(a**n<=x<=b**n,'exact outward root residual')
    return a,b
def strings(pair):return list(map(str,pair))
def coarse(pair,scale=10**12):
    return [str(Q((pair[0]*scale).__floor__(),scale)),str(Q((pair[1]*scale).__ceil__(),scale))]
try:
    pins={}
    for row in read(P/'source_manifest.json')['files']:
        p=Path(row['path']);ck(sha(p)==row['sha256'],'source hash');ck(p.stat().st_size==row['bytes'],'source bytes');pins[str(p)]=sha(p)
    A=Path(reg['author_directory']);man=read(A/'manifest.json')
    ck(sha(A/'manifest.json')==reg['author_manifest_sha256'],'author manifest')
    for row in man['files']:ck(sha(A/row['path'])==row['sha256'] and (A/row['path']).stat().st_size==row['bytes'],'complete author byte whitelist only')
    ar=read(A/'registration.json');mc=read(A/'results.json');mcrc=read(A/'run_receipt.json');rel=read(A/'primary_release.json')
    ck(ar['theorem_sha256']==rel['theorem_sha256']==sha(A/'theorem.md'),'frozen proof')
    ck(ar['source_sha256']==rel['source_sha256']==sha(A/'probe.py') and ar['input_sha256']==rel['input_sha256']==sha(A/'inputs.json'),'preserved MC inputs/source, not executed')
    ck(mc['status']==mcrc['status']=='COMPLETE_MC_DIAGNOSTIC_NOT_EXPECTATION_CERTIFICATE' and mc['n']==8192,'MC remains diagnostic metadata only')
    ck(not mcrc['new_C6_bound'] and not mcrc['certified_EZ_or_global_g0'],'no statistical/Cn endorsement')
    mb=read(Path(reg['main_matlab_results']))
    ck(mb['status']=='PASS_3_SYMBOLIC_NORMALIZATION_IDENTITIES' and mb['scope']=='symbolic normalization only; no expectation or C6 certification','MATLAB symbolic provenance only')
    stage='exact_normalization_and_Arb512'
    q=Q(6,5);g0=Q(1,65000)
    ck(Q(5,2)*q==3 and Q(-1,2)+Q(1,2)/q==Q(-1,12),'shape norm determinant exponent')
    ck(6*2*32==384 and Q(3,16)**6/Q(384)**5==Q(3,2**59),'whole-boundary normalization')
    ck(Q(6,32)==Q(3,16) and 32*Q(3,16)==6 and Q(3,16)/12==Q(1,64),'CB leverage and product exponents')
    ck(32*Q(1,64)==Q(1,2) and 16+Q(1,2)==Q(33,2) and Q(1,2)+Q(1,64)==Q(33,64),'Dirichlet beta parameters')
    Dcoef=Q(2*math.factorial(32),4**16*math.factorial(16)*math.factorial(15))
    ck(Dcoef==Q(300540195,67108864)==Q(mc['D_rational_times_sqrt_pi']),'exact radial Gamma factorial coefficient')
    threshold6=Q(8,9)**6*Q(24576)**5/Q(65000)**6
    ck(threshold6==Q(36893488147419103232,629200415790081024169921875)==Q(2)**83/(Q(3)**7*Q(65000)**6),'complete rational threshold')
    ck(Q(mc['threshold_M_power6_exact'])==threshold6==Q(mb['threshold_sixth_power_exact']) and Q(mb['coefficient_exact'])==Q(3,2**59),'source thresholds and MATLAB exact identities')
    # Only the previously unevaluated explicit analytic scalar baseline.
    import flint
    from flint import arb,fmpq,ctx
    ctx.prec=512
    def B(x):
        x=Q(x);return arb(fmpq(x.numerator,x.denominator))
    def ends(x):
        out=[]
        for side in [x.lower(),x.upper()]:
            m,e=side.man_exp();out.append(Q(int(m))*Q(2)**int(e))
        ck(0<out[0]<=out[1],'positive Arb scalar endpoints')
        return tuple(out)
    sqrtpi=arb.pi().sqrt()
    mg=2/B(32).sqrt()*(B(Q(33,64)).gamma()/sqrtpi)**32
    D=B(Dcoef)*sqrtpi
    mgI,DI=ends(mg),ends(D)
    ck(mgI[1]<DI[0]/32,'Gamma baseline tighter than deterministic D/32')
    crit=root_bounds(threshold6,6)
    cI=Q(3,2**59)
    cU=Q(6)**6*cI/g0**6
    ck(cU==Q(9,8)**6/(Q(384)**5*g0**6),'same-generator cost constant')
    ck(cU*threshold6==Q(64)**5,'exact competition threshold equivalence')
    Ibase=(root_bounds(cI*mgI[0]**6,5)[0],root_bounds(cI*mgI[1]**6,5)[1])
    Ubase=(root_bounds(cU*mgI[0]**6,5)[0],root_bounds(cU*mgI[1]**6,5)[1])
    target=root_bounds(Q(64)**5*(g0/6)**6,5)
    Sbase=(target[0]-Ibase[1],target[1]-Ibase[0])
    critZ=(crit[0]/DI[1],crit[1]/DI[0])
    gammaZ=(mgI[0]/DI[1],mgI[1]/DI[0])
    ck(crit[0]**6<=threshold6<=crit[1]**6,'exact threshold sixth-root enclosure')
    ck(Ibase[0]**5<=cI*mgI[0]**6<=cI*mgI[1]**6<=Ibase[1]**5,'complete analytic I upper coefficient')
    ck(Ubase[0]**5<=cU*mgI[0]**6<=cU*mgI[1]**6<=Ubase[1]**5,'conditional Gamma cost expression')
    if mgI[1]**6<threshold6:comparison='ANALYTIC_GAMMA_BASELINE_FITS_CONDITIONAL_MOMENT_THRESHOLD'
    elif mgI[0]**6>threshold6:comparison='ANALYTIC_GAMMA_BASELINE_STRICTLY_NONCOMPETITIVE'
    else:comparison='UNSEPARATED_AT_REGISTERED_PRECISION_STOP'
    ck(comparison!='UNSEPARATED_AT_REGISTERED_PRECISION_STOP','fixed scalar comparison separated')
    ck((Ubase[0]>64 and Sbase[1]<0) if comparison.endswith('NONCOMPETITIVE') else (Ubase[1]<64 and Sbase[0]>0),'threshold/cost/surplus directions agree')
    for p,digest in pins.items():ck(sha(Path(p))==digest,'source unchanged after analytic scalar run')
    intervals={'M_Gamma_analytic_upper_constant':mgI,'D_radial_factor':DI,'M_critical_at_unproved_g0':crit,'EZ_critical_at_unproved_g0':critZ,'Gamma_EZ_upper_constant':gammaZ,'I_Gamma_complete_upper_constant':Ibase,'U_Gamma_conditional_certificate_expression':Ubase,'S_Gamma_conditional_allocation_expression':Sbase,'critical_I_target_at_unproved_g0':target}
    save('review_results.json',{'status':'PASS_ONE_ANALYTIC_GAMMA_SCALAR_CONFIGURATION','comparison':comparison,'assertions':count,'arb_precision_bits':512,'flint_version':getattr(flint,'__version__','unavailable'),'exact_threshold_power6':str(threshold6),'exact_I_power5_coefficient':'3/2^59','exact_D_coefficient_times_sqrtpi':str(Dcoef),'scalar_intervals':{k:strings(v) for k,v in intervals.items()},'coarse_rational_enclosures':{k:coarse(v) for k,v in intervals.items()},'Gamma_baseline_is_actual_M_I_upper':True,'interval_lower_endpoints_are_not_actual_M_I_lower_bounds':True,'MC_statistical_or_expectation_certification':False,'MC_samples_replayed':0,'matlab_executed':False,'source_paths_pre_post':len(pins),'source_hashes':pins,'old_program_kernel_MC_native_optimization_MATLAB_children':0,'max_science_CPU_cores':1,'competitive_moment_and_globalg0_unpaid':True,'new_Cn':False,'wall_before_result_write':time.monotonic()-start,'CPU_before_result_write':time.process_time()-cpu})
    status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    save('run_receipt.json',{'status':status,'pid':os.getpid(),'assertions':count,'stage':stage,'failure':failure,'attempts':1,'automatic_retries':0,'children':0,'science_CPU_cores':1,'current_user_runtime_ceiling':None,'wall_before_final_receipt_write':time.monotonic()-start,'CPU_since_start_before_final_receipt_write':time.process_time()-cpu,'CPU_process_total_before_final_receipt_write':time.process_time(),'timing_scope':'one scalar process before final receipt/stdout/exit; complete tool roundtrip recorded separately','script_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json')})
    print(json.dumps({'status':status,'assertions':count,'pid':os.getpid()}))
if status!='COMPLETED':raise SystemExit(1)
